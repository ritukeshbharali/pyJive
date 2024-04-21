import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn
from names import PropNames as prn

from models.framemodel import FrameModel
from models.framemodel import MP, DOFTYPES, LINEAR, NONLIN

from utils.node import Node
from utils.xtable import XTable

from scipy.optimize import root_scalar

import sympy as sp


class PlasticFrameModel(FrameModel):
    def take_action(self, action, params, globdat):

        if action == act.GETMATRIX0:
            self._get_matrix(params, globdat)
        elif action == act.COMMIT:
            self._commit()
        elif action == act.GETTABLE:
            if 'stress' in params[pn.TABLENAME]:
                self._get_stress_table(params, globdat)
        else:
            super().take_action(action,params,globdat)

    def configure(self, props, globdat):
        super().configure(props,globdat)

        alpha = sp.symbols('alpha')
        mp  = sp.sympify(str(props.get(MP)))
        der = sp.diff(mp,alpha)

        self._mp = sp.lambdify(alpha,mp,'numpy')
        self._mpder = sp.lambdify(alpha,der,'numpy')

        self._old_alpha = np.zeros(len(self._elems)*self._ipcount)
        self._new_alpha = np.zeros(len(self._elems)*self._ipcount)
        
        self._old_kappap = np.zeros(len(self._elems)*self._ipcount)
        self._new_kappap = np.zeros(len(self._elems)*self._ipcount)

    def _get_matrix(self, params, globdat):
        ipoint = 0

        for elem in self._elems:
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES)
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=0)[:, :]

            d0 = coords[1, :] - coords[0, :]
            phi = np.arctan2(d0[1], d0[0])
            l_0 = np.linalg.norm(d0)
            coords1d = np.array([0, l_0])

            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients([coords1d])
            elmat = np.zeros((6, 6))
            elfor = np.zeros(6)

            ue = [globdat[gn.STATE0][i] for i in idofs]

            if self._subtype == LINEAR: 
                for ip in range(self._ipcount):
                    N = sfuncs[:, ip]
                    dN = grads[:, 0, ip]


                    B = self._get_B_matrix(N=N, dN=dN, omega=phi)

                    evec = np.matmul(B,ue)

                    svec, D = self._update (evec,ipoint)

                    elmat += weights[ip] * np.matmul(B.transpose(), np.matmul(D, B))
                    elfor += weights[ip] * np.matmul(B.transpose(), svec)

                    ipoint += 1

            elif self._subtype == NONLIN:
                if self._shape.node_count() > 2:
                    raise NotImplementedError('nonlinear strain only implemented for 2-node element')


                d = d0 + ue[3:5] - ue[0:2]
                lcps = (d0[0] * d[0] + d0[1] * d[1]) / l_0
                lsps = (d0[0] * d[1] - d0[1] * d[0]) / l_0

                for ip in range(self._ipcount):
                    N = sfuncs[:, ip]
                    dN = grads[:, 0, ip]

                    theta = np.matmul(N, ue[2::3])
                    kappa = np.matmul(dN, ue[2::3])
                    omega = phi + theta
                    gamma = (np.cos(theta) * lsps - np.sin(theta) * lcps) / l_0
                    eps = (np.cos(theta) * lcps + np.sin(theta) * lsps) / l_0 - 1
                    evec = [eps, gamma, kappa]

                    svec, D = self._update (evec,ipoint) 

                    B = self._get_B_matrix(N=N, dN=dN, omega=omega, gamma=gamma, eps=eps)
                    WN = self._get_WN_matrix(N=N, dN=dN, omega=omega, eps=eps)
                    WV = self._get_WV_matrix(N=N, dN=dN, omega=omega, gamma=gamma)
                    Kmat = weights[ip] * np.matmul(B.transpose(), np.matmul(D, B))
                    Kgeo = svec[0] * WN + svec[1] * WV  # *l0/l0

                    elmat += Kmat + Kgeo
                    elfor += weights[ip] * np.matmul(B.transpose(), svec)

                    ipoint += 1

            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elmat
            params[pn.INTFORCE][np.ix_(idofs)] += elfor

    def _get_stress(self, globdat, disp):
        stressmat = np.zeros((len(self._elems), self._strcount, 2))
        # D = self._get_D_matrix()

        ipoint = 0

        for i, elem in enumerate(self._elems):
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES)
            coords = np.stack([globdat[gn.NSET][j].get_coords() for j in inodes], axis=0)[:, :]

            d0 = coords[1, :] - coords[0, :]
            phi = np.arctan2(d0[1], d0[0])
            l_0 = np.linalg.norm(d0)
            coords1d = np.array([0, l_0])

            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients([coords1d])
            m1 = 0
            m2 = 0

            ue = [disp[j] for j in idofs]

            d = d0 + ue[3:5] - ue[0:2]
            lcps = (d0[0] * d[0] + d0[1] * d[1]) / l_0
            lsps = (d0[0] * d[1] - d0[1] * d[0]) / l_0

            assert(self._ipcount == 1) # TODO: this is ugly
            ip = 0

            if self._subtype == LINEAR:
                N = sfuncs[:, ip]
                dN = grads[:, 0, ip]

                B = self._get_B_matrix(N=N, dN=dN, omega=phi)
                evec = np.matmul(B, ue)

                svec, _ = self._update(evec,ipoint)

                m1 += weights[ip] * np.dot(B[:,2], svec)
                m2 += weights[ip] * np.dot(B[:,5], svec)

                ipoint += 1

            elif self._subtype == NONLIN:
                N = sfuncs[:, ip]
                dN = grads[:, 0, ip]

                theta = np.matmul(N, [ue[2], ue[5]])
                kappa = np.matmul(dN, [ue[2], ue[5]])
                gamma = (np.cos(theta) * lsps - np.sin(theta) * lcps) / l_0
                eps = (np.cos(theta) * lcps + np.sin(theta) * lsps) / l_0 - 1
                evec = [eps, gamma, kappa]

                svec, _ = self._update(evec,ipoint)

                g = gamma / 2
                e = (1 + eps) / 2
                t = 1 / l_0

                m1 += weights[ip] * np.matmul([g, -e, -t], svec)
                m2 += weights[ip] * np.matmul([g, -e, t], svec)

                ipoint += 1

            svec[2] = -m1;
            stressmat[i,:,0] = svec;
            svec[2] = m2;
            stressmat[i,:,1] = svec;
        return stressmat

    def _get_stress_table(self, params, globdat):
        table = params[pn.TABLE]

        # Convert the table to an XTable and store the original class
        cls_ = table.__class__
        table.__class__ = XTable

        # Add the columns of all stress components to the table
        jcols = table.add_columns(['N', 'V', 'M'])

        if gn.STATE0 in params:
            disp = params[gn.STATE0]
        else:
            disp = globdat[gn.STATE0]

        smat = self._get_stress(globdat, disp)

        for jcol in jcols:
            table.add_col_values(None, jcol, smat[:,jcol,:].flatten())

        # Convert the table back to the original class
        table.__class__ = cls_

        params[pn.TABLEWEIGHTS] = np.ones(table['N'].shape)


    def _update(self,strain,ipoint):

        D = np.diag([self._EA, self._GAs, self._EI])

        stress = np.matmul(D,strain)

        kappa_elastic = strain[2] - self._old_kappap[ipoint]

        self._m_trial = self._EI * kappa_elastic

        trial_f = np.abs(self._m_trial) - self._mp(self._old_alpha[ipoint])

        tol = 1.e-10

        if trial_f < -tol:
            self._new_kappap[ipoint] = self._old_kappap[ipoint]
            self._new_alpha[ipoint]  = self._old_alpha[ipoint]
            stress[2] = self._m_trial
        else:
            sol = root_scalar(self._yield_function,
                              args = (ipoint,True),
                              fprime = True,
                              rtol = tol,
                              maxiter = 1000,
                              method = 'newton',
                              x0 = 0.0)

            if not sol.converged or sol.root < -tol:
                upper = tol
                while self._yield_function(upper,ipoint)[0] > 0.0:
                    upper *= 2

                print('Trying bisection with upper bound',upper)
                sol = root_scalar(self._yield_function,
                                  args = (ipoint,False),
                                  rtol = tol,
                                  maxiter = 1000,
                                  method = 'bisect',
                                  bracket = [0,upper])

                if not sol.converged:
                    raise RuntimeError('Return mapping failed')

            ddkappa = self._EI**2 / (self._EI + self._mpder(self._new_alpha[ipoint]))

            stress[2] = self._m
            D[2,2]    = self._EI - ddkappa

        return stress, D

    def _yield_function(self,dgam,ipoint,fprime=True):
        n = np.sign(self._m_trial)

        self._new_alpha[ipoint] = self._old_alpha[ipoint] + dgam
        self._new_kappap[ipoint]= self._old_kappap[ipoint] + n * dgam

        self._m = self._m_trial - self._EI * n * dgam

        f = np.abs(self._m) - self._mp(self._new_alpha[ipoint])
        der = -self._EI * n * np.sign(self._m) - self._mpder(self._new_alpha[ipoint])

        if fprime:
            return f, der
        else:
            return f

    def _commit(self):
        self._old_kappap = np.copy(self._new_kappap)
        self._old_alpha  = np.copy(self._new_alpha)

def declare(factory):
    factory.declare_model('PlasticFrame', PlasticFrameModel)
