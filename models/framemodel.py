import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn
from names import PropNames as prn

from models.model import Model

from utils.node import Node
from utils.xtable import XTable

ELEMENTS = 'elements'
SUBTYPE = 'subtype'
LINEAR = 'linear'
NONLIN = 'nonlin'
EA = 'EA'
GAs = 'GAs'
EI = 'EI'
RHOA = "rhoA"
RHOI = "rhoI"
PLASTIC = 'plastic'
MP = 'Mp'
SHAPE = 'shape'
INTSCHEME = 'intScheme'
DOFTYPES = ['dx', 'dy', 'phi']


class FrameModel(Model):
    def take_action(self, action, params, globdat):

        if action == act.GETMATRIX0:
            self._get_matrix(params, globdat)
        elif action == act.GETMATRIX2:
            self._get_mass_matrix(params, globdat)
        elif action == act.GETMATRIXLB:
            self._get_matrix_lb(params, globdat)
        elif action == act.CHECKCOMMIT:
            self._check_commit(params, globdat)
        elif action == act.GETTABLE:
            if 'stress' in params[pn.TABLENAME]:
                self._get_stress_table(params, globdat)

    def configure(self, props, globdat):
        self._subtype = str(props[SUBTYPE])
        self._EA = float(props[EA])
        self._GAs = float(props[GAs])
        self._EI = float(props[EI])
        self._rhoA = float(props.get(RHOA, 0))
        self._rhoI = float(props.get(RHOI, 0))

        self._plastic = bool(eval(props.get(PLASTIC, 'False')))
        if self._plastic:
            self._mp = float(props.get(MP, None))
            globdat[pn.HINGENODES] = []
            globdat[pn.HINGESTEPS] = []
        self._nhinges = 0
        self._hingedofs = []
        self._hingemoments = []
        self._hingenodes = []

        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][prn.TYPE], props[SHAPE][INTSCHEME])
        egroup = globdat[gn.EGROUPS][props[ELEMENTS]]
        self._elems = [globdat[gn.ESET][e] for e in egroup]

        self._ipcount = self._shape.ipoint_count()
        self._dofcount = 3 * self._shape.node_count()
        self._strcount = 3

        nodes = np.unique([node for elem in self._elems for node in elem.get_nodes()])
        for doftype in DOFTYPES:
            globdat[gn.DOFSPACE].add_type(doftype)
            for node in nodes:
                globdat[gn.DOFSPACE].add_dof(node, doftype)

    def _get_matrix(self, params, globdat):
        D = self._get_D_matrix()
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
                    elmat += weights[ip] * np.matmul(B.transpose(), np.matmul(D, B))
                    elfor += np.matmul(elmat, ue)

            elif self._subtype == NONLIN:
                if self._shape.node_count() > 2:
                    raise NotImplementedError('nonlinear strain only implemented for 2-node element')

                # TODO: use shape functions for evaluating psi and l_ locally

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
                    svec = np.matmul(D, evec)

                    B = self._get_B_matrix(N=N, dN=dN, omega=omega, gamma=gamma, eps=eps)
                    WN = self._get_WN_matrix(N=N, dN=dN, omega=omega, eps=eps)
                    WV = self._get_WV_matrix(N=N, dN=dN, omega=omega, gamma=gamma)
                    Kmat = weights[ip] * np.matmul(B.transpose(), np.matmul(D, B))
                    Kgeo = svec[0] * WN + svec[1] * WV  # *l0/l0

                    elmat += Kmat + Kgeo
                    elfor += weights[ip] * np.matmul(B.transpose(), svec)

            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elmat
            params[pn.INTFORCE][np.ix_(idofs)] += elfor

        params[pn.INTFORCE][np.ix_(self._hingedofs)] += self._hingemoments

    def _get_matrix_lb(self, params, globdat):
        D = self._get_D_matrix()
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
            elmatM = np.zeros((6, 6))
            elmatG = np.zeros((6, 6))

            for ip in range(self._ipcount):
                N = sfuncs[:, ip]
                dN = grads[:, 0, ip]

                B = self._get_B_matrix(N=N, dN=dN, omega=phi)

                ue = [globdat[gn.STATE0][i] for i in idofs]
                evec = np.matmul(B, ue)
                svec = np.matmul(D, evec)

                WN = self._get_WN_matrix(N=N, dN=dN, omega=phi, eps=0)
                WV = self._get_WV_matrix(N=N, dN=dN, omega=phi, gamma=0)

                elmatM += weights[ip] * np.matmul(B.transpose(), np.matmul(D, B))
                elmatG += svec[0] * WN + svec[1] * WV  # *l0/l0

            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elmatM
            params[pn.MATRIX1][np.ix_(idofs, idofs)] += elmatG

    def _get_mass_matrix(self, params, globdat):
        M_disp = np.array([[self._rhoA]])
        M_rot = np.array([[self._rhoI]])
        
        for elem in self._elems:
            inodes = elem.get_nodes()
            idofs_dx = globdat[gn.DOFSPACE].get_dofs(inodes, [DOFTYPES[0]])
            idofs_dy = globdat[gn.DOFSPACE].get_dofs(inodes, [DOFTYPES[1]])
            idofs_phi = globdat[gn.DOFSPACE].get_dofs(inodes, [DOFTYPES[2]])
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=0)[:, :]
            
            d0 = coords[1, :] - coords[0, :]
            phi = np.arctan2(d0[1], d0[0])
            l_0 = np.linalg.norm(d0)
            coords1d = np.array([0, l_0])
            
            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords1d)
            elmat_dx = np.zeros((2, 2))
            elmat_dy = np.zeros((2, 2))
            elmat_phi = np.zeros((2, 2))
            
            for ip in range(self._ipcount):
                N = np.array([sfuncs[:, ip]])
                elmat_dx += weights[ip] * np.matmul(N.transpose(), np.matmul(M_disp, N))
                elmat_dy += weights[ip] * np.matmul(N.transpose(), np.matmul(M_disp, N))
                elmat_phi += weights[ip] * np.matmul(N.transpose(), np.matmul(M_rot, N))
                # elmat dx and dy are identical, could remove one
            
            params[pn.MATRIX2][np.ix_(idofs_dx, idofs_dx)] += elmat_dx
            params[pn.MATRIX2][np.ix_(idofs_dy, idofs_dy)] += elmat_dy
            params[pn.MATRIX2][np.ix_(idofs_phi, idofs_phi)] += elmat_phi
            
    def _check_commit(self, params, globdat):
        globdat[gn.ACCEPTED] = True
        if self._plastic:
            if self._mp is None:
                raise RuntimeError('Plastic moment ' + MP + ' has not been defined')

            smat = self._get_stress(globdat, globdat[gn.STATE0])
            moments = smat[:,2,:].reshape(-1,1)

            # Get element and node where highest moment occurs
            maxarg = np.argmax(abs(moments))
            maxratio = np.max(abs(moments)) / self._mp
            hingeelem = maxarg // 2
            enodes = globdat[gn.ESET][hingeelem].get_nodes()
            hingenode = enodes[maxarg % 2]
            if (maxarg%2==0) ^ (moments[maxarg]>0):
                sign = 1
            else:
                sign = -1

            if maxratio > 1.0001:
                globdat[gn.ACCEPTED] = False
                self._add_plastic_hinge(globdat, params, hingenode, hingeelem, sign)

    def _add_plastic_hinge(self, globdat, params, hingenode, hingeelem, sign):
        print('Adding plastic hinge on node %i (in element %i)' % (hingenode, hingeelem))

        # Add node
        oldnode = hingenode
        coords = globdat[gn.NSET][oldnode].get_coords()
        globdat[gn.NSET].append(Node(coords))
        newnode = len(globdat[gn.NSET]) - 1
        globdat[gn.MASTERNODES].append(newnode)

        # Duplicate dofs and add new phi dof
        globdat[gn.DOFSPACE].set_dof(oldnode, newnode, 'dx')
        globdat[gn.DOFSPACE].set_dof(oldnode, newnode, 'dy')
        globdat[gn.DOFSPACE].add_dof(newnode, 'phi')

        # Initialize new dof
        oldphidof = globdat[gn.DOFSPACE].get_dof(oldnode, 'phi')
        newphidof = globdat[gn.DOFSPACE].get_dof(newnode, 'phi')
        state0 = globdat[gn.STATE0]
        oldstate0 = globdat[gn.OLDSTATE0]
        globdat[gn.STATE0] = np.append(state0,state0[oldphidof])
        globdat[gn.OLDSTATE0] = np.append(oldstate0,oldstate0[oldphidof])

        # Add history
        if gn.HISTORY in globdat:
            globdat[gn.HISTORY] = np.column_stack(( globdat[gn.HISTORY], globdat[gn.HISTORY][:,oldphidof]))

        # Update element connectivity
        globdat[gn.ESET][hingeelem].change_node(oldnode, newnode)

        # Modify hinge variables
        self._nhinges += 1
        self._hingemoments.append(-self._mp * sign)
        self._hingemoments.append(self._mp * sign)
        self._hingedofs.append(newphidof)
        self._hingedofs.append(oldphidof)
        globdat[pn.HINGENODES].append(oldnode)
        globdat[pn.HINGESTEPS].append(globdat[gn.TIMESTEP])

    def _get_stress(self, globdat, disp):
        stressmat = np.zeros((len(self._elems), self._strcount, 2))
        D = self._get_D_matrix()
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

            # TODO: make robust implementation for of N and V for more than 1 ip
            assert(self._ipcount == 1)
            ip = 0

            if self._subtype == LINEAR:
                N = sfuncs[:, ip]
                dN = grads[:, 0, ip]

                B = self._get_B_matrix(N=N, dN=dN, omega=phi)
                evec = np.matmul(B, ue)
                svec = np.matmul(D, evec)

                m1 += weights[ip] * np.dot(B[:,2], svec)
                m2 += weights[ip] * np.dot(B[:,5], svec)

            elif self._subtype == NONLIN:
                N = sfuncs[:, ip]
                dN = grads[:, 0, ip]

                theta = np.matmul(N, [ue[2], ue[5]])
                kappa = np.matmul(dN, [ue[2], ue[5]])
                gamma = (np.cos(theta) * lsps - np.sin(theta) * lcps) / l_0
                eps = (np.cos(theta) * lcps + np.sin(theta) * lsps) / l_0 - 1
                evec = [eps, gamma, kappa]
                svec = np.matmul(D, evec)
                g = gamma / 2
                e = (1 + eps) / 2
                t = 1 / l_0

                m1 += weights[ip] * np.matmul([g, -e, -t], svec)
                m2 += weights[ip] * np.matmul([g, -e, t], svec)

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


    def _get_B_matrix(self, N, dN, omega, gamma=0, eps=0):
        B = np.zeros((self._strcount, self._dofcount))
        for inode in range(self._shape.node_count()):
            i = 3 * inode
            c = np.cos(omega) * dN[inode]
            s = np.sin(omega) * dN[inode]
            B[:, i:(i + 3)] = np.array([[c, s, N[inode] * gamma],
                                        [-s, c, -N[inode] * (1 + eps)],
                                        [0, 0, dN[inode]]])
        return B

    def _get_WN_matrix(self, N, dN, omega, eps=0):
        WN = np.zeros((6, 6))
        dn = dN[1]
        l0 = 1 / dn
        c = np.cos(omega)
        s = np.sin(omega)
        n1 = N[0]
        n2 = N[1]

        WN[0, 2] = n1 * s
        WN[0, 5] = n2 * s
        WN[1, 2] = -n1 * c
        WN[1, 5] = -n2 * c
        WN[2, 2] = -n1 ** 2 * l0 * (1 + eps)
        WN[2, 3] = -n1 * s
        WN[2, 4] = n1 * c
        WN[2, 5] = -n1 * n2 * l0 * (1 + eps)
        WN[3, 5] = -n2 * s
        WN[4, 5] = n2 * c
        WN[5, 5] = -n2 ** 2 * l0 * (1 + eps)

        WN = WN + WN.transpose() - np.diag(np.diag(WN))  # Add symmetric values in lower triangle
        return WN

    def _get_WV_matrix(self, N, dN, omega, gamma=0):
        WV = np.zeros((6, 6))
        dn = dN[1]
        l0 = 1 / dn
        c = np.cos(omega)
        s = np.sin(omega)
        n1 = N[0]
        n2 = N[1]

        WV[0, 2] = n1 * c
        WV[0, 5] = n2 * c
        WV[1, 2] = n1 * s
        WV[1, 5] = n2 * s
        WV[2, 2] = -n1 ** 2 * l0 * gamma
        WV[2, 3] = -n1 * c
        WV[2, 4] = -n1 * s
        WV[2, 5] = -n1 * n2 * l0 * gamma
        WV[3, 5] = -n2 * c
        WV[4, 5] = -n2 * s
        WV[5, 5] = -n2 ** 2 * l0 * gamma

        WV = WV + WV.transpose() - np.diag(np.diag(WV))  # Add symmetric values in lower triangle
        return WV

    def _get_D_matrix(self):
        D = np.diag([self._EA, self._GAs, self._EI])
        return D
    
    def get_hinges(self):
        return self._hingenodes, self._hingedofs, self._hingemoments


def declare(factory):
    factory.declare_model('Frame', FrameModel)
