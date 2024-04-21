import warnings
from materials.elasticmaterial import ElasticMaterial
from materials.elasticmaterial import E_PROP, NU_PROP, RHO_PROP, ANMODEL_PROP, PLANE_STRESS, BAR
from utils import proputils as pu

import numpy as np
from math import exp
from math import isnan
import copy

from scipy.optimize import root_scalar

import sympy as sp

YIELD = 'yield'         # yield stress as function of kappa (hardening)
TOLERANCE = 'tolerance' # return mapping tolerance
MAXITER = 'maxIter'     # maximum number of iterations for return mapping scheme

class J2Material(ElasticMaterial):

    def __init__ (self,rank):
        super().__init__(rank)

        self._y = []
        self._tol = 1.0e-6
        self._maxiter = 100

    def configure(self, props, globdat):
        # read basic elastic properties with inherited class

        super().configure(props,globdat)

        # read new plasticity properties

        self._yfunc = str(props.get(YIELD))
        self._tol = float(props.get(TOLERANCE,self._tol))
        self._maxiter = int(props.get(MAXITER,self._maxiter))

        self._globdat = globdat

        globdat['maxkappa'] = []

    def update(self,strain,ipoint):
        # get elastic strains from converged plastic strain values

        eps_elastic = self._y[ipoint].elastic_strain(strain)

        # compute trial stress using parent class (elastic)

        trial_stiff, trial_stress = super().update(eps_elastic,ipoint)
        
        # check if yield criterium is being violated

        trial_f = self._y[ipoint].yield_trial(trial_stress) 

        # if the criterium is satisfied, the update is elastic
        # if not, perform a return mapping correction

        if trial_f < -self._tol:
            return trial_stiff, trial_stress
        else:

            sol = root_scalar(self._y[ipoint].yield_function,
                              args = (True),
                              fprime=True,
                              rtol=self._tol,
                              maxiter=self._maxiter,
                              method='newton',
                              x0=0.0)

            # if a Newton solver fails, try bisection with increasingly large bounds

            if not sol.converged or sol.root < -self._tol:
                warnings.warn('Newton return mapping failed.')

                upper = self._tol
                while self._y[ipoint].yield_function(upper)[0] > 0.0:
                    upper *= 2

                print('Trying bisection with upper bound',upper)
                sol = root_scalar(self._y[ipoint].yield_function,
                                  args = (False),
                                  rtol=self._tol,
                                  maxiter=self._maxiter,
                                  method='bisect',
                                  bracket=[0,upper])

                if not sol.converged:
                    raise RuntimeError('Return mapping failed')

            # get updated stress and tangent stiffness from the auxiliary class

            stress = self._y[ipoint].updated_stress()
            stiff  = self._y[ipoint].tangent_stiffness()

            return stiff, stress

    def commit(self, ipoint=None):
        # store converged plastic strain and kappa values for the next step

        for y in self._y:
            y.commit()

        maxkappa = 0

        for y in self._y:
            if y._old_kappa > maxkappa:
                maxkappa = y._old_kappa

        self._globdat['maxkappa'].append(maxkappa)

    def create_material_points(self, npoints):
        for i in range(npoints):
            if self._anmodel == BAR:
                self._y.append(self._Yield1D(self._E,0.0,self._yfunc,1))
            elif self._anmodel == PLANE_STRESS:
                self._y.append(self._YieldPlaneStress(self._E,self._nu,self._yfunc,3))
            else:
                self._y.append(self._Yield3D(self._E,self._nu,self._yfunc,self._strcount))

    '''
    Abstract parent class to handle updating plasticity variables
    The idea is to use an optimizer to find the root of the yield
    criterium. This is done with the 'yield_function' member, which
    should return f(dgam) and dfddgam(dgam). After convergence,
    member functions 'updated_stress' and 'tangent_stiffness' can 
    be called.
    '''

    class _Yield:
        def __init__(self,young,poisson,yfunc,strcount):
            self._E = young
            self._nu = poisson
            self._strcount = strcount

            self._G = self._E / 2.0 / (1.0 + self._nu)
            self._K = self._E / 3. / (1. - 2. * self._nu)

            self._old_kappa = 0.0
            self._new_kappa = 0.0

            kappa = sp.symbols('kappa')
            sigy = sp.sympify(yfunc)
            deriv = sp.diff(sigy,kappa)

            self._sigy = sp.lambdify(kappa,sigy,'numpy')
            self._dsigy = sp.lambdify(kappa,deriv,'numpy')

        def elastic_strain(self,strain):
            pass

        def yield_trial(self,sigtr):
            pass

        def yield_function(self,dgam,fprime=True): 
            pass

        def updated_stress(self):
            pass

        def tangent_stiffness(self):
            pass

        def commit(self):
            self._old_kappa = self._new_kappa
            self._old_epsp = np.copy(self._new_epsp)

    '''
    Class that handles 1D plane stress plasticity
    '''

    class _Yield1D (_Yield):
        def __init__(self,young,poisson,yfunc,strcount):
            super().__init__(young,poisson,yfunc,strcount)

            self._old_epsp = np.zeros(1)
            self._new_epsp = np.zeros(1)

            self._stress = np.zeros(1)

        def elastic_strain(self,strain):
            return strain - self._old_epsp

        def yield_trial(self,sigtr):
            self._sigtr = sigtr

            self._new_kappa = self._old_kappa
            self._new_epsp  = self._old_epsp

            return np.abs(self._sigtr) - self._sigy(self._old_kappa)

        def yield_function(self,dgam,fprime=True): 
            n = np.sign(self._sigtr[0])
            
            self._new_kappa = self._old_kappa + dgam
            self._new_epsp  = self._old_epsp + n * dgam

            self._stress = self._sigtr - self._E * n * dgam

            f = np.abs(self._stress) - self._sigy(self._new_kappa)
            der = -self._E * n * np.sign(self._stress) - self._dsigy(self._new_kappa)

            if fprime:
                return f, der
            else:
                return f

        def updated_stress(self):
            return self._stress

        def tangent_stiffness(self):
            h = self._dsigy(self._new_kappa)

            der_depsp = self._E * self._E / (self._E + h)

            return (self._E - der_depsp) * np.ones((1,1))

    '''
    Class that handles 2D plane stress plasticity.
    '''

    class _YieldPlaneStress (_Yield):
        def __init__(self, young, poisson,yfunc,strcount):
            super().__init__(young,poisson,yfunc,strcount)

            self._P = np.zeros((3, 3))

            self._P[0, 0] = self._P[1, 1] = 2. / 3.
            self._P[0, 1] = self._P[1, 0] = -1. / 3.
            self._P[2, 2] = 2.

            self._stress = np.zeros(3)

            self._old_epsp = np.zeros(strcount)
            self._new_epsp = np.zeros(strcount)

        def elastic_strain(self,strain):
            return strain - self._old_epsp

        def yield_trial(self,sigtr):
            self._sigtr = sigtr

            self._a1tr = (sigtr[0] + sigtr[1])**2.
            self._a2tr = (sigtr[1] - sigtr[0])**2.
            self._a3tr = sigtr[2]**2.

            self._ksitr = self._a1tr / 6.0 + 0.5 * self._a2tr + 2.0 * self._a3tr

            self._new_kappa = self._old_kappa
            self._new_epsp  = self._old_epsp

            return 0.5 * self._ksitr - self._sigy(self._old_kappa)**2. / 3.0

        def yield_function(self,dgam,fprime=True): 
            f1 = 1. + self._E * dgam / 3. / (1. - self._nu)
            f2 = 1. + 2. * self._G * dgam

            f12 = 6.*f1**2. 
            f13 = 9.*f1**3.
            f22 = f2**2.
            f23 = f2**3.

            ksi = self._a1tr / f12 + (0.5 * self._a2tr + 2. * self._a3tr) / f22
            ksider = -self._a1tr * self._E/(1.-self._nu) / f13 - 2. * self._G * (self._a2tr + 4.*self._a3tr) / f23

            self._new_kappa = self._old_kappa + dgam * (2./3.*ksi)**0.5

            sy = self._sigy(self._new_kappa)
            h = self._dsigy(self._new_kappa)

            hbar = 2.*(2./3.)**0.5 * sy*h*(ksi**0.5 + dgam*ksider/2./(ksi**0.5))

            f   = 0.5 * ksi - sy**2.0 / 3.0
            der = 0.5 * ksider - hbar / 3.0

            a1 = 3. * (1.-self._nu) / (3. * (1. - self._nu) + self._E * dgam)
            a2 = 1. / f2

            a_mat = np.zeros((3,3))

            a_mat[[0,1],[0,1]] = 0.5*(a1+a2)
            a_mat[[0,1],[1,0]] = 0.5*(a1-a2)
            a_mat[2,2] = a2

            self._stress = np.matmul(a_mat,self._sigtr)

            self._dgam = dgam
            self._new_epsp = self._old_epsp + self._dgam * np.matmul(self._P,self._stress)

            if fprime:
                return f, der
            else:
                return f

        def updated_stress(self):
            return self._stress

        def tangent_stiffness(self):
            h = self._dsigy(self._new_kappa)

            psig = np.matmul(self._P,self._stress)

            k = np.dot(psig,self._stress)

            e_mat = np.zeros((3,3))

            e1 = 3. * self._E / (3. * (1. - self._nu) + self._E * self._dgam)
            e2 = 2. * self._G / (1. + 2. * self._G * self._dgam)

            e_mat[[0,1],[0,1]] = .5 * (e1 + e2)
            e_mat[[0,1],[1,0]] = .5 * (e1 - e2)
            e_mat[2, 2] = .5 * e2

            n = np.matmul(e_mat,psig)

            alpha = 1. / (np.dot(psig,n) + 2.*k*h/(3.-2.*h*self._dgam))

            return e_mat - alpha*np.outer(n,n)
    '''
    Class that handles plane strain or 3D plasticity.
    The idea is to bring everything to 3D, perform the plasticity
    computations and go back to the original space
    '''

    class _Yield3D (_Yield):
        def __init__(self,young,poisson,yfunc,strcount):
            super().__init__(young,poisson,yfunc,strcount)

            self._P = np.zeros((6,6))

            self._P[[0,1,2],[0,1,2]] = 2./3.
            self._P[[0,0,1,1,2,2],[1,2,0,2,0,1]] = -1./3.
            self._P[[3,4,5],[3,4,5]] = 1./2.

            self._efac = self._E * (1.-self._nu) / (1.+self._nu) / (1.-2.*self._nu)
            self._nfac = self._E * self._nu / (1.+self._nu) / (1.-2.*self._nu)

            self._stress = np.zeros(6)

            self._old_epsp = np.zeros(6)
            self._new_epsp = np.zeros(6)

        def elastic_strain(self,strain):
            self._epse = strain - self._shrink_vector(self._old_epsp)
            return self._epse

        def yield_trial(self,sigtr):
            sigtr6 = self._expand_stress(sigtr)

            i1tr = np.sum(sigtr6[:3])
            i2tr = (i1tr**2. - np.sum(sigtr6[:3]**2.)) / 2.

            self._ptr = i1tr/3.
            self._devtr = sigtr6 - self._ptr*np.array([1,1,1,0,0,0])

            j2tr = i1tr**2.0 / 3. - i2tr

            self._ksitr = (3.*j2tr)**0.5

            self._new_kappa = self._old_kappa
            self._new_epsp  = self._old_epsp

            return self._ksitr - self._sigy(self._old_kappa)

        def yield_function(self,dgam,fprime=True): 
            self._new_kappa = self._old_kappa + dgam

            ksi = self._ksitr - 3. * self._G * dgam

            f = ksi - self._sigy(self._new_kappa)
            der = -3. * self._G - self._dsigy(self._new_kappa)

            frac = 1. - dgam * 3. * self._G / self._ksitr

            dev = frac * self._devtr

            self._stress = dev + self._ptr * np.array([1,1,1,0,0,0])

            dev_norm = (np.dot(dev[:3],dev[:3]) + 2.*np.dot(dev[3:],dev[3:]))**0.5

            depsp = dgam * (3./2.)**0.5 * dev / dev_norm
            depsp[3:] *= 2.

            self._dgam = dgam
            self._new_epsp = self._old_epsp + depsp

            if fprime:
                return f, der
            else:
                return f

        def updated_stress(self):
            return self._shrink_vector(self._stress)

        def tangent_stiffness(self):
            devtr_norm = (np.dot(self._devtr[:3],self._devtr[:3]) + 2.*np.dot(self._devtr[3:],self._devtr[3:]))**0.5

            n = self._devtr / devtr_norm

            fac1 = 2.*self._G * (1.-self._dgam*3.*self._G/self._ksitr)
            fac2 = 6.*self._G**2. * (self._dgam/self._ksitr - 1./(3.*self._G + self._dsigy(self._new_kappa)))

            n_mat = np.outer(n,n)
            i_mat = np.block([[np.ones((3,3)), np.zeros((3,3))], [np.zeros((3,3)), np.zeros((3,3))]])

            tangent = fac1*self._P + fac2*n_mat + self._K*i_mat

            return self._shrink_matrix(tangent)
        
        def _shrink_vector(self,vec):
            ret = np.zeros(self._strcount)

            if self._strcount == 3:
                ret[:2] = vec[:2]
                ret[2] = vec[3]
            else:
                ret = vec

            return ret

        def _expand_stress(self,vec):
            ret = np.zeros(6)

            if self._strcount == 3: 
                ret[:2] = vec[:2]
                ret[3]  = vec[2]

                # Make PLANE_STRAIN adjustments

                ret[0] -= self._old_epsp[2] * self._nfac
                ret[1] -= self._old_epsp[2] * self._nfac
                ret[2]  = np.sum(self._epse[:2]) * self._nfac - self._old_epsp[2] * self._efac
            else:
                ret = vec

            return ret

        def _shrink_matrix(self,mat):
            ret = np.zeros((self._strcount,self._strcount))

            if self._strcount == 3:
                ret[:2,:2] = mat[:2,:2]
                ret[:2,2]  = mat[:2,3]
                ret[2,:2]  = mat[3,:2]
                ret[2,2]   = mat[3,3]
            else:
                ret = mat

            return ret


