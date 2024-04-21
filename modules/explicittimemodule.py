import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

from modules.module import Module
from modules.controlmodule import ControlModule
from utils.constrainer import Constrainer
from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

STOREMATRIX = 'storeMatrix'
STOREMASSMATRIX = 'storeMassMatrix'
DELTATIME = 'deltaTime'

class ExplicitTimeModule(ControlModule):
    def init(self, props, globdat):
        super().init(props, globdat)
        myprops = props[self._name]
        self._dtime = float(myprops[DELTATIME])
        self._store_matrix = bool(eval(myprops.get(STOREMATRIX, 'False')))
        self._store_mass_matrix = bool(eval(myprops.get(STOREMASSMATRIX,'False')))
        self._c1 = 1 / (self._dtime**2)
        self._c2 = 1 / (2 * self._dtime)
    
    def run(self, globdat):
        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[gn.MODEL]
        
        # Damping currently not used
        K = np.zeros((dc, dc))
        C = np.zeros((dc, dc))
        M = np.zeros((dc, dc))
        f_ext = np.zeros(dc)
        f_int = np.zeros(dc)
        c = Constrainer()

        params = {pn.MATRIX0: K, pn.MATRIX1: C, pn.MATRIX2: M, pn.EXTFORCE: f_ext, 
                  pn.INTFORCE: f_int, pn.CONSTRAINTS: c}

        # Get previous timesteps, set if zero, assemble fi
        if self._step == 0:
            a_n = np.zeros(dc)
            a_min1 = np.zeros(dc)
            fi = np.zeros(dc)
            globdat[gn.STATE0] = np.zeros(dc)
        
        elif self._step == 1:
            a_n = globdat[gn.STATE0]
            a_min1 = np.zeros(dc)
            model.take_action(act.GETINTFORCE, params, globdat)
        
        else:
            a_n = globdat[gn.STATE0]
            a_min1 = globdat[gn.OLDSTATE0]
            model.take_action(act.GETINTFORCE, params, globdat)
        
        # Advance time step
        super().advance(globdat)
        model.take_action(act.ADVANCE, params, globdat)
        
        # Assemble K
        model.take_action(act.GETMATRIX0, params, globdat)
        
        # Assemble C
        # -
        
        # Assemble M
        if self._step == 0:
            model.take_action(act.GETMATRIX2, params, globdat)
            self._Ml = np.sum(M, axis=1)
            self._Cl = np.sum(C, axis=1)
            self._Mhat = self._c1 * self._Ml  + self._c2 * self._Cl

        # Assemble fe
        model.take_action(act.GETEXTFORCE, params, globdat)
        
        # Construct fhat
        fhat = f_ext + self._c1 * self._Ml * (2 * a_n - a_min1) + self._c2 * self._Cl * a_min1 - f_int
            
        # Get constraints
        model.take_action(act.GETCONSTRAINTS, params, globdat)

        # Constrain M_hat and f_hat
        Mc, fc = c.constraindiag(self._Mhat, fhat)
        
        # Solve
        a_plus1 = fc / Mc;
        
        # Store solution in Globdat, move old solution
        globdat[gn.STATE0] = a_plus1
        globdat[gn.OLDSTATE0] = a_n

        # Store velocities and accelerations (of the previous step)
        globdat[gn.STATE1] = (a_plus1 - a_min1) / 2. / self._dtime
        globdat[gn.STATE2] = (a_min1 - 2. * a_n + a_plus1) / self._dtime / self._dtime

        # Optionally store stiffness matrix in Globdat
        if self._store_matrix:
            globdat[gn.MATRIX0] = K
        
        # Optionally get the mass matrix
        if self._store_mass_matrix:
            M = np.zeros((dc, dc))
            globdat[gn.MATRIX2] = M

        return super().run(globdat)
    
    def shutdown(self, globdat):
        pass
    
    def __solve(self, globdat):
        pass

def declare(factory):
    factory.declare_module('ExplicitTime', ExplicitTimeModule)
