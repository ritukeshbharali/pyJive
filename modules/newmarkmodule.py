import warnings
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
GETMASSMATRIX = 'getMassMatrix'
DELTATIME = 'deltaTime'
GAMMA = 'gamma'
BETA = 'beta'

class NewmarkModule(ControlModule):
    def init(self, props, globdat):
        super().init(props, globdat)
        myprops = props[self._name]
        self._dtime = float(myprops[DELTATIME])
        self._store_matrix = bool(eval(myprops.get(STOREMATRIX, 'False')))
        self._get_mass_matrix = bool(eval(myprops.get(GETMASSMATRIX,'False')))
        self._gamma = float(myprops.get(GAMMA, 0.5))
        self._beta = float(myprops.get(BETA, 0.25))
        self._c0 = 1 / (self._beta * self._dtime**2)
        self._c1 = 1 / (self._beta * self._dtime)
        self._c2 = 1 / (2 * self._beta) - 1
        self._c3 = (1 - self._gamma) * self._dtime
        self._c4 = self._gamma * self._dtime

    def run(self, globdat):
        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[gn.MODEL]
        
        K = np.zeros((dc, dc))
        C = np.zeros((dc, dc))
        M = np.zeros((dc, dc))
        fe = np.zeros(dc)
        c = Constrainer()

        params = {pn.MATRIX0: K, pn.MATRIX1: C, pn.MATRIX2: M, pn.EXTFORCE: fe, 
                  pn.INTFORCE: np.zeros(dc), pn.CONSTRAINTS: c}

        if self._step == 0:
            a_n = np.zeros(dc)
            a_n_dot = np.zeros(dc)
            a_n_ddot = np.zeros(dc)
        else:
            a_n = globdat[gn.STATE0]
            a_n_dot = globdat[gn.STATE1]
            a_n_ddot = globdat[gn.STATE2]

        # Advance time step
        super().advance(globdat)
        model.take_action(act.ADVANCE, params, globdat)
        
        # Assemble K
        model.take_action(act.GETMATRIX0, params, globdat)
        
        # Assemble C
        # -
        
        # Assemble M
        model.take_action(act.GETMATRIX2, params, globdat)

        # Assemble fe
        model.take_action(act.GETEXTFORCE, params, globdat)

        fhat = M @ (self._c0 * a_n + self._c1 * a_n_dot + self._c2 * a_n_ddot) + fe
        Khat = self._c0 * M + K

        # Get constraints
        model.take_action(act.GETCONSTRAINTS, params, globdat)

        # Constrain M_hat and f_hat
        Kc, fc = c.constrain(Khat, fhat)
        
        # Sparsify and solve
        smat = sparse.csr_matrix(Kc)
        a_plus1 = linalg.spsolve(smat, fc)

        globdat[gn.STATE0] = a_plus1
        globdat[gn.STATE2] = self._c0 * (a_plus1 - a_n) - self._c1 * a_n_dot - self._c2 * a_n_ddot
        globdat[gn.STATE1] = a_n_dot + self._c3	 * a_n_ddot + self._c4 * globdat[gn.STATE2]            

        # Optionally store stiffness matrix in Globdat
        if self._store_matrix:
            globdat[gn.MATRIX0] = K
        
        # Optionally get the mass matrix
        if self._get_mass_matrix:
            M = np.zeros((dc, dc))
            params[pn.MATRIX2] = M
            globdat[gn.MATRIX2] = M
            model.take_action(act.GETMATRIX2, params, globdat)
            
        return super().run(globdat)
    
    def shutdown(self, globdat):
        pass
    
    def __solve(self, globdat):
        pass

def declare(factory):
    factory.declare_module('Newmark', NewmarkModule)
