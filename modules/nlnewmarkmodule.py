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

ITERMAX = 'itermax'
TOLERANCE = 'tolerance'
LENIENT = 'lenient'

class NLNewmarkModule(ControlModule):
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

        self._itermax = int(myprops.get(ITERMAX, 100))
        self._tolerance = float(myprops.get(TOLERANCE, 1e-6))
        self._lenient = bool(eval(myprops.get(LENIENT,'True')))

        dc = globdat[gn.DOFSPACE].dof_count()
        globdat[gn.STATE1] = np.zeros(dc)
        globdat[gn.STATE2] = np.zeros(dc)
        
    def run(self, globdat):
        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[gn.MODEL]
        
        C = np.zeros((dc, dc))
        M = np.zeros((dc, dc))

        fext = np.zeros(dc)
        c = Constrainer(globdat[gn.STATE0])

        params = {pn.MATRIX1: C, pn.MATRIX2: M, pn.EXTFORCE: fext, pn.CONSTRAINTS: c}

        a_n = np.copy(globdat[gn.STATE0])
        a_n_dot = np.copy(globdat[gn.STATE1])
        a_n_ddot = np.copy(globdat[gn.STATE2])

        # Advance time step
        super().advance(globdat)
        model.take_action(act.ADVANCE, params, globdat)

        # Iteration count
        iteration = 0
        
        # Assemble K and fint
        params[pn.MATRIX0] = np.zeros((dc,dc))
        params[pn.INTFORCE] = np.zeros(dc)
        model.take_action(act.GETMATRIX0, params, globdat)
        
        # Assemble C
        # -
        
        # Assemble M
        model.take_action(act.GETMATRIX2, params, globdat)

        # Assemble fe
        model.take_action(act.GETEXTFORCE, params, globdat)

        Khat = self._c0 * M + params[pn.MATRIX0]
        fhat = fext + M @ (self._c1 * a_n_dot + self._c2 * a_n_ddot) - params[pn.INTFORCE]

        # Get constraints
        model.take_action(act.GETCONSTRAINTS, params, globdat)
        cdofs, cvals = c.get_constraints()
        fdofs = [i for i in range(dc) if i not in cdofs]

        # Constrain M_hat and f_hat
        Kc, fc = c.constrain(Khat, fhat)

        # Sparsify and solve
        smat = sparse.csr_matrix(Kc)
        u = linalg.spsolve(smat, fc)

        globdat[gn.STATE0] += u

        # Calculate first residual
        r = fhat

        rel = 1
        ref = max(np.linalg.norm(r), np.linalg.norm(fext),1)

        while rel > self._tolerance and iteration < self._itermax:
            iteration += 1

            params[pn.MATRIX0] = np.zeros((dc,dc))
            params[pn.INTFORCE] = np.zeros(dc)
            model.take_action(act.GETMATRIX0,params,globdat)

            globdat[gn.STATE2] = self._c0 * (globdat[gn.STATE0] - a_n) - self._c1 * a_n_dot - self._c2 * a_n_ddot
            r = fext - M @ globdat[gn.STATE2] - params[pn.INTFORCE]

            Khat = self._c0 * M + params[pn.MATRIX0]

            c.set_zero()
            Kc, rc = c.constrain(Khat, r)
            smat = sparse.csr_matrix(Kc)
            du = linalg.spsolve(smat,rc)

            globdat[gn.STATE0] += du
            
            rel = np.linalg.norm(r[np.ix_(fdofs)]) / ref

            print('Iteration %i, relative residual norm: %.4e' %(iteration,rel))

            if np.isnan(rel):
                raise RuntimeError ('NaN residual in time step %i' % self._step)

        # Alert if not convergence
        if rel > self._tolerance:
            if rel > 1:
                raise RuntimeError('Divergence in time step %i' % self._step)
            elif not self._lenient:
                raise RuntimeError('No convergence in time step %i' % self._step)
            else:
                warnings.warn('No convergence in time step %i' % self._step)

        # Check commit
        model.take_action(act.CHECKCOMMIT, params, globdat)

        if globdat[gn.ACCEPTED]:
            print(f"Converged after {iteration} iterations\n")
            model.take_action(act.COMMIT, params, globdat)

        # Update accelerations
        globdat[gn.STATE1] = a_n_dot + self._c3	 * a_n_ddot + self._c4 * globdat[gn.STATE2]            

        # Optionally store stiffness matrix in Globdat
        if self._store_matrix:
            globdat[gn.MATRIX0] = params[pn.MATRIX0]
        
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
    factory.declare_module('NLNewmark', NLNewmarkModule)
