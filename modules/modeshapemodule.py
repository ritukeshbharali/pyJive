import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg
import scipy
from warnings import warn

from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

from modules.module import Module

from utils.constrainer import Constrainer

class ModeShapeModule(Module):
    
    def init(self, props, globdat):
        myprops = props[self._name]
    
    def run(self, globdat):
        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[gn.MODEL]
        rigidBody = False
        
        print('ModeShapeModule: running eigenvalue problem...')
        
        K = np.zeros((dc, dc))
        M = np.zeros((dc, dc))
        f_int = np.zeros(dc)
        globdat[gn.STATE0] = np.zeros(dc)
        c = Constrainer()
        
        params = {pn.MATRIX0: K, pn.MATRIX2: M, pn.CONSTRAINTS: c, pn.INTFORCE: f_int}
        
        model.take_action(act.GETMATRIX0, params, globdat)

        model.take_action(act.GETMATRIX2, params, globdat)

        model.take_action(act.GETCONSTRAINTS, params, globdat)
        
        cdofs, cvals = c.get_constraints()
        fdofs = [i for i in range(dc) if i not in cdofs]
        if len(cvals) > 0:
            assert max(max(cvals), -min(cvals)) < 1.e-10, 'ModeShapeModule does not work with nonzero Dirichlet BCs'
            if len(cvals) < 3: 
                print("Less than three Dirichlet BCs, at least one rigid body mode will be found")
                rigidBody = True
        else:
            print("No Dirichlet BCs, three rigid body modes will be found")
            rigidBody = True
        
        nfreqs_sq, modes = sparse.linalg.eigs(K[np.ix_(fdofs, fdofs)], len(fdofs), M[np.ix_(fdofs, fdofs)])
        
        nfreqs = np.sqrt(nfreqs_sq)
        
        for idx in np.argsort(cdofs):
            modes = np.insert(modes, cdofs[idx], cvals[idx], axis=0)
        
        z = list(zip(nfreqs, modes.transpose()))
        zs = sorted(z, key=lambda f: abs(f[0]))
        nfreqs_sorted, modes_sorted = list(zip(*zs))
        
        globdat[gn.EIGENFREQS] = np.real_if_close(nfreqs_sorted)
        globdat[gn.MODALSHAPES] = modes_sorted
        globdat[gn.HISTORY] = np.asarray(modes_sorted)
        
        print(f'ModeShapeModule: smallest natural frequency  {np.real_if_close(nfreqs_sorted[0]):.4e} rad / s') 
        
        return 'exit'
    
    def shutdown(self, globdat):
        pass

    def __solve(self, globdat):
        pass

def declare(factory):
    factory.declare_module('ModeShape', ModeShapeModule)
