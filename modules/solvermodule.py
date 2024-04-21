import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

from modules.module import Module
from modules.controlmodule import ControlModule
from utils.constrainer import Constrainer

STOREMATRIX = 'storeMatrix'
GETMASSMATRIX = 'storeMassMatrix'
STORECONSTRAINTS = 'storeConstraints'


class SolverModule(ControlModule):

    def init(self, props, globdat):
        super().init(props, globdat)
        myprops = props[self._name]
        self._store_matrix = bool(eval(myprops.get(STOREMATRIX, 'False')))
        self._store_mass_matrix = bool(eval(myprops.get(GETMASSMATRIX,'False')))
        self._store_constraints = bool(eval(myprops.get(STORECONSTRAINTS, 'False')))

    def run(self, globdat):
        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[gn.MODEL]

        K = np.zeros((dc, dc))
        f = np.zeros(dc)
        f_int = np.zeros(dc)
        c = Constrainer()

        params = {pn.MATRIX0: K, pn.EXTFORCE: f, pn.INTFORCE: f_int, pn.CONSTRAINTS: c}

        # Advance time step
        super().advance(globdat)
        model.take_action(act.ADVANCE, params, globdat)

        # Assemble K
        model.take_action(act.GETMATRIX0, params, globdat)

        # Assemble f
        model.take_action(act.GETEXTFORCE, params, globdat)

        # Get constraints
        model.take_action(act.GETCONSTRAINTS, params, globdat)

        # Constrain K and f
        Kc, fc = c.constrain(K, f)

        # Sparsify and solve
        smat = sparse.csr_matrix(Kc)
        u = linalg.spsolve(smat, fc)

        # Store rhs and solution in Globdat
        globdat[gn.EXTFORCE] = f
        globdat[gn.STATE0] = u

        # Optionally store stiffness matrix in Globdat
        if self._store_matrix:
            globdat[gn.MATRIX0] = K

        # Optionally store the mass matrix
        if self._store_mass_matrix:
            M = np.zeros((dc, dc))
            params[pn.MATRIX2] = M
            globdat[gn.MATRIX2] = M
            model.take_action(act.GETMATRIX2, params, globdat)

        # Optionally store the constrainer in Globdat
        if self._store_constraints:
            globdat[gn.CONSTRAINTS] = c

        return super().run(globdat)

    def shutdown(self, globdat):
        pass

    def __solve(self, globdat):
        pass


def declare(factory):
    factory.declare_module('Solver', SolverModule)
