import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg
import scipy

from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

from modules.module import Module

from utils.constrainer import Constrainer


class LinBuckModule(Module):

    def init(self, props, globdat):
        myprops = props[self._name]
        globdat[gn.ACCEPTED] = True

    def run(self, globdat):
        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[gn.MODEL]

        print('LinBuckModule: running unit load analysis...')

        K = np.zeros((dc, dc))
        f = np.zeros(dc)
        f_int = np.zeros(dc)
        globdat[gn.STATE0] = np.zeros(dc)
        c = Constrainer()

        params = {pn.MATRIX0: K, pn.CONSTRAINTS: c, pn.EXTFORCE: f, pn.INTFORCE: f_int}

        model.take_action(act.GETMATRIX0, params, globdat)

        model.take_action(act.GETEXTFORCE, params, globdat)

        model.take_action(act.GETCONSTRAINTS, params, globdat)

        Kc, fc = c.constrain(K, f)

        smat = sparse.csr_matrix(Kc)
        u = linalg.spsolve(smat, fc)

        globdat[gn.STATE0] = u

        print('LinBuckModule: running eigenvalue problem...')

        KM = np.zeros((dc, dc))
        KG = np.zeros((dc, dc))

        params[pn.MATRIX0] = KM
        params[pn.MATRIX1] = KG

        model.take_action(act.GETMATRIXLB, params, globdat)

        cdofs, cvals = c.get_constraints()
        fdofs = [i for i in range(dc) if i not in cdofs]
        assert max(max(cvals), -min(cvals)) < 1.e-10, 'LinBuckModule does not work with nonzero Dirichlet BCs'

        ls, vs = scipy.linalg.eig(KM[np.ix_(fdofs, fdofs)], -KG[np.ix_(fdofs, fdofs)])

        for idx in np.argsort(cdofs):
            vs = np.insert(vs, cdofs[idx], cvals[idx], axis=0)

        z = list(zip(ls,vs.transpose()))
        zs = sorted(z, key=lambda f: abs(f[0]))
        lss, vss = list(zip(*zs))
 
        globdat[gn.LBFACTORS] = np.real_if_close(lss)
        globdat[gn.HISTORY] = np.asarray(vss)
        print('LinBuckModule: critical load factor:  %8.3e' % np.real_if_close(lss[0]))

        return 'exit'

    def shutdown(self, globdat):
        pass

    def __solve(self, globdat):
        pass


def declare(factory):
    factory.declare_module('LinBuck', LinBuckModule)
