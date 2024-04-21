import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn
from names import PropNames as prn

from models.model import Model

ELEMENTS = 'elements'
EI = 'EI'
GAs = 'GAs'
SHAPE = 'shape'
INTSCHEME = 'intScheme'
DOFTYPES = ['phi', 'dy']


class TimoshenkoModel(Model):
    def take_action(self, action, params, globdat):
        print('Model taking action', action)

        if action == act.GETMATRIX0:
            self.__stiffness(params, globdat)

    def configure(self, props, globdat):
        self._EI = float(props[EI])
        self._GAs = float(props[GAs])
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][prn.TYPE], props[SHAPE][INTSCHEME])
        egroup = globdat[gn.EGROUPS][props[ELEMENTS]]
        self._elems = [globdat[gn.ESET][e] for e in egroup]

        self._ipcount = self._shape.ipoint_count()
        self._dofcount = 2 * self._shape.node_count()

        nodes = np.unique([node for elem in self._elems for node in elem.get_nodes()])
        for doftype in DOFTYPES:
            globdat[gn.DOFSPACE].add_type(doftype)
            for node in nodes:
                globdat[gn.DOFSPACE].add_dof(node, doftype)

    def __stiffness(self, params, globdat):
        for elem in self._elems:
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES)
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:1, :]
            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords)

            elmat = np.zeros((4, 4))
            for ip in range(self._ipcount):
                B_theta = np.zeros((1, 4))
                N_theta = np.zeros((1, 4))
                B_v = np.zeros((1, 4))
                N_v = np.zeros((1, 4))
                B_theta[:, 0::2] = grads[:, :, ip].transpose()
                B_v[:, 1::2] = grads[:, :, ip].transpose()
                N_theta[:, 0::2] = sfuncs[:, ip].transpose()
                N_v[:, 1::2] = sfuncs[:, ip].transpose()

                elmat += weights[ip] * (np.matmul(B_theta.transpose() * self._EI, B_theta) + \
                                        np.matmul(N_theta.transpose() * self._GAs, N_theta) - \
                                        np.matmul(N_theta.transpose() * self._GAs, B_v) - \
                                        np.matmul(B_v.transpose() * self._GAs, N_theta) + \
                                        np.matmul(B_v.transpose() * self._GAs, B_v))

            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elmat


def declare(factory):
    factory.declare_model('Timoshenko', TimoshenkoModel)
