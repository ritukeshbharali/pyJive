import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn

from models.model import Model

ELEMENTS = 'elements'
EA = 'EA'
SHAPE = 'shape'
TYPE = 'type'
INTSCHEME = 'intScheme'
DOFTYPES = ['dx']


class BarModel(Model):
    def take_action(self, action, params, globdat):
        if action == act.GETMATRIX0:
            self._get_matrix(params, globdat)

    def configure(self, props, globdat):
        # Get basic parameter values
        self._EA = float(props[EA])

        # Get shape and element info
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][TYPE], props[SHAPE][INTSCHEME])
        egroup = globdat[gn.EGROUPS][props[ELEMENTS]]
        self._elems = [globdat[gn.ESET][e] for e in egroup]

        # Get basic dimensionality info
        self._rank = self._shape.global_rank()
        self._ipcount = self._shape.ipoint_count()
        self._nodecount = self._shape.node_count()
        self._dofcount = self._rank * self._nodecount
        self._strcount = self._rank * (self._rank + 1) // 2   # 1-->1, 2-->3, 3-->6

        # Create a new dof for every node and dof type
        nodes = np.unique([node for elem in self._elems for node in elem.get_nodes()])
        for doftype in DOFTYPES[0:self._rank]:
            globdat[gn.DOFSPACE].add_type(doftype)
            for node in nodes:
                globdat[gn.DOFSPACE].add_dof(node, doftype)

    def _get_matrix(self, params, globdat):
        D = np.array([[self._EA]])

        for elem in self._elems:
            # Get the nodal coordinates of each element
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:self._rank, :]

            # Get the shape functions, gradients and weights of each integration point
            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords)

            # Reset the element stiffness matrix
            elmat = np.zeros((self._dofcount, self._dofcount))

            for ip in range(self._ipcount):
                # Get the B matrix for each integration point
                B = np.zeros((1, self._nodecount))
                B[0, :] = grads[:, :, ip].transpose()

                # Compute the element stiffness matrix
                elmat += weights[ip] * np.matmul(np.transpose(B), np.matmul(D, B))

            # Add the element stiffness matrix to the global stiffness matrix
            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elmat

def declare(factory):
    factory.declare_model('Bar', BarModel)
