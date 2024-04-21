import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn
from models.model import Model
from materials.material import new_material
from utils.xtable import XTable
from utils import proputils as pu

ELEMENTS = 'elements'
SHAPE = 'shape'
TYPE = 'type'
INTSCHEME = 'intScheme'
DOFTYPES = ['dx', 'dy', 'dz']
MATERIAL = 'material'
THICKNESS_PROP = 'thickness'
AREA_PROP = 'area'
GRAVITY = 'gravity'


class SolidModel(Model):
    def take_action(self, action, params, globdat):

        if action == act.GETMATRIX0:
            self._get_matrix(params, globdat)
        elif action == act.GETEXTFORCE:
            self._get_body_force(params, globdat)
        elif action == act.COMMIT:
            self._commit(params, globdat)
        elif action == act.CHECKCOMMIT:
            self._check_commit(params, globdat)
        elif action == act.GETTABLE:
            if 'stress' in params[pn.TABLENAME]:
                self._get_stresses(params, globdat)
            if 'strain' in params[pn.TABLENAME]:
                self._get_strains(params, globdat)
            if 'stiffness' in params[pn.TABLENAME]:
                self._get_stiffness(params, globdat)

    def configure(self, props, globdat):

        # Configure the material
        matprops = props[MATERIAL]
        self._mat = new_material(matprops)
        self._mat.configure(matprops, globdat)

        # Get shape and element info
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][TYPE], props[SHAPE][INTSCHEME])
        egroup = globdat[gn.EGROUPS][props[ELEMENTS]]
        self._elems = [globdat[gn.ESET][e] for e in egroup]

        # Make sure the shape rank and mesh rank are identitcal
        if self._shape.global_rank() != globdat[gn.MESHRANK]:
            raise RuntimeError('ElasticModel: Shape rank must agree with mesh rank')

        # Get basic dimensionality info
        self._rank = self._shape.global_rank()
        self._ipcount = self._shape.ipoint_count()
        self._dofcount = self._rank * self._shape.node_count()
        self._strcount = self._rank * (self._rank + 1) // 2  # 1-->1, 2-->3, 3-->6

        self._rankfac = 1.0

        if self._rank == 1:
            self._rankfac = props.get(AREA_PROP, self._rankfac)
            self._rankfac = pu.soft_cast(self._rankfac, float)
        elif self._rank == 2:
            self._rankfac = props.get(THICKNESS_PROP, self._rankfac)
            self._rankfac = pu.soft_cast(self._rankfac, float)

        self._gravity = bool(eval(props.get(GRAVITY, 'False')))

        # Initialize integration point
        total_ipcount =  self._ipcount * len(self._elems)
        self._mat.create_material_points( total_ipcount )

        # Create a new dof for every node and dof type
        nodes = np.unique([node for elem in self._elems for node in elem.get_nodes()])
        for doftype in DOFTYPES[0:self._rank]:
            globdat[gn.DOFSPACE].add_type(doftype)
            for node in nodes:
                globdat[gn.DOFSPACE].add_dof(node, doftype)

    def _commit(self, params, globdat):
        self._mat.commit()

    def _check_commit(self, params, globdat):
        self._mat.check_commit(params, globdat)

    def _get_matrix(self, params, globdat):
        # _get_matrix Directly stores the internal force as well to avoid the material being computed twice

        ipoint = 0

        # Get the STATE0 vector
        disp = globdat[gn.STATE0]

        for elem in self._elems:
            # Get the nodal coordinates of each element
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:self._rank, :]

            # Get the gradients, weights and coordinates of each integration point
            grads, weights = self._shape.get_shape_gradients(coords)
            if self._rank < 3:
                weights *= self._rankfac

            # Reset the element stiffness matrix
            elmat = np.zeros((self._dofcount, self._dofcount))

            # Get the displacements at the element nodes.
            eldisp = disp[idofs]

            # Reset the element force vector
            elfor = np.zeros(self._dofcount)

            for ip in range(self._ipcount):
                # Get the B matrix for this integration point
                B = self._get_B_matrix(grads[:, :, ip])

                # Compute the strain vector of this integration point
                strain = np.matmul(B, eldisp)

                # Update the material model and obtain the stiffness given the current strain
                D, stress = self._mat.update(strain, ipoint)

                # Compute the element stiffness matrix and force
                elmat += weights[ip] * np.matmul(np.transpose(B), np.matmul(D, B))
                elfor += weights[ip] * np.matmul(np.transpose(B), stress).reshape(-1)

                ipoint += 1

            # Add the element stiffness matrix to the global stiffness matrix
            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elmat

            # Add the element internal force vector to the global internal force vector
            params[pn.INTFORCE][idofs] += elfor

    def _get_body_force(self, params, globdat):

        if not self._gravity:
            return

        if self._rank == 1:
            gravity = np.array([1])
        elif self._rank == 2:
            gravity = np.array([0, -1])
        else:
            return

        for elem in self._elems:
            # Get the nodal coordinates of each element
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:self._rank, :]

            # Get the shape functions, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            weights = self._shape.get_integration_weights(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            if self._rank < 3:
                weights *= self._rankfac

            # Reset the element force vector
            elfor = np.zeros(self._dofcount)

            for ip in range(self._ipcount):
                # Get the N matrix and b vector for each integration point
                N = self._get_N_matrix(sfuncs[:, ip])
                b = self._mat.self_weight_at_point(gravity, ipcoords[:, ip])

                # Compute the element force vector
                elfor += weights[ip] * np.matmul(np.transpose(N), b)

            # Add the element force vector to the global force vector
            params[pn.EXTFORCE][idofs] += elfor

    def _get_strains(self, params, globdat):
        xtable = params[pn.TABLE]
        tbwts = params[pn.TABLEWEIGHTS]

        # Convert the table to an XTable and store the original class
        cls_ = xtable.__class__
        xtable.__class__ = XTable

        # Get the STATE0 vector if no custom displacement field is provided
        if pn.SOLUTION in params:
            disp = params[pn.SOLUTION]
        else:
            disp = globdat[gn.STATE0]

        # Add the columns of all stress components to the table
        if self._rank == 1:
            jcols = xtable.add_columns(['xx'])
        elif self._rank == 2:
            jcols = xtable.add_columns(['xx', 'yy', 'xy'])
        elif self._rank == 3:
            jcols = xtable.add_columns(['xx', 'yy', 'zz', 'xy', 'yz', 'zx'])

        for elem in self._elems:
            # Get the nodal coordinates of each element
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:self._rank, :]

            # Get the shape functions, gradients, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            # Get the nodal displacements
            eldisp = disp[idofs]

            # Reset the element stress matrix and weights
            eleps = np.zeros((self._shape.node_count(), self._strcount))
            elwts = np.zeros(self._shape.node_count())

            for ip in range(self._ipcount):
                # Get the B matrix for each integration point
                B = self._get_B_matrix(grads[:, :, ip])

                # Get the strain of the element in the integration point
                strain = np.matmul(B, eldisp)

                # Compute the element strain and weights
                eleps += np.outer(sfuncs[:, ip], strain)
                elwts += sfuncs[:, ip].flatten()

            # Add the element weights to the global weights
            tbwts[inodes] += elwts

            # Add the element stresses to the global stresses
            xtable.add_block(inodes, jcols, eleps)

        # Convert the table back to the original class
        xtable.__class__ = cls_

    def _get_stresses(self, params, globdat):
        xtable = params[pn.TABLE]
        tbwts = params[pn.TABLEWEIGHTS]

        # Convert the table to an XTable and store the original class
        cls_ = xtable.__class__
        xtable.__class__ = XTable

        # Get the STATE0 vector if no custom displacement field is provided
        if pn.SOLUTION in params:
            disp = params[pn.SOLUTION]
        else:
            disp = globdat[gn.STATE0]

        # Add the columns of all stress components to the table
        if self._rank == 1:
            jcols = xtable.add_columns(['xx'])
        elif self._rank == 2:
            jcols = xtable.add_columns(['xx', 'yy', 'xy'])
        elif self._rank == 3:
            jcols = xtable.add_columns(['xx', 'yy', 'zz', 'xy', 'yz', 'zx'])

        for elem in self._elems:
            # Get the nodal coordinates of each element
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:self._rank, :]

            # Get the shape functions, gradients, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            # Get the nodal displacements
            eldisp = disp[idofs]

            # Reset the element stress matrix and weights
            elsig = np.zeros((self._shape.node_count(), self._strcount))
            elwts = np.zeros(self._shape.node_count())

            for ip in range(self._ipcount):
                # Get the B matrix for each integration point
                B = self._get_B_matrix(grads[:, :, ip])

                # Get the strain and stress of the element in the integration point
                strain = np.matmul(B, eldisp)
                stress = self._mat.stress_at_point(strain, ipcoords[:, ip])

                # Compute the element stress and weights
                elsig += np.outer(sfuncs[:, ip], stress)
                elwts += sfuncs[:, ip].flatten()

            # Add the element weights to the global weights
            tbwts[inodes] += elwts

            # Add the element stresses to the global stresses
            xtable.add_block(inodes, jcols, elsig)

        # Convert the table back to the original class
        xtable.__class__ = cls_

    def _get_stiffness(self, params, globdat):
        xtable = params[pn.TABLE]
        tbwts = params[pn.TABLEWEIGHTS]

        # Convert the table to an XTable and store the original class
        cls_ = xtable.__class__
        xtable.__class__ = XTable

        # Add the column of the Young's modulus to the table
        jcol = xtable.add_column('')

        for elem in self._elems:
            # Get the nodal coordinates of each element
            inodes = elem.get_nodes()
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:self._rank, :]

            # Get the shape functions, gradients, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            ipcoords = self._shape.get_global_integration_points(coords)

            # Reset the element stress matrix and weights
            elyoung = np.zeros((self._shape.node_count()))
            elwts = np.zeros(self._shape.node_count())

            for ip in range(self._ipcount):
                # Get the stiffness in the integration point
                E = self._mat._get_E(ipcoords[:, ip])

                # Compute the element stiffness and weights
                elyoung += E * sfuncs[:, ip]
                elwts += sfuncs[:, ip].flatten()

            # Add the element weights to the global weights
            tbwts[inodes] += elwts

            # Add the element stiffness to the global stiffness
            xtable.add_col_values(inodes, jcol, elyoung)

        # Divide the stresses by the shape function weights
        values = xtable.get_col_values(None, jcol)
        xtable.set_col_values(None, jcol, values / tbwts)

        # Convert the table back to the original class
        xtable.__class__ = cls_

    def _get_N_matrix(self, sfuncs):
        N = np.zeros((self._rank, self._dofcount))
        for i in range(self._rank):
            N[i, i::self._rank] = sfuncs.transpose()
        return N

    def _get_B_matrix(self, grads):
        B = np.zeros((self._strcount, self._dofcount))
        if self._rank == 1:
            B = grads.transpose()
        elif self._rank == 2:
            for inode in range(self._shape.node_count()):
                i = 2 * inode
                gi = grads[inode, :]
                B[0:3, i:(i + 2)] = [
                    [gi[0], 0.],
                    [0., gi[1]],
                    [gi[1], gi[0]]
                ]
        elif self._rank == 3:
            B = np.zeros((6, self._dofcount))
            for inode in range(self._shape.node_count()):
                i = 3 * inode
                gi = grads[inode, :]
                B[0:6, i:(i + 3)] = [
                    [gi[0], 0., 0.],
                    [0., gi[1], 0.],
                    [0., 0., gi[2]],
                    [gi[1], gi[0], 0.],
                    [0., gi[2], gi[1]],
                    [gi[2], 0., gi[0]]
                ]
        return B


def declare(factory):
    factory.declare_model('Solid', SolidModel)
