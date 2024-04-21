import numpy as np

from utils.shape import Shape

class Tri3Shape(Shape):
    def __init__(self, intscheme):
        print('Creating Tri3Shape...')

        # Set the nodecount and rank of the elements
        self._ncount = 3
        self._rank = 2

        # Refer to the Shape class to handle the rest of the initialization
        super().__init__(intscheme)

    def get_local_node_coords(self):
        # Return the standard triangle with nodes at (0,0), (1,0) and (0,1)
        loc_coords = np.zeros((self._rank, self._ncount))

        loc_coords[0, 0] = 0.0
        loc_coords[0, 1] = 1.0
        loc_coords[0, 2] = 0.0
        loc_coords[1, 0] = 0.0
        loc_coords[1, 1] = 0.0
        loc_coords[1, 2] = 1.0

        return loc_coords

    def eval_shape_functions(self, loc_point):
        # Evalulate the shape functions in the local coordinate system
        sfuncs = np.zeros(self._ncount)

        sfuncs[0] = 1.0 - loc_point[0] - loc_point[1]
        sfuncs[1] = loc_point[0]
        sfuncs[2] = loc_point[1]

        return sfuncs

    def eval_shape_gradients(self, loc_point):
        # Evaluate the shape gradients in the local coordinate system
        # Note that no weights are applied!
        sgrads = np.zeros((self._ncount, self._rank))

        sgrads[0, 0] = -1.0
        sgrads[0, 1] = -1.0
        sgrads[1, 0] = 1.0
        sgrads[1, 1] = 0.0
        sgrads[2, 0] = 0.0
        sgrads[2, 1] = 1.0

        return sgrads

    def get_local_point(self, glob_point, glob_coords):
        # Return the local coordinates corresponding to the given global point
        Ax = glob_coords[0,0]
        Ay = glob_coords[1,0]
        Bx = glob_coords[0,1]
        By = glob_coords[1,1]
        Cx = glob_coords[0,2]
        Cy = glob_coords[1,2]

        mat = np.zeros((self._rank, self._rank))
        rhs = np.zeros(self._rank)

        mat[0,0] = Bx - Ax
        mat[0,1] = Cx - Ax
        mat[1,0] = By - Ay
        mat[1,1] = Cy - Ay

        rhs[0] = glob_point[0] - Ax
        rhs[1] = glob_point[1] - Ay

        return np.linalg.solve(mat, rhs)

    def contains_local_point(self, loc_point, tol=0.0):
        # Return whether or not the local point falls inside or on the element boundaries
        if loc_point[0] + loc_point[1] > 1 + tol:
            return False
        elif loc_point[0] < 0 - tol:
            return False
        elif loc_point[1] < 0 - tol:
            return False
        else:
            return True

class Tri6Shape(Shape):
    def __init__(self, intscheme):
        print('Creating Tri6Shape...')

        # Set the nodecount and rank of the elements
        self._ncount = 6
        self._rank = 2

        # Refer to the Shape class to handle the rest of the initialization
        super().__init__(intscheme)

    def get_local_node_coords(self):
        # Return the standard triangle with nodes at (0,0), (1,0) and (0,1)
        loc_coords = np.zeros((self._rank, self._ncount))

        loc_coords[0, 0] = 0.0
        loc_coords[0, 1] = 1.0
        loc_coords[0, 2] = 0.0
        loc_coords[0, 3] = 0.5
        loc_coords[0, 4] = 0.5
        loc_coords[0, 5] = 0.0
        loc_coords[1, 0] = 0.0
        loc_coords[1, 1] = 0.0
        loc_coords[1, 2] = 1.0
        loc_coords[1, 3] = 0.0
        loc_coords[1, 4] = 0.5
        loc_coords[1, 5] = 0.5

        return loc_coords

    def eval_shape_functions(self, loc_point):
        # Evalulate the shape functions in the local coordinate system
        sfuncs = np.zeros(self._ncount)

        sfuncs[0] = 2 * (0.5 - loc_point[0] - loc_point[1]) * (1 - loc_point[0] - loc_point[1])
        sfuncs[1] = -2 * loc_point[0] * (0.5 - loc_point[0])
        sfuncs[2] = -2 * loc_point[1] * (0.5 - loc_point[1])
        sfuncs[3] = 4 * loc_point[0] * (1 - loc_point[0] - loc_point[1])
        sfuncs[4] = 4 * loc_point[0] * loc_point[1]
        sfuncs[5] = 4 * loc_point[1] * (1 - loc_point[0] - loc_point[1])

        return sfuncs

    def eval_shape_gradients(self, loc_point):
        # Evaluate the shape gradients in the local coordinate system
        # Note that no weights are applied!
        sgrads = np.zeros((self._ncount, self._rank))

        sgrads[0, 0] = -3 + 4 * loc_point[0] + 4 * loc_point[1]
        sgrads[0, 1] = -3 + 4 * loc_point[0] + 4 * loc_point[1]
        sgrads[1, 0] = -1 + 4 * loc_point[0]
        sgrads[1, 1] = 0.0
        sgrads[2, 0] = 0.0
        sgrads[2, 1] = -1 + 4 * loc_point[1]
        sgrads[3, 0] = 4 - 8 * loc_point[0] - 4 * loc_point[1]
        sgrads[3, 1] = -4 * loc_point[0]
        sgrads[4, 0] = 4 * loc_point[1]
        sgrads[4, 1] = 4 * loc_point[0]
        sgrads[5, 0] = -4 * loc_point[1]
        sgrads[5, 1] = 4 - 4 * loc_point[0] - 8 * loc_point[1]

        return sgrads

    def contains_local_point(self, loc_point, tol=0.0):
        # Return whether or not the local point falls inside or on the element boundaries
        if loc_point[0] + loc_point[1] > 1 + tol:
            return False
        elif loc_point[0] < 0 - tol:
            return False
        elif loc_point[1] < 0 - tol:
            return False
        else:
            return True

class Quad4Shape(Shape):
    def __init__(self, intscheme):
        print('Creating Quad4Shape...')

        # Set the nodecount and rank of the elements
        self._ncount = 4
        self._rank = 2

        # Refer to the Shape class to handle the rest of the initialization
        super().__init__(intscheme)

    def get_local_node_coords(self):
        # Return the standard triangle with nodes at (0,0), (1,0) and (0,1)
        loc_coords = np.zeros((self._rank, self._ncount))

        loc_coords[0, 0] = -1.0
        loc_coords[0, 1] = 1.0
        loc_coords[0, 2] = 1.0
        loc_coords[0, 3] = -1.0
        loc_coords[1, 0] = -1.0
        loc_coords[1, 1] = -1.0
        loc_coords[1, 2] = 1.0
        loc_coords[1, 3] = 1.0

        return loc_coords

    def eval_shape_functions(self, loc_point):
        # Evalulate the shape functions in the local coordinate system
        sfuncs = np.zeros(self._ncount)

        x = loc_point[0]
        y = loc_point[1]

        sfuncs[0] = 0.25 * (1 - x) * (1 - y)
        sfuncs[1] = 0.25 * (1 + x) * (1 - y)
        sfuncs[2] = 0.25 * (1 + x) * (1 + y)
        sfuncs[3] = 0.25 * (1 - x) * (1 + y)

        return sfuncs

    def eval_shape_gradients(self, loc_point):
        # Evaluate the shape gradients in the local coordinate system
        # Note that no weights are applied!
        sgrads = np.zeros((self._ncount, self._rank))

        x = loc_point[0]
        y = loc_point[1]

        sgrads[0, 0] = -0.25 * (1 - y)
        sgrads[0, 1] = -0.25 * (1 - x)
        sgrads[1, 0] = 0.25 * (1 - y)
        sgrads[1, 1] = -0.25 * (1 + x)
        sgrads[2, 0] = 0.25 * (1 + y)
        sgrads[2, 1] = 0.25 * (1 + x)
        sgrads[3, 0] = -0.25 * (1 + y)
        sgrads[3, 1] = 0.25 * (1 - x)

        return sgrads

    def contains_local_point(self, loc_point, tol=0.0):
        # Return whether or not the local point falls inside or on the element boundaries
        if loc_point[0] < -1 - tol:
            return False
        elif loc_point[0] > 1 + tol:
            return False
        elif loc_point[1] < -1 - tol:
            return False
        elif loc_point[1] > 1 + tol:
            return False
        else:
            return True

class Line2Shape(Shape):
    def __init__(self, intscheme):
        print('Creating Line2Shape...')

        # Set the nodecount and rank of the elements
        self._ncount = 2
        self._rank = 1

        # Refer to the Shape class to handle the rest of the initialization
        super().__init__(intscheme)

    def get_local_node_coords(self):
        # Return the standard line with nodes at (-1) and (1)
        loc_coords = np.zeros((self._rank, self._ncount))

        loc_coords[0, 0] = -1.0
        loc_coords[0, 1] = 1.0

        return loc_coords

    def eval_shape_functions(self, loc_point):
        # Evalulate the shape functions in the local coordinate system
        sfuncs = np.zeros(self._ncount)

        sfuncs[0] = 0.5 - 0.5 * loc_point[0]
        sfuncs[1] = 0.5 + 0.5 * loc_point[0]

        return sfuncs

    def eval_shape_gradients(self, loc_point):
        # Evaluate the shape gradients in the local coordinate system
        # Note that no weights are applied!
        sgrads = np.zeros((self._ncount, self._rank))

        sgrads[0, 0] = -0.5
        sgrads[1, 0] = 0.5

        return sgrads

    def get_local_point(self, glob_point, glob_coords):
        # Return the local coordinates corresponding to the given global point
        loc_point = np.zeros(self._rank)

        A = glob_coords[0,0]
        B = glob_coords[0,1]
        X = glob_point[0]

        loc_point[0] = (A + B - 2*X)/(A - B)

        return loc_point

    def contains_local_point(self, loc_point, tol=0.0):
        # Return whether or not the local point falls inside or on the element boundaries
        if loc_point[0] > 1 + tol:
            return False
        elif loc_point[0] < -1 - tol:
            return False
        else:
            return True


class Line3Shape(Shape):
    def __init__(self, intscheme):
        print('Creating Line2Shape...')

        # Set the nodecount and rank of the elements
        self._ncount = 3
        self._rank = 1

        # Refer to the Shape class to handle the rest of the initialization
        super().__init__(intscheme)

    def get_local_node_coords(self):
        # Return the standard line with nodes at (-1) and (1)
        loc_coords = np.zeros((self._rank, self._ncount))

        loc_coords[0, 0] = -1.0
        loc_coords[0, 1] = 0.0
        loc_coords[0, 2] = 1.0

        return loc_coords

    def eval_shape_functions(self, loc_point):
        # Evalulate the shape functions in the local coordinate system
        sfuncs = np.zeros(self._ncount)

        sfuncs[0] = 0.5 * loc_point[0] * (loc_point[0] - 1)
        sfuncs[1] = 1 - loc_point[0]**2
        sfuncs[2] = 0.5 * loc_point[0] * (loc_point[0] + 1)

        return sfuncs

    def eval_shape_gradients(self, loc_point):
        # Evaluate the shape gradients in the local coordinate system
        # Note that no weights are applied!
        sgrads = np.zeros((self._ncount, self._rank))

        sgrads[0, 0] = loc_point[0] - 0.5
        sgrads[1, 0] = -2 * loc_point[0]
        sgrads[2, 0] = loc_point[0] + 0.5

        return sgrads

    def contains_local_point(self, loc_point, tol=0.0):
        # Return whether or not the local point falls inside or on the element boundaries
        if loc_point[0] > 1 + tol:
            return False
        elif loc_point[0] < -1 - tol:
            return False
        else:
            return True


def declare(factory):
    factory.declare_shape('Triangle3', Tri3Shape)
    factory.declare_shape('Triangle6', Tri6Shape)
    factory.declare_shape('Quad4', Quad4Shape)
    factory.declare_shape('Line2', Line2Shape)
    factory.declare_shape('Line3', Line3Shape)
