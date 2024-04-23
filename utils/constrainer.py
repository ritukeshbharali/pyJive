import numpy as np
import copy


class Constrainer:
    def __init__(self, state0=None):
        self._dofs = []
        self._vals = []
        self._state0 = state0

    def add_constraint(self, dof, val):
        self._dofs.append(dof)

        if self._state0 is None:
            self._vals.append(val)
        else:
            self._vals.append(val - self._state0[dof])

    def constrain(self, k, f):

        for dof, val in zip(self._dofs, self._vals):
            for i in range(k.get_shape()[0]):
                if i == dof:
                    f[i] = val
                else:
                    f[i] -= k[i, dof] * val

            k[:, dof] = k[dof, :] = 0.0
            k[dof, dof] = 1.0

    def constraindiag(self, k, f):
        kc = np.copy(k)
        fc = np.copy(f)

        for dof, val in zip(self._dofs, self._vals):
            fc[dof] = val
            kc[dof] = 1.0

        return kc, fc
   
    def set_zero(self):
            self._vals = np.zeros(len(self._dofs))

    def set_zero(self):
        self._vals = np.zeros(len(self._dofs))

    def get_constraints(self):
        return self._dofs, self._vals
