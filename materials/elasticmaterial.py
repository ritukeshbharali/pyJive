from materials.material import Material
import numpy as np

E_PROP = 'E'
NU_PROP = 'nu'
RHO_PROP = 'rho'
ANMODEL_PROP = 'anmodel'
BAR = 'bar'
PLANE_STRESS = 'plane_stress'
PLANE_STRAIN = 'plane_strain'
SOLID = 'solid'


class ElasticMaterial(Material):

    def __init__(self, rank):
        super().__init__(rank)

        self._rank = rank
        self._strcount = self._rank * (self._rank + 1) // 2

        self._anmodel = ''

        self._E = 1.0
        self._nu = 0.0
        self._rho = 0.0

        self._stiff_matrix = np.zeros((self._strcount, self._strcount))
        self._mass_matrix = np.zeros((self._rank, self._rank))

    def configure(self, props, globdat):

        self._anmodel = props.get(ANMODEL_PROP, self._anmodel)
        assert self._is_valid_anmodel(self._anmodel), 'Analysis model ' + self._anmodel + ' not valid for rank ' + str(
            self._rank)

        self._E = float(props.get(E_PROP, self._E))
        self._nu = float(props.get(NU_PROP, self._nu))
        self._rho = float(props.get(RHO_PROP, self._rho))

        self._stiff_matrix = self._compute_stiff_matrix()
        self._mass_matrix = self._compute_mass_matrix()

    def get_config(self):
        config = {
            ANMODEL_PROP: self._anmodel,
            E_PROP: self._E,
            NU_PROP: self._nu,
            RHO_PROP: self._rho
        }

        return config

    def stress_at_point(self, strain, ipoint=None):
        stiff = self.stiff_at_point(ipoint)
        stress = stiff @ strain

        return stress

    def stiff_at_point(self, ipoint=None):
        return self._stiff_matrix

    def self_weight_at_point(self, gravity, ipoint=None):
        rho = self._get_rho(ipoint)
        weight = rho * gravity

        return weight

    def mass_at_point(self, ipoint=None):
        return self._mass_matrix

    def update(self, strain, ipoint=None):
        stiff = self.stiff_at_point(ipoint)
        stress = self.stress_at_point(strain, ipoint)
        return stiff, stress

    def _compute_stiff_matrix(self, ipoint=None):
        E = self._get_E(ipoint)
        nu = self._get_nu(ipoint)
        stiff = np.zeros((self._strcount, self._strcount))

        if self._anmodel == BAR:
            stiff[0, 0] = E 

        elif self._anmodel == PLANE_STRESS:
            stiff[[0, 1], [0, 1]] = E / (1 - nu ** 2)
            stiff[[0, 1], [1, 0]] = (nu * E) / (1 - nu ** 2)
            stiff[2, 2] = 0.5 * E / (1 + nu)

        elif self._anmodel == PLANE_STRAIN:
            d = (1 + nu) * (1 - 2 * nu)
            stiff[[0, 1], [0, 1]] = E * (1 - nu) / d
            stiff[[0, 1], [1, 0]] = E * nu / d
            stiff[2, 2] = 0.5 * E / (1 + nu)

        elif self._anmodel == SOLID:
            d = (1 + nu) * (1 - 2 * nu)
            stiff[[0, 1, 2], [0, 1, 2]] = E * (1 - nu) / d
            stiff[[0, 0, 1, 1, 2, 2], [1, 2, 0, 2, 0, 1]] = E * nu / d
            stiff[[3, 4, 5], [3, 4, 5]] = 0.5 * E / (1 + nu)

        return stiff

    def _compute_mass_matrix(self, ipoint=None):
        rho = self._get_rho(ipoint)
        mass = np.zeros((self._rank, self._rank))

        if self._rank == 1:
            mass[0, 0] = rho
        elif self._rank == 2:
            mass[[0, 1], [0, 1]] = rho
        elif self._rank == 3:
            mass[[0, 1, 2], [0, 1, 2]] = rho

        return mass

    def _is_valid_anmodel(self, anmodel):
        valid = False

        if self._rank == 1 and self._anmodel == BAR:
            valid = True
        elif self._rank == 2 and (self._anmodel == PLANE_STRESS or self._anmodel == PLANE_STRAIN):
            valid = True
        elif self._rank == 3 and self._anmodel == SOLID:
            valid = True

        return valid

    def _get_E(self, ipoint=None):
        return self._E

    def _get_nu(self, ipoint=None):
        return self._nu

    def _get_rho(self, ipoint=None):
        return self._rho
