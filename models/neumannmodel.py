import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn

from models.model import Model

from utils import proputils as pu

GROUPS = 'groups'
DOFS = 'dofs'
VALS = 'values'
INCR = 'loadIncr'
SIGNAL = 'timeSignal'
TIMEINCR = 'deltaTime'


class NeumannModel(Model):
    def take_action(self, action, params, globdat):
        if action == act.GETEXTFORCE:
            self._get_ext_force(params, globdat)
        if action == act.GETUNITFORCE:
            self._get_unit_force(params, globdat)
        if action == act.ADVANCE:
            self._advance_step(params, globdat)

    def configure(self, props, globdat):
        self._groups = pu.parse_list(props[GROUPS])
        self._dofs = pu.parse_list(props[DOFS])
        self._vals = pu.parse_list(props[VALS], float)
        self._init_load = self._vals

        self._sigtable = False

        self._dt = float(props.get(TIMEINCR,1.0))
        self._signal = str(props.get(SIGNAL,''))

        try:
            self._table = np.loadtxt(self._signal)
            self._sigtable = True
        except:
            pass

        if INCR in props:
            self._loadIncr = pu.parse_list(props[INCR], float)
        else: 
            self._loadIncr = np.zeros(len(self._vals))

    def _get_ext_force(self, params, globdat):
        ds = globdat[gn.DOFSPACE]
        for group, dof, val in zip(self._groups, self._dofs, self._vals):
            for node in globdat[gn.NGROUPS][group]:
                idof = ds.get_dof(node, dof)
                params[pn.EXTFORCE][idof] += val

    def _get_unit_force(self, params, globdat):
        ds = globdat[gn.DOFSPACE]
        for group, dof, incr in zip(self._groups, self._dofs, self._loadIncr):
            for node in globdat[gn.NGROUPS][group]:
                idof = ds.get_dof(node, dof)
                params[pn.UNITFORCE][idof] += incr

    def _advance_step(self, params, globdat):
        globdat[gn.TIME] = (globdat[gn.TIMESTEP]+1) * self._dt

        if not self._signal:
            self._vals = np.array(self._init_load) + globdat[gn.TIMESTEP] * np.array(self._loadIncr)
        elif self._sigtable:
            self._vals = np.array(self._init_load) * self._table[globdat[gn.TIMESTEP]]
        else:
            self._vals = np.array(self._init_load) * eval(self._signal,{"t": globdat[gn.TIME],"np":np})


def declare(factory):
    factory.declare_model('Neumann', NeumannModel)
