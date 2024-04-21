import numpy as np

from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

from modules.module import Module

from utils import proputils as pu

GROUPS = 'groups'
DISP = 'disp'
LOAD = 'load'

class LoadDispModule(Module):

    def init(self, props, globdat):

        if self._name in props:
            myprops = props.get(self._name)

            if GROUPS in myprops:
                self._groups = pu.parse_list(myprops[GROUPS])

        mydata = {}

        for group in self._groups:
            mydata[group] = {}
            mydata[group][LOAD] = {}
            mydata[group][DISP] = {}
            for typ in globdat[gn.DOFSPACE].get_types():
                mydata[group][LOAD][typ] = []
                mydata[group][DISP][typ] = []

        globdat[self._name] = mydata

    def run(self, globdat):

        model = globdat[gn.MODEL]
        nodes = globdat[gn.NSET]
        disp = globdat[gn.STATE0]
        dofs = globdat[gn.DOFSPACE]
        types = dofs.get_types()
        dc = dofs.dof_count()
        fint = np.zeros(dc)
        mydata = globdat[self._name]

        params = {}
        params[pn.INTFORCE] = fint
        params[pn.MATRIX0] = np.zeros((dc, dc))

        model.take_action(act.GETMATRIX0, params, globdat)

        if globdat[gn.ACCEPTED]:
            for group in self._groups:
                for typ in types:
                    idofs = dofs.get_dofs(globdat[gn.NGROUPS][group], [typ])
                    mydata[group][DISP][typ].append( np.mean(disp[idofs]) )
                    mydata[group][LOAD][typ].append( np.sum(fint[idofs]) )

        return 'ok'

    def shutdown(self, globdat):
        for group in self._groups:
            groupData = globdat[self._name][group]
            for typ in globdat[gn.DOFSPACE].get_types():
                groupData[DISP][typ] = np.array( groupData[DISP][typ] )
                groupData[LOAD][typ] = np.array( groupData[LOAD][typ] )

def declare(factory):
    factory.declare_module('LoadDisp', LoadDispModule)
