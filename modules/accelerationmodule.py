import numpy as np

from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

from modules.module import Module

from utils import proputils as pu

GROUPS = 'groups'
DISP = 'disp'
VELO = 'velo'
ACCEL = 'accel'

class AccelerationModule(Module):

    def init(self, props, globdat):

        if self._name in props:
            myprops = props.get(self._name)

            if GROUPS in myprops:
                self._groups = pu.parse_list(myprops[GROUPS])

        mydata = {}

        for group in self._groups:
            mydata[group] = {}
            mydata[group][DISP] = {}
            mydata[group][VELO] = {}
            mydata[group][ACCEL] = {}
            for typ in globdat[gn.DOFSPACE].get_types():
                mydata[group][DISP][typ] = []
                mydata[group][VELO][typ] = []
                mydata[group][ACCEL][typ] = []

        globdat[self._name] = mydata

    def run(self, globdat):

        model = globdat[gn.MODEL]
        nodes = globdat[gn.NSET]
        disp  = globdat[gn.STATE0]
        velo  = globdat[gn.STATE1]
        accel = globdat[gn.STATE2]
        dofs = globdat[gn.DOFSPACE]
        types = dofs.get_types()

        mydata = globdat[self._name]

        if globdat[gn.ACCEPTED]:
            for group in self._groups:
                for typ in types:
                    idofs = dofs.get_dofs(globdat[gn.NGROUPS][group], [typ])
                    mydata[group][DISP][typ].append( np.mean(disp[idofs]) )
                    mydata[group][VELO][typ].append( np.mean(velo[idofs]) )
                    mydata[group][ACCEL][typ].append( np.mean(accel[idofs]) )

        return 'ok'

    def shutdown(self, globdat):
        for group in self._groups:
            groupData = globdat[self._name][group]
            for typ in globdat[gn.DOFSPACE].get_types():
                groupData[DISP][typ] = np.array( groupData[DISP][typ] )
                groupData[VELO][typ] = np.array( groupData[VELO][typ] )
                groupData[ACCEL][typ] = np.array( groupData[ACCEL][typ] )

def declare(factory):
    factory.declare_module('Acceleration', AccelerationModule)
