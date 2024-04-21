import numpy as np

from modules.module import Module

from names import GlobNames as gn

NSTEPS = 'nsteps'
STOREHISTORY = 'storeHistory'

class ControlModule(Module):

    def init(self, props, globdat):
        myprops = props[self._name]
        self._step = 0
        self._nsteps = int(myprops.get(NSTEPS, 1))
        self._store_history = bool(eval(myprops.get(STOREHISTORY, 'True')))
        globdat[gn.ACCEPTED] = True

        dc = globdat[gn.DOFSPACE].dof_count()

        if self._step == 0:
            globdat[gn.STATE0] = np.zeros(dc)
            globdat[gn.OLDSTATE0] = np.zeros(dc)

    def advance(self, globdat):
        globdat[gn.TIMESTEP] = self._step

        if self._nsteps >= 100:
            if  self._step % ( self._nsteps // 100 ) == 0 :
                print('Running time step', self._step)
        else:
            print('Running time step', self._step)

        if globdat[gn.ACCEPTED]:
            globdat[gn.OLDSTATE0] = np.copy(globdat[gn.STATE0])

    def run(self, globdat):
        if globdat[gn.ACCEPTED]:
            self._step += 1
            if self._store_history:
                if gn.HISTORY not in globdat:
                    globdat[gn.HISTORY] = np.array([globdat[gn.STATE0]])
                else:
                    globdat[gn.HISTORY] = np.vstack((globdat[gn.HISTORY], globdat[gn.STATE0]))

        if self._step >= self._nsteps:
            return 'exit'
        else:
            return 'ok'

        return 'ok'

