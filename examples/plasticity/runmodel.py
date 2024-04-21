import sys

sys.path.append('../../')

import numpy as np
import main

from utils import proputils as pu

props = pu.parse_file('simple.pro')

globdat = main.jive(props)
u = globdat['state0']

disps = globdat['lodi']['right']['disp']['dx']
loads = globdat['lodi']['right']['load']['dx']

with open ('loaddisp.dat','w') as lodi:
    for load, disp in zip(loads,disps):
        lodi.write(str(disp) + ' ' + str(load) + '\n')

print('Load-displacement data written to loaddisp.dat')
