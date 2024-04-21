import sys
sys.path.append('../../')

from math import exp
import matplotlib.pyplot as plt
import numpy as np
import main
from utils import proputils as pu

def mesher_lin(L, n):
    dx = L / n
    with open('bar.mesh', 'w') as fmesh:
        fmesh.write('nodes (ID, x, [y], [z])\n')
        for i in range(n + 1):
            fmesh.write('%d %f\n' % (i, i * dx))
        fmesh.write('elements (node#1, node#2, [node#3, ...])\n')
        for i in range(n):
            fmesh.write('%d %d\n' % (i, i + 1))

props = pu.parse_file('bar.pro')

P = 1
L = 10

mesher_lin(L,4)

globdat = main.jive(props)

K = globdat['matrix0']
u = globdat['state0']

print('Solution',u)
print('Stiffness matrix\n',K)
