import numpy as np
import matplotlib.pyplot as plt

from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

from modules.module import Module

from utils import proputils as pu

XDATA = 'xData'
YDATA = 'yData'
LEGEND = 'legend'
XLABEL = 'xlabel'
YLABEL = 'ylabel'

class GraphModule(Module):

    def init(self, props, globdat):

        if self._name in props:
            myprops = props.get(self._name)

            self._xdata = pu.parse_list(myprops[XDATA])
            self._ydata = pu.parse_list(myprops[YDATA])
            self._legend = []
            self._xlabel = ''
            self._ylabel = ''
            if LEGEND in myprops:
                self._legend = pu.parse_list(myprops[LEGEND])
            if XLABEL in myprops:
                self._xlabel = myprops.get(XLABEL)
            if YLABEL in myprops:
                self._ylabel = myprops.get(YLABEL)

    def run(self, globdat):
        return 'ok'
                
    def shutdown(self, globdat):
        self._make_graph(globdat)

    def _make_graph(self, globdat):
        fig = plt.figure()
        for xdat, ydat in zip(self._xdata, self._ydata):
            [module, group, ld, typ] = xdat.split('.')
            x = globdat[module][group][ld][typ]
            [module, group, ld, typ] = ydat.split('.')
            y = globdat[module][group][ld][typ]
            plt.plot(x, y, '.-')
            if len(self._legend) > 0:
                plt.legend(self._legend)
            if len(self._xlabel) > 0:
                plt.xlabel(self._xlabel)
            if len(self._ylabel) > 0:
                plt.ylabel(self._ylabel)
        plt.show()

def declare(factory):
    factory.declare_module('Graph', GraphModule)
