import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.widgets import Slider
import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.text as text

from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

from modules.module import Module

from utils.table import Table
from utils.xtable import XTable

LINEWIDTH ='linewidth'
PLOT = 'plot'
NCOLORS = 'ncolors'
NTICKS = 'nticks'
DEFORM = 'deform'
COLORMAP = 'colorMap'
INTERACTIVE = 'interactive'
MAXSTEP = 'maxStep'
DEFAULTSTEP = 'step0'
LABEL = 'label'
CONSTANTLEVELS = 'constantLevels'
CONSTANTTICKS = 'constantTicks'

class ViewModule(Module):

    def init(self, props, globdat):
        if self._name in props:
            myprops = props.get(self._name)

            self._interactive = bool(eval(myprops.get(INTERACTIVE, 'True')))
            self._maxStep = int(myprops.get(MAXSTEP, -1))
            self._defaultStep = int(myprops.get(DEFAULTSTEP, max(self._maxStep,0)))
            if not gn.HISTORY in globdat:
                self._defaultStep = 0
            elif self._defaultStep >= len(globdat[gn.HISTORY]):
                print('using last step ',len(globdat[gn.HISTORY])-1)
                self._defaultStep = len(globdat[gn.HISTORY])-1
            self._label = myprops.get(LABEL, 'Step')
            self._constantLevels = bool(eval(myprops.get(CONSTANTLEVELS, 'False')))
            self._constantTicks = bool(eval(myprops.get(CONSTANTTICKS, 'True')))

            self._linewidth = float(myprops.get(LINEWIDTH, 0.2))
            self._pname = myprops.get(PLOT, "")
            self._scale = float(myprops.get(DEFORM, 0.0))
            self._ncolors = int(myprops.get(NCOLORS, 100))
            self._nticks = int(myprops.get(NTICKS, 5))
            self._colormap = myprops.get(COLORMAP, 'viridis')

        self._modelname = myprops.get(gn.MODEL, gn.MODEL)
    
    def run(self, globdat):
        return 'ok'

    def shutdown(self, globdat):
        
        assert globdat[gn.MESHSHAPE] == 'Triangle3', 'ViewModule only supports triangles for now'
        
        if self._maxStep != -1:
            self._nsteps_plot = self._maxStep
            assert gn.HISTORY in globdat, "No history stored!"
        elif gn.HISTORY in globdat:
            self._nsteps_plot = len(globdat[gn.HISTORY])
        else:
            self._nsteps_plot = 1
        
        z, levels, ticks, triang = self.fill_data(globdat)

        fig = plt.figure()
        contour_ax = plt.gca()
        plt.ion()
        plt.cla()
        plt.axis('equal')
        plt.axis('off')
        
        def update(val):
            # Reset plot
            contour_ax.clear()
            contour_ax.axis('equal')
            contour_ax.axis('off')
            if self._interactive:
                step = int(s_step.val)
            else:
                step = self._defaultStep
            
            if self._pname != "":
                self._tricontour = contour_ax.tricontourf(triang[step], z[step], levels=levels[step], cmap=self._colormap)
                if self._cbar:
                    self._cbar.remove()
                self._cbar = plt.colorbar(self._tricontour, ticks=ticks[step], ax=contour_ax)
            
            contour_ax.triplot(triang[step],'k-',linewidth=self._linewidth, figure=fig)

            plt.show(block=False)
        
        if not self._interactive:
            self._cbar =  None
            update(self._defaultStep)
            
        elif self._interactive:
            if gn.HISTORY not in globdat:
                raise RuntimeError('FrameViewModule:' + gn.HISTORY + ' has not been defined')
            
            if self._defaultStep < 0 or self._defaultStep >= len(globdat[gn.HISTORY]):
                print('using last step ',len(globdat[gn.HISTORY])-1)
                self._defaultStep = len(globdat[gn.HISTORY])-1
            
            # Slider axes
            axcolor = 'lightgoldenrodyellow'
            axstep = plt.axes([0.1, 0.1, 0.8, 0.03], facecolor=axcolor)

            # Create slider object
            if  self._maxStep < 0 :
                self._maxStep = len(globdat[gn.HISTORY])

            s_step = Slider(ax=axstep, valmin=0., label=self._label,
                            valmax=self._maxStep-1, valinit=self._defaultStep, valstep=1.)

            self._cbar =  None
            update(self._defaultStep)
            
            s_step.on_changed(update)

            if gn.SLIDERS not in globdat:
                globdat[gn.SLIDERS] = []

            globdat[gn.SLIDERS].append(s_step)

        return 'ok'

    def fill_data(self, globdat):
        nodes = globdat[gn.NSET]
        elems = globdat[gn.ESET]
        if gn.HISTORY in globdat:
            disp = globdat[gn.HISTORY]
        else:
            assert not self._interactive, "Cannot use interactive plot without stored history"
            disp = np.array([globdat[gn.STATE0]]).reshape(1, len(globdat[gn.STATE0]))
        dofs = globdat[gn.DOFSPACE]
        types = dofs.get_types()

        x = np.zeros((len(nodes)))
        y = np.zeros(len(nodes))
        el= np.zeros((len(elems),3),dtype=int)

        for n, node in enumerate(nodes):
            coords = node.get_coords()

            x[n] = coords[0]
            y[n] = coords[1]

        for e, elem in enumerate(elems):
            inodes = elem.get_nodes()

            el[e,:] = inodes
        
        plot_dof = False # plot one of the dof values
        plot_other = False # plot another quantitiy
        
        z = np.zeros((self._nsteps_plot, len(nodes))) # Backup in case z gets used but not filled
        
        if self._pname != '':
            if 'solution' in self._pname:
                comp = self._pname.split('[')[1].split(']')[0]
                assert comp in types, 'Invalid DOF type: %s' % comp
                plot_dof = True
            
            else: #TODO: Poisson model does not have GETTABLE, Elastic model does but only for STATE0 -> always takes final value
                # assert not self._interactive, "Table results are only supported in non-interactive mode"
                name = self._pname.split('[')[0]
                comp = self._pname.split('[')[1].split(']')[0]
                self._write_table(name, globdat)
                table = globdat[gn.TABLES][name]
                assert comp in table, 'Invalid component: %s' % comp
                plot_other = True

        dx = np.ones((self._nsteps_plot, len(nodes))) * x
        dy = np.ones((self._nsteps_plot, len(nodes))) * y
        triang = []
        
        levels = np.zeros((self._nsteps_plot, self._ncolors + 1))
        ticks = np.zeros((self._nsteps_plot, self._nticks))

        mesh_displaced = 'dx' in types and 'dy' in types

        for s in range(self._nsteps_plot):
            disp_s = disp[s]
            
            if mesh_displaced:
                for n in range(len(nodes)):
                    idofs = dofs.get_dofs([n], ['dx', 'dy'])
                    du = disp_s[idofs]

                    dx[s, n] += self._scale * du[0]
                    dy[s, n] += self._scale * du[1]

            if plot_dof:
                for n in range(len(nodes)):
                    z[s, n] = disp_s[dofs.get_dof(n,comp)]
            
            elif plot_other:
                for n in range(len(nodes)):
                    z[s, n] = table[comp][n]
    
            triang.append(tri.Triangulation(dx[s], dy[s], el))

            # check if field is zero or uniform per step, fill levels and ticks
            if np.max(np.abs(z[s])) < 1e-6:
                levels[s] = np.linspace(-1e-6, 1e-6, self._ncolors + 1)
                ticks[s] = np.linspace(-1e-6, 1e-6, self._nticks)
            elif np.abs(z[s].max() - z[s].min()) < 1e-6:
                levels[s] = np.linspace(0.9*z[s].max(), 1.1*z[s].max(), self._ncolors + 1)
                ticks[s] = np.linspace(0.9*z[s].max(), 1.1*z[s].max(), self._nticks)
            else:
                levels[s] = np.linspace(z[s].min(), z[s].max(), self._ncolors + 1)
                ticks[s] = np.linspace(z[s].min(), z[s].max(), self._nticks)
        
        if self._constantLevels:
            conslevels = np.linspace(levels.min(), levels.max(), self._ncolors + 1)
            levels = np.ones_like(levels) * conslevels
            if self._constantTicks:
                consticks = np.linspace(levels.min(), levels.max(), self._nticks)
                ticks = np.ones_like(ticks) * consticks
        
        return z, levels, ticks, triang
    
    def plot(self, globdat):
        self.shutdown(globdat)

    def _write_table(self, name, globdat):
        nodes = globdat[gn.NSET]
        model = globdat[self._modelname]

        globdat[gn.TABLES] = {}

        params = {}
        params[pn.TABLE] = Table()
        params[pn.TABLENAME] = name
        params[pn.TABLEWEIGHTS] = np.zeros(len(nodes))

        model.take_action(act.GETTABLE, params, globdat)

        table = params[pn.TABLE]

        # convert table from Table to xTable
        tblcls = table.__class__
        table.__class__ = XTable

        # divide cols by node weights
        for icol in range(table.column_count()):
            table.set_col_values(None, icol,
                                 table.get_col_values(None, icol)/params[pn.TABLEWEIGHTS])

        # convert table back
        table.__class__ = tblcls

        # store table
        globdat[gn.TABLES][name] = params[pn.TABLE]

def declare(factory):
    factory.declare_module('View', ViewModule)
