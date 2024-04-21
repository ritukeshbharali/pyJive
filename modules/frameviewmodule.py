import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.text as text
from matplotlib.widgets import Slider
import warnings

from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

from modules.module import Module

from utils.constrainer import Constrainer
from utils.table import Table
from utils.xtable import XTable

LINEWIDTH = 'linewidth'
DEFORM = 'deform'
INTERACTIVE = 'interactive'
LABEL = 'label'
MAXSTEP = 'maxStep'
DEFAULTSTEP = 'step0'
MARKERSIZE = 'markerSize'
VECTORSIZE = 'vectorSize'
PLOTSTRESS = 'plotStress'
PLOTUNDEFORMED = 'plotUndeformed'
PLOTDIRICHLET = 'plotDirichlet'
PLOTNEUMANN = 'plotNeumann'
PLOTPHINGES = 'plotPlasticHinges'
STRESS = 'stress'
MODEL = 'model'
FRAME = 'Frame'
PLASTIC = "plastic"
DIRICHLET = 'Dirichlet'
NEUMANN = 'Neumann'
TYPE = 'type'

class FrameViewModule(Module):

    def init(self, props, globdat):
        myprops = props[self._name]
        self._lw = float(myprops.get(LINEWIDTH, 0.5))
        self._scale = float(myprops.get(DEFORM, 1.0))
        self._interactive = bool(eval(myprops.get(INTERACTIVE, 'True')))
        self._step = 0
        self._label = myprops.get(LABEL, 'Step')
        self._maxStep = int(myprops.get(MAXSTEP, -1))
        self._defaultStep = int(myprops.get(DEFAULTSTEP, max(self._maxStep,0)))
        self._msize = float(myprops.get(MARKERSIZE, 4))
        self._vsize = float(myprops.get(VECTORSIZE, -1))
        self._plotStress = myprops.get(PLOTSTRESS,'').upper()
        self._plotUndeformed = bool(eval(myprops.get(PLOTUNDEFORMED, 'True')))
        self._plotNeumann = bool(eval(myprops.get(PLOTNEUMANN, 'True')))
        self._plotDirichlet = bool(eval(myprops.get(PLOTDIRICHLET, 'True')))
        self._plotPlasticHinges = bool(eval(myprops.get(PLOTPHINGES, 'True')))
        self._diriprops = {}
        self._neumprops = {}
        
        for k in props[MODEL]:
            modelprops = props[MODEL][k]
            if ( isinstance(modelprops, dict) )  and ( TYPE in modelprops ):
                if modelprops[TYPE] == DIRICHLET:
                    if self._diriprops:
                        warnings.warn('FrameViewModule: ignoring second or later Dirichlet model')
                    else:
                        self._diriprops = modelprops
                if modelprops[TYPE] == NEUMANN:
                    if self._neumprops:
                        warnings.warn('FrameViewModule: ignoring second or later Neumann model')
                    else:
                        self._neumprops = modelprops
                if modelprops[TYPE] == FRAME:
                    if not PLASTIC in modelprops:
                        self._plotPlasticHinges = False

        if len(self._plotStress) > 0:
            if len(self._plotStress) > 1 or self._plotStress not in 'NVM':
                raise SyntaxError(PLOTSTRESS + "should be 'N', 'V' or 'M'")

        self._modelname = myprops.get(gn.MODEL, gn.MODEL)

    def run(self, globdat):
        return 'ok'

    def shutdown(self, globdat):
        nodes = globdat[gn.NSET]
        elems = globdat[gn.ESET]
        dofs = globdat[gn.DOFSPACE]
    
        if gn.HISTORY in globdat:
            disp = globdat[gn.HISTORY]
            if self._defaultStep >= len(globdat[gn.HISTORY]):
                warnings.warn(f'using last step as default: {len(globdat[gn.HISTORY])-1}')
                self._defaultStep = len(globdat[gn.HISTORY])-1
        else:
            assert not self._interactive, "Cannot use interactive plot without stored history"
            disp = np.array([globdat[gn.STATE0]]).reshape(1, len(globdat[gn.STATE0]))
            self._defaultStep = 0

        if self._maxStep != -1:
            self._nsteps_plot = self._maxStep
            assert gn.HISTORY in globdat, "No history stored!"
        elif gn.HISTORY in globdat:
            self._nsteps_plot = len(globdat[gn.HISTORY])
        else:
            self._nsteps_plot = 1
        
        x0 = np.zeros(len(nodes))
        y0 = np.zeros(len(nodes))

        for n, node in enumerate(nodes):
            coords = node.get_coords()
            x0[n] = coords[0]
            y0[n] = coords[1]

        xu = np.zeros((self._nsteps_plot, len(nodes)))
        yu = np.zeros((self._nsteps_plot, len(nodes)))

        xp = []
        yp = []
        
        for s in range(self._nsteps_plot):
            disp_s = disp[s]
            xu[s] = np.copy(x0)
            yu[s] = np.copy(y0)
            for n in range(len(nodes)):
                xu[s, n] += self._scale * disp_s[dofs.get_dof(n, 'dx')]
                yu[s, n] += self._scale * disp_s[dofs.get_dof(n, 'dy')]

            for i, elem in enumerate(elems):
                inodes = elem.get_nodes()
                if i == 0:
                    x = xu[s, inodes]
                    y = yu[s, inodes]
                else:
                    x = np.hstack((x, xu[s, inodes]))
                    y = np.hstack((y, yu[s, inodes]))
                x = np.hstack((x,np.nan))
                y = np.hstack((y,np.nan))
            
            xp.append(x)
            yp.append(y)
        
        if not self._interactive:

            disp_s = disp[self._defaultStep]

            fig, ax = plt.subplots()
            plt.subplots_adjust(left=0.1, bottom=0., right=0.9, top=1.)
            plt.ion()
            plt.cla()
            plt.axis('equal')
            plt.axis('off')

            plt.plot(xp[self._defaultStep], yp[self._defaultStep], 'k-o', linewidth=self._lw, markersize=self._msize, zorder=3)
            
            if self._vsize == -1:
                xmin, xmax, ymin, ymax = plt.axis()
                plotsize = np.min(np.abs([xmax - xmin, ymax - ymin]))
                self._vsize = 0.25 * plotsize

            if self._plotUndeformed:
                self._plot_undeformed(globdat, elems, x0, y0)
            
            if self._plotDirichlet or self._plotNeumann:
                self._plot_boundaries(globdat, xu, yu, step=self._defaultStep, init=True)
            
            if self._plotPlasticHinges:
                self._plot_hinges(globdat, xu, yu, step=self._defaultStep, init=True)
            
            if len(self._plotStress) > 0:
                for istep in range(len(globdat[gn.HISTORY])):
                    self._write_table(istep, globdat)
                self.plot_stress(self._defaultStep, globdat, xu, yu, ax)
            
            plt.show(block=False)

        elif self._interactive:
            if gn.HISTORY not in globdat:
                raise RuntimeError('FrameViewModule:' + gn.HISTORY + ' has not been defined')

            if len(self._plotStress) > 0:
                for istep in range(len(globdat[gn.HISTORY])):
                    self._write_table(istep, globdat)

            if self._defaultStep < 0 or self._defaultStep >= len(globdat[gn.HISTORY]):
                print('using last step ',len(globdat[gn.HISTORY])-1)
                self._defaultStep = len(globdat[gn.HISTORY])-1

            fig, ax = plt.subplots()
            plt.subplots_adjust(left=0.1, bottom=0., right=0.9, top=1.)
            plt.ion()
            plt.cla()
            plt.axis('equal')
            plt.axis('off')
            
            line, = plt.plot(xp[self._defaultStep], yp[self._defaultStep], 'k-o', linewidth=self._lw, markersize=self._msize, zorder=3)

            if self._vsize == -1.0:
                xmin, xmax, ymin, ymax = plt.axis()
                plotsize = np.max(np.abs([xmax - xmin, ymax - ymin]))
                self._vsize = 0.25 * plotsize
            
            if self._plotDirichlet or self._plotNeumann:
                self._plot_boundaries(globdat, xu, yu, step=self._defaultStep, init=True)
            
            if self._plotUndeformed:
                self._plot_undeformed(globdat, elems, x0, y0)

            if self._plotPlasticHinges:
                self._plot_hinges(globdat, xu, yu, step=self._defaultStep, init=True)
            
            if len(self._plotStress) > 0:
                self.plot_stress(self._defaultStep, globdat, xu, yu, ax)

            # Slider axes
            axcolor = 'lightgoldenrodyellow'
            axstep = plt.axes([0.1, 0.1, 0.8, 0.03], facecolor=axcolor)

            # Create slider object
            if  self._maxStep < 0 :
                self._maxStep = len(globdat[gn.HISTORY])

            s_step = Slider(ax=axstep, label=self._label, valmin=0.,
                            valmax=self._maxStep-1, valinit=self._defaultStep, valstep=1.)

            # Slider function. Updates drawing with new x,y coordinates
            def update(val):
                # Remove old patches
                for p in list(ax.patches):
                    if not isinstance(p, matplotlib.patches.FancyArrow):
                        p.set_visible(False)
                        p.remove()
                for child in ax.get_children():
                    if isinstance(child, text.Annotation):
                        child.set_visible(False)
                        child.remove()

                step = int(s_step.val)

                if self._plotDirichlet or self._plotNeumann:
                    self._plot_boundaries(globdat, xu, yu, step)
                
                line.set_xdata(xp[step])
                
                line.set_ydata(yp[step])
                
                if self._plotPlasticHinges:
                    self._plot_hinges(globdat, xu, yu, step)

                if len(self._plotStress) > 0:
                    self.plot_stress(step, globdat, xu, yu, ax)
                
                fig.canvas.draw_idle()

            s_step.on_changed(update)
            plt.show()

            if gn.SLIDERS not in globdat:
                globdat[gn.SLIDERS] = []

            globdat[gn.SLIDERS].append(s_step)

    def _plot_undeformed(self, globdat, elems, x0, y0):
        for i, elem in enumerate(elems):
            inodes = elem.get_nodes()
            if i == 0:
                x = x0[inodes]
                y = y0[inodes]
            else:
                x = np.hstack((x, x0[inodes]))
                y = np.hstack((y, y0[inodes]))
            x = np.hstack((x,np.nan))
            y = np.hstack((y,np.nan))

            line, = plt.plot(x, y, color='lightgray', marker="o", linewidth=self._lw, markersize=self._msize, zorder=1)
    
    def _plot_hinges(self, globdat, xu, yu, step, init=False):
        if pn.HINGENODES in globdat.keys():
            if init:
                self._hplotdict = {}
                self._hplotdict["phinge_nodes"] = globdat[pn.HINGENODES]
                self._hplotdict["phinge_steps"] = globdat[pn.HINGESTEPS]
                self._hplotdict["phinge_x"] = []
                self._hplotdict["phinge_y"] = []
            
                for s in range(self._nsteps_plot): 
                    hx = []
                    hy = []
                    for n, hnode in enumerate(self._hplotdict["phinge_nodes"]):
                        if self._hplotdict["phinge_steps"][n] + 1 <= step:
                            hx.append(xu[s, hnode])
                            hy.append(yu[s, hnode])
                    self._hplotdict["phinge_x"].append(hx)
                    self._hplotdict["phinge_y"].append(hy)
            
                self.hmarkers, = plt.plot(self._hplotdict["phinge_x"][step], self._hplotdict["phinge_y"][step],
                                         "o", markeredgecolor="black", markeredgewidth=0.5*self._msize,
                                          markerfacecolor="white", markersize=self._msize*1.5, zorder=5)
            
            else:
                self.hmarkers.set_xdata(self._hplotdict["phinge_x"][step])
                self.hmarkers.set_ydata(self._hplotdict["phinge_y"][step])
        
    def _plot_boundaries(self, globdat, xu, yu, step, init=False):
        nodegroups = globdat[gn.NGROUPS]

        if step == 'end':
            disp_s = globdat[gn.STATE0]
            step = globdat[gn.TIMESTEP]
        else:
            disp_s = globdat[gn.HISTORY][step - 1, :]
        
        if step == 0:
            disp_s = np.zeros(len(disp_s))
        
        # Dirichlet
        if self._diriprops and self._plotDirichlet:
            if init:
                self._dplotdict = {"dx": [], "dy": [], "phi": []}
                dirigroups = self._diriprops['groups'][1:-1].split(",") # list from str
                diridofs = self._diriprops['dofs'][1:-1].split(",")
                for i, dgroup in enumerate(dirigroups):
                    dnodes = nodegroups[dgroup]
                    ddof = diridofs[i]
                    for dnode in dnodes:
                        self._dplotdict[ddof].append(dnode)
                
                for dof in ["dx", "dy", "phi"]:
                    self._dplotdict[f"{dof}_x"] = np.zeros((self._nsteps_plot, len(self._dplotdict[f"{dof}"])))
                    self._dplotdict[f"{dof}_y"] = np.zeros((self._nsteps_plot, len(self._dplotdict[f"{dof}"])))

                    for s in range(self._nsteps_plot):
                        for n, dnode in enumerate(self._dplotdict[f"{dof}"]): # n = index from 0 to no. of nodes with diri x (- 1), dxnode is no. of node in mesh
                            self._dplotdict[f"{dof}_x"][s, n] = xu[s, dnode]
                            self._dplotdict[f"{dof}_y"][s, n] = yu[s, dnode]
            
                # could also loop over [dx, dy, phi] here but would require lists for markers, sizes, less room for adjustment
                self._dplotdict["dx_markers"], = plt.plot(self._dplotdict["dx_x"][step], self._dplotdict["dx_y"][step], 
                                            marker=5, linestyle="None", color="turquoise", markersize=self._msize*3, zorder=2)
                self._dplotdict["dy_markers"], = plt.plot(self._dplotdict["dy_x"][step], self._dplotdict["dy_y"][step], 
                                            marker=6, linestyle="None", color="turquoise", markersize=self._msize*3, zorder=2)
                self._dplotdict["phi_markers"], = plt.plot(self._dplotdict["phi_x"][step], self._dplotdict["phi_y"][step], 
                                            marker="s", linestyle="None", color="turquoise", markersize=self._msize*2, zorder=2)
                
            else:
                for dof in ["dx", "dy", "phi"]:
                    self._dplotdict[f"{dof}_markers"].set_xdata(self._dplotdict[f"{dof}_x"][step])
                    self._dplotdict[f"{dof}_markers"].set_ydata(self._dplotdict[f"{dof}_y"][step])

        # Neumann
        if self._neumprops and self._plotNeumann:
            if init:
                self._nplotdict = {"dx_nodes": [], "dy_nodes": [], "phi_nodes": [],
                                   "dx_val": [], "dy_val": [], "phi_val": [],
                                   "dx_inc": [], "dy_inc": [], "phi_inc": [],
                                   "dx_markers": [], "dy_markers": [], "phi_markers": [],
                                   "dx_scale": [], "dy_scale": [], "phi_scale": []}
                neumgroups = self._neumprops['groups'][1:-1].split(",") # list from str
                neumdofs = self._neumprops['dofs'][1:-1].split(",")
                neumvals = self._neumprops['values'][1:-1].split(",")
                if "loadIncr" in self._neumprops:
                    neumincrs = self._neumprops['loadIncr'][1:-1].split(",")
                else:
                    neumincrs = np.zeros(len(neumvals))
                
                neumvals = [float(nv) for nv in neumvals]
                neumincrs = [float(ni) for ni in neumincrs]
                
                for i, ngroup in enumerate(neumgroups):
                    nnodes = nodegroups[ngroup]
                    dof = neumdofs[i]
                    for nnode in nnodes:
                        self._nplotdict[f"{dof}_nodes"].append(nnode)
                        self._nplotdict[f"{dof}_val"].append(float(neumvals[i]))
                        self._nplotdict[f"{dof}_inc"].append(float(neumincrs[i]))
                
                initial_forces = np.array(self._nplotdict["dx_val"] + self._nplotdict["dy_val"])
                increment_forces = np.array(self._nplotdict["dx_inc"] + self._nplotdict["dy_inc"])
                initial_moments = np.array(self._nplotdict["phi_val"])
                increment_moments = np.array(self._nplotdict["phi_inc"])
                
                for dof in ["dx", "dy", "phi"]:
                    self._nplotdict[f"{dof}_x"] = np.zeros((self._nsteps_plot, len(self._nplotdict[f"{dof}_nodes"])))
                    self._nplotdict[f"{dof}_y"] = np.zeros((self._nsteps_plot, len(self._nplotdict[f"{dof}_nodes"])))
                    self._nplotdict[f"{dof}_scale"] = np.zeros((self._nsteps_plot, len(self._nplotdict[f"{dof}_nodes"])))

                for s in range(self._nsteps_plot):
                    for dof in neumdofs:
                        for n, dnode in enumerate(self._nplotdict[f"{dof}_nodes"]):
                            self._nplotdict[f"{dof}_x"][s, n] = xu[s, dnode]
                            self._nplotdict[f"{dof}_y"][s, n] = yu[s, dnode]
            
                    if len(initial_forces) != 0:
                            total_forces = initial_forces + s * increment_forces
                            relative_forces = total_forces / np.max(np.abs(total_forces))
                            relative_forces = np.nan_to_num(relative_forces)
                            self._nplotdict["dx_scale"][s] = self._vsize * relative_forces[0:len(self._nplotdict["dx_nodes"])]
                            self._nplotdict["dy_scale"][s] = self._vsize * relative_forces[len(self._nplotdict["dx_nodes"]):]
        
                    if len(initial_moments) != 0:
                        total_moments = initial_moments + s * increment_moments
                        relative_moments =  total_moments / np.max(np.abs(total_moments))
                        relative_moments = np.nan_to_num(relative_moments)
                        self._nplotdict["phi_scale"][s] = self._vsize * relative_moments
            
                for i, nxnode in enumerate(self._nplotdict["dx_nodes"]):
                    length = self._nplotdict["dx_scale"][step, i]
                    arrow = plt.arrow(x=self._nplotdict["dx_x"][step, i] - length, y=self._nplotdict["dx_y"][step, i], 
                                      dx=length, dy=0, color="firebrick",
                                      width=1/100 * self._vsize, head_width=1/15 * self._vsize,
                                      length_includes_head=True, zorder=10)
                    self._nplotdict["dx_markers"].append(arrow)
                for i, nynode in enumerate(self._nplotdict["dy_nodes"]):
                    length = self._nplotdict["dy_scale"][step, i]
                    arrow = plt.arrow(x=self._nplotdict["dy_x"][step, i], y=self._nplotdict["dy_y"][step, i] - length, 
                                      dx=0, dy=length, color="firebrick",
                                      width=1/100 * self._vsize, head_width=1/15 * self._vsize,
                                      length_includes_head=True, zorder=10)
                    self._nplotdict["dy_markers"].append(arrow)
                theta = np.linspace(0.25*np.pi, 1.75*np.pi, 100)
                xc = np.cos(theta)
                yc = np.sin(theta)
                for i, nphinode in enumerate(self._nplotdict["phi_nodes"]):
                    radius = 0.25 * self._nplotdict["phi_scale"][step, i]
                    theta = np.linspace(0.25*np.pi, 1.75*np.pi, 100)
                    x_mom = self._nplotdict["phi_x"][step, i] + radius * xc
                    y_mom = self._nplotdict["phi_y"][step, i] + radius * yc
                    circle, = plt.plot(x_mom, y_mom, color="firebrick")
                    if self._nplotdict["phi_val"][i] + step * self._nplotdict["phi_inc"][i] <= 0:
                        arrow = plt.arrow(x_mom[0], y_mom[0], -0.0001, 0.0001, 
                                          head_width=1/30 * self._vsize, color="firebrick", zorder=10)
                    else:
                        arrow = plt.arrow(x_mom[-1], y_mom[-1], 0.0001, 0.0001, 
                                          head_width=1/30 * self._vsize, color="firebrick", zorder=10)
                    if np.abs(radius) < 1e-11:
                        arrow.set_visible(False)
                    else:
                        arrow.set_visible(True)
                    self._nplotdict["phi_markers"].append([circle, arrow])
          
            else: # not init
                for i, arrow in enumerate(self._nplotdict["dx_markers"]):
                    length = self._nplotdict["dx_scale"][step, i]
                    arrow.set_data(x=self._nplotdict["dx_x"][step, i] - length, y=self._nplotdict["dx_y"][step, i], dx=length)
                for i, arrow in enumerate(self._nplotdict["dy_markers"]):
                    length = self._nplotdict["dy_scale"][step, i]
                    arrow.set_data(x=self._nplotdict["dy_x"][step, i], y=self._nplotdict["dy_y"][step, i] - length, dy=length)
                theta = np.linspace(0.25*np.pi, 1.75*np.pi, 100)
                xc = np.cos(theta)
                yc = np.sin(theta)
                for i, marker in enumerate(self._nplotdict["phi_markers"]):
                    radius = 0.25 * self._nplotdict["phi_scale"][step, i]
                    x_mom = self._nplotdict["phi_x"][step, i] + radius * xc
                    y_mom = self._nplotdict["phi_y"][step, i] + radius * yc
                    marker[0].set_xdata(x_mom)
                    marker[0].set_ydata(y_mom)
                    if self._nplotdict["phi_val"][i] + step * self._nplotdict["phi_inc"][i] <= 0:
                        marker[1].set_data(x=x_mom[0], y=y_mom[0], dx=-0.0001, dy=0.0001)
                    else:
                        marker[1].set_data(x=x_mom[-1], y=y_mom[-1], dx=0.0001, dy=0.0001)
                    if np.abs(radius) < 1e-11:
                        marker[1].set_visible(False)
                    else:
                        marker[1].set_visible(True)
    
    def plot_stress(self, step, globdat, xu, yu, ax):
        elems = globdat[gn.ESET]
        data = globdat[gn.TABLES][STRESS][step][self._plotStress]
        eps = 1.e-8
        data[np.abs(data) < eps] = 0
        dxall = max(xu[step]) - min(xu[step])
        dyall = max(yu[step]) - min(yu[step])
        smax = max(abs(data))
        scalef = -0.1*max(dxall,dyall)/smax
        XY = np.zeros((len(elems), 4, 2))

        for i, elem in enumerate(elems):
            inodes = elem.get_nodes()
            p0 = np.array([xu[step, inodes[0]],yu[step, inodes[0]]])
            p1 = np.array([xu[step, inodes[1]],yu[step, inodes[1]]])
            dx = p1-p0
            dn = np.array([-dx[1],dx[0]])/np.linalg.norm(dx)
            p2 = p1 + scalef*dn*data[i*2+1]
            p3 = p0 + scalef*dn*data[i*2]

            XY[i,:,:] = np.array([p0,p1,p2,p3])

            compound = path.Path.make_compound_path_from_polys(XY)
            pathpatch = patches.PathPatch(compound)
            patch = ax.add_patch(pathpatch)
            patch.set(edgecolor=patch.get_facecolor())

            if inodes[0] in globdat[gn.MASTERNODES]:
                pmid = 0.5*(p0+p3)
                ax.annotate("{:.4f}".format(data[i*2]),xy=pmid)
            if inodes[1] in globdat[gn.MASTERNODES]:
                pmid = 0.5*(p1+p2)
                ax.annotate("{:.4f}".format(data[i*2+1]),xy=pmid)

    def _write_table(self, istep, globdat):
        # Element table does not require weights
        model = globdat[self._modelname]
        name = STRESS

        if gn.TABLES not in globdat:
            globdat[gn.TABLES] = {}

        if name not in globdat[gn.TABLES]:
            globdat[gn.TABLES][name] = []

        params = {}
        params[pn.TABLE] = Table()
        params[pn.TABLENAME] = name
        params[gn.STATE0] = globdat[gn.HISTORY][istep]

        model.take_action(act.GETTABLE, params, globdat)

        table = params[pn.TABLE]
        tbwts = params[pn.TABLEWEIGHTS]

        # convert table from Table to XTable
        tblcls = table.__class__
        table.__class__ = XTable

        # divide cols by node weights
        for icol in range(table.column_count()):
            table.set_col_values(None, icol,
                                 table.get_col_values(None, icol)/params[pn.TABLEWEIGHTS])

        # convert table back
        table.__class__ = tblcls

        if len(globdat[gn.TABLES][name]) == istep:
            globdat[gn.TABLES][name] = np.append(globdat[gn.TABLES][name], table)
        elif len(globdat[gn.TABLES][name]) > istep:
            globdat[gn.TABLES][name][istep] = table
        else:
            raise RuntimeError('cannot store table for time step %i' % istep)

def declare(factory):
    factory.declare_module('FrameView', FrameViewModule)
