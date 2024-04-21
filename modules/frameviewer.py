import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.text as text
from matplotlib.widgets import Slider
import warnings

from modules.frameviewmodule import FrameViewModule
from names import GlobNames as gn

class FrameViewer(FrameViewModule):
    def __init__(self, props, globdat):
        self._name = "frameview" #TODO: find name from props so it becomes case-insensitive
        self._globdat = globdat
        self._props = props
        super().init(self._props, self._globdat)
    
    def _update_props(self, new_props):
        self._props[self._name].update(new_props)
        super().init(self._props, self._globdat)
    
    def _set_title(self, title=None, track=None):
        def update_title(val):
            self._title = title + f"{self._globdat[gn.HISTORY][int(val -1)][self._trackdof]:.3e}" 
            plt.suptitle(self._title)
        
        if track:
            if self._storeHistory:
                self._track = track
                self._trackdof = self._globdat[gn.DOFSPACE].get_dof(*self._track)
                self._title = title + f"{self._globdat[gn.HISTORY][int(self._defaultStep)][self._trackdof]:.3e}"
                plt.suptitle(self._title)
                self._slider.on_changed(update_title)
            
            else:
                self._track = track
                self._trackdof = self._globdat[gn.DOFSPACE].get_dof(*self._track)
                self._title = title + f"{self._globdat[gn.STATE0][self._trackdof]:.3e}"
                plt.suptitle(self._title)
            
        else:
            self._title = title
            plt.suptitle(self._title)
    
    def save(self, path, step='current'):
        if type(step) in [int, float]:
            self._slider.set_val(int(step))
        plt.draw()
        self._fig.canvas.draw()
        self._fig.canvas.flush_events()
        plt.savefig(path)
    
    @staticmethod
    def plot(props, globdat, new_props=None, **kwargs):
        # get kwargs
        title = kwargs.get('title', None)
        track = kwargs.get('track', None)
        save = kwargs.get('save', None)

        FV = FrameViewer(props, globdat)
        if type(new_props) == dict:
            FV._update_props(new_props)
        FV.shutdown(globdat)

        FV._fig = plt.gcf()
        FV._slider = FV._globdat[gn.SLIDERS][-1]

        if title:
            FV._set_title(title, track)
        
        if save:
            FV.save(save, FV._defaultStep)
        
        return FV

        # currently not used, but might be useful for other implementations
        # fig = plt.gcf()
        # plotax, sliderax = fig.get_axes()
        
        # fig.canvas.draw()
        # fig.canvas.flush_events()
