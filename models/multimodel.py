import numpy as np

from models.model import Model
from utils import proputils as pu

MODELS = 'models'
TYPE = 'type'

class MultiModel(Model):
    def take_action(self, action, params, globdat):
        for model in self._models:
            model.take_action(action, params, globdat)

    def configure(self, props, globdat):
        models = pu.parse_list(props[MODELS])

        mfac = globdat['modelFactory']

        self._models = []

        for m in models:
            modelprops = props[m]

            model = mfac.get_model(modelprops[TYPE], m)
            model.configure(modelprops, globdat)

            self._models.append(model)


def declare(factory):
    factory.declare_model('Multi', MultiModel)
