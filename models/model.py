class ModelFactory:
    def __init__(self):
        self._creators = {}

    def declare_model(self, typ, creator):
        self._creators[typ] = creator

    def get_model(self, typ, name):
        creator = self._creators.get(typ)
        if not creator:
            raise ValueError(typ)
        return creator(name)

    def is_model(self, typ):
        return typ in self._creators


class Model:
    def __init__(self, name):
        self._name = name

    def take_action(self, action, params, globdat):
        print('Empty model takeAction')

    def configure(self, props, globdat):
        print('Empty model configure')
