from names import GlobNames as gn

from modules.module import Module

class OutputModule(Module):

    def init(self, props, globdat):
        pass

    def run(self, globdat):
        # Temporary strat
        fname = 'step' + str(globdat[gn.TIMESTEP]) + '.disp'
        u = globdat[gn.STATE0]
        print(u)

        with open(fname, 'w') as out:
            for val in u:
                out.write(str(val) + '\n')

        return 'ok'

    def shutdown(self, globdat):
        pass


def declare(factory):
    factory.declare_module('Output', OutputModule)
