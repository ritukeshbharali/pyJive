import declare
from names import GlobNames as gn

def jive(props, globdat=None):
    # Initialize global database, declare models and modules

    if not globdat:
        globdat = {}
        declare.declare_models(globdat)
        declare.declare_modules(globdat)
        declare.declare_shapes(globdat)

    # Build main Module chain
    print('Initializing module chain...')
    modulefac = globdat[gn.MODULEFACTORY]
    modelfac = globdat[gn.MODELFACTORY]

    chain = []

    for name in props:
        # Get the name of each item in the property file
        if 'type' in props[name]:
            typ = props[name]['type']
        else:
            typ = name.title()

        # If it refers to a module (and not to a model), add it to the chain
        if modulefac.is_module(typ):
            chain.append(modulefac.get_module(typ, name))
        elif not modelfac.is_model(typ):
            raise ValueError('%s is neither a module nor a model' % typ)

    # Initialize chain
    for module in chain:
        module.init(props, globdat)

    # Run chain until one of the modules ends the computation
    print('Running chain...')

    keep_going = True

    while keep_going:
        for module in chain:
            if 'exit' in module.run(globdat):
                keep_going = False

    # Run postprocessing routines
    for module in chain:
        module.shutdown(globdat)

    print('End of execution')

    return globdat
