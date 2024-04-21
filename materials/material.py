
MATERIAL = 'material'
TYPE = 'type'
RANK = 'rank'


def new_material(props):
    typ = props[TYPE]
    rank = int(props[RANK])

    if typ == 'Elastic':
        from materials.elasticmaterial import ElasticMaterial
        mat = ElasticMaterial(rank)
    elif typ == 'J2':
        from materials.j2material import J2Material
        mat = J2Material(rank)
    else:
        raise ValueError(typ + ' is not a valid material')

    return mat


class Material:

    def __init__(self, rank):
        pass

    def configure(self, props, globdat):
        pass

    def get_config(self):
        pass

    def commit(self, ipoint=None):
        pass

    def check_commit(self, params, globdat):
        pass

    def create_material_points(self, npoints):
        pass
