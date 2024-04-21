import numpy as np

class Element:
    def __init__(self, nodes):
        self._nodes = np.array(nodes, dtype=int)

    def get_node_count(self):
        return len(self._nodes)

    def get_nodes(self):
        return self._nodes

    def set_nodes(self, nodes):
        if len(nodes) != self.get_node_count():
            raise ValueError('set_nodes cannot change the element node count')
        self._nodes = np.array(nodes, dtype=int)

    def change_node(self, oldnode, newnode):
        self._nodes[self._nodes == oldnode] = newnode
