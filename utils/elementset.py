import numpy as np

from utils.itemset import ItemSet, XItemSet
from utils.element import Element

class ElementSet(ItemSet):

    def __init__(self, nodes, elems=None):
        super().__init__(elems)
        self._nodes = nodes
        self._maxNodeCount = 0

    def find_element(self, elem_id):
        return self.find_item(elem_id)

    def get_elem_id(self, ielem):
        return self.get_item_id(ielem)

    def max_elem_node_count(self):
        return self._maxNodeCount

    def max_elem_node_count_of(self, ielems):
        maxNodeCount = 0
        for ielem in ielems:
            maxNodeCount = max(maxNodeCount, self.get_elem_node_count(ielem))
        return maxNodeCount

    def get_elem_node_count(self, ielem):
        return self._data[ielem].get_node_count()

    def get_elem_nodes(self, ielem):
        return self._data[ielem].get_nodes()

    def get_nodes(self):
        return self._nodes

    def get_some_elem_nodes(self, index, inode):
        return np.array(self.get_elem_nodes(inode)[index], dtype=int)

    def get_nodes_of(self, ielems):
        inodes = []
        for ielem in ielems:
            for inode in self.get_elem_nodes(ielem):
                inodes.append(inode)
        return np.array(inodes, dtype=int)

    def get_unique_nodes_of(self, ielems):
        return np.unique(self.get_nodes_of(ielems))


class XElementSet(ElementSet, XItemSet):

    def add_element(self, inodes, elem_id=None):
        elem = Element(inodes)
        nodeCount = elem.get_node_count()

        if self._maxNodeCount < nodeCount:
            self._maxNodeCount = nodeCount

        self.add_item(elem, elem_id)

    def erase_element(self, ielem):
        nodeCount = self.get_elem_node_count(ielem)

        self.erase_item(ielem)

        if nodeCount == self.max_elem_node_count():
            for elem in self._data:
                if elem.get_node_count() == nodeCount:
                    break
            else:
                self._maxNodeCount = self.max_elem_node_count_of(range(self.size()))

    def set_elem_nodes(self, ielem, nodes):
        self._data[ielem].set_nodes(nodes)

    def to_elementset(self):
        self.__class__ = ElementSet
        return self

def to_xelementset(elems):
    elems.__class__ = XElementSet
    return elems
