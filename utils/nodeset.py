import numpy as np

from utils.itemset import ItemSet, XItemSet
from utils.node import Node

class NodeSet(ItemSet):

    def __init__(self, nodes=None):
        super().__init__(nodes)
        self._rank = 0

    def rank(self):
        return self._rank

    def find_node(self, node_id):
        return self.find_item(node_id)

    def find_nodes(self, node_ids):
        return self.find_items(node_ids)

    def get_node_id(self, inode):
        return self.get_item_id(inode)

    def get_node_ids(self, inodes):
        return self.get_item_ids(inodes)

    def get_node_coords(self, inode):
        return self._data[inode].get_coords()

    def get_coords(self):
        coords = []
        for inode in range(self.size()):
            coords.append(self.get_node_coords(inode))
        return np.array(coords).T

    def get_some_coords(self, inodes):
        coords = []
        for inode in inodes:
            coords.append(self.get_node_coords(inode))
        return np.array(coords).T


class XNodeSet(NodeSet, XItemSet):

    def add_node(self, coords, node_id=None):
        node = Node(coords)
        if self.size() == 0:
            self._rank = node.rank()
        else:
            assert self._rank == node.rank()
        self.add_item(node, node_id)

    def erase_node(self, inode):
        self.erase_item(inode)

    def set_node_coords(self, inode, coords):
        self._data[inode].set_coords(coords)

    def set_coords(self, coords):
        if coords.shape[0] != self.size():
            raise ValueError('first dimension of coords does not match the number of nodes')
        for inode in range(self.size()):
            self.set_node_coords(inode, coords[inode])

    def set_some_coords(self, inodes, coords):
        if coords.shape[0] != self.size():
            raise ValueError('first dimension of coords does not match the size of inodes')
        for i, inode in enumerate(inodes):
            self.set_node_coords(inode, coords[i])

    def to_nodeset(self):
        self.__class__ = NodeSet
        return self

def to_xnodeset(nodes):
    nodes.__class__ = XNodeSet
    return nodes
