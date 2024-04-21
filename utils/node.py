import numpy as np

class Node:
    def __init__(self, coords):
        self._coords = np.array(coords, dtype=float)

    def rank(self):
        return len(self._coords)

    def get_coords(self):
        return self._coords

    def set_coords(self, coords):
        if len(coords) != self.rank():
            raise ValueError('set_coords cannot change the node rank')
        self._coords = np.array(coords, dtype=float)
