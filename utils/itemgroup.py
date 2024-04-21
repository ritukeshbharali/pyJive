import numpy as np

from utils.itemset import ItemSet

class ItemGroup():

    def __init__(self, items, data=None):
        if data is None:
            self._data = []
        else:
            self._data = data

        self._items = items

        assert isinstance(self._items, ItemSet)

    def __len__(self):
        return self.size()

    def __iter__(self):
        return iter(self._data)

    def __next__(self):
        return next(self._data)

    def __contains__(self, iitem):
        return self.contains(iitem)

    def size(self):
        return len(self._data)

    def get_ids(self, iitems):
        self._items.get_item_ids(self.get_indices())

    def get_indices(self):
        return np.array(self._data, dtype=int)

    def get_items(self):
        return self._items

    def contains(self, iitem):
        return iitem in self.get_indices()

    def find_members(self, iitems):
        jitems = []
        for iitem in iitems:
            if self.contains(iitem):
                jitems.append(iitem)
        return np.array(jitems, dtype=int)

    def find_non_members(self, iitems):
        jitems = []
        for iitem in iitems:
            if not self.contains(iitem):
                jitems.append(iitem)
        return np.array(jitems, dtype=int)


class XItemGroup(ItemGroup):

    def clear(self):
        self._data = []

    def add_item(self, iitem):
        self._data.append(iitem)

    def add_items(self, iitems):
        for iitem in iitems:
            self.add_item(iitem)

    def erase_item(self, iitem):
        self._data.pop(iitem)

    def erase_items(self, iitems):
        for iitem in iitems:
            self.erase_items(iitem)
