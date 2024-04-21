import numpy as np

from utils.table import Table

class XTable(Table):

    def __init__(self, tbl=None):
        super().__init__()

        if not tbl is None:
            self._data = tbl._data
            self._header = tbl._header

    def clear_data(self):
        self._data = np.zeros((0, self._header.size()))

    def clear_all(self):
        self._header = np.zeros(0)
        self.clear_data()

    def reserve(self, rowcount):
        if self.row_count() < rowcount:
            tmp = np.zeros((rowcount, self.column_count()))
            tmp[:self.row_count(), :self.column_count()] = self._data
            self._data = tmp

    def add_column(self, name):
        if self.find_column(name) < 0:
            self._header = np.append(self._header, name)
            self._data.resize((self.row_count(), self.column_count() + 1))
        return self.find_column(name)

    def add_columns(self, names):
        for name in names:
            self.add_column(name)
        return self.find_columns(names)

    def set_value(self, irow, jcol, value):
        self.reserve(irow+1)
        self._data[irow, jcol] = value

    def add_value(self, irow, jcol, value):
        self.reserve(irow+1)
        self._data[irow, jcol] += value

    def set_block(self, irows, jcols, block):
        self.reserve(max(irows)+1)
        self._data[np.ix_(irows, jcols)] = block

    def add_block(self, irows, jcols, block):
        self.reserve(max(irows)+1)
        self._data[np.ix_(irows, jcols)] += block

    def set_row_values(self, irow, jcols, values):
        self.reserve(irow+1)
        if jcols is None:
            self._data[irow, :] = values
        else:
            self._data[irow, jcols] = values

    def add_row_values(self, irow, jcols, values):
        self.reserve(irow+1)
        if jcols is None:
            self._data[irow, :] += values
        else:
            self._data[irow, jcols] += values

    def set_col_values(self, irows, jcol, values):
        if irows is None:
            self.reserve(len(values))
            self._data[:, jcol] = values
        else:
            self.reserve(max(irows)+1)
            self._data[irows, jcol] = values

    def add_col_values(self, irows, jcol, values):
        if irows is None:
            self.reserve(len(values))
            self._data[:, jcol] += values
        else:
            self.reserve(max(irows)+1)
            self._data[irows, jcol] += values

    def to_table(self):
        self.__class__ = Table
