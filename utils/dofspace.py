class DofSpace:
    def __init__(self):
        self._dofs = {}
        self._count = 0

    def add_type(self, name):
        if name not in self._dofs:
            self._dofs[name] = {}

    def add_dof(self, inode, tname):
        if inode not in self._dofs[tname]:
            self._dofs[tname][inode] = self._count
            self._count += 1

    def type_count(self):
        return len(self._dofs)

    def dof_count(self):
        return self._count

    def get_types(self):
        return self._dofs.keys()

    def find_dof(self, inode, typ):
        return self._dofs[typ].get(inode)

    def get_dof(self, inode, typ):
        idof = self.find_dof(inode, typ)
        if idof is None:
            raise RuntimeError('DofSpace: Non-existent DOF')
        return idof

    def get_dofs(self, inodes, types):
        idofs = []

        for i in inodes:
            for t in types:
                idof = self._dofs[t].get(i)
                if idof is not None:
                    idofs.append(idof)

        return idofs

    def set_dof(self, oldnode, newnode, typ):
        dof = self.get_dof(oldnode, typ)
        self._dofs[typ][newnode] = dof
