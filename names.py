class GlobNames:
    MODEL = 'model'
    DOFSPACE = 'dofSpace'
    MESHSHAPE = 'meshShape'
    MESHRANK = 'meshRank'
    ESET = 'elemSet'
    NSET = 'nodeSet'
    NGROUPS = 'nodeGroups'
    EGROUPS = 'elemGroups'
    MODELFACTORY = 'modelFactory'
    MODULEFACTORY = 'moduleFactory'
    SHAPEFACTORY = 'shapeFactory'
    TIMESTEP = 'timeStep'
    TIME     = 'time'
    STATE0 = 'state0'
    STATE1 = 'state1'
    STATE2 = 'state2'
    OLDSTATE0 = 'oldstate0'
    OLDSTATE1 = 'oldstate2'
    OLDSTATE2 = 'oldstate2'
    OLDERSTATE0 = 'olderstate0'
    BACKUPSTATE0 = 'backupstate0'
    HISTORY = 'history'
    MATRIX0 = 'matrix0'
    MATRIX2 = 'matrix2'
    EXTFORCE = 'extForce'
    CONSTRAINTS = 'constraints'
    TABLES = 'tables'
    LAMBDA = 'lambda'
    MODULEDATA = 'module'
    SLIDERS = 'sliders'
    ACCEPTED = 'accepted'
    LBFACTORS = 'lbFactors'
    MASTERNODES = 'masterNodes'
    EIGENFREQS = "eigenFrequencies"
    MODALSHAPES = "modalShapes"
    MASTERNODES = 'masterNodes'
    HINGENODES = 'hingeNodes'
    HINGESTEPS = "hingeSteps"


class PropNames:
    TYPE = 'type'


class Actions:
    GETMATRIX0 = 'getMatrix0'
    GETMATRIX2 = 'getMatrix2'
    GETMATRIXLB = 'getMatrixLB'
    GETEXTFORCE = 'getExtForce'
    GETINTFORCE = 'getIntForce'
    GETUNITFORCE = 'getUnitForce'
    GETINTFORCE = 'getIntForce'
    COMMIT = 'commit'
    CHECKCOMMIT = 'checkCommit'
    GETCONSTRAINTS = 'getConstraints'
    GETTABLE = 'getTable'
    ADVANCE = 'advance'


class ParamNames:
    MATRIX0 = 'matrix0'
    MATRIX1 = 'matrix1'
    MATRIX2 = 'matrix2'
    INTFORCE = 'intForce'
    EXTFORCE = 'extForce'
    UNITFORCE = 'unitForce'
    CONSTRAINTS = 'constraints'
    TABLE = 'table'
    TABLENAME = 'tableName'
    TABLEWEIGHTS = 'tableWeights'
    HINGENODES = 'hingeNodes'
    HINGESTEPS = "hingeSteps"
    SOLUTION = 'solution'
