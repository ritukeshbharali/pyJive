# pyJive
A python Finite Element library inspired by jive

## Contributors
- Iuri Rocha (i.rocha@tudelft.nl)
- Frans van der Meer (f.p.vandermeer@tudelft.nl)
- Anne Poot
- Jelle Knibbe
- Andres Martinez Colan

## Philosophy
The library works with *modules* and *models*. Modules define the flow of the program, e.g. the `solvermodule` defines that a matrix and RHS vector need to be assembled after which a linear system of equation is solved. Models define how the main task are implemented, e.g. the `elasticmodel` that assembles the stiffness matrix for an elasticity problem.

Modules and models interact through the `take_action` function. 

## Input
Problem-specific input is defined in a properties file (with .pro extension). The repository contains sample input files. 

## Important classes
- `dofspace`: degree of freedom management to provide mapping between dof index and node number and dof type 
- `shape`: abstract base class that defines the element type 
- `constrainer`: class that handles the application of constrains to the as-assembled stiffness matrix and force vector, and returns constrained versions of them ready for solving 

## Implemented shapes
Generates shape functions, gradients and weights. The following classes are implemented in `paramshapes.py`:
- `Tri3Shape`: 3-node linear triangular shape 
- `Tri6Shape`: 6-node quadratic triangular shape
- `Quad4Shape`: 4-node linear quadrilateral shape
- `Line2Shape`: 2-node linear line shape
- `Line3Shape`: 3-node quadratic line shape

## Utility functions
- `declare.py`: this is where the available models and modules are defined (this is needed to be able to construct a problem dependent module-model scheme)
- `runPro.py`: simple master script that parses an input .pro file through `proputils.py` and calls `main.py`
- `main.py`: defines chain of modules

## Modules
Preprocessing
- `initmodule`: initializes `globdat` object with global data including mesh and nodegroups. Accepts mesh files from `gmsh` and `meshio` as well as manually generated mesh files (see syntaxis on `initmodule.py`)

Main solver routines
- `linbuckmodule`: solves eigenvalue problem for linear buckling analysis
- `solvermodule`: assembles and solves system of equations on a given number of steps
- `nonlinmodule`: newton-raphson solver for nonlinear system of equations
- `modeshapemodule`: solves eigenvalue problem for free vibration

Postprocessing
- `vtkoutmodule`: writes output to vtk file (compatible with paraview)
- `viewmodule`: visualization of full field data
- `frameviewmodule`: dedicated postprocessor for frame structures
- `loaddispmodule`: stores load displacement data for specified nodegroups in globdat
- `graphmodule`: plots data for instance data stored by loaddispmodule
- `outputmodule`: writes primary unknown to file

## Models
- `multimodel`: provides a fork from a single model into a list of models
- `elasticmodel`: assembles matrix for elasticity problem
- `barmodel`: assembles matrix for a 1D bar problem
- `timoshenkomodel`: assembles matrix for a Timoshenko beam problem
- `framemodel`: assembles matrix for frame problem with extensible Timoshenko elements
- `dirimodel`: implements Dirichlet boundary conditions
- `neumannmodel`: implements Neumann boundary conditions
