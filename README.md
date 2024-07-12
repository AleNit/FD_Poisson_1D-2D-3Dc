# FD_Poisson_1D-2D-3Dc
Contains several scripts for solving the Poisson equation by a finite difference approximation in 1D, 2D or 3D in a cylindrical coordinate system over non-uniform structured grids.

Derivatives are approximated by a 2-nd order accurate centered difference over a cell-centered grid. Grid stretching can be prescribed. In case of uniform grid arrangement, a one directional DFT (Discrete Fourier Transform) is applied to accelerate the problem solution.
Dirichelet/Neumann/periodic boundary conditions are prescribed by a ghost cell method. The generic boundary condition is formulated as a Robin condition, see personal notes (.pdf file) for details.

The discretization arrangement and the DFT application rely on the research papers:
- Farnell, L. "Solution of Poisson equations on a nonuniform grid." Journal of Computational Physics 35.3 (1980): 408-425.
