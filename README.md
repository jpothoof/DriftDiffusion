# DriftDiffusion

Note: These functions were built using the US NIST FiPy 3.4 software.
https://www.ctcms.nist.gov/fipy/

There are two main functions included, ``dd_charge`` and ``dd_relax``. Both of these call a separate function ``mesh_1d`` to construct a 1D mesh cell arangement with length (in meters) and a number of cells, cell_num. These functions currently only work in one-dimension.

``dd_charge`` is designed to be used before ``dd_relax`` as the outputted data can then be used as an input. The parameters to be worked with for these functions are steps, timestep, length, cell_num, left_bias, right_bias, Pini, Nini, mu_p_value, mu_n_value, epsilion_r, T, and k_rec.

- steps: Number of timsteps to take or how many iterations the function will calculate.
- timestep: Length in time between each step/calculation in seconds. May need to optimize to ensure stability of calculation.
- length: Length of mesh in meters.
- cell_num: Number of cells to include in mesh. This value may also need to be optimized for stability sake.
- left_bias: The potential in volts that will be applied to the left-most mesh cell. The first cell in the mesh is constrained to this value using Direchlet boundary conditions.
- right_bias: The potential in volts that will be applied to the right-most mesh cell. The first cell in the mesh is constrained to this value using Direchlet boundary conditions.
- Pini: Initial positive-species charge density. Input can be an integer, float, or array of same size as mesh, but results may vary. Explained later.
- Nini: Initial negative-species charge density. Input can be an integer, float, or array of same size as mesh, but results may vary. Explained later.
- mu_p_value: Mobility of positive-species to be defined at every cell face. Input can be an integer, float, or array of same size as mesh, but results may vary. Explained Later.
- mu_n_value: Mobility of negative-species to be defined at every cell face. Input can be an integer, float, or array of same size as mesh, but results may vary. Explained Later.
- epsilon_r: Relative permattivity of system.
- T: Temperature of system, defaulted to 298 K.
- k_rec: Recombination term, defaulted to 0.

### Units for inputs
Pini and Nini should both be inputted in units of $m^{-3}$