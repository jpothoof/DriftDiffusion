# DriftDiffusion

Note: These functions were built using the US NIST FiPy 3.4 software.
https://www.ctcms.nist.gov/fipy/

There are two main functions included, ``dd_charge`` and ``dd_relax``. Both of these call a separate function ``mesh_1d`` to construct a 1D mesh cell arangement with length (in meters) and a number of cells, cell_num. These functions currently only work in one-dimension.

These functions solved coupled partial differential equations that describe drift and diffusion of charged particles. Initially, it was designed to look at a system of positively and negatively charged ions, but can be used for electronic carriers as well depending on inputted charge densities and mobilities.

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

As the data for each function is saved as a text file, the final timestep from ``dd_charge`` can be imported and used to define the initial charge densities, Pini and Nini, for ``dd_relax``. Additionally, ``dd_relax`` defaults left_bias and right_bias to 0 to simulate the electric bias being removed and the system grounded.

### Units for inputs
length should be inputted in units of m
left_bias and right_bias should be inputted in units of V
Pini and Nini should both be inputted in units of m<sup>-3</sup>
mu_p_value and mu_n_value should both be inputted in units of m<sup>2</sup>/V s
epsilon_r is unitless
T should be inputted in units of K

### Varying Inputs
As previously stated, Pini, Nini, mu_p_value, and mu_n_value can be inputted with various data types: integer, float, and array. If inputted as single variables (integer/float), ``dd_charge`` automattically applies an insulated layer to the first and last 2.5% of the mesh. So, the positive- and negative-species charge densities and positive and negative mobilities are set to zero in those cells. The internal cells and faces are set to a uniform/flat distribution of charged-species and mobilites, respectively.

Therefore, the user should be careful about using mixed data type inputs (i.e. single variables for the charge-densities and arrays for the mobilities). The function will automatically shut down if mixed data types are used for Pini and Nini or for mu_p_value and mu_n_value.

It is recommended for the user to either use single-variable inputs for all parameters, or array inputs for all parameters. Otherwise, be mindful of the 2.5% rule on either side of the mesh and construct the array accordingly.

### Function Outputs
As either function runs, a plot will be displayed showing the evolution of the potential, positive- and negative-species charge densities, and the electric field over time. After completion, ``dd_charge`` will save four txt files to the user's current directory with data structure **cell_num X steps** for "Positive Ion Density_Charge.txt", "Negative Ion Density_Charge.txt", "Potential_Charge.txt", and data structure **cell_num + 1 X steps** for "Electric Field_Charge.txt". ``dd_relax`` does the same, but with the filenames containing "_Relax" instead of "_Charge".
