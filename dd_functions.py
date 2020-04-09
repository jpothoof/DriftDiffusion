import numpy as np

import matplotlib

matplotlib.use('TkAgg', force=True)
from matplotlib import pyplot as plt

from fipy import Variable, FaceVariable, CellVariable, Grid1D, ExplicitDiffusionTerm, TransientTerm, DiffusionTerm, \
    Viewer, ImplicitSourceTerm, ConvectionTerm


def mesh_1d(length, cell_num):
    '''
    This function creates a 1D mesh from an inputted length and number of cells. This mesh is used in other FiPy
    functions to establish and calculate different variables.
    Inputs:
        length - Length of system in meters
        cell_num - Number of cells in mesh
    Outputs:
        mesh
    '''
    dx = length / cell_num  # Width of each mesh cell
    mesh = Grid1D(dx=dx, nx=cell_num)
    return mesh


def dd_charge(steps, timestep, length, cell_num, left_bias, right_bias, Pini, Nini, mu_p_value, mu_n_value, epsilon_r,
              T=298, k_rec=0):
    '''
    This functions solves drift-diffusion equations for a lateral bi-electrode system. This function is  designed for
    applying an external bias to the system to drive charge carriers to either electrode. From a set of inputs, it will
    display a figure showing the evolution of positive and negative charge carriers over time and changes in the
    electric field and potential curve. It is recommended to use the same data type for Pini, Nini, mu_p_value, and
    mu_n_value (i.e. float for all values, or array for all values). WARNING: Does not currently function in Jupyter
    Notebooks. Tested and works in PyCharm and Spyder.
    Inputs:
        steps - Number of iterations to be performed
        timestep - Time in seconds between each iteration. If this value is too large, the simulation may crash.
        length - Length of system in meters
        cell_num - Number of cells in mesh
        left_bias - Electrical bias in Volts to be applied to the left electrode.
        right_bias - Electrical bias in Volts to be applied to the right electrode. right_bias=0 is grounded.
        Pini - Single value input = Uniform initial positive charge density distribution, defaults to 2.5% of mesh on
        either side being insulated. An array of the same shape as the mesh will allow a non-uniform initial
        distribution.
        Nini - Single value input = Uniform initial positive charge density distribution, defaults to 2.5% of mesh on
        either side being insulated. An array of the same shape as the mesh will allow a non-uniform initial
        distribution.
        mu_p_value - Single value input = Uniform positive charge mobility at every face, defaults to 2.5% of mesh on
        either side being insulated. An array of the shape cell_num + 1 will allow varying mobilities at each cell face.
        mu_n_value - Single value input = Uniform positive charge mobility at every face, defaults to 2.5% of mesh on
        either side being insulated. An array of the shape cell_num + 1 will allow varying mobilities at each cell face.
        epsilon_r - Relative permittivity of system.
        T - Temperature is defaulted to 298 K (Room Temperature).
        k_rec - Recombination Term. Defaulted to be 'off'. Will implement later.
    '''

    q = 1.602e-19  # C
    kb = 1.38e-23  # J/K
    epsilon = epsilon_r * 8.854e-12  # Permittivity C/V*m
    f = kb * T / q  # Einstein Relation V

    mesh = mesh_1d(length, cell_num)  # Creates 1D mesh of described length and number of mesh cells
    mesh_insulated = np.int(np.round(cell_num*0.025, decimals=0))
    face_insulated = mesh_insulated + 1
    x = np.linspace(0, length, cell_num)  # Numpy array of same length and mesh for plotting
    x2 = np.linspace(0, length, cell_num+1)

    if isinstance(Pini, (int, float)) == True and isinstance(Nini, (int, float)) == True:
        p_initial = np.zeros(int(cell_num))
        n_initial = np.zeros(int(cell_num))
        p_initial[mesh_insulated:-mesh_insulated] = Pini
        n_initial[mesh_insulated:-mesh_insulated] = Nini

        if isinstance(mu_p_value, (int, float)) == True and isinstance(mu_n_value, (int,float)) == True:
            mu_n_initial = np.zeros(int(cell_num+1))
            mu_p_initial = np.zeros(int(cell_num+1))
            mu_n_initial[face_insulated:-face_insulated] = mu_n_value
            mu_p_initial[face_insulated :-face_insulated] = mu_p_value

        elif isinstance(mu_p_value, np.ndarray) == True and isinstance(mu_n_value, np.ndarray) == True:
            mu_n_initial = mu_n_value
            mu_p_initial = mu_p_value

        elif (isinstance(mu_p_value, (int, float)) == True and isinstance(mu_n_value, (int, float)) == False) or \
                (isinstance(mu_p_value, (int, float)) == False and isinstance(mu_n_value, (int, float)) == True):
            print("Error: Pini and Nini must both be an integer/float, or both must be an array value")
            raise SystemExit(0)

        else:
            pass

    elif (isinstance(Pini, (int, float)) == True and isinstance(Nini, (int, float)) == False) or \
            (isinstance(Pini, (int, float)) == False and isinstance(Nini, (int, float)) == True):
        print("Error: Pini and Nini must both be an integer/float, or both must be an array value")
        raise SystemExit(0)

    elif isinstance(Pini, np.ndarray) == True and isinstance(Nini, np.ndarray) == True:
        p_initial = Pini
        n_initial = Nini

        if isinstance(mu_p_value, (int, float)) == True and isinstance(mu_n_value, (int,float)) == True:
            mu_n_initial = np.zeros(int(cell_num+1))
            mu_p_initial = np.zeros(int(cell_num+1))
            mu_n_initial[face_insulated:-face_insulated] = mu_n_value
            mu_p_initial[face_insulated :-face_insulated] = mu_p_value

        elif isinstance(mu_p_value, np.ndarray) == True and isinstance(mu_n_value, np.ndarray) == True:
            mu_n_initial = mu_n_value
            mu_p_initial = mu_p_value

        elif (isinstance(mu_p_value, (int, float)) == True and isinstance(mu_n_value, (int, float)) == False) or \
                (isinstance(mu_p_value, (int, float)) == False and isinstance(mu_n_value, (int, float)) == True):
            print("Error: Pini and Nini must both be an integer/float, or both must be an array value")
            raise SystemExit(0)

        else:
            pass

    # VARIABLE SETUP:
    # We establish each of our variable as FiPy CellVariable for ones establish at mesh cell centers: positive ions
    # negative ions, and potential. Variables that describe the flux of particles are described as FaceVariables, which
    # establishes them at mesh cell interfaces. We can also set initial values here for positive and negative ion
    # densities.
    Pion = CellVariable(mesh=mesh, name='Positive ion Charge Density', value=p_initial)
    Nion = CellVariable(mesh=mesh, name='Negative ion Charge Density', value=n_initial)
    mu_n = FaceVariable(mesh=mesh, name='Negative ion mobility', value=mu_n_initial)
    mu_p = FaceVariable(mesh=mesh, name='Positive ion mobility', value=mu_p_initial)
    Dn = FaceVariable(mesh=mesh, name='Negative ion diffusivity', value=f * mu_n)
    Dp = FaceVariable(mesh=mesh, name='Positive ion diffusivity', value=f * mu_p)
    potential = CellVariable(mesh=mesh, name='Potential')

    # EQUATION SETUP:
    # The differential equations are defined by different terms. Pertinent for Drift-Diffusion equations are
    # TransientTerm = dVar/dt
    # ConvectionTerm = dVar/dx
    # DiffusionTerm = d^2Var/dx^2

    # In English:  dPion/dt = -1/q * divergence.Jp(x,t) - k_rec * Nion(x,t) * Pion(x,t)     where
    #             Jp = q * mu_p * E(x,t) * Pion(x,t) - q * Dp * grad.Pion(x,t)   and   E(x,t) = -grad.potential(x,t)

    Pion_equation = TransientTerm(coeff=1., var=Pion) == ConvectionTerm(coeff=mu_p * potential.faceGrad, var=Pion) +  \
                    DiffusionTerm(coeff=Dp, var=Pion) - k_rec * Pion * Nion

    # In English:  dNion/dt = 1/q * divergence.Jn(x,t) - k_rec * Nion(x,t) * Pion(x,t)   where
    #             Jn = q * mu_n * E(x,t) * Nion(x,t) - q * Dn * grad.Nion(x,t)   and   E(x,t) = -grad.potential(x,t)

    Nion_equation = TransientTerm(coeff=1., var=Nion) == ConvectionTerm(coeff=-mu_n * potential.faceGrad, var=Nion) +  \
                    DiffusionTerm(coeff=Dn, var=Nion) - k_rec * Pion * Nion

    # In English:  d^2potential/dx^2 = -q/epsilon * Charge_Density   and   Charge Density= Pion - Nion
    # Poisson's Equation

    potential_equation = DiffusionTerm(coeff=1., var=potential) == -(q / epsilon) * (Pion - Nion)

    # BOUNDARY CONDITIONS
    # Boundary conditions can be applied only to CellVariables. Direchlet Boundaries can easily be applied by using
    # the CellVariable.constrain method with a value and a definition of where. mesh.facesLeft or Right apply the
    # boundary condition to the furthest left and furthest right cell. One can also use a Boolean array where True
    # values define where to apply the boundary conditon.
    potential.constrain(left_bias, where=mesh.facesLeft)
    potential.constrain(right_bias, where=mesh.facesRight)

    # Coupling equations to be solved together
    eq = Pion_equation & Nion_equation & potential_equation

    # PLOT FIGURE SET UP
    # Potential and Electric Field are plotted together in the first subplot. Positive and Negative charge densities
    # are plotted together in the second subplot.
    plt.ion()
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    plt.subplots_adjust(wspace=0.3)
    ax1.axhline(y=0)
    ax1.set_ylabel('Potential (V)', c='tab:green')
    ax1.tick_params(axis='y', labelcolor='tab:green')
    ax1.set_ylim([np.min([left_bias, right_bias]), np.max([left_bias, right_bias])+1])
    ax1.set_xticks(np.linspace(0, length, 8))
    ax1.ticklabel_format(axis='x', style='sci', scilimits=(0, 1))
    ax1.set_xlabel('Distance (m)')
    ax2.set_ylabel('Charge Density (cm^-3)')
    if isinstance(Pini, (int, float)) == True:
        ax2.set_ylim([0, Pini * 20])
    else:
        ax2.set_ylim([0, max(Pini) * 20])
    ax2.set_xticks(np.linspace(0, length, 8))
    ax2.ticklabel_format(axis='x', style='sci', scilimits=(0, 1))
    ax2.set_xlabel('Distance (m)')
    ax3.set_ylabel('Electric Field (V/m)', c='grey')
    ax3.tick_params(axis='y', labelcolor='grey')
    ax3.ticklabel_format(axis='y', style='sci', scilimits=(0, 1))
    ax3.set_xticks(np.linspace(0, length, 8))
    ax3.ticklabel_format(axis='x', style='sci', scilimits=(0, 1))
    ax3.set_xlabel('Distance (m)')

    # Establish arrays to save the data to.
    Efield_save = np.empty([cell_num+1, steps])
    potential_save = np.empty([cell_num, steps])
    Pion_save = np.empty_like(potential_save)
    Nion_save = np.empty_like(potential_save)

    # Iterate through a number of steps to solve equations with a set timestep
    for i in range(steps):
        eq.solve(dt=timestep)  # Solves the coupled differential equations
        # print(potential.grad())
        Efield = np.array(-potential.faceGrad[0])

        l4 = ax3.plot(x2, Efield, c='grey')
        l1 = ax1.plot(x, potential.value, c='tab:green')
        l2 = ax2.plot(x, Pion.value, label='Positive Density', c='tab:blue')
        l3 = ax2.plot(x, Nion.value, label='Negative Density', c='tab:red')

        #Legend only needs to be set once
        if i == 0:
            ax2.legend(loc="upper right")
            ax3.set_ylim([0, max(Efield) * 10])

        else:
            pass

        # Fill data into array for each iteration
        Efield_save[:, i] = Efield
        potential_save[:, i] = potential.value
        Pion_save[:, i] = Pion.value
        Nion_save[:, i] = Nion.value

        fig.suptitle('Time: ' + str(i * timestep)[:4])
        plt.pause(0.00000001)
        plt.draw()

        # Need to remove the data from the lines every step, otherwise it plots them all and looks messy.
        for l in [l1, l2, l3, l4]:
            l[0].remove()

    # Saves data for each desired variable to a txt file to your current directory.
    np.savetxt('Positive Ion Density_Charge.txt', Pion_save, delimiter='\t')
    np.savetxt('Negative Ion Density_Charge.txt', Nion_save, delimiter='\t')
    np.savetxt('Electric Field_Charge.txt',Efield_save, delimiter='\t')
    np.savetxt('Potential_Charge.txt', potential_save, delimiter='\t')

def dd_relax(steps, timestep, length, cell_num, Pini, Nini, mu_p_value, mu_n_value,
              epsilon_r, left_bias=0, right_bias=0, T=298, k_rec=0.):
    '''
    This functions solves drift-diffusion equations for a lateral bi-electrode system. This function is designed to
    simulate how the charge carriers relax after the charging process. Either use the last timestep of data for positive
    ions and negative ions as Pini and Nini (from the saved txt files), or experimental data for the initial values.
    From a set of inputs, it will display a figure showing the evolution of positive and negative charge carriers over
    time and changes in the electric field and potential curve. It is recommended to use the same data type for Pini,
    Nini, mu_p_value, and mu_n_value (i.e. float for all values, or array for all values). WARNING: Does not currently
    function in Jupyter Notebooks. Tested and works in PyCharm and Spyder.
    Inputs:
        steps - Number of iterations to be performed
        timestep - Time in seconds between each iteration. If this value is too large, the simulation may crash.
        length - Length of system in meters
        cell_num - Number of cells in mesh
        Pini - Single value input = Uniform initial positive charge density distribution, defaults to 2.5% of mesh on
        either side being insulated. An array of the same shape as the mesh will allow a non-uniform initial
        distribution.
        Nini - Single value input = Uniform initial positive charge density distribution, defaults to 2.5% of mesh on
        either side being insulated. An array of the same shape as the mesh will allow a non-uniform initial
        distribution.
        mu_p_value - Single value input = Uniform positive charge mobility at every face, defaults to 2.5% of mesh on
        either side being insulated. An array of the shape cell_num + 1 will allow varying mobilities at each cell face.
        mu_n_value - Single value input = Uniform positive charge mobility at every face, defaults to 2.5% of mesh on
        either side being insulated. An array of the shape cell_num + 1 will allow varying mobilities at each cell face.
        left_bias - Electrical bias in Volts to be applied to the left electrode. Defaulted to 0 for grounding.
        right_bias - Electrical bias in Volts to be applied to the right electrode. Defaulted to 0 for grounding
        T - Temperature is defaulted to 298 K (Room Temperature).
        k_rec - Recombination Term. Defaulted to be 'off'. Will implement later.
    '''

    q = 1.602e-19  # C
    kb = 1.38e-23  # J/K
    epsilon = epsilon_r * 8.854e-12  # Permittivity C/V*m
    f = kb * T / q  # Einstein Relation V

    mesh = mesh_1d(length, cell_num)  # Creates 1D mesh of described length and number of mesh cells
    mesh_insulated = np.int(np.round(cell_num*0.025, decimals=0))
    face_insulated = mesh_insulated + 1
    x = np.linspace(0, length, cell_num)  # Numpy array of same length and mesh for plotting
    x2 = np.linspace(0, length, cell_num+1)

    if isinstance(Pini, (int, float)) == True and isinstance(Nini, (int, float)) == True:
        p_initial = np.zeros(int(cell_num))
        n_initial = np.zeros(int(cell_num))
        p_initial[mesh_insulated:-mesh_insulated] = Pini
        n_initial[mesh_insulated:-mesh_insulated] = Nini

        if isinstance(mu_p_value, (int, float)) == True and isinstance(mu_n_value, (int,float)) == True:
            mu_n_initial = np.zeros(int(cell_num+1))
            mu_p_initial = np.zeros(int(cell_num+1))
            mu_n_initial[face_insulated:-face_insulated] = mu_n_value
            mu_p_initial[face_insulated :-face_insulated] = mu_p_value

        elif isinstance(mu_p_value, np.ndarray) == True and isinstance(mu_n_value, np.ndarray) == True:
            mu_n_initial = mu_n_value
            mu_p_initial = mu_p_value

        elif (isinstance(mu_p_value, (int, float)) == True and isinstance(mu_n_value, (int, float)) == False) or \
                (isinstance(mu_p_value, (int, float)) == False and isinstance(mu_n_value, (int, float)) == True):
            print("Error: Pini and Nini must both be an integer/float, or both must be an array value")
            raise SystemExit(0)

        else:
            pass

    elif (isinstance(Pini, (int, float)) == True and isinstance(Nini, (int, float)) == False) or \
            (isinstance(Pini, (int, float)) == False and isinstance(Nini, (int, float)) == True):
        print("Error: Pini and Nini must both be an integer/float, or both must be an array value")
        raise SystemExit(0)

    elif isinstance(Pini, np.ndarray) == True and isinstance(Nini, np.ndarray) == True:
        p_initial = Pini
        n_initial = Nini

        if isinstance(mu_p_value, (int, float)) == True and isinstance(mu_n_value, (int,float)) == True:
            mu_n_initial = np.zeros(int(cell_num+1))
            mu_p_initial = np.zeros(int(cell_num+1))
            mu_n_initial[face_insulated:-face_insulated] = mu_n_value
            mu_p_initial[face_insulated:-face_insulated] = mu_p_value

        elif isinstance(mu_p_value, np.ndarray) == True and isinstance(mu_n_value, np.ndarray) == True:
            mu_n_initial = mu_n_value
            mu_p_initial = mu_p_value

        elif (isinstance(mu_p_value, (int, float)) == True and isinstance(mu_n_value, (int, float)) == False) or \
                (isinstance(mu_p_value, (int, float)) == False and isinstance(mu_n_value, (int, float)) == True):
            print("Error: Pini and Nini must both be an integer/float, or both must be an array value")
            raise SystemExit(0)

        else:
            pass

    # VARIABLE SETUP:
    # We establish each of our variable as FiPy CellVariable for ones establish at mesh cell centers: positive ions
    # negative ions, and potential. Variables that describe the flux of particles are described as FaceVariables, which
    # establishes them at mesh cell interfaces. We can also set initial values here for positive and negative ion
    # densities.
    Pion = CellVariable(mesh=mesh, name='Positive ion Charge Density', value=p_initial)
    Nion = CellVariable(mesh=mesh, name='Negative ion Charge Density', value=n_initial)
    mu_n = FaceVariable(mesh=mesh, name='Negative ion mobility', value=mu_n_initial)
    mu_p = FaceVariable(mesh=mesh, name='Positive ion mobility', value=mu_p_initial)
    Dn = FaceVariable(mesh=mesh, name='Negative ion diffusivity', value=f * mu_n)
    Dp = FaceVariable(mesh=mesh, name='Positive ion diffusivity', value=f * mu_p)
    potential = CellVariable(mesh=mesh, name='Potential')

    # EQUATION SETUP:
    # The differential equations are defined by different terms. Pertinent for Drift-Diffusion equations are
    # TransientTerm = dVar/dt
    # ConvectionTerm = dVar/dx
    # DiffusionTerm = d^2Var/dx^2

    # In English:  dPion/dt = -1/q * divergence.Jp(x,t) - k_rec * Nion(x,t) * Pion(x,t)     where
    #             Jp = q * mu_p * E(x,t) * Pion(x,t) - q * Dp * grad.Pion(x,t)   and   E(x,t) = -grad.potential(x,t)

    Pion_equation = TransientTerm(coeff=1., var=Pion) == ConvectionTerm(coeff=mu_p * potential.faceGrad, var=Pion) + \
                    DiffusionTerm(coeff=Dp, var=Pion) - k_rec * Pion * Nion

    # In English:  dNion/dt = 1/q * divergence.Jn(x,t) - k_rec * Nion(x,t) * Pion(x,t)   where
    #             Jn = q * mu_n * E(x,t) * Nion(x,t) - q * Dn * grad.Nion(x,t)   and   E(x,t) = -grad.potential(x,t)

    Nion_equation = TransientTerm(coeff=1., var=Nion) == ConvectionTerm(coeff=-mu_n * potential.faceGrad, var=Nion) + \
                    DiffusionTerm(coeff=Dn, var=Nion) - k_rec * Pion * Nion

    # In English:  d^2potential/dx^2 = -q/epsilon * Charge_Density   and   Charge Density= Pion - Nion
    # Poisson's Equation

    potential_equation = DiffusionTerm(coeff=1., var=potential) == -(q / epsilon) * (Pion - Nion)

    # BOUNDARY CONDITIONS
    # Boundary conditions can be applied only to CellVariables. Direchlet Boundaries can easily be applied by using
    # the CellVariable.constrain method with a value and a definition of where. mesh.facesLeft or Right apply the
    # boundary condition to the furthest left and furthest right cell. One can also use a Boolean array where True
    # values define where to apply the boundary conditon.
    potential.constrain(left_bias, where=mesh.facesLeft)
    potential.constrain(right_bias, where=mesh.facesRight)

    # Coupling equations to be solved together
    eq = Pion_equation & Nion_equation & potential_equation

    # PLOT FIGURE SET UP
    # Potential and Electric Field are plotted together in the first subplot. Positive and Negative charge densities
    # are plotted together in the second subplot.
    plt.ion()
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    plt.subplots_adjust(wspace=0.3)

    ax1.axhline(y=0)
    ax1.set_ylabel('Potential (V)', c='tab:green')
    ax1.tick_params(axis='y', labelcolor='tab:green')
    ax1.set_xticks(np.linspace(0, length, 8))
    ax1.ticklabel_format(axis='x', style='sci', scilimits=(0, 1))
    ax1.set_xlabel('Distance (m)')
    ax2.set_ylabel('Charge Density (cm^-3)')
    if isinstance(Pini, (int, float)) == True:
        ax2.set_ylim([0, Pini * 20])
    else:
        ax2.set_ylim([0, max(Pini)])
    ax2.set_xticks(np.linspace(0, length, 8))
    ax2.ticklabel_format(axis='x', style='sci', scilimits=(0, 1))
    ax2.set_xlabel('Distance (m)')
    ax3.set_ylabel('Electric Field (V/m)', c='grey')
    ax3.tick_params(axis='y', labelcolor='grey')
    ax3.ticklabel_format(axis='y', style='sci', scilimits=(0, 1))
    ax3.set_xticks(np.linspace(0, length, 8))
    ax3.ticklabel_format(axis='x', style='sci', scilimits=(0, 1))
    ax3.set_xlabel('Distance (m)')

    # Establish arrays to save the data to.
    Efield_save = np.empty([cell_num+1, steps])
    potential_save = np.empty([cell_num, steps])
    Pion_save = np.empty_like(potential_save)
    Nion_save = np.empty_like(potential_save)

    # Iterate through a number of steps to solve equations with a set timestep
    for i in range(steps):
        eq.solve(dt=timestep)  # Solves the coupled differential equations
        Efield = np.array(-potential.faceGrad[0])
        print(Nion.value[0], sum(Nion.value))

        l4 = ax3.plot(x2, Efield, c='grey')
        l1 = ax1.plot(x, potential.value, c='tab:green')
        l2 = ax2.plot(x, Pion.value, label='Positive Density', c='tab:blue')
        l3 = ax2.plot(x, Nion.value, label='Negative Density', c='tab:red')

        # Legend only needs to be set once
        if i == 0:
            ax3.set_ylim([-max(Efield) * 0.5, max(Efield) * 1.5])
            ax2.legend(loc="upper right")

        else:
            pass
        ax1.set_ylim([min(np.array(potential.value)), max(np.array(potential.value))])

        # Fill data into array for each iteration
        Efield_save[:, i] = Efield
        potential_save[:, i] = potential.value
        Pion_save[:, i] = Pion.value
        Nion_save[:, i] = Nion.value

        fig.suptitle('Time: ' + str(i * timestep)[:4])
        plt.pause(0.00000001)
        plt.draw()

        # Need to remove the data from the lines every step, otherwise it plots them all and looks messy.
        for l in [l1, l2, l3, l4]:
            l[0].remove()

    # Saves data for each desired variable to a txt file to your current directory.
    np.savetxt('Positive Ion Density_Relax.txt', Pion_save, delimiter='\t')
    np.savetxt('Negative Ion Density_Relax.txt', Nion_save, delimiter='\t')
    np.savetxt('Electric Field_Relax.txt', Efield_save, delimiter='\t')
    np.savetxt('Potential_Relax.txt', potential_save, delimiter='\t')
