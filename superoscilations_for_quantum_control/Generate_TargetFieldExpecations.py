#####################################################################
# This file takes an analytical eq. for a driving laser pulse and   #
# and uses the Hubbard Model to generate expectations based on the  #
# parameters you provide.                                           #
# If you are trying to generate expectations from external field    #
# or need to double check tracking results, please use              #
# Generate_ImportedFieldExpectations.py instead                     #
# Based on the Quspin package.                                      #
# For details, see:  http://weinbe58.github.io/QuSpin/index.html    #
#####################################################################
from __future__ import print_function, division
import os
import sys

"""Open MP and MKL should speed up the time required to run these simulations!"""
# threads = sys.argv[1]
threads = 16
os.environ['NUMEXPR_MAX_THREADS'] = '{}'.format(threads)
os.environ['NUMEXPR_NUM_THREADS'] = '{}'.format(threads)
os.environ['OMP_NUM_THREADS'] = '{}'.format(threads)
os.environ['MKL_NUM_THREADS'] = '{}'.format(threads)
# line 4 and line 5 below are for development purposes and can be remove
from quspin.operators import hamiltonian, exp_op, quantum_operator  # operators
from quspin.basis import spinful_fermion_basis_1d  # Hilbert space basis
from quspin.tools.measurements import obs_vs_time  # calculating dynamics
from quspin.tools.evolution import evolve  # evolving system
import numpy as np  # general math functions
from time import time  # tool for calculating computation time
from tqdm import tqdm
import matplotlib.pyplot as plt  # plotting library
from scipy.signal.windows import blackman
from scipy import fftpack
from quspin.tools.measurements import project_op
from tools import HubbardModel as fhmodel, InitializeArchive
import psutil



# note cpu_count for logical=False returns the wrong number for multi-socket CPUs.
print("logical cores available {}".format(psutil.cpu_count(logical=True)))
t_init = time()
np.__config__.show()


########################################################################################################################
# Declare Parameters
########################################################################################################################

"""Laser Pulse parameters"""
field = 32.9  # field angular frequency THz
F0 = 4.5  # Field amplitude MV/cm
a = 4   # Lattice constant Angstroms

"""Parameters for a target or reference field"""
# Hubbard model
L = 10             # system size
N_up = L // 2 + L % 2   # number of fermions with spin up
N_down = L // 2     # number of fermions with spin down
N = N_up + N_down   # number of particles
t0 = 0.52       # hopping strength
U = 0.0 * t0    # interaction strength
pbc = True

# Parameters for evolving the system
cycles = 10     # time in cycles of field frequency
n_steps = 10000  # Number of steps for time resolution

# Bundle parameters to pass to Hubbard Model class for unit conversion
tgt_params = dict(
    nx=L,
    hopping=t0,
    interaction=U,
    n_up=N_up,
    n_down=N_down,
    angular_frequency=field,
    lattice_constant=a,
    field_amplitude=F0,
    chem_potential=0,
    cycles=cycles,
    n_steps=n_steps,
    ny=0,       # 1D simulations do not use y-axis
    soc=0,      # No spin orbit coupling
    gamma=0     #
)

# get the converted units for creating a target field
tgt = fhmodel(**tgt_params)
start = 0.0
stop = tgt.stop
times, dt = np.linspace(start, stop, num=n_steps, endpoint=True, retstep=True)

########################################################################################################################
# Get the directory and set tag for data saving
########################################################################################################################

# Is this your first time running any program? If so, make sure your directories are set up but setting to true
make_archive = False
# Reference: 0 for Targets, 1 for Tracking, 2 for superoscillations, 3 for importing anything else
sim_type = 0

# Build out the directory
directory = InitializeArchive(directory_number=sim_type)
if make_archive:
    directory.build_archive()

lib = directory.get_dir()
data_path = lib['data path']
plot_path = lib['plots path']
plot_name = lib['plot name']

# Get the unique tag for these simulations. Done externally in order to ensure consistency
filetag = lib['tag'] + tgt.tag
########################################################################################################################
# Define the relevant functions
########################################################################################################################


def phi(current_time):
    """
    Defines Phi(t) for an individual time step
    :param current_time: Current time
    :return: Phi(current time)
    """
    return (tgt.a * tgt.F0 / tgt.omega) * (np.sin(tgt.omega * current_time / (2. * cycles)) ** 2.) * np.sin(
        tgt.omega * current_time)


def expiphi(current_time):
    """
    Exponentiation of Phi(t): e^(i*phi(t))
    :param current_time: time
    :return: e^(1j*phi(current time))
    """
    return -tgt.t0 * np.exp(1j * phi(current_time))


def expiphiconj(current_time):
    """
    Complex conjugate of exponentiation of Phi(t): e^-(i*phi)
    :param current_time: time
    :return: e^(-1j*phi(current_time)
    """
    return -tgt.t0 * np.exp(-1j * phi(current_time))


########################################################################################################################
# Define the Hamiltonian
# H = -t0 sum(e^(-i phi(t)) c^*_{i,s} c_{i+1,s} - e^(i phi(t)) c_{i,s} c^*_{i+1,s}) + U sum(n_{i,up} n_i,down})
# or H = sum(-t0 * expiphiconj * hop_right + t0 * expiphi * hop_left) + sum(U * interactions)
########################################################################################################################
# Specify the basis
basis = spinful_fermion_basis_1d(L=L, Nf=(N_up, N_down))
no_checks = dict(check_pcon=False, check_symm=False, check_herm=False)

# Define onsite coupling lists, U * sum(n_{i, up}, n_{i, down})
interactions = [[tgt.U, _, _] for _ in range(L)]
# For hopping, the previous code sometimes uses -t0 for hop_left, and sometimes -t0 for hop right instead
# To be clear, the following code is correct.
hop_left = [[1, _, (_ + 1)] for _ in range(L - 1)]         # Allow OBC, but then add PBC
hop_right = [[-1, _, (_ + 1)] for _ in range(L - 1)]       # Allow OBC, but then add PBC

if pbc:
    hop_left.append([1, L-1, 0])
    hop_right.append([-1, L-1, 0])

# Build the Hamiltonian
static = [
    ['n|n', interactions]   # up-down interactions
]

# Because the hop-left and -right operators are time dependent, they are placed here
"""
dynamic_arguments is meant to be a list of arguments passed to the function in the dynamic Hamiltonian
However, we do not have a list of times to pass to the function, so right now it is blank
This may be a source of problems.
"""
dynamic_arguments = []
dynamic = [
    ["+-|", hop_left, expiphiconj, dynamic_arguments],          # Spin-up hop left
    ["|+-", hop_left, expiphiconj, dynamic_arguments],          # Spin-down hop left
    ["-+|", hop_right, expiphi, dynamic_arguments],     # Spin-up hop right
    ["|-+", hop_right, expiphi, dynamic_arguments]      # Spin-down hop right
]

# Actually build the Hamiltonian
# For finding ground-state using Imaginary Time
gs_ham = hamiltonian(static + [_[:2] for _ in dynamic], [], basis=basis, **no_checks)
# Hamiltonian for propagation
ham = hamiltonian(static, dynamic, basis=basis, **no_checks)

########################################################################################################################
# Define the operator dictionary. This can be used to get expectations, such as current
########################################################################################################################

op_dict = dict(H=ham)

""" We are doing this to eventually get expectation of the current operator
^J(t) = -i a t0 sum(e^(-i phi(t) c+_{j,s} c_{j+1,s} - H.c.)
or J(t) = -1j * a * t0 * sum(expiphiconj * hop_right - Hermitian conjugate)
"""

# Neighbor interactions, terms are listed in static hamiltonian as there is no time dependence
op_dict['neighbor'] = hamiltonian([["+-|", hop_left], ["|+-", hop_left]], [], basis=basis, **no_checks)

# Hop Right terms, Empty lists for static Hamiltonian section.
# Again, an empty list is used for the dynamic arguments. This could be an error source
op_dict['left_hop_up'] = hamiltonian([], [["-+|", hop_left, expiphiconj, dynamic_arguments]],
                                      basis=basis, **no_checks)
op_dict['left_hop_down'] = hamiltonian([], [["|-+", hop_left, expiphiconj, dynamic_arguments]],
                                        basis=basis, **no_checks)

# Next, we need to obtain spin densities per site
# This is an inefficient method, but at the current time, it works, so it is low priority on the changelist
for _ in range(L):
    # spin up densities for each site
    op_dict['n_up' + str(_)] = hamiltonian([["n|", [[1.0, _]]]], [], basis=basis, **no_checks)
    # spin down
    op_dict['n_down' + str(_)] = hamiltonian([["|n", [[1.0, _]]]], [], basis=basis, **no_checks)
    # Doublon densities
    op_dict["doublons" + str(_)] = hamiltonian([["n|n", [[1.0, _, _]]]], [], basis=basis, **no_checks)

# Determine Hopping Operators for Later calculation
hop_left_operator = hamiltonian([["+-|", hop_left], ["|+-", hop_left]], [], basis=basis, **no_checks)
hop_right_operator = hamiltonian([["-+|", hop_right], ["|-+", hop_right]], [], basis=basis, **no_checks)

########################################################################################################################
# Get the Ground State, Evolve the Hamiltonian, And Calculate Expectations
########################################################################################################################
# Get the ground state
it_groundstate = True
if it_groundstate:
    E, psi_0 = tgt.get_groundstate(H=gs_ham, imaginary_time=True, n_steps=10)
else:
    E, psi_0 = tgt.get_groundstate(H=ham, imaginary_time=False, n_steps=0)

# Evolve the Hamiltonian
print('Now evolving system')
t_evolve_start = time()
# psi = ham.evolve(v0=psi_0, t0=start, times=times, iterate=True, verbose=False)
# Obs_vs_Time is designed to work with iterate=True, but the generator object is unrealistically difficult to work with
# Until this is fixed, we can simply squeeze the result and directly return the array
psi = np.squeeze(ham.evolve(v0=psi_0, t0=start, times=times, verbose=False))
t_evolve_mid = time()
print('System evolution complete. It took {:.3f} seconds.'.format(t_evolve_mid-t_evolve_start))

print('Now Calculating Expectations')
expectations = obs_vs_time(psi_t=psi, times=times, Obs_dict=op_dict)


current_sum_term = (expectations['left_hop_up'] + expectations['left_hop_down'])
current = 1j * tgt.a * (current_sum_term - current_sum_term.conjugate())

# Store the field, current, times, and dt
results = directory.bundle_results(expectations, phi(times), current, times, dt, psi_0, E)


print('Checking for Solution Uniqueness')
x_den = 2 * tgt.a * tgt.t0 * np.abs(hop_left_operator.expt_value(psi))
x_t = current / x_den
if all(_ < 0.9999999 for _ in np.abs(x_t)):
    print('Solution is unique. Max value of X(t) = {:.7f}'.format(np.abs(x_t).max()))
else:
    print('Solution is not guaranteed to be unique. Max value of X(t) = {:.7f}'.format(np.abs(x_t).max()))

print('Expectations Calculated. This took an additional {:.3f} seconds.'.format(time()-t_evolve_mid))

########################################################################################################################
# Check linearity of results
########################################################################################################################


def plot_spectrum(f, t):
    """
    Plot the High Harmonic Generation spectrum
    """
    # Power spectrum emitted is calculated using the Larmor formula
    #   (https://en.wikipedia.org/wiki/Larmor_formula)
    # which says that the power emitted is proportional to the square of the acceleration
    # i.e., the RHS of the second Ehrenfest theorem

    N = len(f)
    k = np.arange(N)

    # frequency range
    omegas = (k - N / 2) * np.pi / (0.5 * t.max())

    # spectra of the
    spectrum = np.abs(
        # used windows fourier transform to calculate the spectra
        # rhttp://docs.scipy.org/doc/scipy/reference/tutorial/fftpack.html
        fftpack.fft((-1) ** k * blackman(N) * f)
    ) ** 2
    spectrum /= spectrum.max()
    plt.semilogy(omegas, spectrum)
    plt.ylabel('spectrum (arbitrary units)')
    plt.xlabel('frequency')
    plt.ylim([1e-15, 1.])


########################################################################################################################
# Save Results and the Plot
########################################################################################################################
"""Save Results"""
outfile = data_path + filetag + '.npz'
print('Saving results here: {}'.format(outfile))
np.savez(outfile, **results)

fignum = 1
plt_params = {'legend.fontsize': 'x-large',
              'figure.figsize': (8, 7),
              'axes.labelsize': 'x-large',
              'axes.titlesize': 'xx-large',
              'xtick.labelsize': 'x-large',
              'ytick.labelsize': 'x-large'
}
plt.rcParams.update(plt_params)

""" Plot the Results """
plt.figure(fignum)
fignum += 1
plt.plot(times, results['phi']/(0.5*np.pi))
plt.xlabel('Time (t)')
plt.ylabel('$\\Phi$ in units of $\\pi /2$')
plt.tight_layout()
plt.savefig(plot_path + 'Field_' + filetag + '.pdf')

plt.figure(fignum)
fignum += 1
plt.plot(times, results['current'])
plt.xlabel('Time (t)')
plt.ylabel('Current')
plt.tight_layout()
plt.savefig(plot_path + 'Current_' + filetag + '.pdf')

plt.figure(fignum)
fignum += 1
plt.plot(times, np.abs(x_t))
plt.hlines(1.0, xmin=start, xmax=stop, colors='r', linestyles='--')
plt.xlabel('Time (t)')
plt.ylabel('|X(t)|')
plt.tight_layout()
plt.savefig(plot_path + 'X(t)_' + filetag + '.pdf')

plt.show()
