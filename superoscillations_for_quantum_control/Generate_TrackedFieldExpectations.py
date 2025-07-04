#####################################################################
# This file loads the expectations generated by some target field   #
# and uses tracking control to match the optical response (current) #
# to create a new field of given parameters.                        #
# To double-check results, consider importing the tracked field to  #
# Generate_ImportedFieldExpectations.py instead                     #
# For Further details see Phys. Rev. Lett. 124, 183201              #
# Based on the Quspin package.                                      #
# For more info, see:  http://weinbe58.github.io/QuSpin/index.html  #
#####################################################################
from __future__ import print_function, division
import os
import sys
from quspin.operators import hamiltonian, exp_op, quantum_operator  # operators
from quspin.basis import spinful_fermion_basis_1d  # Hilbert space basis
from quspin.tools.measurements import obs_vs_time  # calculating dynamics
from quspin.tools.evolution import evolve   # ODE evolve tool
import numpy as np  # general math functions
from scipy.sparse.linalg import eigsh
from time import time  # tool for calculating computation time
from scipy.interpolate import UnivariateSpline
from scipy.signal.windows import blackman
from scipy import fftpack
import matplotlib.pyplot as plt  # plotting library
from tqdm import tqdm
from tools import HubbardModel as fhmodel, InitializeArchive
import psutil
from tools import parameter_instantiate as hhg  # Used for scaling units.
"""Multi-processing"""
# threads = sys.argv[1]
threads = 16
os.environ['NUMEXPR_MAX_THREADS'] = '{}'.format(threads)
os.environ['NUMEXPR_NUM_THREADS'] = '{}'.format(threads)
os.environ['OMP_NUM_THREADS'] = '{}'.format(threads)
os.environ['MKL_NUM_THREADS'] = '{}'.format(threads)

# note cpu_count for logical=False returns the wrong number for multi-socket CPUs.
print("logical cores available {}".format(psutil.cpu_count(logical=True)))
sim_time_start = time()
np.__config__.show()

########################################################################################################################
# Declare parameters of the target system such that they may be loaded
########################################################################################################################
"""Laser Pulse parameters. These should remain constant across fields (with the exception of a)"""
field = 32.9  # field angular frequency THz
F0 = 10  # Field amplitude MV/cm
a_target = 4   # Lattice constant Angstroms

"""Parameters for a target or reference field"""
# Hubbard model
L_Target = 10             # system size
N_up_target = L_Target // 2 + L_Target % 2   # number of fermions with spin up
N_down_target = L_Target // 2           # number of fermions with spin down
N_target = N_up_target + N_down_target  # number of particles
t0_target = 0.52       # hopping strength
U_target = 0.0 * t0_target    # interaction strength
pbc = True

# Parameters for evolving the system
cycles = 10     # time in cycles of field frequency
n_steps_target = 10000  # Number of steps for time resolution

# Bundle parameters to pass to Hubbard Model class for unit conversion
tgt_params = dict(nx=L_Target, hopping=t0_target, interaction=U_target, n_up=N_up_target, n_down=N_down_target,
                  angular_frequency=field, lattice_constant=a_target, field_amplitude=F0, chem_potential=0,
                  cycles=cycles, n_steps=n_steps_target, ny=0, soc=0, gamma=0)

# Convert units using built-in functions
tgt = fhmodel(**tgt_params)
# Timing parameters must remain the same for both systems
start = 0.0
stop = tgt.stop
times_target, dt_target = np.linspace(start, stop, num=n_steps_target, endpoint=True, retstep=True)

########################################################################################################################
# Load the Target system:
# 0 for a field and Target current generated by Generate_TargetFieldExpectations.py
# 1 for another Tracked field generated by this file
# 2 for a Superoscillating field generated externally
# 3 for Importing anything else or for testing
########################################################################################################################
"""Declare the simulation type you wish to track to using the above as reference"""
sim_type_to_be_loaded = 0
# Build out the directory
loading_lib = InitializeArchive(directory_number=sim_type_to_be_loaded).get_dir()
load_path = loading_lib['data path']
loadfile = load_path + loading_lib['tag'] + tgt.tag + '.npz'
tgt_expectations = dict(np.load(loadfile))
tgt_current = tgt_expectations['current']
tgt_phi = tgt_expectations['phi']
tgt_neighbor = tgt_expectations['neighbor']
tgt_psi_0 = tgt_expectations['psi_0']

# Interpolate the current and field to ensure consistency with time scale
tgt_current_interp = UnivariateSpline(times_target, tgt_current, k=3, s=0)
tgt_phi_interp = UnivariateSpline(times_target, tgt_phi, k=3, s=0)

########################################################################################################################
# Declare the parameters of the system which will be tracked to a target
########################################################################################################################

"""Parameters for the tracking system"""
# Hubbard model
L = L_Target                # number of sites should be the same as the target field
N_up = N_up_target          # number of fermions with spin up should be same as target field
N_down = N_down_target      # number of fermions with spin down should be same as target field
N = N_up + N_down           # number of particles (should be same as target field)
t0 = t0_target              # hopping strength, (Usually same as target field)
U = 1.0 * t0                # interaction strength

# Scaling terms
scaling = True
a_scaling = 5              # Used if you need to scale the lattice constant
j_scaling = 1               # Used if you are scaling the current up instead of the lattice constant
a = a_target * a_scaling    # sets the lattice constant for the current system
n_steps = 10000

track_params = dict(nx=L, hopping=t0, interaction=U, n_up=N_up, n_down=N_down, angular_frequency=field,
                    lattice_constant=a, field_amplitude=F0, chem_potential=0, cycles=cycles, n_steps=n_steps,
                    ny=0, soc=0, gamma=0, tracking=True, int_track=U_target, a_scale=scaling, lat_track=a_target)
track = fhmodel(**track_params)
# Timing parameters must remain the same for both systems
start = 0.0
stop = track.stop
times, dt = np.linspace(start, stop, num=n_steps, endpoint=True, retstep=True)

########################################################################################################################
# Get the directory and set tag for data saving
########################################################################################################################
archive = InitializeArchive(directory_number=1)
lib = archive.get_dir()
data_path = lib['data path']
plot_path = lib['plots path']
plot_name = lib['plot name']

# Get the unique tag for these simulations. Done externally in order to ensure consistency
filetag = lib['tag'] + track.tag


########################################################################################################################
# Define the Hamiltonian
# H = -t0 sum(e^(-i phi(t)) c^*_{i,s} c_{i+1,s} - e^(i phi(t)) c_{i,s} c^*_{i+1,s}) + U sum(n_{i,up} n_i,down})
# or H = sum(-t0 * expiphiconj * hop_right + t0 * expiphi * hop_left) + sum(U * interactions)
########################################################################################################################


# Specify the basis
basis = spinful_fermion_basis_1d(L=track.L, Nf=(track.N_up, track.N_down))
no_checks = dict(check_pcon=False, check_symm=False, check_herm=False)

# Define onsite coupling lists, U * sum(n_{i, up}, n_{i, down})
interactions = [[track.U, _, _] for _ in range(L)]
hop_left = [[1, _, (_ + 1)] for _ in range(L - 1)]         # Allow OBC, but then add PBC
hop_right = [[-1, _, (_ + 1)] for _ in range(L - 1)]       # Allow OBC, but then add PBC

if pbc:
    hop_left.append([1, L-1, 0])
    hop_right.append([-1, L-1, 0])

# Build the Hamiltonian
static = [
    ['n|n', interactions],  # up-down interactions
]

# On-site Hamiltonian for tracking
ham_onsite = hamiltonian(static, [], basis=basis, **no_checks)

# set up operators for the ground state Hamiltonian
hop_left_operator = hamiltonian([["+-|", hop_left], ["|+-", hop_left]], [], basis=basis, **no_checks)
hop_right_operator = hamiltonian([["-+|", hop_right], ["|-+", hop_right]], [], basis=basis, **no_checks)

# Hamiltonian for the groundstate
gs_ham = ham_onsite - track.t0 * (hop_right_operator + hop_left_operator)


# Get the ground state
#psi_0 = tgt_psi_0
it_groundstate = True

# Gets ground state based on eigsh or imaginary time
if it_groundstate:
    E, psi_0 = tgt.get_groundstate(H=gs_ham, imaginary_time=True, n_steps=10)
else:
    E, psi_0 = tgt.get_groundstate(H=gs_ham, imaginary_time=False, n_steps=0)

########################################################################################################################
# Declare functions for tracking and operator dictionaries
########################################################################################################################


def wrapper(new_phi, prev_phi):
    """
    Eliminates shifts of n2pi from phi such that any shift outside of +-pi is kept inside pi
    :param new_phi: The current value of phi
    :param prev_phi: The previous phi which we want to be within +-pi of the new phi
    :return: a new value for phi with phases of 2pi added or subtracted until it is within +-pi of prev phi
    """
    val = new_phi - prev_phi
    while val > np.pi:
        new_phi -= 2 * np.pi
        val = new_phi - prev_phi
    while val < -np.pi:
        new_phi += 2 * np.pi
        val = new_phi - prev_phi
    if val > np.pi:
        new_phi -= 2 * np.pi
        val = new_phi - prev_phi
    elif val < -np.pi:
        new_phi += 2 * np.pi
        val = new_phi - prev_phi
    return new_phi


def nl_evolution(t, psi):
    """
    Solves the non-linear evolution equation: dPsi/dt = H(J(t), Psi)|Psi>
    :param t: Current time
    :param psi: Psi as it is evolved
    :return: dPsi/dt
    """
    nn_expect = hop_left_operator.expt_value(psi)
    phi_t = track.tracking1D(tgt_current_interp, t, nn_expect)[0]
    dPsi_dt = -1j * (ham_onsite.dot(psi))
    dPsi_dt += 1j * track.t0 * np.exp(-1j * phi_t) * hop_left_operator.dot(psi)
    dPsi_dt += 1j * track.t0 * np.exp(1j * phi_t) * hop_right_operator.dot(psi)
    return dPsi_dt

def nl_evolution_(t, psi, phi_t):
    """
    Solves the non-linear evolution equation: dPsi/dt = H(J(t), Psi)|Psi>
    :param t: Current time
    :param psi: Psi as it is evolved
    :return: dPsi/dt
    """
    dPsi_dt = -1j * (ham_onsite.dot(psi))
    dPsi_dt += 1j * track.t0 * np.exp(-1j * phi_t) * hop_left_operator.dot(psi)
    dPsi_dt += 1j * track.t0 * np.exp(1j * phi_t) * hop_right_operator.dot(psi)
    return dPsi_dt

op_dict = dict(H_onsite=ham_onsite)

# Creat lists for tracking
nn_expectations = []
X_track = []
R_track = []
theta_track = []
phi_track = []
current_track = []
psi_t = np.squeeze(psi_0)


# Get the initial conditions for all tracked quantities
nn_expectations.append(hop_left_operator.expt_value(psi_t))
phi_0, theta_0, X_0, R_0 = track.tracking1D(tgt_current_interp, times[0], nn_expectations[0])
phi_track.append(phi_0)
R_track.append(R_0)
theta_track.append(theta_0)
X_track.append(X_0)
init_current_partial = -track.t0 * np.exp(-1j * phi_track[0]) * nn_expectations[0]
current_track.append(1j * track.a * (init_current_partial - init_current_partial.conjugate()))

# Make a finer time time grid

# Implement tqdm as a progress bar
for _ in tqdm(times[:-1]):
    psi_t = np.squeeze(evolve(psi_t, _, np.array([_ + dt]), nl_evolution))
    nn_expectations.append(hop_left_operator.expt_value(psi_t))
    next_phi, next_theta, next_X, next_R = track.tracking1D(tgt_current_interp, _ + dt, nn_expectations[-1])
    next_phi = wrapper(next_phi, phi_track[-1])
    X_track.append(next_X)
    R_track.append(next_R)
    phi_track.append(next_phi)
    theta_track.append(next_theta)
    current_partial = -track.t0 * np.exp(-1j * phi_track[-1]) * nn_expectations[-1]
    current_track.append(1j * track.a * (current_partial - current_partial.conjugate()))

# Implement tqdm as a progress bar
"""
for _ in tqdm(times[:-1]):
    psi_t = np.squeeze(evolve(psi_t, _, np.array([_ + dt]), nl_evolution_, f_params=(phi_track[-1],)))
    nn_expectations.append(hop_left_operator.expt_value(psi_t))
    next_phi, next_theta, next_X, next_R = track.tracking1D(tgt_current_interp, _ + dt, nn_expectations[-1])
    next_phi = wrapper(next_phi, phi_track[-1])
    X_track.append(next_X)
    R_track.append(next_R)
    phi_track.append(next_phi)
    theta_track.append(next_theta)
    current_partial = -track.t0 * np.exp(-1j * phi_track[-1]) * nn_expectations[-1]
    current_track.append(1j * track.a * (current_partial - current_partial.conjugate()))
"""
########################################################################################################################
# Propagate the full Phi to verify the tracking was successful
########################################################################################################################

# Interpolate the Tracking Field to get elements per time
tracking_phi = UnivariateSpline(times, phi_track, k=3, s=0)

def expiphi(t):
    return -track.t0 * np.exp(1j * tracking_phi(t))    # * -track.t0


def expiphiconj(t):
    return -track.t0 * np.exp(-1j * tracking_phi(t))    # * -track.t0

########################################################################################################################
# Add dynamic elements to Hamiltonian so that it may be done through previous method
########################################################################################################################

dynamic_arguments = []
dynamic = [
    ["+-|", hop_left, expiphiconj, dynamic_arguments],          # Spin-up hop left
    ["|+-", hop_left, expiphiconj, dynamic_arguments],          # Spin-down hop left
    ["-+|", hop_right, expiphi, dynamic_arguments],     # Spin-up hop right
    ["|-+", hop_right, expiphi, dynamic_arguments]      # Spin-down hop right
]

# Build Hamiltonian with dynamic elements
# For finding ground-state using Imaginary Time
res_gs_ham = hamiltonian(static + [_[:2] for _ in dynamic], [], basis=basis, **no_checks)
# Hamiltonian for propagation
res_ham = hamiltonian(static, dynamic, basis=basis, **no_checks)

########################################################################################################################
# Define the operator dictionary. This can be used to get expectations, such as current
########################################################################################################################

op_dict['H_full'] = res_ham
# Neighbor interactions, terms are listed in static hamiltonian as there is no time dependence
op_dict['neighbor'] = hamiltonian([["+-|", hop_left], ["|+-", hop_left]], [], basis=basis, **no_checks)
op_dict['left_hop_up'] = hamiltonian([], [["+-|", hop_left, expiphiconj, dynamic_arguments]],
                                      basis=basis, **no_checks)
op_dict['left_hop_down'] = hamiltonian([], [["|+-", hop_left, expiphiconj, dynamic_arguments]],
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


########################################################################################################################
# Get the Ground State, Evolve the Hamiltonian, And Calculate Expectations
########################################################################################################################
# Get the ground state again so that it can be verified

"""if it_groundstate:
    res_E, res_psi_0 = track.get_groundstate(H=res_gs_ham, imaginary_time=True, n_steps=10)
else:
    res_E, res_psi_0 = track.get_groundstate(H=res_ham, imaginary_time=False, n_steps=0)"""

res_E = E
res_psi_0 = psi_0
# Evolve the Hamiltonian
print('Now evolving system')
t_evolve_start = time()
# psi = ham.evolve(v0=psi_0, t0=start, times=times, iterate=True, verbose=False)
# Obs_vs_Time is designed to work with iterate=True, but the generator object is unrealistically difficult to work with
# Until this is fixed, we can simply squeeze the result and directly return the array
psi_res = np.squeeze(res_ham.evolve(v0=res_psi_0, t0=start, times=times, verbose=False))
t_evolve_mid = time()
print('System evolution complete. It took {:.3f} seconds.'.format(t_evolve_mid-t_evolve_start))

print('Now Calculating Expectations')
expectations = obs_vs_time(psi_t=psi_res, times=times, Obs_dict=op_dict)


current_sum_term = (expectations['left_hop_up'] + expectations['left_hop_down'])
tracked_current = 1j * track.a * (current_sum_term - current_sum_term.conjugate())

results = archive.bundle_results(expectations, tracking_phi(times), tracked_current, times, dt, res_psi_0)

x_t = tracked_current / 2 * track.a * np.abs(hop_left_operator.expt_value(psi_res))
print('Expectations Calculated. This took an additional {:.3f} seconds.'.format(time()-t_evolve_mid))

########################################################################################################################
# Save Results and the Plot
########################################################################################################################
"""Save Results"""
outfile = data_path + filetag + '.npz'
print('Saving expectations here: {}'.format(outfile))
np.savez(outfile, **results)

def plot_spectrum(f, t, sim, **kwargs):
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
    plt.semilogy(omegas / sim.omega, spectrum, **kwargs, label=get)
    plt.ylabel('spectrum (arbitrary units)')
    plt.xlabel(r'frequency / $\omega$')
    plt.xlim([0, 20])
    plt.ylim([1e-15, 1.])

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
plt.plot(times_target, tgt_phi_interp(times_target)*(1/(0.5 * np.pi)), label='Original $\\Phi (t)$')
plt.plot(times, results['phi']/(0.5*np.pi), linestyle='--', label='$\\Phi (t)$ for Verification')
plt.plot(times, np.array(phi_track)/(0.5*np.pi), alpha=0.8, label='$\\Phi (t)$ from Single Step Evolution')
plt.xlabel('Time (t)')
plt.ylabel('$\\Phi$ in units of $\\pi /2$')
plt.legend()
plt.tight_layout()
plt.savefig(plot_path + 'Field_' + filetag + '.pdf')

plt.figure(fignum)
fignum += 1
plt.plot(times, tgt_current_interp(times), label='Target Current')
plt.plot(times, current_track, linestyle='-.', label='Current from Single Step Evolution')
plt.plot(times, results['current'], linestyle='--', label='Verified Current from Tracking')
plt.xlabel('Time (a.u.)')
plt.ylabel('Current (a.u.)')
plt.legend()
plt.tight_layout()
plt.savefig(plot_path + 'Current_' + filetag + '.pdf')

plt.figure(fignum)
fignum += 1
plt.plot(times, np.abs(x_t), label='Verified |X(t)| from Tracking')
plt.plot(times, np.abs(X_track), linestyle='--', label='|X(t)| from Single Step Evolution')
plt.hlines(1.0, xmin=start, xmax=stop, colors='r', linestyles='--')
plt.xlabel('Time (a.u.)')
plt.ylabel('|X(t)|')
plt.legend()
plt.tight_layout()
plt.savefig(plot_path + 'X(t)_' + filetag + '.pdf')

plt.figure(fignum)
fignum += 1
pm = np.array(phi_track)-np.array(theta_track)
plt.plot(times, pm)
plt.xlabel('Time (a.u.)')
plt.ylabel('Phi - Theta')


plt.figure(fignum)
fignum += 1
get = 'Transform Limited Pulse'
plot_spectrum(tgt_phi_interp(times), times, tgt)
get = 'Tracking Field'
plot_spectrum(results['phi'], times, track)
get = 'Current'
plot_spectrum(results['current'], times, track, linestyle='--')
plt.legend()
plt.tight_layout()
plt.savefig(plot_path + 'SpectrumAnalysis' + filetag + '.pdf')
plt.show()

print('Done')
