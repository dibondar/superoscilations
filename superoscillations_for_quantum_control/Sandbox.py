#####################################################################
# This file is for doing whatever and probably shouldn't have been  #
# committed to Github, but there's an error in removing it so we    #
# have to live with it now, I guess.                                #
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
from scipy.interpolate import UnivariateSpline
from quspin.tools.measurements import project_op
from tools import HubbardModel as fhmodel, InitializeArchive
import psutil
import pandas as pd
import seaborn as sns



# note cpu_count for logical=False returns the wrong number for multi-socket CPUs.
print("logical cores available {}".format(psutil.cpu_count(logical=True)))
t_init = time()
np.__config__.show()


########################################################################################################################
# Declare Parameters of system to be loaded
########################################################################################################################

"""Laser Pulse parameters"""
field = 32.9  # field angular frequency THz
F0 = 10 # Field amplitude MV/cm
a = 4   # Lattice constant Angstroms

"""Parameters for a target or reference field"""
# Hubbard model
L = 10             # system size
N_up = L // 2 + L % 2   # number of fermions with spin up
N_down = L // 2     # number of fermions with spin down
N = N_up + N_down   # number of particles
t0 = 0.52       # hopping strength
U = 0.00 * t0    # interaction strength
pbc = True

# Parameters for evolving the system
cycles = 10     # time in cycles of field frequency
n_steps = 10000  # Number of steps for time resolution

########################################################################################################################
# Do you need to compare these results to something else
########################################################################################################################
compare = True

if compare:
    sim_type_to_compare = 0
    comp = fhmodel(
        nx=L,                       # L for system you want to compare to (should be the same)
        hopping=t0,                 # t0 of system you want to compare to (should be the same)
        interaction=0.0,            # U of system you want to compare to
        n_up=N_up,                  # (should be the same)
        n_down=N_down,              # (should be the same)
        angular_frequency=field,    # (should be the same)
        lattice_constant=a,         # a of system you want to compare to
        field_amplitude=4.0,         # (should be the same)
        chem_potential=0,           # (should be the same)
        cycles=cycles,              # (should be the same)
        n_steps=10000,            # number of steps for system you want to compare to
        ny=0,                       # 1D simulations do not use y-axis
        soc=0,                      # No spin orbit coupling
        gamma=0,                    # No gamma
        tracking=False,             # Are you comparing to a field for tracking
        int_track=1                 # If so, you need to list the U for the system you are tracking to
    )
    comp_lib = InitializeArchive(directory_number=sim_type_to_compare).get_dir()
    comp_path = comp_lib['data path']
    comp_tag = comp_lib['tag'] + comp.tag
    comp_file = np.load(comp_path + comp_tag + '.npz')
    comp_phi = UnivariateSpline(comp_file['times'], comp_file['phi'], k=3, s=0)
    comp_J = UnivariateSpline(comp_file['times'], comp_file['current'], k=3, s=0)

########################################################################################################################
# Load the relevant file
########################################################################################################################
# Reference: 0 for Targets, 1 for Tracking, 2 for superoscillations, 3 for importing anything else
sim_type_to_be_loaded = 1

# Target the filed to be loaded
loader = InitializeArchive(directory_number=sim_type_to_be_loaded).get_dir()
load_path = loader['data path']
load_tag = loader['tag']

# If pasting a direct location, use this
#loadfile = './Data/BestFit_SCF/SCFparams_10sites-0,049TrackedTo0,000U-0,52t0-4,5F0-4TrackedTo4a-10cycles-4000steps_N133pulses.npz'
# Otherwise have the class file find it for you

modU = np.array([0.025, 0.05, 0.075, 0.0875, 0.09375, 0.1])
modU *= t0
pulses = np.array([25, 50, 100, 200, 400, 800])
loadingreaping = np.array([True, True, True, True, True, False])
loadingwithoutreaping = np.array([False, False, False, False, False, False])
counter = np.array([0, 1, 2, 3, 4, 5])

# Now build directory for saving files
directory = InitializeArchive(directory_number=3)
lib = directory.get_dir()

data_path = lib['data path']
plot_path = lib['plots path']
plot_name = lib['plot name']
filetag = lib['tag'] + comp.tag
#filetag = lib['tag'] + 'params_10sites-0,049TrackedTo0,000U-0,52t0-4,5F0-4TrackedTo4a-10cycles-4000steps_N133pulses'

params = dict(
    nx=L,
    hopping=t0,
    n_up=N_up,
    n_down=N_down,
    angular_frequency=field,
    lattice_constant=a,
    field_amplitude=F0,
    chem_potential=0,
    cycles=cycles,
    n_steps=n_steps,
    ny=0,  # 1D simulations do not use y-axis
    soc=0,  # No spin orbit coupling
    gamma=0,  #
    tracking=True,  # Are you loading a field for tracking
    int_track=0,  # If so, you need to list the U for
    scf=True,
)
reapingloadfile = []
noreapingloadfile = []
U_tracker = []
pulses_tracker = []
for _ in modU:
    params['interaction'] = _
    for c in counter:
        # Bundle parameters to pass to Hubbard Model class for unit conversion
        params['pulses'] = pulses[c]
        params['reaping'] = loadingreaping[c]
        here = fhmodel(**params)
        reapingloadfile.append(load_path + load_tag + fhmodel(**params).tag + '.npz')
        U_tracker.append(_)
        pulses_tracker.append(pulses[c])
        print('Gained loading location: {}'.format(reapingloadfile[c]))
        params['reaping'] = loadingwithoutreaping[c]
        noreapingloadfile.append(load_path + load_tag + fhmodel(**params).tag + '.npz')
        print('Gained loading location: {}'.format(noreapingloadfile[c]))

# get the converted units for creating a target field
start = 0.0
stop = here.stop
times, dt = np.linspace(start, stop, num=n_steps, endpoint=True, retstep=True)

phi_holder = []
J_holder = []
error_tracker = []

for _ in noreapingloadfile:
    load = np.load(_)
    phi_interp = UnivariateSpline(load['times'], load['phi'], k=3, s=0)
    J_interp = UnivariateSpline(load['times'], load['current'], k=3, s=0)
    phi_holder.append([phi_interp(load['times'])])
    J_holder.append([J_interp(load['times'])])
    error_tracker.append(float(load['error']))



########################################################################################################################
# Plot results?
########################################################################################################################


fignum = 1

plt.figure(fignum)
fignum += 1
error_tracker = np.array(error_tracker)
s = 10 * error_tracker / error_tracker.min()
plt.scatter((np.array(U_tracker)*1/t0), np.log10(np.array(error_tracker)), s=pulses_tracker, c=pulses_tracker,
            cmap='spring', alpha=0.6)
plt.colorbar(cmap='spring')
plt.xlabel('U ($t_0$)')
plt.ylabel('log$\\int Error^2 dt$')
plt.tight_layout()
plt.show()


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
    plt.semilogy(omegas / sim.omega, spectrum, **kwargs)
    plt.ylabel('spectrum (arbitrary units)')
    plt.xlabel(r'frequency / $\omega$')
    plt.xlim([0, 20])
    plt.ylim([1e-15, 1.])


"""colouring = list(
        plt.cm.get_cmap('plasma')(np.linspace(0., 0.8, 4))
    )

plt.figure(fignum)
fignum += 1
plt.plot(times_target, tgt_current_interp(times_target), color=colouring[0])
plt.plot(times_target, tgt_phi_interp(times_target), color=colouring[1])
plt.plot(times, track_current_interp(times), color=colouring[2])
plt.plot(times, track_phi_interp(times), color=colouring[3])
plt.legend(['Target J', 'Target $\\Phi$', 'Track J', 'Track $\\Phi$'])
plt.xlabel('Time')



plt.figure(fignum)
fignum += 1
plot_spectrum(tgt_current, times_target, tgt, color=colouring[0])
plot_spectrum(tgt_phi, times_target, tgt, color=colouring[1])
plot_spectrum(track_current, times, track, color=colouring[2])
plot_spectrum(track_phi, times, track, color=colouring[3])
plt.legend(['Target J', 'Target $\\Phi$', 'Track J', 'Track $\\Phi$'])
plt.show()"""

