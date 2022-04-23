import sys
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from matplotlib import cm as cm
from scipy.signal import stft
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from scipy import signal

sys.path.append('../')  # Go up in folder to get to tools
from tools import parameter_instantiate as hhg  # Used for scaling units.

"""Hubbard model parameters for the reference field (what you are comparing to) """
L_ref = 10  # system size
N_up = L_ref // 2 + L_ref % 2  # number of fermions with spin up
N_down = L_ref // 2  # number of fermions with spin down
N = N_up + N_down  # number of particles
t0_ref = 0.52  # hopping strength
U_ref = 2.0 * t0_ref      # Original Field interaction strength
a_scale_ref = 60
J_scale = 1
pbc = True

"""Laser pulse parameters"""
field = 32.9  # field angular frequency THz
F0 = 10  # Field amplitude MV/cm
a_ref = a_scale_ref * 4  # Lattice constant Angstroms
lat = hhg(field=field, nup=N_up, ndown=N_down, nx=L_ref, ny=0, U=U_ref, t=t0_ref, F0=F0, a=a_ref, pbc=pbc)

"""Hubbard model parameters for the field you are testing """
L_comp = L_ref  # system size
N_up_comp = L_comp // 2 + L_comp % 2  # number of fermions with spin up
N_down_comp = L_comp // 2  # number of fermions with spin down
N_comp = N_up_comp + N_down_comp  # number of particles
t0_comp = 0.52  # hopping strength
U_comp = 2.0 * t0_comp  # Field interaction strength
a_scale_comp = 60
J_scale_comp = 1
pbc = True
a_comp = a_scale_comp * 4  # Lattice constant Angstroms


"""System Evolution Time"""
cycles = 10  # time in cycles of field frequency
n_steps = 2000
start = 0
# real time
# stop = cycles / lat.freq
# scaling time to frequency
stop = cycles
times, delta = np.linspace(start, stop, num=n_steps, endpoint=True, retstep=True)


"""load the expectations of the two fields you want to compare"""
############################################################################
# Directory:                                                               #
#   1. Reference or Target field from Generate_TargetFieldExpectations.py  #
#   2. Tracked control field from Generate_TrackedFieldExpectations.py     #
#   3. Fitted superoscillating control field from FittingSCF.ipynb         #
#   4. An imported field from Generate_ImportedFieldExpectations.py        #
############################################################################

directory = {
    'location': {
        1: './Data/Expectations_TargetsForTracking/',
        2: './Data/Expectations_TrackingResults/',
        3: './Data/BestFit_SCF/',
        4: './Data/Expectations_ImportedField/'
    },
    'tag': {
        1: 'TGT',
        2: 'Tracked',
        3: 'SCF',
        4: 'ImportField'
    },
    'field': {
        1: 'phi',
        2: 'tracking_phi',
        3: 'field',
        4: 'phi'
    },
    'current': {
        1: 'current',
        2: 'tracking_current',
        3: 'current',
        4: 'current'
    },
    'plot names': {
        1: 'Target',
        2: 'Tracked',
        3: 'Superoscillating',
        4: 'Interpolated'
    }
}

# Use the above directory to input number of corresponding expectations to be loaded
dir_num_ref = 2
dir_num_comp = 4
ref_loc = directory['location'][dir_num_ref]
ref_tag = directory['tag'][dir_num_ref] + 'params_{}sites-{}U-{}t0-{}a-{}cycles-{}steps-{}pbc'.format(
    L_ref, U_ref, t0_ref, a_ref, cycles, n_steps, pbc)
out_reference = ref_loc + ref_tag + '.npz'
comp_loc = directory['location'][dir_num_comp]
comp_tag = directory['tag'][dir_num_comp] + 'params_{}sites-{}U-{}t0-{}a-{}cycles-{}steps-{}pbc'.format(
    L_comp, U_comp, t0_comp, a_comp, cycles, n_steps, pbc)
out_compare = comp_loc + comp_tag + '.npz'
exp_ref = np.load(out_reference)
exp_comp = np.load(out_compare)

"""show the expectations available here"""
print('expectations available: {}'.format(exp_ref.files))
sites = [j for j in range(L_ref)]
sites_scf = sites
"""Begin Comparison"""
figparams = directory['plot names'][dir_num_ref] + 'params_{}sites-{}U-{}t0-{}a-{}cycles-{}steps-{}pbc'.format(
    L_ref, U_ref, t0_ref, a_ref, cycles, n_steps, pbc) + '_vs_' + directory['plot names'][dir_num_comp] + \
            'params_{}sites-{}U-{}t0-{}a-{}cycles-{}steps-{}pbc'.format(
                L_comp, U_comp, t0_comp, a_comp, cycles, n_steps, pbc)


J_field_ref = exp_ref[directory['current'][dir_num_ref]] # / a_scale_ref
J_field_comp = exp_comp[directory['current'][dir_num_comp]] # / a_scale_comp
E_ref = exp_ref['H']
E_comp = exp_comp['H']
phi_ref = exp_ref[directory['field'][dir_num_ref]]
phi_comp = exp_comp[directory['field'][dir_num_comp]]

plt.xlabel("Time (cycles)")
plt.ylabel("$J(t)$")
plt.grid(True)
# plt.plot(times, J_field, label='Original Field')
plt.plot(times, J_field_ref, label=directory['plot names'][dir_num_ref] + ' Current')
plt.plot(times, J_field_comp, label=directory['plot names'][dir_num_comp] + ' Current', linestyle='-.')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig('./Plots/Comparison_plots/current-' + figparams + '.pdf', bbox_inches='tight')
plt.show()

plt.xlabel("Time (cycles)")
plt.ylabel("$\\Phi(t)$")
plt.grid(True)
plt.plot(times, phi_ref, label=directory['plot names'][dir_num_ref] + ' Field')
plt.plot(times, phi_comp, linestyle='--', label=directory['plot names'][dir_num_comp] + ' Field')
plt.legend()
plt.tight_layout()
plt.savefig('./Plots/Comparison_plots/phi-' + figparams + '.pdf', bbox_inches='tight')
plt.show()
