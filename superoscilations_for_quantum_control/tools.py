"""Contains tools for generating expectations, tracking control, and fitting or interpolating fields"""

import os
import sys
from typing import Dict, Union, List

import psutil

threads = 16
os.environ['NUMEXPR_MAX_THREADS'] = '{}'.format(threads)
os.environ['NUMEXPR_NUM_THREADS'] = '{}'.format(threads)
os.environ['OMP_NUM_THREADS'] = '{}'.format(threads)
os.environ['MKL_NUM_THREADS'] = '{}'.format(threads)
from quspin.operators import hamiltonian, exp_op, quantum_operator  # operators
from quspin.basis import spinful_fermion_basis_1d  # Hilbert space basis
from quspin.tools.measurements import obs_vs_time, project_op  # calculating dynamics
from quspin.tools.evolution import evolve  # evolving system
import numpy as np  # general math functions
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.sparse.linalg import eigsh
from scipy.optimize import curve_fit, nnls
from scipy.constants import hbar, pi
from time import time  # tool for calculating computation time
from tqdm import tqdm
from itertools import product
import matplotlib.pyplot as plt  # plotting library


class InitializeArchive:
    def __init__(self, directory_number):
        self.directory = {
            'data_loc': [
                './Data/TargetsForTracking/',
                './Data/TrackingResults/',
                './Data/BestFit_SCF/',
                './Data/ImportedFieldResults/'
            ],
            'plots_loc': [
                './Plots/TargetsForTracking/',
                './Plots/TrackingResults/',
                './Plots/BestFit_SCF/',
                './Plots/ImportedFieldResults/'
            ],
            'tag': [
               'TGT',
               'Tracked',
               'SCF',
               'ImportField'
            ],
            'current': [
                'current',
                'tracking_current',
                'current',
                'current'
            ],
            'plot name': [
                'Target',
                'Tracked',
                'Superoscillating',
                'Interpolated'
            ]
        }
        self.directory_number = directory_number

    def build_archive(self):
        if os.path.exists('./Data'):
            print('Data folder exists. The location is {}/Data'.format(os.getcwd()))
        else:
            os.mkdir('./Data')
        if os.path.exists('./Plots'):
            print('Plots folder exists. The location is {}/Plots'.format(os.getcwd()))
        else:
            os.mkdir('.Plots')
        print('Proceeding to make individual directories')
        # Make the directories for unique types
        for i in range(4):
            if os.path.exists(self.directory['data_loc'][i]):
                print('Data directory appears to be set up. Trying individual folders: {}/8'.format(2*i + 1))
            else:
                os.mkdir(self.directory['data_loc'][i])
            if os.path.exists(self.directory['plots_loc'][i]):
                print('Plots directory appears to be set up. Trying individual folders: {}/8'.format(2*i+2))
            else:
                os.mkdir(self.directory['plots_loc'][i])
        return self.directory

    def get_dir(self):
        return {
            'data path': self.directory['data_loc'][self.directory_number],
            'plots path': self.directory['plots_loc'][self.directory_number],
            'tag': self.directory['tag'][self.directory_number],
            'plot name': self.directory['plot name'][self.directory_number],
            'current': self.directory['current'][self.directory_number],
        }

    def bundle_results(self, expectations, phi, current, times, dt, psi_0=['NotGiven'], gs_energy=['NotGiven']):
        results = expectations
        results['phi'] = phi
        results['current'] = current
        results['times'] = times
        results['dt'] = dt
        results['psi_0'] = psi_0
        results['gs_energy'] = gs_energy
        return results


class HubbardModel:
    def __init__(self, hopping, interaction, n_up, n_down, angular_frequency, lattice_constant, field_amplitude,
                 chem_potential=0, nx=4, ny=0, cycles=10, n_steps=2000, soc=0, gamma=0, tracking=False, scf=False,
                 reaping=False, int_track=0, lat_track=4, pulses=0):
        """
        A class for Establishing the Hubbard Model. Adapted to
        input units: THz (field), eV (t, U), MV/cm (peak amplitude), Angstroms (lattice cst)
        converts to a'.u, which are atomic units but with energy normalised to t, so
        Note, hbar=e=m_e=1/4pi*ep_0=1, and c=1/alpha=137
        :param hopping: (eV) Hopping strength term, t0
        :param interaction: (eV) Interaction strength term, U
        :param n_up: Number of spin-up Fermions
        :param n_down: Number of spin-down Fermions
        :param angular_frequency: (THz) Laser field angular frequency
        :param lattice_constant: (Angstroms) Lattice constant
        :param field_amplitude: (MV/cm) Laser field amplitude
        :param chem_potential: (eV) Chemical potential
        :param nx: Number of sites on the x-axis, L
        :param ny: Number of sites on the y-axis (usually 0)
        :param cycles: Time in units of field frequency
        :param n_steps: Number of steps for time resolution
        :param soc: Spin-Orbit Coupling
        :param gamma:
        :param tracking: Specify True/False for if this model is for tracking (matters for saving results)
        :param int_track: If tracking, this is the Interaction Strength U of the system you are tracking to
        :param lat_track: If tracking, this is the lattice constant of the system you are tracking to
        """
        # Hubbard Model Parameters
        self.t0_units = hopping
        self.energy_conversion = 27.2113961 / self.t0_units   # Converts energy (27.2113961) eV/Hartree to units of t0
        self.L = nx
        self.t0 = 1.
        self.U = interaction/self.t0_units
        self.mu = chem_potential
        self.N_up = n_up
        self.N_down = n_down
        self.N = self.N_up + self.N_down
        self.SO = soc
        self.gamma = gamma/np.sqrt(self.energy_conversion)
        # Laser Parameters
        self.omega = angular_frequency * self.energy_conversion * 0.0001519828442
        self.frequency = self.omega / (2 * np.pi)
        self.a = (lattice_constant * 1.889726125) / self.energy_conversion
        self.F0 = field_amplitude * 1.944689151e-4 * (self.energy_conversion ** 2)
        # Time Parameters
        self.cycles = cycles
        self.n_steps = n_steps
        self.stop = self.cycles / self.frequency
        if tracking:
            tag = 'params_{}sites-{:.3f}TrackedTo{:.3f}U-{:.2f}t0-{}F0-{}TrackedTo{}a-{}cycles-{}steps'.format(
                self.L, interaction, int_track, hopping, field_amplitude, lattice_constant, lat_track, self.cycles,
                self.n_steps
            )
        else:
            tag = 'params_{}sites-{:.3f}U-{:.2f}t0-{}F0-{}a-{}cycles-{}steps'.format(
                self.L, interaction, hopping, field_amplitude, lattice_constant, self.cycles, self.n_steps
            )
        if scf:
            if reaping:
                tag += '-{}reaped_pulses'.format(pulses)
            else:
                tag += '-{}pulses'.format(pulses)
        self.tag = tag.replace('.', ',')
        print("Angular frequency= %.3f" % self.omega)
        print("Frequency= %.3f" % self.frequency)
        print("Lattice constant= %.3f" % self.a)
        print("Field Amplitude= %.3f" % self.F0)


    def get_groundstate(self, H, imaginary_time=True, n_steps=5):
        """
        Gets the ground state Energy and Psi(0) using Imaginary time or Eigsh
        :param H: The constructed hamiltonian
        :param imaginary_time: True for using imaginary time method, otherwise uses Eigsh
        :param n_steps: number of steps to use for the imaginary time method
        :return:
        """
        start_time = time()
        if imaginary_time:
            print('Calculating ground state using imaginary time method with Quspin Evolve')
            it_steps = np.array([n_steps])          # Number of steps for imaginary time propagation
            psi_0 = H.evolve(np.ones([H.Ns, 1]) + 0j, 0, it_steps, imag_time=True, verbose=False).reshape(-1)
            psi_0 /= np.linalg.norm(psi_0)
            gs_energy = H.expt_value(psi_0)
        else:
            print('Calculating ground state using Quspin Eigsh')
            gs_energy, psi_0 = H.eigsh(k=1, which='SA')
            gs_energy = gs_energy[-1]
            psi_0 /= np.linalg.norm(psi_0)
        gs_energy = np.array([np.real(gs_energy)])[0]
        print('Ground state calculation complete. This took {:.1f} seconds'.format(time()-start_time))
        print('Normalization of Psi_0: {:.4f}'.format(np.linalg.norm(psi_0.reshape(-1))))
        print('Ground State Energy: {:.4f}'.format(gs_energy))
        return gs_energy, psi_0

    def tracking1D(self, target_current, t, nearest_neighbor_expectations):
        current = target_current(t)
        R = np.abs(nearest_neighbor_expectations)
        theta = np.angle(nearest_neighbor_expectations)
        X = current / (2 * self.a * self.t0 * R)
        phi_track = np.real(np.arcsin(-current / (2 * self.a * self.t0 * R) + 0j) + theta)
        return phi_track, theta, X, R

########################################################################################################################
# Old method of scaling units
########################################################################################################################

class parameter_instantiate:
    def __init__(self, field=0, nup=1, ndown=1, nx=2, ny=0, U=0.52, t=0.52, SO=0, F0=10., a=4., pbc=True, gamma=0,
                 mu=0):
        self.nx = nx
        self.pbc = pbc
        if pbc:
            print("Periodic Boundary conditions")
        else:
            print("Open Boundary conditions")
        self.nup = nup
        print("%s up electrons" % self.nup)
        self.ndown = ndown
        print("%s down electrons" % self.nup)
        self.ne = nup + ndown
        # input units: THz (field), eV (t, U), MV/cm (peak amplitude), Angstroms (lattice cst)
        # converts to a'.u, which are atomic units but with energy normalised to t, so
        # Note, hbar=e=m_e=1/4pi*ep_0=1, and c=1/alpha=137
        print("Scaling units to energy of t_0")
        factor = 1. / (t * 0.036749323)
        # factor=1
        self.factor = factor
        # self.factor=1
        self.mu = mu
        print('mu={:.3f}'.format(self.mu))
        self.U = U / t
        self.SO = SO / t
        self.gamma = gamma / np.sqrt(factor)
        # self.gamma=gamma/t
        print('gamma= {:.3f} (t_0)^0.5'.format(self.gamma))
        if type(self.U) is float:
            print("U= %.3f t_0" % self.U)
        else:
            print('onsite potential U list:')
            print(self.U)
        print("SO= %.3f t_0" % self.SO)
        # self.U=U
        self.t = 1.
        print("t_0 = %.3f" % self.t)
        # self.t=t
        # field is the angular frequency, and freq the frequency = field/2pi
        self.field = field * factor * 0.0001519828442
        print("angular frequency= %.3f" % self.field)
        self.freq = self.field / (2. * 3.14159265359)
        print("frequency= %.3f" % self.freq)
        self.a = (a * 1.889726125) / factor
        print("lattice constant= %.3f" % self.a)
        self.F0 = F0 * 1.944689151e-4 * (factor ** 2)
        print("Field Amplitude= %.3f" % self.F0)
        assert self.nup <= self.nx, 'Too many ups!'
        assert self.ndown <= self.nx, 'Too many downs!'
