"""
Module for analyzing superoscillations in time-domain
"""
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.linalg import eigh
from scipy.integrate import simps
import re

def get_overlap(pulses:pd.DataFrame, times:np.ndarray):
    """
    Return the overlap matrix using the Simpson integration rule
    :param pulses: pandas dataframe containing time sampled pulses
    :param times: 1D numpy array of time intervals
    """
    return np.array(
        [
            [simps(pulses[Ei] * pulses[Ej], times) for Ej in pulses.columns]
                for Ei in pulses.columns
        ]
    )

class CTSuperoscillations(object):
    """
    Class for analyzing superoscillations in time-domain using the definition of superoscillations as eigenvectors.
    """
    def __init__(self, pulses:pd.DataFrame, times:np.ndarray):
        """
        Loading data
        :param pulses: pandas dataframe containing time sampled pulses
        :param times: 1D numpy array of time intervals
        """
        self.pulses = pulses
        self.times = times

        # extract numbers from the column names as str
        self.omega = [
            re.search('[0-9]+.[0-9]+', _).group() for _ in self.pulses.columns
        ]

        # extract the column corresponding to the largest frequency
        self.largest_freq = pulses.columns[
            np.array(self.omega, dtype=np.float).argmax()
        ]

        # Calculate the overlap matrix
        # $$
        #     S_{ij}(-\infty, \infty) = \int_{-\infty}^{\infty} E_i(t) E_j(t) dt
        # $$
        self.S_infty = get_overlap(pulses, times)

    def get_pulse_time_shifts(self):
        """
        Calculate the time shift of pulses such that their peaks are aligned,
        and irregular pre-pulses are ignored. Save the results of calculations as a dictionary self.pulse_shifts.

        Also, this method calculates the observational window (self.observational_window),
        i.e., the period of the largest frequency.

        :return: None
        """

        # dictionary to save save the determined pulse shifts
        self.pulse_shifts = dict()

        # save the positions of peaks for each pulse
        peaks_pulses = dict()

        for colname in self.pulses.columns:
            field = self.pulses[colname]

            # find the positions of peaks
            peaks_indx = find_peaks(field.abs())[0]

            peaks_pulses[colname] = peaks_indx

            # ignore all the peaks before the maximum
            peaks_indx = peaks_indx[peaks_indx >= field.abs().argmax()]

            # find unique spacing between peaks
            spacing_peaks, indx, counts = np.unique(
                np.diff(peaks_indx), return_counts=True, return_index=True
            )

            # find the most common spacing between peaks
            most_freq_spacing = spacing_peaks[counts.argmax()]

            # allow for +/- dt mistake in the peak position determination
            _ = np.argwhere(
                (most_freq_spacing - 1 <= spacing_peaks) &
                (spacing_peaks <= most_freq_spacing + 1)
            )

            self.pulse_shifts[colname] = peaks_indx[indx[_].min()]

        ############################################################################################################

        # Find the observational window for the largest frequency field
        peaks = peaks_pulses[self.largest_freq]
        shift = self.pulse_shifts[self.largest_freq]

        peak_before_shift = peaks[peaks < shift].max()
        peak_after_shift = peaks[peaks > shift].min()

        self.observational_window = peak_after_shift - peak_before_shift

        # Additional offset the pulses in time such that the observational time windows starts at time index 0
        add_offset = shift - peak_before_shift

        for colname in self.pulse_shifts:
            self.pulse_shifts[colname] -= add_offset

    def time_align_pulses(self, pulse_shifts=None):
        """
        Time shift of pulses such that their peaks are aligned, and irregular pre-pulses are ignored.

        :param pulse_shifts: dict. If not specified, the method self.get_pulse_time_shifts() will be called
                            to get self.pulse_shifts
        :return: None
        """

        if pulse_shifts is None:
            self.get_pulse_time_shifts()
            pulse_shifts = self.pulse_shifts

        # temporally shift all the pulses
        for colname, shift in pulse_shifts.items():
            self.pulses[colname] = np.roll(self.pulses[colname], -shift)
