#!/usr/bin/env python

"""
Data structures and utilities related to spectra
"""

import logging
from pyvalise.proteomics import peptides
import pyximport;
pyximport.install()
from pyvalise.proteomics import binning
from pyvalise.util import charts
import math
import matplotlib.pyplot as plt


__author__ = "Damon May"
__copyright__ = "Copyright (c) 2012-2014 Fred Hutchinson Cancer Research Center"
__license__ = ""
__version__ = ""

log = logging.getLogger(__name__)


def sqrt_normalize_intensities(intensity_values):
    result = [math.sqrt(x) for x in intensity_values]
    int_sum = sum(intensity_values)
    result = [x / int_sum for x in result]
    return result


def plot_two_spectra(ms2_spectrum1, ms2_spectrum2,
                     peptide_sequence1=None, peptide_sequence2=None,
                     title=None):
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    if title:
        ax.set_title(title)
    normalized_intensities_1 = sqrt_normalize_intensities(ms2_spectrum1.intensity_array)
    normalized_intensities_2 = sqrt_normalize_intensities(ms2_spectrum2.intensity_array)
    add_spectrumpeaks_to_ax(ax, ms2_spectrum1.mz_array, normalized_intensities_1)
    aa_mods = [peptides.MODIFICATION_IODOACETAMIDE_STATIC]
    if peptide_sequence1:
        theoretical_peak_mzs = peptides.calc_theoretical_peak_mzs(peptide_sequence1, [1, 2],
                                                                  aa_mods, 0, 5000)
#        print("theo")
#        print(theoretical_peak_mzs)
        match_idxs = binning.bin_compare_two_spectra(ms2_spectrum1.mz_array, theoretical_peak_mzs,
                                                     0, 5000, binning.DEFAULT_BIN_SIZE)
#        print("match")
#        print(match_idxs)
        add_spectrumpeaks_to_ax(ax, [ms2_spectrum1.mz_array[i] for i in match_idxs],
                                [normalized_intensities_1[i] for i in match_idxs], color='red',
                                should_invert=False)
    add_spectrumpeaks_to_ax(ax, ms2_spectrum2.mz_array, normalized_intensities_2, should_invert=True)
    if peptide_sequence2:
        theoretical_peak_mzs = peptides.calc_theoretical_peak_mzs(peptide_sequence2, [1, 2],
                                                                  aa_mods, 0, 5000)
        match_idxs = binning.bin_compare_two_spectra(ms2_spectrum2.mz_array, theoretical_peak_mzs,
                                                     0, 5000, binning.DEFAULT_BIN_SIZE)
        add_spectrumpeaks_to_ax(ax, [ms2_spectrum2.mz_array[i] for i in match_idxs],
                                [normalized_intensities_2[i] for i in match_idxs], color='red',
                                should_invert=True)
    min_mz = min(min(ms2_spectrum1.mz_array), min(ms2_spectrum2.mz_array))
    max_mz = max(max(ms2_spectrum1.mz_array), max(ms2_spectrum2.mz_array))
    ax.plot([min_mz, max_mz], [0.0, 0.0], color='black')

    return figure


def add_spectrumpeaks_to_ax(ax, mz_values, intensity_values, should_invert=False,
                            color='black'):
    if should_invert:
        intensity_values = [-x for x in intensity_values]
    charts.lineplot_peaks(mz_values, intensity_values, ax=ax, linewidth=0.3,
                          color=color)

class MSSpectrum(object):
    """
    represents a single MS spectrum
    """
    def __init__(self, scan_number, retention_time, level, mz_array, intensity_array):
        assert(len(mz_array) == len(intensity_array))
        self.scan_number = scan_number
        self.retention_time = retention_time
        self.mz_array = mz_array
        self.intensity_array = intensity_array
        self.level = level


class MS2Spectrum(MSSpectrum):
    """
    represents a single MS/MS spectrum
    """
    def __init__(self, scan_number, retention_time, mz_array, intensity_array,
                 precursor_mz, charge):
        MSSpectrum.__init__(self, scan_number, retention_time, 2, mz_array, intensity_array)
        self.precursor_mz = precursor_mz
        self.charge = charge

    def generate_mz_intensity_pairs(self):
        for i in xrange(0, len(self.mz_array)):
            yield(self.mz_array[i], self.intensity_array[i])



