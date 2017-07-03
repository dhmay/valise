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

# defaults for min and max m/z when generating theoretical spectra
DEFAULT_MIN_THEORETICAL_FRAGMENT_MZ = 0.0
DEFAULT_MAX_THEORETICAL_FRAGMENT_MZ = 5000.0

# value to use for all theoretical peak intensities. Shouldn't matter unless
# the result is passed to something that does filtering.
THEORETICAL_PEAK_INTENSITY = 100.0

log = logging.getLogger(__name__)


def make_peptide_theoretical_spectra(peptide_sequence, aa_modifications,
                                     precursor_charges=(2, 3),
                                     fragment_charges=(1,),
                                     min_fragment_mz=DEFAULT_MIN_THEORETICAL_FRAGMENT_MZ,
                                     max_fragment_mz=DEFAULT_MAX_THEORETICAL_FRAGMENT_MZ):
    """
    given a peptide sequence and a list of modifications, make all the theoretical
    spectra for all modified versions of the peptide. Return a map from ModifiedPeptide to MS2Spectrum
    :param peptide_sequences: 
    :param aa_modifications: 
    :param precursor_charges: charge states for precursor
    :param fragment_charges: 
    :param min_fragment_mz: 
    :param max_fragment_mz: 
    :return: 
    """
    modpep_spectrum_map = {}
    modified_peptides_this_seq = peptides.calc_uniquemass_modpeps(peptide_sequence, aa_modifications)
    for modpep in modified_peptides_this_seq:
        # assign scan numbers sequentially. Dummy RTs of 0.0
        for charge in precursor_charges:
            modpep_spectrum_map[modpep] = make_theoretical_spectrum(modpep, charge, fragment_charges,
                                                                    min_fragment_mz, max_fragment_mz,
                                                                    len(modpep_spectrum_map) + 1,
                                                                    0.0)
    return modpep_spectrum_map


def make_theoretical_spectrum(modified_peptide, charge,
                              fragment_charges=(1,),
                              min_fragment_mz=DEFAULT_MIN_THEORETICAL_FRAGMENT_MZ,
                              max_fragment_mz=DEFAULT_MAX_THEORETICAL_FRAGMENT_MZ,
                              scan_number=0, retention_time=0.0):
    """
    Make a theoretical spectrum for a theoretical peptide with completely-defined modifications
    :param peptide_sequence: 
    :param position_massdiff_list: 
    :param fragment_charges: 
    :param n_terminal_massdiff: 
    :param c_terminal_massdiff: 
    :param min_fragment_mz: 
    :param max_fragment_mz: 
    :param scan_number: 
    :param retention_time: 
    :return: 
    """
    peptide_sequence = modified_peptide.sequence
    position_massdiff_list, nterm_deltamass, cterm_deltamass = modified_peptide.make_posmassdifflist_ntermdiff_ctermdiff()
    theoretical_peak_mzs = peptides.calc_theoretical_peak_mzs(peptide_sequence, fragment_charges,
                                                              position_massdiff_list,
                                                              min_fragment_mz, max_fragment_mz,
                                                              nterm_deltamass=nterm_deltamass,
                                                              cterm_deltamass=cterm_deltamass)
    theoretical_peak_intensities = [THEORETICAL_PEAK_INTENSITY] * len(theoretical_peak_mzs)
    precursor_mz = peptides.calc_mz_from_mplush_charge(modified_peptide.mass, charge)

    return MS2Spectrum(scan_number, retention_time, theoretical_peak_mzs, theoretical_peak_intensities,
                       precursor_mz, charge)


def sqrt_normalize_intensities(intensity_values):
    result = [math.sqrt(x) for x in intensity_values]
    int_sum = sum(intensity_values)
    result = [x / int_sum for x in result]
    return result


def make_binidx_matchcount_map(mzs, fragment_min_mz, fragment_max_mz, bin_size):
    """
    Utility function for calc_proportion_fragments_incommon.
    Notionally bin the mzs from a list of mzs using the specified bin size 
    (don't actually build the binned array), with specified
    lower and upper limits. Return a map from bin indexes to a count of fragments
    :param mzs: 
    :param fragment_min_mz:
    :param fragment_max_mz:
    :param bin_size:
    :return:  a map from bin index to a count of fragments in that bin
    """
    binidx_matchcount_map = {}
    for mz in mzs:
        if mz < fragment_min_mz or mz > fragment_max_mz:
            continue
        bin_idx = int((mz - fragment_min_mz) / bin_size)
        if bin_idx not in binidx_matchcount_map:
            binidx_matchcount_map[bin_idx] = 0
        binidx_matchcount_map[bin_idx] += 1
    return binidx_matchcount_map


def bin_compare_two_spectra(mzs_1, mzs_2, fragment_min_mz, fragment_max_mz,
                            bin_size=binning.DEFAULT_BIN_SIZE):
    """
    :param mzs_1: 
    :param mzs_2: 
    :param fragment_min_mz: 
    :param fragment_max_mz: 
    :param bin_size: 
    :return: indexes into mzs_1 indicating overlap with mzs_2.
    """
    bin_peakcounts_2 = make_binidx_matchcount_map(mzs_2,fragment_min_mz, fragment_max_mz, bin_size=bin_size)
    result_idxs = []
    for i in xrange(0, len(mzs_1)):
        bin_idx = int((mzs_1[i] - fragment_min_mz) / bin_size)
        if bin_idx in bin_peakcounts_2:
            result_idxs.append(i)
    return result_idxs


def plot_ms2_spectrum(ms2_spectrum, peptide_sequence=None, title=None):
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    if title:
        ax.set_title(title)
    normalized_intensities_1 = sqrt_normalize_intensities(ms2_spectrum.intensity_array)
    add_spectrumpeaks_to_ax(ax, ms2_spectrum.mz_array, normalized_intensities_1)

    aa_mods = [peptides.MODIFICATION_IODOACETAMIDE_STATIC]
    if peptide_sequence:
        massdeltas_1, ntermdiff, ctermdiff = peptides.apply_modifications_to_sequence(peptide_sequence, aa_mods)
        theoretical_peak_mzs = peptides.calc_theoretical_peak_mzs(peptide_sequence, [1, 2],
                                                                  massdeltas_1, 0, 5000,
                                                                  nterm_deltamass=ntermdiff,
                                                                  cterm_deltamass=ctermdiff)
        match_idxs = bin_compare_two_spectra(ms2_spectrum.mz_array, theoretical_peak_mzs,
                                                     0, 5000, binning.DEFAULT_BIN_SIZE)
        add_spectrumpeaks_to_ax(ax, [ms2_spectrum.mz_array[i] for i in match_idxs],
                                [normalized_intensities_1[i] for i in match_idxs], color='red',
                                should_invert=False)
    min_mz = min(ms2_spectrum.mz_array)
    max_mz = max(ms2_spectrum.mz_array)
    ax.plot([min_mz, max_mz], [0.0, 0.0], color='black')

    return figure


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
        massdeltas_1, ntermdiff, ctermdiff = peptides.apply_modifications_to_sequence(peptide_sequence1, aa_mods)
        theoretical_peak_mzs = peptides.calc_theoretical_peak_mzs(peptide_sequence1, [1, 2],
                                                                  massdeltas_1, 0, 5000,
                                                                  nterm_deltamass=ntermdiff,
                                                                  cterm_deltamass=ctermdiff)
#        print("theo")
#        print(theoretical_peak_mzs)
        match_idxs = bin_compare_two_spectra(ms2_spectrum1.mz_array, theoretical_peak_mzs,
                                                     0, 5000, binning.DEFAULT_BIN_SIZE)
#        print("match")
#        print(match_idxs)
        add_spectrumpeaks_to_ax(ax, [ms2_spectrum1.mz_array[i] for i in match_idxs],
                                [normalized_intensities_1[i] for i in match_idxs], color='red',
                                should_invert=False)
    add_spectrumpeaks_to_ax(ax, ms2_spectrum2.mz_array, normalized_intensities_2, should_invert=True)
    if peptide_sequence2:
        massdeltas_2, ntermdiff, ctermdiff = peptides.apply_modifications_to_sequence(peptide_sequence1, aa_mods)

        theoretical_peak_mzs = peptides.calc_theoretical_peak_mzs(peptide_sequence2, [1, 2],
                                                                  massdeltas_2, 0, 5000,
                                                                  nterm_deltamass=ntermdiff,
                                                                  cterm_deltamass=ctermdiff)
        match_idxs = bin_compare_two_spectra(ms2_spectrum2.mz_array, theoretical_peak_mzs,
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
        self.info_name_value_dict = None


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



