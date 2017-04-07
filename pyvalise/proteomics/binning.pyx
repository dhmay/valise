#!/usr/bin/env python
"""
A shell of a Python script
"""

import logging
from scipy import sparse
cimport cython
from cpython cimport array
import numpy as np
cimport numpy as np

NDARRAY_DTYPE = np.float32
ctypedef np.float32_t NDARRAY_DTYPE_t

DEFAULT_BIN_SIZE = 1.000495

# default window around each precursor m/z to exclude signal from, when calculating TIC
DEFAULT_PRECURSOR_MZ_WINDOW_EXCLUDE_UP_DOWN = DEFAULT_BIN_SIZE * 2

logger = logging.getLogger(__name__)


def calc_nbins(fragment_min_mz, fragment_max_mz, bin_size):
    """
    Calculate the number of bins needed to cover the specified range.
    This is exposed as a convenience method.
    :param fragment_min_mz:
    :param fragment_max_mz:
    :param bin_size:
    :return:
    """
    return int(float(fragment_max_mz - fragment_min_mz) / float(bin_size)) + 1



def bin_spectra(spectra, fragment_min_mz, fragment_max_mz, bin_size=DEFAULT_BIN_SIZE):
    """
    Bin the MS spectra, processing and yielding one at a time.
    :param spectra: a generator or list. type is spectra.MSSpectrum
    :param fragment_min_mz:
    :param fragment_max_mz:
    :param bin_size:
    :return:
    """
    for spectrum in spectra:
        yield bin_spectrum(spectrum.mz_array, spectrum.intensity_array,
                           fragment_min_mz, fragment_max_mz, bin_size)


def bin_compare_two_spectra(mzs_1, mzs_2, fragment_min_mz, fragment_max_mz, bin_size):
    """
    :param mzs_1:
    :param mzs_2:
    :param fragment_min_mz:
    :param fragment_max_mz:
    :return: the indexes into mzs_1 that match peaks in mzs_2
    """
    mzs_2_bin_set = set()
    for mz in mzs_2:
        if mz < fragment_min_mz or mz > fragment_max_mz:
            continue
        bin_idx = int((mz - fragment_min_mz) / bin_size)
        mzs_2_bin_set.add(bin_idx)
    result = []
    for i in xrange(0, len(mzs_1)):
        mz = mzs_1[i]
        if mz < fragment_min_mz or mz > fragment_max_mz:
            continue
        bin_idx = int((mz - fragment_min_mz) / bin_size)
        if bin_idx in mzs_2_bin_set:
            result.append(i)
    return result




def bin_spectrum(mz_array, intensity_array,
                 float fragment_min_mz, float fragment_max_mz, float bin_size,
                 float precursor_mz, should_normalize, should_normalize_exclude_precursor,
                 float window_exclude_precursor_signal=DEFAULT_PRECURSOR_MZ_WINDOW_EXCLUDE_UP_DOWN):
    """
    Given an array of m/z values and an array of intensity values for fragment peaks
    from one spectrum, produce a binned representation of the spectrum.

    Values for each bin represent the intensity of the most-intense peak falling into the bin.

    This code has been somewhat optimized for speed. It is far faster than a full-Python implementation.

    :param mz_array:
    :param intensity_array:
    :param spectra: a generator or list. type is spectra.MS2Spectrum
    :param fragment_min_mz: low end of the m/z range to represent. Values below this limit ignored.
    :param fragment_max_mz: high end of the m/z range to represent. Values below this limit ignored.
    :param bin_size:
    :param precursor_mz:
    :param should_normalize:
    :param should_normalize_exclude_precursor:
    :param window_exclude_precursor_signal:
    :return: an ndarray
    """
    cdef int bin_idx, peak_idx
    cdef int i
    cdef float mz
    cdef float intensity
    cdef int nbins = calc_nbins(fragment_min_mz, fragment_max_mz, bin_size)
    cdef np.ndarray[NDARRAY_DTYPE_t, ndim=2] scan_matrix = np.zeros((1, nbins), dtype=NDARRAY_DTYPE)
    cdef float frag_intensity_sum = 0.0


    for peak_idx in xrange(0, len(mz_array)):
        mz = mz_array[peak_idx]
        intensity = intensity_array[peak_idx]
        if mz < fragment_min_mz or mz > fragment_max_mz:
            continue
        bin_idx = int((mz - fragment_min_mz) / bin_size)
        if bin_idx < 0 or bin_idx > nbins - 1:
            continue
        # if we're not excluding the precursor from TIC, or the peak is sufficiently far from the precursor,
        # add it to the spectrum intensity sum
        if (not should_normalize_exclude_precursor) or \
                (abs(precursor_mz - mz) > window_exclude_precursor_signal):
            frag_intensity_sum += intensity
        scan_matrix[0, bin_idx] = max(scan_matrix[0, bin_idx], intensity)
    if should_normalize:
        if frag_intensity_sum > 0:
            # divide intensities by sum
            scan_matrix /= frag_intensity_sum
        else:
            # this can happen if the precursor is the only signal in the spectrum!
            logger.debug("0-intensity spectrum (excluding precursor)! Precursor: %d" % precursor_mz)
    return scan_matrix


