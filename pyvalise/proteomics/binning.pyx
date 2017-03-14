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


def bin_spectra(spectra, fragment_min_mz, fragment_max_mz, bin_size):
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



def bin_spectrum(mz_array, intensity_array,
                 float fragment_min_mz, float fragment_max_mz, float bin_size):
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
    :return:
    """
    cdef int bin_idx, peak_idx
    cdef int i
    cdef float mz
    cdef float intensity
    cdef int nbins = calc_nbins(fragment_min_mz, fragment_max_mz, bin_size)
    cdef np.ndarray[NDARRAY_DTYPE_t, ndim=2] scan_matrix = np.zeros((1, nbins), dtype=NDARRAY_DTYPE)

    for peak_idx in xrange(0, len(mz_array)):
        mz = mz_array[peak_idx]
        intensity = intensity_array[peak_idx]
        if mz < fragment_min_mz or mz > fragment_max_mz:
            continue
        bin_idx = int((mz - fragment_min_mz) / bin_size)
        if bin_idx < 0 or bin_idx > nbins - 1:
            continue
        scan_matrix[0, bin_idx] = max(scan_matrix[0, bin_idx],  intensity)
    return scan_matrix



