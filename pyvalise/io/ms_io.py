#!/usr/bin/env python

"""
Reading and writing of MS-related files (abstracted from ms2_io, mzml_io)
"""

import logging
import ms2_io
import mzml_io
import gzip

__author__ = "Damon May"
__copyright__ = "Copyright (c) 2012-2014 Fred Hutchinson Cancer Research Center"
__license__ = ""
__version__ = ""

logger = logging.getLogger(__name__)


def read_spectra(spectra_file, scan_numbers_to_keep=None, level=None):
    """
    use the appropriate module to read the spectra from a file, lazily
    :param spectra_file:
    :param scan_numbers_to_keep: if specified, a list of scan numbers to keep
    :param level: level of scans to return. Ignored if scan_numbers_to_keep provided
    :return:
    """
    handle = spectra_file
    if spectra_file.name.endswith('.gz'):
        handle = gzip.open(spectra_file.name)
    if '.ms2' in spectra_file.name:
        io_module = ms2_io
    elif '.mzML' in spectra_file.name:
        io_module = mzml_io
    else:
        raise ValueError('yield_kept_spectra, can\'t determine file type from name. Name=%s' % spectra_file.name)
    if scan_numbers_to_keep is not None:
        for spectrum in io_module.retrieve_scans(handle, scan_numbers_to_keep):
            if not level or level == spectrum.level:
                yield(spectrum)
    else:
        for spectrum in io_module.read_ms2_scans(handle):
            if scan_numbers_to_keep is not None and spectrum.scan_number not in scan_numbers_to_keep:
                continue
            yield spectrum

