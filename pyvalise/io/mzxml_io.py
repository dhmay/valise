#!/usr/bin/env python

"""
mzML reading
"""

import logging
from pyvalise.proteomics import spectra
from pyteomics import mzxml


__author__ = "Damon May"
__copyright__ = "Copyright (c) 2012-2014 Fred Hutchinson Cancer Research Center"
__license__ = ""
__version__ = ""

PEPXML_NS_URL = "http://regis-web.systemsbiology.net/pepXML"
PEPXML_NS = "{" + PEPXML_NS_URL + "}"

# minimum ratio, for pegging 0-ratio values
MIN_RATIO = 0.001

log = logging.getLogger(__name__)


def retrieve_scans(mzxml_file, scan_numbers):
    """
    retrieve scans from a preindexed mzML file.
    :param mzxml_file:
    :param scan_numbers
    :return:
    """
    with mzxml.MzXML(mzxml_file) as reader:
        for scan_number in scan_numbers:
            spectrum = read_scan(reader.get_by_id("controllerType=0 controllerNumber=1 scan=%d" % scan_number))
            yield spectrum


def read_ms2_scans(mzxml_file):
    return read_scans(mzxml_file, ms_levels=[2])


def read_scans(mzxml_file, ms_levels=(1, 2)):
    """
    yields all spectra from an mzML file with level in ms_levels, or
    all processable scans if ms_levels not specified
    :param mzxml_file:
    :param ms_levels:
    :param min_pprophet:
    :return:
    """
    with mzxml.MzXML(mzxml_file) as reader:
        for scan in reader:
            if int(scan['msLevel']) in ms_levels:
                yield read_scan(scan)


def read_scan(scan):
    """

    :param scan:
    :return:
    """
    # see below for the byzantine intricacies of the scan object
    ms_level = int(scan['msLevel'])
    scan_number = int(scan['id'])
    mz_array = scan['m/z array']
    intensity_array = scan['intensity array']
    retention_time = float(scan['retentionTime'])
    if ms_level == 1:
        return spectra.MSSpectrum(scan_number, retention_time, ms_level,
                                  mz_array,
                                  intensity_array)
    elif ms_level == 2:
        print(scan)
        precursor_selected_ion_map = scan['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]
        precursor_mz = precursor_selected_ion_map['selected ion m/z']
        if 'charge state' in precursor_selected_ion_map:
            charge = precursor_selected_ion_map['charge state']
        elif 'possible charge state' in precursor_selected_ion_map:
            charge = precursor_selected_ion_map['possible charge state']
        else:
            raise ValueError("Could not find charge for scan %d" % scan_number)

        return spectra.MS2Spectrum(scan_number, retention_time,
                                   mz_array,
                                   intensity_array,
                                   precursor_mz, charge)
    else:
        raise ValueError("Unhandleable scan level %d" % ms_level)


# scan keys:
# ['polarity', 'basePeakIntensity', 'scanType', 'intensity array',
# 'collisionEnergy', 'retentionTime', 'basePeakMz', 'peaksCount',
# 'msLevel', 'lowMz', 'num', 'm/z array', 'totIonCurrent', 'filterLine',
# 'highMz', 'id', 'precursorMz']

# scan:
# {'polarity': '+',
# 'basePeakIntensity': 3566.73,
# ro'scanType': 'Full',
# 'intensity array': array([ 1644.8190918 ,  1435.2199707 ,  1511.34741211,  1842.39599609,
#        1631.34814453,  1739.75708008,  1924.64733887,  2468.96166992,
#        2293.41625977,  2619.66186523,  2719.72924805,  2413.22680664,
#        2607.99462891,  2888.88549805,  2599.54492188,  3323.59326172,
#        3566.73071289], dtype=float32), 'collisionEnergy': 25.0,
# 'retentionTime': 1205.99,
# 'basePeakMz': 1726.38,
# 'peaksCount': 17,
# 'msLevel': '2', 'lowMz': 172.891,
# 'num': '14',
# 'm/z array': array([  172.89068604,   179.21194458,   212.1013031 ,   222.95936584,
#         229.01470947,   246.40187073,   262.47314453,   332.75473022,
#         446.77319336,   470.28076172,   473.73666382,   668.23144531,
#         883.44836426,  1028.12353516,  1282.01599121,  1608.28796387,
#        1726.38085938], dtype=float32),
# 'totIonCurrent': 39231.3,
# 'filterLine': 'FTMS + p NSI d Full ms2 678.34@hcd25.00 [140.00-2100.00]',
#  'highMz': 1726.38, 'id': '14',
# 'precursorMz': [{'precursorIntensity': 6760.85, 'activationMethod': 'HCD', 'precursorMz': 678.0021461}]}


