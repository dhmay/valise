#!/usr/bin/env python

"""
.ms2 reading.
Some bits grabbed from Alex Hu
See comments at the bottom for example file contents
"""

import logging
from pyvalise.proteomics import spectra

__author__ = "Damon May"
__copyright__ = "Copyright (c) 2012-2014 Fred Hutchinson Cancer Research Center"
__license__ = ""
__version__ = ""

logger = logging.getLogger(__name__)


def read_ms2(ms2_filepath):
    """
    yield all scans in the file at ms2_filepath
    :param ms2_filepath:
    :return:
    """
    ms2file = open(ms2_filepath)

    # store the values we care about for the current scan
    precursor_mz = None
    scan_number = None
    retention_time = None
    charge = None
    fragment_mzs = []
    fragment_intensities = []

    line_number = 0

    n_yielded = 0

    in_preamble = True
    for line in ms2file:
        line_number += 1
        line = line.rstrip()
        chunks = line.split()

        # ignore these types of lines
        if ((chunks[0] == "H") or
            (chunks[0] == "D")):
            continue

        elif chunks[0] == "S":
            logger.debug("S line")
            if len(chunks) != 4:
                raise ValueError("Misformatted line %d.\n%s\n" % (line_number, line))
            # this begins a new scan, so start by writing the old one
            if scan_number and retention_time and fragment_mzs and fragment_intensities and precursor_mz and charge:
                yield spectra.MS2Spectrum(scan_number,
                                          retention_time,
                                          fragment_mzs,
                                          fragment_intensities,
                                          precursor_mz, charge)
                # zero out everything not on this line so that we know if we got it for the next scan
                charge = None
                fragment_mzs = []
                fragment_intensities = []

                n_yielded += 1
                logger.debug("Yielded #%d" % n_yielded)
            else:
                # missing that stuff is only OK if this is the first scan
                if not in_preamble:
                    raise ValueError("Line %d. Tried to write scan with not all values collected. Values: %s" %
                                 (line_number, (scan_number, retention_time, fragment_mzs,
                                                fragment_intensities, precursor_mz, charge)))

            in_preamble = False
            precursor_mz = float(chunks[3])
            scan_number = int(chunks[1])

        elif chunks[0] == "I":
            if chunks[1] == "RTime" or chunks[1] == "RetTime":
                retention_time = float(chunks[2])

        elif chunks[0] == "Z":
            logger.debug("Z line")
            if len(chunks) != 3:
                raise ValueError("Misformatted Z line:\n%s\n" % line)
            charge = int(chunks[1])
        # must be a peak line or junk
        elif len(chunks) == 4:
            fragment_mzs.append(float(chunks[0]))
            fragment_intensities.append(float(chunks[1]))
        # not a recognized line type. Barf.
        else:
            raise ValueError("len(chunks) == %d\n" % len(chunks))

    if scan_number and retention_time and fragment_mzs and fragment_intensities and precursor_mz and charge:
        yield spectra.MS2Spectrum(scan_number,
                                  retention_time,
                                  fragment_mzs,
                                  fragment_intensities,
                                  precursor_mz, charge)
        n_yielded += 1
    else:
        raise ValueError("Tried to write scan with not all values collected")
    logger.debug("Returned %d spectra" % n_yielded)


# here's the preamble and first scan (just the first few peaks) of an .ms2 file

# H       Creation Date   6/22/2016 4:53:49 PM
# H       Extractor       RawConverter
# H       ExtractorVersion        1.0.0.x
# H       Comments        RawConverter written by Lin He, 2014
# H       Comments        RawConverter modified by Yen-Yin Chu, 2015
# H       ExtractorOptions        MSn
# H       AcquisitionMethod       Data-Dependent
# H       InstrumentType  FTMS
# H       DataType        Centroid
# H       ScanType        MS2
# H       Resolution
# H       IsolationWindow
# H       FirstScan       1
# H       LastScan        105807
# S       000514  000514  400.72803
# I       RetTime 4.68
# I       IonInjectionTime        250
# I       ActivationType  HCD
# I       InstrumentType  FTMS
# I       TemperatureFTAnalyzer   -1
# I       Filter  FTMS + p NSI d Full ms2 400.7280@hcd35.00 [110.0000-1613.0000]
# I       PrecursorScan   513
# I       PrecursorInt    56751.9
# I       Predicted Precursor: 400.72803 x 4
# Z       4       1599.89027
# 113.9151 383.3 0 14400
# 129.1021 2006.7 0 20000
# 130.0864 954.2 0 16500
# 131.8942 351.4 0 12900
# 134.075 385.3 0 14400

