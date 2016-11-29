#!/usr/bin/env python

"""
.ms2 reading.
Some bits grabbed from Alex Hu
See comments at the bottom for example file contents
"""

import logging
from pyvalise.proteomics import spectra, peptides
from pyvalise.util import charts
import numpy as np
import datetime

__author__ = "Damon May"
__copyright__ = "Copyright (c) 2012-2014 Fred Hutchinson Cancer Research Center"
__license__ = ""
__version__ = ""

logger = logging.getLogger(__name__)


def retrieve_scans(ms2_file, scan_numbers, precursor_from_zline=True, should_calc_zs_mz_diffs=False):
    """
    retrieve only the scans in the scan_number list
    :param ms2_file:
    :param scan_numbers:
    :param precursor_from_zline:
    :param should_calc_zs_mz_diffs:
    :return:
    """
    for scan in read_scans(ms2_file, precursor_from_zline=precursor_from_zline,
                           should_calc_zs_mz_diffs=should_calc_zs_mz_diffs):
        if scan.scan_number in scan_numbers:
            yield scan


def read_ms2_scans(ms2_file, precursor_from_zline=True, should_calc_zs_mz_diffs=False):
    """
    Silly cover method, because ms2 files only contain ms2 scans. For consistency with mzml_io
    :param ms2_file:
    :param precursor_from_zline:
    :param should_calc_zs_mz_diffs:
    :return:
    """
    return read_scans(ms2_file, precursor_from_zline=precursor_from_zline,
                      should_calc_zs_mz_diffs=should_calc_zs_mz_diffs)


def read_scans(ms2_file, precursor_from_zline=True, should_calc_zs_mz_diffs=False,
               require_rt=False):
    """
    yield all scans in the file at ms2_filepath
    :param ms2_file:
    :param precursor_from_zline:
    :param should_calc_zs_mz_diffs:
    :param require_rt:
    :return:
    """

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

    # differences in m/z between that on the s-line and that calculated from the z-line
    zline_sline_precursor_deltas = []
    zline_sline_masses = []

    for line in ms2_file:
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
            if scan_number and (retention_time or not require_rt) and fragment_mzs and fragment_intensities \
                    and precursor_mz and charge:
                if not retention_time:
                    # if we get here, rt is allowed to be missing, so set it to be 0.0
                    retention_time = 0.0
                # sometimes, a Z line will have a 0 charge. Punt on those
                if charge is not None and charge > 0:
                    logger.debug("0 charge!")
                    yield spectra.MS2Spectrum(scan_number,
                                              retention_time,
                                              fragment_mzs,
                                              fragment_intensities,
                                              precursor_mz, charge)
                    n_yielded += 1
                    logger.debug("Yielded #%d" % n_yielded)
                # zero out everything not on this line so that we know if we got it for the next scan
                charge = None

            else:
                logger.debug("Incomplete scan!")

            # new scan, so reset mz and intensity lists
            fragment_mzs = []
            fragment_intensities = []

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
            z_precursor_mplush = float(chunks[2])
            
            if charge != 0:
                zline_precursor_mz = peptides.calc_mz_from_mplush_charge(z_precursor_mplush, charge)
            if should_calc_zs_mz_diffs:
                if abs(zline_precursor_mz - precursor_mz) > 0.5:
                    pass
                    #print("scan %s, z=%d, S* %f Z* %f, Z_mz* %f" % (scan_number, charge, precursor_mz, z_precursor_mplush, zline_precursor_mz))

                    if (abs((zline_precursor_mz - precursor_mz) * charge)) > 0.5:
                        print("%d    %f    %f" % (scan_number, (zline_precursor_mz - precursor_mz) * charge,
                                            ((zline_precursor_mz - precursor_mz) * charge) % peptides.HYDROGEN_MASS))
                diff_mod = abs((zline_precursor_mz - precursor_mz) * charge) % 1.000495
                while diff_mod > 1.000495 / 2:
                    diff_mod -= 1.000495
                if abs(diff_mod) > 0.2:
                    print("HUH?! %f" % diff_mod)
                #zline_sline_precursor_deltas.append(diff_mod)
                zline_sline_precursor_deltas.append((zline_precursor_mz - precursor_mz) * charge)
                zline_sline_masses.append(precursor_mz * charge)
            if precursor_from_zline:
                precursor_mz = zline_precursor_mz
        # must be a peak line or junk
        elif len(chunks) == 4 or len(chunks) == 2:
            fragment_mzs.append(float(chunks[0]))
            fragment_intensities.append(float(chunks[1]))
        # not a recognized line type. Barf.
        else:
            print("Bad line:\n*\n%s\n*" % line)
            raise ValueError("len(chunks) == %d\n" % len(chunks))

    if scan_number and (retention_time or not require_rt) and fragment_mzs and \
            fragment_intensities and precursor_mz and charge:
        yield spectra.MS2Spectrum(scan_number,
                                  retention_time,
                                  fragment_mzs,
                                  fragment_intensities,
                                  precursor_mz, charge)
        n_yielded += 1
    else:
        logger.debug("Tried to write scan with not all values collected")
    logger.debug("Returned %d spectra" % n_yielded)

    # this is for figuring out jus what is up with the Z-line and S-line precursor m/z / mass values
    if should_calc_zs_mz_diffs:
        mycharts = [charts.hist(zline_sline_precursor_deltas, title='zline-sline precursor deltas')]
        mycharts.append(charts.scatterplot(zline_sline_masses, zline_sline_precursor_deltas, title='zline-sline vs. mass',
                                     xlabel='mass', ylabel='delta'))
        print("Median z-s precursor delta: %f" % np.median(zline_sline_precursor_deltas))
        print("Mean z-s precursor delta: %f" % np.mean(zline_sline_precursor_deltas))
        print("Min z-s precursor delta: %f" % np.min(zline_sline_precursor_deltas))
        print("Max z-s precursor delta: %f" % np.max(zline_sline_precursor_deltas))
        chartfile = open('deltas.pdf','w')
        charts.write_pdf(mycharts, chartfile)
        chartfile.close()


def write_ms2(scans, outfile):
    """

    :param scans:
    :param outfile:
    :return:
    """
    write_header(outfile)
    for scan in scans:
        write_scan(scan, outfile)


def write_header(outfile):
    """

    :param outfile:
    :return:
    """
    outfile.write("H\tCreationDate %s\n" % datetime.date.today())
    outfile.write("H\tExtractor\tDamon Homebrew\n")
    outfile.write("H\tExtractor version\t0.0\n")


def write_scan(ms2_spectrum, outfile):
    """

    :param ms2_spectrum:
    :param outfile:
    :return:
    """
    outfile.write("S\t%d\t%d\t%f\n" % (ms2_spectrum.scan_number, ms2_spectrum.scan_number,
                                     ms2_spectrum.precursor_mz))
    outfile.write("I\tRTime\t%f\n" % ms2_spectrum.retention_time)
    zline_mplush = peptides.calc_mplush_from_mz_charge(ms2_spectrum.precursor_mz,  ms2_spectrum.charge)
    outfile.write("Z\t%d\t%f\n" % (ms2_spectrum.charge, zline_mplush))
    for i in xrange(0, len(ms2_spectrum.intensity_array)):
        outfile.write("%s %s\n" % (ms2_spectrum.mz_array[i], ms2_spectrum.intensity_array[i]))


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

