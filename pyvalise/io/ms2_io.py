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
import StringIO

__author__ = "Damon May"
__copyright__ = "Copyright (c) 2012-2014 Fred Hutchinson Cancer Research Center"
__license__ = ""
__version__ = ""

logger = logging.getLogger(__name__)


def retrieve_scans(ms2_file, scan_numbers, precursor_from_zline=True, should_calc_zs_mz_diffs=False):
    """
    retrieve only the scans in the scan_number list.
    For a huge performance boost if we're reading a relatively small portion of the file,
    I'm reading the lines containing just the scans we care about into a buffer, then
    passing that buffer to read_ms2_scans.
    That could be made far more memory-efficient by parsing the scans one at a time.
    :param ms2_file:
    :param scan_numbers:
    :param precursor_from_zline:
    :param should_calc_zs_mz_diffs:
    :return:
    """
    logger.debug("Building buffer with just kept scans...")
    buffer_with_scans = StringIO.StringIO()
    is_in_header = True
    is_in_kept_scan = False
    n_kept_in_buffer = 0
    for line in ms2_file:
        if line.startswith("S"):
            is_in_header = False
            chunks = line.split()
            scan_number = int(chunks[1])
            if scan_number in scan_numbers:
                is_in_kept_scan = True
                n_kept_in_buffer += 1
            else:
                is_in_kept_scan = False
        if is_in_header or is_in_kept_scan:
            buffer_with_scans.write(line)
    logger.debug("Built buffer. It has %d scans in it" % n_kept_in_buffer)

    for scan in read_scans(StringIO.StringIO(buffer_with_scans.getvalue()), precursor_from_zline=precursor_from_zline,
                           should_calc_zs_mz_diffs=should_calc_zs_mz_diffs):
        if scan.scan_number in scan_numbers:
            yield scan


def read_header_namevalues(ms2_file):
    """
    Read a dict of name-value pairs from the header of a .ms2 file
    :param ms2_file:
    :return:
    """
    result = {}
    for line in ms2_file:
        chunks = line.rstrip().split('\t')
        if chunks[0] == "H":
            if chunks[1].startswith("@") and "=" in chunks[1]:
                name, value = chunks[1][1:].split("=")
                result[name] = value
        else:
            # after the header lines, quit
            break
    return result


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


def prepare_spectrum(scan_number, precursor_mz, retention_time, fragment_mzs, fragment_intensities,
                        all_zline_mplush, all_charges, info_name_value_dict, require_rt):
    """
    Try to prepare a spectrum from the information given. If can't, raise ValueError
    :param scan_number:
    :param precursor_mz:
    :param retention_time:
    :param fragment_mzs:
    :param fragment_intensities:
    :param all_zline_mplush:
    :param all_charges:
    :param info_name_value_dict:
    :param require_rt:
    :return:
    """
    should_write_scan = False
    if scan_number and (retention_time or not require_rt) and fragment_mzs and fragment_intensities \
            and precursor_mz is not None and all_charges:
        if not retention_time:
            # if we get here, rt is allowed to be missing, so set it to be 0.0
            retention_time = 0.0

        # use the first z-line for charge determination. This is totally arbitrary!
        # Must make use of all_charges.
        charge = all_charges[0]

        # try to figure out mz, charges and masses
        if precursor_mz == 0:
            # we don't have a precursor m/z. Try to infer it from z-lines
            if len(all_charges) == len(set(all_charges)):
                # only one line per charge. Good! Pick the first for precursor m/z determination.
                zline_mass = all_zline_mplush[0]
                precursor_mz = peptides.calc_mz_from_mplush_charge(zline_mass, charge)
                #if scan_number in [15904, 16011]:
                #    print("Calced mz: mass={}, charge={}, mz={}".format(zline_mass, charge, precursor_mz))
            else:
                # cannot determine precursor m/z unambiguously, so don't write this scan
                pass
        else:
            # we do have a precursor m/z, so we don't care if there's multiple reporting of the same charge.
            # In that case, collapse them.
            all_charges = list(set(all_charges))
        if precursor_mz > 0:
            should_write_scan = True
    if should_write_scan:
        spectrum = spectra.MS2Spectrum(scan_number, retention_time, fragment_mzs, fragment_intensities,
                                       precursor_mz, charge)
        spectrum.all_charges = all_charges
        info_name_value_dict['all_charges'] = ','.join([str(x) for x in all_charges])
        spectrum.info_name_value_dict = info_name_value_dict
        return spectrum
    else:
        raise ValueError("Cannot write spectrum. Data: scan_number: %s. rt: %s, precursor_mz: %s, fragment_mzs: %s, fragment_intensities: %s, charges: %s" %
                             (str(scan_number), str(retention_time),
                              str(precursor_mz), len(fragment_mzs) > 0, len(fragment_intensities) > 0, str(all_charges)))


def read_scans(ms2_file, precursor_from_zline=True, should_calc_zs_mz_diffs=False,
               require_rt=False, should_renumber_if_needed=False):
    """
    yield all scans in the file at ms2_filepath.
    If the scan has multiple Z lines, keep the first one.
    :param ms2_file:
    :param precursor_from_zline:
    :param should_calc_zs_mz_diffs:
    :param require_rt:
    :param should_renumber_if_needed: Should we renumber the scans if they're all numbered 0?
    Alternatively, refuse to yield scans with scan number 0.
    :return:
    """

    # store the values we care about for the current scan
    precursor_mz = None
    scan_number = 0
    retention_time = None
    info_name_value_dict = {}
    fragment_mzs = []
    fragment_intensities = []
    all_charges = list()
    all_zline_mplush = list()

    line_number = 0

    n_yielded = 0

    # did we actually renumber any scans?
    renumbered_scans = False

    is_before_first_scan = True
    for line in ms2_file:
        line_number += 1
        line = line.rstrip()
        chunks = line.split()

        zline_precursor_mz = None

        # ignore these types of lines
        if ((chunks[0] == "H") or
            (chunks[0] == "D")):
            continue

        elif chunks[0] == "S":
            logger.debug("S line")
            if len(chunks) != 4:
                raise ValueError("Misformatted line %d.\n%s\n" % (line_number, line))
            # this begins a new scan, so start by writing the old one

            try:
                spectrum = prepare_spectrum(scan_number, precursor_mz, retention_time, fragment_mzs, fragment_intensities,
                        all_zline_mplush, all_charges, info_name_value_dict, require_rt)
                yield spectrum
                n_yielded += 1
                logger.debug("Yielded #%d" % n_yielded)
            except Exception as e:
                if not is_before_first_scan:
                    logger.debug("Incomplete scan!")
                    logger.debug(str(e))

            # done with writing old scan.
            # new scan, so reset mz and intensity lists and other things
            retention_time = None
            all_charges = list()
            all_zline_mplush = list()
            info_name_value_dict = {}
            fragment_mzs = []
            fragment_intensities = []

            precursor_mz = float(chunks[3])
            old_scan_number = scan_number
            scan_number = int(chunks[1])
            if scan_number == 0 and should_renumber_if_needed:
                renumbered_scans = True
                scan_number = old_scan_number + 1
            is_before_first_scan = False

        elif chunks[0] == "I":
            logger.debug("I line: %s" % line)
            if chunks[1] == "RTime" or chunks[1] == "RetTime":
                retention_time = float(chunks[2])
            elif chunks[1].startswith("@") and "=" in line:
                name, value = chunks[1].split("=")
                name = name[1:]
                info_name_value_dict[name] = value

        elif chunks[0] == "Z":
            logger.debug("Z line: %s" % line)
            if len(chunks) != 3:
                raise ValueError("Misformatted Z line:\n%s\n" % line)
            thisline_charge = int(chunks[1])
            if thisline_charge is None or thisline_charge == 0:
                continue
            if thisline_charge > 0:
                all_charges.append(thisline_charge)
                all_zline_mplush.append(float(chunks[2]))

        # must be a peak line or junk
        elif len(chunks) == 4 or len(chunks) == 2:
            fragment_mzs.append(float(chunks[0]))
            fragment_intensities.append(float(chunks[1]))
        # not a recognized line type. Barf.
        else:
            print("Bad line:\n*\n%s\n*" % line)
            raise ValueError("len(chunks) == %d\n" % len(chunks))

    try:
        spectrum = prepare_spectrum(scan_number, precursor_mz, retention_time, fragment_mzs, fragment_intensities,
                all_zline_mplush, all_charges, info_name_value_dict, require_rt)
        yield spectrum
        n_yielded += 1
        logger.debug("Yielded #%d" % n_yielded)
    except Exception as e:
        logger.debug("Incomplete scan!")
        logger.debug(str(e))
    if renumbered_scans:
        logger.debug("Renumbered one or more scans.")
    logger.debug("Returned %d spectra" % n_yielded)

#    # this is for figuring out jus what is up with the Z-line and S-line precursor m/z / mass values
#    if should_calc_zs_mz_diffs:
#        mycharts = [charts.hist(zline_sline_precursor_deltas, title='zline-sline precursor deltas')]
#        mycharts.append(charts.scatterplot(zline_sline_masses, zline_sline_precursor_deltas, title='zline-sline vs. mass',
#                                     xlabel='mass', ylabel='delta'))
#        print("Median z-s precursor delta: %f" % np.median(zline_sline_precursor_deltas))
#        print("Mean z-s precursor delta: %f" % np.mean(zline_sline_precursor_deltas))
#        print("Min z-s precursor delta: %f" % np.min(zline_sline_precursor_deltas))
#        print("Max z-s precursor delta: %f" % np.max(zline_sline_precursor_deltas))
#        chartfile = open('deltas.pdf','w')
#        charts.write_pdf(mycharts, chartfile)
#        chartfile.close()


def write_ms2(scans, outfile, header_name_value_dict=None):
    """

    :param scans:
    :param outfile:
    :return:
    """
    write_header(outfile, header_name_value_dict)
    for scan in scans:
        write_scan(scan, outfile)


def write_header(outfile, name_value_dict=None):
    """

    :param outfile:
    :return:
    """
    outfile.write("H\tCreationDate %s\n" % datetime.date.today())
    outfile.write("H\tExtractor\tDamon Homebrew\n")
    outfile.write("H\tExtractor version\t0.0\n")
    if name_value_dict:
        for name in name_value_dict:
            outfile.write("H\t@%s=%s\n" % (name, name_value_dict[name]))


def write_scan(ms2_spectrum, outfile):
    """

    :param ms2_spectrum:
    :param outfile:
    :return:
    """
    outfile.write("S\t%d\t%d\t%f\n" % (ms2_spectrum.scan_number, ms2_spectrum.scan_number,
                                     ms2_spectrum.precursor_mz))
    if ms2_spectrum.retention_time:
        outfile.write("I\tRTime\t%f\n" % ms2_spectrum.retention_time)
    if ms2_spectrum.info_name_value_dict:
        for name in ms2_spectrum.info_name_value_dict:
            outfile.write("I\t@%s=%s\n" % (name, ms2_spectrum.info_name_value_dict[name]))
    for i in xrange(0, len(ms2_spectrum.all_charges)):
        zline_mplush = peptides.calc_mplush_from_mz_charge(ms2_spectrum.precursor_mz,  ms2_spectrum.all_charges[i])
        outfile.write("Z\t%d\t%f\n" % (ms2_spectrum.all_charges[i], zline_mplush))
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

