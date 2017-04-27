#!/usr/bin/env python

"""
file reading/writing
"""

import logging
import csv
import gzip


log = logging.getLogger(__name__)


def read_strings_from_file(infile):
    """read a file that consists of floats, one per line, into a list."""
    lines = infile.read().splitlines()
    return [line.strip() for line in lines]


def read_floats_from_file(infile):
    """read a file that consists of floats, one per line, into a list."""
    return [float(line) for line in read_strings_from_file(infile)]


def write_stringlist_to_file(stringlist, outfile):
    """write a list of strings to a file, one per line"""
    for stri in stringlist:
        outfile.write(stri + "\n")
    outfile.close()


def open_maybe_gzipped(filepath, mode='r'):
    """
    If filepath ends in .gz, open with gzip. Otherwise, open regularly
    :param filepath: 
    :param mode: 
    :return: 
    """
    if filepath.endswith('.gz'):
        return gzip.open(filepath, mode)
    else:
        return open(filepath, mode)


def read_namevalue_dict_file(infile, sep='='):
    result = {}
    for line in read_strings_from_file(infile):
        # comment line
        if line.startswith('#'):
            continue
        if sep not in line:
            raise ValueError("Separator character not found on line")
        chunks3 = line.partition(sep)
        name = chunks3[0].strip()
        value = chunks3[2].strip()
        result[name] = value
    return result


def read_spreadsheet_columns_withheader(infile, sep='\t'):
    """
    read a csv spreadsheet into a map from column names to lists of values.
    Values will be strings.
    """
    # TODO: handle non-strings in some fun way.

    dictreader = csv.DictReader(BlankCommentCSVFile(infile), delimiter=sep)

    result = {}
    for fieldname in dictreader.fieldnames:
        result[fieldname] = []
    for row in dictreader:
        for fieldname in result:
            result[fieldname].append(row[fieldname])
    return result


def read_spreadsheet_rowmaps_withheader(infile, sep='\t'):
    """
    read a csv spreadsheet into a list of mapss from column names to values.
    Values will be strings.
    """
    # TODO: handle non-strings in some fun way.

    dictreader = csv.DictReader(BlankCommentCSVFile(infile), delimiter=sep)
    log.debug("read_spreadsheet_rowmaps_withheader: loading rows from file %s. Columns:" % infile.name)
    log.debug(dictreader.fieldnames)
    return list(dictreader)


class BlankCommentCSVFile:
    """Wrapper to ignore comment or blank lines when reading a csv file. From here:
    http://bytes.com/topic/python/answers/513222-csv-comments"""
    def __init__(self, fp):
        self.fp = fp

    def __iter__(self):
        return self

    def next(self):
        line = self.fp.next()
        if not line.strip() or line[0] == "#":
            return self.next()
        return line
