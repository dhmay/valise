#!/usr/bin/env python

"""
CSV helpers. Depends on csv
"""

import csv

import logging


log = logging.getLogger(__name__)

def load_csv_rowmap_by_col(csv_file, colname, delimiter='\t'):
    """
    Load rows from a CSV file into maps, one per row. Build a map from values in column
    colname to those maps. Non-unique key values currently handled haphazardly -- only one row
    is returned.
    throw an error if column colname is missing from the file.
    Ignore rows where colname is missing.
    :param csv_file:
    :return: a map from values in column colname to rows
    """
    result = {}
    reader = csv.DictReader(csv_file, delimiter=delimiter)
    for rowmap in reader:
        if colname not in rowmap:
            raise ValueError("File %s does not have column %s" % (csv_file.name, colname))
        if rowmap[colname] is not None:
            result[rowmap[colname]] = rowmap
    return result
