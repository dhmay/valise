#!/usr/bin/env python

"""
For writing xlsx
"""

import xlsxwriter
import logging


logger = logging.getLogger(__name__)


def isfloat(value):
    """
    Utility method to check for floatness of values
    :param value:
    :return:
    """
    try:
        float(value)
        return True
    except ValueError:
        return False


def write_xlsx(fieldnames, rows, outfilename, should_format_numbers=True):
    """
    Read in all the data from a .tsv file, store it column by column, and then write those columns
    to an xlsx file. This is somewhat wasteful, because double-storing the spreadsheet, but that's necessary
    for formatting the numbers, because apparently xlsxwriter doesn't allow you to actually access the
    data stored in a spreadsheet (?!).

    Only tricky thing here is determining which columns have all float values and then casting them,
    if should_format_numbers is True.
    :param fieldnames:
    :param rows:
    :param outfilename:
    :return:
    """
    workbook = xlsxwriter.Workbook(outfilename)
    worksheet = workbook.add_worksheet()
    rowidx = 0
    colidxs_without_allfloats = set()
    # write header
    for i in xrange(0, len(fieldnames)):
        worksheet.write(rowidx, i, fieldnames[i])
    rowidx += 1
    values_bycol = []
    for _ in xrange(0, len(fieldnames)):
        values_bycol.append([])
    logger.debug("reading .tsv file...")
    for row in rows:
        for i in xrange(0, len(fieldnames)):
            strval = row[fieldnames[i]]
            if not isfloat(strval):
                colidxs_without_allfloats.add(i)
            values_bycol[i].append(strval)
        rowidx += 1
    logger.debug("Writing xlsx...")
    for i in xrange(0, len(fieldnames)):
        if i in colidxs_without_allfloats or not should_format_numbers:
            logger.debug("Writing column %d in string format." % i)
            worksheet.write_column(1, i, values_bycol[i])
        else:
            logger.debug("Writing column %d in number format." % i)
            worksheet.write_column(1, i, [float(strval) for strval in values_bycol[i]])
    workbook.close()

