#!/usr/bin/env python

"""
For writing xlsx
"""

import xlsxwriter
import logging


log = logging.getLogger(__name__)

def write_xlsx(fieldnames, rows, outfilename):
    workbook = xlsxwriter.Workbook(outfilename)
    worksheet = workbook.add_worksheet()
    rowidx = 0
    for i in xrange(0, len(fieldnames)):
        worksheet.write(rowidx, i, fieldnames[i])
    rowidx += 1
    for row in rows:
        for i in xrange(0, len(fieldnames)):
            worksheet.write(rowidx, i, row[fieldnames[i]])
        rowidx += 1
    workbook.close()

