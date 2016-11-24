#!/usr/bin/env python
"""
Utilities for writing stuff to PDFs. Depends on reportlab
"""

import logging
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, inch
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle

log = logging.getLogger(__name__)

DEFAULT_ROWHEIGHT = 0.5


def create_pdf(filepath):
    return SimpleDocTemplate(filepath, pagesize=letter)


def make_table(data, columnwidths=None, rowheights=None,
               tableheight=None, tablewidth=6.5):
    """
    :param data: a list of lists. Doesn't need to be strings -- str() will get called on each cell
    :param columnwidths: widths of each column. If not provided, will divide tablewidth evenly
    :param rowheights: heights of each row. If not provided, will divide tableheight evenly
    :param tableheight:
    :param tablewidth:
    :return:
    """
    data_as_strings = []
    if not columnwidths:
        n_columns = len(data[0])
        columnwidths = n_columns * [tablewidth / n_columns]
    if not rowheights:
        n_rows = len(data)
        if tableheight:
            rowheights = n_rows * [tableheight / n_rows]
        else:
            rowheights = n_rows * [DEFAULT_ROWHEIGHT]
    for row in data:
        data_as_strings.append([str(x) for x in row])
    t = Table(data, [x * inch for x in columnwidths],
              [x * inch for x in rowheights])
    t.setStyle(TableStyle([
        ('ALIGN', (1, 1), (-2, -2), 'RIGHT'),
        #                           ('TEXTCOLOR', (1, 1), (-2, -2), colors.red),
        #                           ('VALIGN', (0, 0), (0, -1), 'TOP'),
        #                           ('TEXTCOLOR', (0, 0), (0, -1), colors.blue),
        #                           ('ALIGN', (0, -1), (-1, -1), 'CENTER'),
        #                           ('VALIGN', (0, -1), (-1, -1), 'MIDDLE'),
        #                           ('TEXTCOLOR', (0, -1), (-1, -1), colors.green),
        ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black),
        ('BOX', (0, 0), (-1, -1), 0.25, colors.black),
    ]))
    return t


def build_pdf(pdf_doc, elements):
    pdf_doc.build(elements)

