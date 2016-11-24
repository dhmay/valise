#!/usr/bin/env python
"""
Utilities for writing stuff to PDFs. Depends on reportlab
"""

import logging
import cStringIO
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, inch
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Flowable
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT
from reportlab.lib.utils import ImageReader

log = logging.getLogger(__name__)

DEFAULT_ROWHEIGHT = 0.5


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


class PdfImage(Flowable):

    def __init__(self, img_data, width_inches=5, height_inches=4):
        self.img_width = width_inches * inch
        self.img_height = height_inches * inch
        self.img_data = img_data

    def wrap(self, width, height):
        return self.img_width, self.img_height

    def drawOn(self, canv, x, y, _sW=0):
        if _sW > 0 and hasattr(self, 'hAlign'):
            a = self.hAlign
            if a in ('CENTER', 'CENTRE', TA_CENTER):
                x += 0.5*_sW
            elif a in ('RIGHT', TA_RIGHT):
                x += _sW
            elif a not in ('LEFT', TA_LEFT):
                raise ValueError("Bad hAlign value " + str(a))
        canv.saveState()
        img = self.img_data
        # for later, after I get pdfrw working
#        if isinstance(img, PdfDict):
#            xscale = self.img_width / img.BBox[2]
#            yscale = self.img_height / img.BBox[3]
#            canv.translate(x, y)
#            canv.scale(xscale, yscale)
#            canv.doForm(makerl(canv, img))
#        else:
        canv.drawImage(img, x, y, self.img_width, self.img_height)
        canv.restoreState()


def make_element_from_figure(figure, width_inches=7, height_inches=5):
    # for pdfrw
    #fig.savefig(imgdata, format='pdf' if use_pdfrw else 'png')
    imgdata = cStringIO.StringIO()
    figure.savefig(imgdata, format='png')
    imgdata.seek(0)
    # for pdfr2
    #reader = form_xo_reader if use_pdfrw else ImageReader
    reader = ImageReader
    image = reader(imgdata)
    img = PdfImage(image, width_inches=width_inches, height_inches=height_inches)
    return img


def build_pdf(elements, filepath):
    pdf_doc = SimpleDocTemplate(filepath, pagesize=letter)
    pdf_doc.build(elements)

