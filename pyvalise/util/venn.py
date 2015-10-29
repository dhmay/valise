#!/usr/bin/env python

"""Charting utilities. Consolidated here so that the rest of the codebase
isn't dependent on charting implementation (currently matplotlib)."""

import logging
import locale
import matplotlib
# so that I can run without X
matplotlib.use('agg')
import pylab as plt
from matplotlib_venn import venn2, venn3
logger = logging.getLogger(__name__)


__author__ = "Damon May"
__copyright__ = "Copyright (c) 2012-2014 Damon May"
__license__ = ""
__version__ = ""

log = logging.getLogger(__name__)


def format_int(int_to_format):
    locale.setlocale(locale.LC_ALL, 'en_US')
    return locale.format("%d", int_to_format, grouping=True)


def proportional_venn_from_counts_2d(n_total1, n_total2, n_intersect, title=None, labels=None):
    """
    Create a proportional Venn diagram from counts of two sets' total values and intersection
    :param n_total1:
    :param n_total2:
    :param n_intersect:
    :return:
    """
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    if not labels:
        labels = [ "Set 1", "Set 2"]

    n_1only = n_total1 - n_intersect
    n_2only = n_total2 - n_intersect
    c = venn2([n_1only, n_2only, n_intersect], set_labels=(labels[0], labels[1]), ax=ax)
    # hardcode better colors than default
    c.get_patch_by_id('10').set_color('red')
    c.get_label_by_id('10').set_text(format_int(n_1only))
    c.get_patch_by_id('01').set_color('blue')
    c.get_label_by_id('01').set_text(format_int(n_2only))
    c.get_patch_by_id('10').set_edgecolor('none')
    c.get_patch_by_id('01').set_edgecolor('none')
    c.get_patch_by_id('10').set_alpha(0.4)
    c.get_patch_by_id('01').set_alpha(0.4)
    c.get_patch_by_id('11').set_color('magenta')
    c.get_patch_by_id('11').set_edgecolor('none')
    c.get_label_by_id('11').set_text(format_int(n_intersect))
    c.get_patch_by_id('11').set_alpha(0.4)
    if title:
        ax.set_title(title)
    return figure


def proportional_venn_from_lists(lists, title=None, labels=None):
    """Proportional Venn from two or three lists of items"""
    if len(lists) not in [2, 3]:
        raise Exception("venn() called with %d lists to compare" % len(lists))
    sets = []
    if not labels:
        labels = []
    for i in xrange(0, len(lists)):
        sets.append(set(lists[i]))
        labels.append("Set" + str(i + 1))
    figure = None
    if len(lists) == 2:
        intersection_size = len(sets[0].intersection(sets[1]))
        logger.debug("Intersection size: %d. Size 1: %d. Size2 2: %d" % (intersection_size, len(sets[0]), len(sets[1])))
        figure = proportional_venn_from_counts_2d(len(sets[0]), len(sets[1]),
                                                  intersection_size, labels=(labels[0], labels[1]))
    elif len(lists) == 3:
        figure = plt.figure()
        ax = figure.add_subplot(1, 1, 1)
        n_inter_01 = len(sets[0].intersection(sets[1]))
        n_inter_12 = len(sets[1].intersection(sets[2]))
        n_inter_02 = len(sets[0].intersection(sets[2]))
        n_inter_012 = len(sets[0].intersection(sets[1]).intersection(sets[2]))

        #order: (100, 010, 110, 001, 101, 011, 111)       
        subset_sizes = (len(sets[0]) - n_inter_01 - n_inter_02 + n_inter_012,
                        len(sets[1]) - n_inter_01 - n_inter_12 + n_inter_012,
                        n_inter_01 - n_inter_012,
                        len(sets[2]) - n_inter_12 - n_inter_02 + n_inter_012,
                        n_inter_02 - n_inter_012,
                        n_inter_12 - n_inter_012,
                        n_inter_012 )
        venn3(subsets=subset_sizes,
              set_labels=(labels[0], labels[1], labels[2]), ax=ax)
        if title:
            ax.set_title(title)

    return figure


