#!/usr/bin/env python

"""Charting utilities. Consolidated here so that the rest of the codebase
isn't dependent on charting implementation (currently matplotlib)."""

import logging
import pylab as plt
from matplotlib_venn import venn2, venn3


__author__ = "Damon May"
__copyright__ = "Copyright (c) 2012-2014 Damon May"
__license__ = ""
__version__ = ""

log = logging.getLogger(__name__)


def venn(lists, title=None, labels=None):
    """Proportional Venn"""
    if len(lists) not in [2, 3]:
        raise Exception("venn() called with %d lists to compare" % len(lists))
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    sets = []
    if not labels:
        labels = []
    for i in xrange(0, len(lists)):
        sets.append(set(lists[i]))
        labels.append("Set" + str(i + 1))
    if len(lists) == 2:
        intersection_size = len(sets[0].intersection(sets[1]))
        venn2([len(sets[0]) - intersection_size, len(sets[1]) - intersection_size,
               intersection_size], set_labels=(labels[0], labels[1]), ax=ax)
    elif len(lists) == 3:
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


