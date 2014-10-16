#!/usr/bin/env python

"""Charting utilities. Consolidated here so that the rest of the codebase
isn't dependent on charting implementation (currently matplotlib)."""

import logging
from matplotlib.backends.backend_pdf import PdfPages
import pylab as plt
from statsmodels.graphics import gofplots
from pyptides.util import stats as pyptide_stats
from matplotlib_venn import venn2, venn3
from scipy.stats import gaussian_kde
from numpy import arange


__author__ = "Damon May"
__copyright__ = "Copyright (c) 2012-2014 Damon May"
__license__ = ""
__version__ = ""

log = logging.getLogger(__name__)

# default color sequence for lines and bars. If need more colors, add more colors
COLORS = ['#0000ff',  # blue
          '#ff0000',  # red
          '#00ff00',  # green
          '#ff00ff',  # purple
          '#888888']  # grey

DEFAULT_HIST_BINS = 40

# hack for printing numeric values on piechart segments
piechart_valuesum = 0


def write_pdf(figures, pdf_file):
    """write an iterable of figures to a PDF"""
    pdf = PdfPages(pdf_file)
    for plot in figures:
        plot.savefig(pdf, format='pdf')
    pdf.close()


def hist(values, title=None, bins=DEFAULT_HIST_BINS, color=None):
    """trivial histogram.
    Stupidly, the default range doesn't go from 0."""
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    if not color:
        color = 'black'

    myrange = [min(values), max(values)]
    rangesize = myrange[1] - myrange[0]
    myrange = [myrange[0] - rangesize / 20, myrange[1] + rangesize / 20]
    ax.hist(values,  # histtype = "stepfilled",
            color=color, range=myrange, bins=bins)
    if title:
        ax.set_title(title)
    return figure


def multiboxplot(valueses, title=None, labels=None):
    """multiple boxplots side by side"""
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    ax.boxplot(valueses)
    tick_locs = list()
    for i in xrange(0, len(valueses)):
        tick_locs.append(i + 1)

    if labels:
        plt.xticks(tick_locs, labels)

    if title:
        ax.set_title(title)
    return figure


def multiviolin(valueses, title=None, labels=None):
    """multiple violin plots side by side"""
    figure = plt.figure()

    ax = figure.add_subplot(1, 1, 1)

    violin_plot(ax, valueses, range(len(valueses)))

    tick_locs = list()
    for i in xrange(0, len(valueses)):
        tick_locs.append(i + 1)

    if labels:
        plt.xticks(tick_locs, labels)

    if title:
        ax.set_title(title)
    return figure


def violin_plot(ax, data, pos, bp=True, color='y'):
    """
    create violin plots on an axis.
    Borrowed from here:
    http://pyinsci.blogspot.com/2009/09/violin-plot-with-matplotlib.html
    """
    dist = max(pos) - min(pos)
    w = min(0.15 * max(dist, 1.0), 0.5)
    for d, p in zip(data, pos):
        k = gaussian_kde(d)  # calculates the kernel density
        m = k.dataset.min()  # lower bound of violin
        M = k.dataset.max()  # upper bound of violin
        x = arange(m, M, (M - m) / 100.)  # support for violin
        v = k.evaluate(x)  # violin profile (density curve)
        v = v / v.max() * w  # scaling the violin to the available space
        ax.fill_betweenx(x, p, v + p, facecolor=color, alpha=0.3)
        ax.fill_betweenx(x, p, -v + p, facecolor=color, alpha=0.3)
    if bp:
        ax.boxplot(data, notch=1, positions=pos, vert=1)


def multihist(valueses, title=None, bins=DEFAULT_HIST_BINS, colors=None):
    """histogram of multiple datasets.
    valueses: a list of lists of values"""
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    mymin = float("inf")
    mymax = 0
    for values in valueses:
        print(len(values))
        mymax = max(mymax, max(values))
        mymin = min(mymin, min(values))
    myrange = [mymin, mymax]
    rangesize = mymax - mymin
    myrange = [myrange[0] - rangesize / 20, myrange[1] + rangesize / 20]
    if not colors:
        colors = COLORS[0:len(valueses)]
    colorind = -1
    for values in valueses:
        colorind += 1
        ax.hist(values,  # histtype = "stepfilled",
                color=colors[colorind], alpha=0.5, range=myrange, bins=bins)

    if title:
        ax.set_title(title)
    return figure


def multibar(valueses, labels, title=None, colors=None):
    """barchart of multiple datasets.
    valueses: a list of lists of values. Should have same cardinalities"""
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    if not colors:
        colors = COLORS[0:len(valueses)]
    ind = -1
    for values in valueses:
        ind += 1
        xvals = [0] * len(values)
        for xind in xrange(0, len(xvals)):
            xvals[xind] = 2 * xind * len(valueses) + ind
        ax.bar(xvals, values, color=colors[ind])
    tick_xs = [0] * len(labels)
    for ind in xrange(0, len(tick_xs)):
        tick_xs[ind] = 2 * ind * len(valueses)
    plt.xticks(tick_xs, labels)
    if title:
        ax.set_title(title)
    return figure


def line_plot(x_values, y_values, title=None):
    """trivial line plot"""
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    ax.plot(x_values, y_values)
    if title:
        ax.set_title(title)
    return figure


def multiline(x_valueses, y_valueses, labels, title=None, colors=None,
              legend_loc='upper center'):
    """line plot with multiple lines and a legend"""
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    if not colors:
        colors = COLORS[0:len(x_valueses)]
    for i in xrange(0, len(x_valueses)):
        x_values = x_valueses[i]
        y_values = y_valueses[i]
        color = colors[i]
        ax.plot(x_values, y_values, color=color, label=labels[i])
    if title:
        ax.set_title(title)
    # Now add the legend with some customizations.
    legend = ax.legend(loc=legend_loc, shadow=True)
    # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
    frame = legend.get_frame()
    frame.set_facecolor('0.90')

    # Set the fontsize
    for label in legend.get_texts():
        label.set_fontsize('large')

    for label in legend.get_lines():
        label.set_linewidth(1.5)  # the legend line width
    return figure


def scatterplot(x_values, y_values, title=None):
    """trivial scatterplot"""
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    ax.scatter(x_values, y_values)
    if title:
        ax.set_title(title)
    return figure


def qqplot_2samples(vals1, vals2, title=None,
                    xlabel='sample 1', ylabel='sample 2'):
    """qqplot of the distributions of two samples. The first is taken as
    the reference and the second is plotted against it"""
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    gofplots.qqplot(pyptide_stats.list_to_array(vals2),
                    pyptide_stats.get_scipy_distribution(pyptide_stats.list_to_array(vals1)),
                    line='45', ax=ax)
    if title:
        ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return figure


# def png_to_plot(png_filename):
#     """convert a png file to a Matplotlib plot"""
#     figure = plt.figure()
#     figure.add_subplot(1,1,1)
#     img = PIL.Image.open(png_filename)
#     
#     plt.imshow(numpy.asarray(img), shape = img.size)
#     return figure

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


def pie(values, title=None, labels=None,
        colors=('b', 'g', 'r', 'c', 'm', 'y', 'k', 'w')):
    """trivial pie chart"""
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    ax.pie(values, autopct='%.2f%%', labels=labels, colors=colors)
    if title:
        ax.set_title(title)
    return figure