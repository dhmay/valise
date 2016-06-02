#!/usr/bin/env python

"""Charting utilities. Consolidated here so that the rest of the codebase
isn't dependent on charting implementation (currently matplotlib)."""

import matplotlib
# so that I can run without X
matplotlib.use('Agg')
import logging
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from statsmodels.graphics import gofplots
import statsmodels.api as sm
from pyvalise.util import stats as pyptide_stats
from scipy.stats import gaussian_kde
from numpy import arange
from matplotlib import colors
import numpy as np



COLORMAP_REDBLUE = plt.get_cmap('bwr')

DEFAULT_LOG_BASE = 10
ALPHA_FOR_MULTISCATTER = 0.85

DEFAULT_POINTSIZE = 1

__author__ = "Damon May"
__copyright__ = "Copyright (c) 2012-2014 Damon May"
__license__ = ""
__version__ = ""

logger = logging.getLogger(__name__)

# default color sequence for lines and bars. If need more colors, add more colors
BASE_COLORS = ['#0000ff',  # blue
               '#ff0000',  # red
               '#00ff00',  # green
               '#ff00ff',  # purple
               '#888888',  # grey
               '#FFFF00',  # brown?
               '#222222',  # dark grey
               ]
COLORS = []
for i in xrange(0, 20):
    COLORS.extend(BASE_COLORS)

COLORMAP_BASECOLORS = colors.ListedColormap(BASE_COLORS)


DEFAULT_HIST_BINS = 40

# hack for printing numeric values on piechart segments
piechart_valuesum = 0


def write_pdf(figures, pdf_file):
    """write an iterable of figures to a PDF"""
    with PdfPages(pdf_file) as pdf:
        for plot in figures:
            plot.savefig(pdf, format='pdf')


def hist(values, title=None, bins=DEFAULT_HIST_BINS, color=None,
         should_logx=False, should_logy=False, log_base=DEFAULT_LOG_BASE):
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
    if should_logy:
        plt.yscale('log', basey=log_base)
    if should_logx:
        plt.xscale('log', basex=log_base)
    return figure

def multiboxplot(valueses, title=None, labels=None):
    """multiple boxplots side by side"""
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    ax.boxplot(valueses)
    tick_locs = list()
    for i in xrange(0, len(valueses)):
        tick_locs.append(i)

    if labels:
        plt.xticks(tick_locs, labels)

    if title:
        ax.set_title(title)
    return figure

def multiviolin_fromxy(xvals, yvals, title=None, minforplot=2):
    """Turn xy-pair data into violin plots"""
    xval_yvals_map = {}
    for i in xrange(0, len(xvals)):
        if not xvals[i] in xval_yvals_map:
            xval_yvals_map[xvals[i]] = []
        xval_yvals_map[xvals[i]].append(yvals[i])
    valueses = []
    labels = []

    for xval in sorted(xval_yvals_map.keys()):
        if len(xval_yvals_map[xval]) >= minforplot:
            labels.append(xval)
            valueses.append(xval_yvals_map[xval])
    return multiviolin(valueses, title, labels)


def multiviolin(valueses, title=None, labels=None):
    """multiple violin plots side by side"""
    figure = plt.figure()

    ax = figure.add_subplot(1, 1, 1)

    violin_plot(ax, valueses, range(len(valueses)))

    tick_locs = list()
    for i in xrange(0, len(valueses)):
        tick_locs.append(i)

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


def multihist(valueses, title=None, bins=DEFAULT_HIST_BINS, colors=None,
              legend_labels=None, legend_on_chart=False, histtype='bar',
              should_normalize=False):
    """histogram of multiple datasets.
    valueses: a list of lists of values"""
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    mymin = float("inf")
    mymax = 0
    for values in valueses:
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
        weights = None
        if should_normalize:
            weights = np.ones_like(values)/float(len(values))
        ax.hist(values, color=colors[colorind], alpha=0.5, range=myrange, bins=bins, histtype=histtype,
                weights=weights)
    if legend_labels:
        add_legend_to_chart(ax, legend_on_chart=legend_on_chart, labels=legend_labels)
    if title:
        ax.set_title(title)
    return figure


def multihist_skinnybar(valueses, title=None, bins=DEFAULT_HIST_BINS, colors=None,
                        legend_labels=None, legend_on_chart=False, normed=False):
    """histogram of multiple datasets.
    valueses: a list of lists of values"""
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    if not colors:
        colors = COLORS[0:len(valueses)]

    ax.hist(valueses, bins=bins, color=colors, normed=normed)
    if legend_labels:
        add_legend_to_chart(ax, legend_on_chart=legend_on_chart, labels=legend_labels)
    if title:
        ax.set_title(title)
    return figure


def bar(values, labels, title=None, colors=None, rotate_labels=False):
    return multibar([values], labels, title=title, colors=colors, legend_labels=None, rotate_labels=rotate_labels)


def multibar(valueses, labels, title='', colors=None,
             legend_labels=None, legend_on_chart=True, rotate_labels=False):
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
    locs, ticklabels = plt.xticks(tick_xs, labels)
    if rotate_labels:
        logger.debug("Rotating tick labels 90 degrees.")
        plt.setp(ticklabels, rotation=90)
        # also, shrink the chart vertically to make room
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.3, box.width, box.height * 0.6])
    if legend_labels:
        add_legend_to_chart(ax, legend_on_chart=legend_on_chart, labels=legend_labels, rotate_labels=rotate_labels)
    ax.set_title(title)
    return figure


def line_plot(x_values, y_values, title=None, lowess=False,
              xlabel=None, ylabel=None,
              should_logx=False, should_logy=False, log_base=DEFAULT_LOG_BASE):
    """trivial line plot"""
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    ax.plot(x_values, y_values)
    if lowess:
        lowess = sm.nonparametric.lowess(y_values, x_values, frac=0.1)
        ax.plot(lowess[:, 0], lowess[:, 1])
    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if should_logy:
        plt.yscale('log', basey=log_base)
    if should_logx:
        plt.xscale('log', basex=log_base)
    return figure


def multiline(x_valueses, y_valueses, labels=None, title=None, colors=None,
              linestyles=None, legend_on_chart=False, xlabel=None, ylabel=None,
              should_logx=False, should_logy=False, log_base=DEFAULT_LOG_BASE,
              diff_yaxis_scales=False):
    """

    :param x_valueses:
    :param y_valueses:
    :param labels:
    :param title:
    :param colors:
    :param linestyles:
    :param legend_on_chart:  legend on the chart, rather than off to right
    :param xlabel:
    :param ylabel:
    :param should_logx:
    :param should_logy:
    :param log_base:
    :param diff_yaxis_scales: use different y-axis scales. Only valid for len(y_valueses)==2
    :return:
    """
    if diff_yaxis_scales and len(y_valueses) != 2:
        raise Exception("multiline: diff_yaxis_scales only meaningful if two sets of values.")
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    axes = [ax]
    if not colors:
        colors = COLORS[0:len(x_valueses)]
    for i in xrange(0, len(x_valueses)):
        x_values = x_valueses[i]
        y_values = y_valueses[i]
        color = colors[i]
        if diff_yaxis_scales and i == 1:
            ax = ax.twinx()
            axes.append(ax)
        if labels:
            if linestyles:
                ax.plot(x_values, y_values, color=color, label=labels[i], linestyle=linestyles[i])
            else:
                ax.plot(x_values, y_values, color=color, label=labels[i])
        else:
            if linestyles:
                ax.plot(x_values, y_values, color=color, linestyle=linestyles[i])
            else:
                ax.plot(x_values, y_values, color=color)
    if title:
        ax.set_title(title)
    # Now add the legend with some customizations.
    if labels:
        for ax in axes:
            add_legend_to_chart(ax, legend_on_chart=legend_on_chart)
    if should_logy:
        plt.yscale('log', basey=log_base)
    if should_logx:
        plt.xscale('log', basex=log_base)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    for ax in axes:
        ax.set_aspect(1./ax.get_data_ratio())
    return figure





def multiscatter(x_valueses, y_valueses, title=None,
                 xlabel='', ylabel='', pointsize=DEFAULT_POINTSIZE, labels=None,
                 should_logx=False, should_logy=False, log_base=DEFAULT_LOG_BASE,
                 legend_on_chart=False, colors=COLORS, rotate_labels=False,
                 draw_1to1=False):
    """
    Scatterplot multiple sets of values in different colors
    :param x_valueses:
    :param y_valueses:
    :param title:
    :param xlabel:
    :param ylabel:
    :param pointsize:
    :param colors:
    :param cmap:
    :param show_colorbar:
    :param should_logx:
    :param should_logy:
    :param log_base:
    :return:
    """
    assert(len(x_valueses) == len(y_valueses))

    for i in xrange(0, len(x_valueses)):
       assert(len(x_valueses[i]) == len(y_valueses[i]))
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    if title:
        ax.set_title(title)
    for i in xrange(0, len(x_valueses)):
        logger.debug("multiscatter set %d, color %s" % (i, colors[i]))
        ax.scatter(x_valueses[i], y_valueses[i], s=pointsize, facecolors=colors[i], edgecolors='none',
                   alpha=ALPHA_FOR_MULTISCATTER)
    if labels:
        add_legend_to_chart(ax, legend_on_chart=legend_on_chart, labels=labels, rotate_labels=rotate_labels)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if should_logy:
        plt.yscale('log', basey=log_base)
    if should_logx:
        plt.xscale('log', basex=log_base)
    if draw_1to1:
        lims = [
            np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
            np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
        ]
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
    ax.set_aspect(1./ax.get_data_ratio())
    return figure


def hexbin(x_values, y_values, title=None,
           xlabel='', ylabel='', gridsize=100, should_log_color=False,
           cmap=None, show_colorbar=True,
           should_logx=False, should_logy=False, log_base=DEFAULT_LOG_BASE):
    """
    hexbin plot
    :param x_values:
    :param y_values:
    :param title:
    :param xlabel:
    :param ylabel:
    :param gridsize:
    :param should_log_color: log-transform the color scale?
    :param cmap:
    :param show_colorbar:
    :return:
    """
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    bins = None
    if should_log_color:
        bins = 'log'
    myhexbin = ax.hexbin(x_values, y_values, gridsize=gridsize, bins=bins,
                         cmap=cmap)
    if title:
        ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if should_logy:
        plt.yscale('log', basey=log_base)
    if should_logx:
        plt.xscale('log', basex=log_base)
    if show_colorbar:
        plt.colorbar(myhexbin)
    return figure


def scatterplot(x_values, y_values, title=None, lowess=False,
                xlabel='', ylabel='', pointsize=1, draw_1to1 = False,
                colors=None, cmap=None, show_colorbar=False,
                should_logx=False, should_logy=False, log_base=DEFAULT_LOG_BASE,
                alpha=0.5):
    """
    scatter plot
    :param x_values:
    :param y_values:
    :param title:
    :param lowess:
    :param xlabel:
    :param ylabel:
    :param pointsize:
    :param draw_1to1:
    :param colors:
    :param cmap:
    :param show_colorbar:
    :return:
    """
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    if colors is not None:
        myscatter = ax.scatter(x_values, y_values, s=pointsize, c=colors, cmap=cmap, edgecolors='none', alpha=alpha)
    else:
        myscatter = ax.scatter(x_values, y_values, s=pointsize, cmap=cmap, edgecolors='none', alpha=alpha)
    if draw_1to1:
        lims = [
            np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
            np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
        ]
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
    if title:
        ax.set_title(title)
    if lowess:
        lowess = sm.nonparametric.lowess(y_values, x_values, frac=0.1)
        ax.plot(lowess[:, 0], lowess[:, 1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if should_logy:
        plt.yscale('log', basey=log_base)
    if should_logx:
        plt.xscale('log', basex=log_base)
    if show_colorbar:
        plt.colorbar(myscatter)
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


def pie(values, title=None, labels=None,
        colors=('b', 'g', 'r', 'c', 'm', 'y', 'k', 'w')):
    """trivial pie chart"""
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    ax.pie(values, autopct='%.2f%%', labels=labels, colors=colors)
    if title:
        ax.set_title(title)
    return figure


def heatmap(values_ndarray, xtick_positions=None, xlabels=None,
            ytick_positions=None, ylabels=None, colormap='gray',
            xlabel=None, ylabel=None, title=None, show_colorbar=True,
            width_proportion=1.0, height_proportion=1.0):
    """heatmap. You want axis labels, you provide 'em."""
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    cax = ax.pcolor(values_ndarray, cmap=colormap, vmin=values_ndarray.min(), vmax=values_ndarray.max())
    if xtick_positions:
        assert(xlabels is not None)
        ax.set_xticks(xtick_positions)
        ax.set_xticklabels(xlabels)
    if ytick_positions:
        assert(ylabels is not None)
        ax.set_yticks(ytick_positions)
        ax.set_yticklabels(ylabels)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if show_colorbar:
        figure.colorbar(cax)
    # scale the heatmap's size
    box = ax.get_position()
    ax.set_position([box.x0 + box.width * (1.0 - width_proportion), box.y0 + box.height * (1.0 - height_proportion),
                     box.width * width_proportion, box.height * height_proportion])

    return figure


def surface(x_values, y_values, z_values, title=None,
            xlabel=None, ylabel=None,
            should_logx=False, should_logy=False,
            log_base=DEFAULT_LOG_BASE):
    """3D surface plot"""
    figure = plt.figure()
    ax = figure.gca(projection='3d')
    ax.plot_surface(x_values, y_values, z_values)
    if should_logx:
        plt.xscale('log', basex=log_base)
    if should_logy:
        plt.yscale('log', basey=log_base)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    return figure


def scatter3d(x_values, y_values, z_values, title=None,
                xlabel='', ylabel='', zlabel='', pointsize=1, alpha=0.5,
                colors=None, should_logx=False, should_logy=False, cmap=None,
              log_base=DEFAULT_LOG_BASE, show_colorbar=False):
    """3D scatterplot"""
    figure = plt.figure()
    # ax = Axes3D(figure)
    ax = figure.add_subplot(111, projection='3d')
    myscatter = None
    if colors:
        myscatter = ax.scatter(x_values, y_values, z_values, s=pointsize, c=colors, alpha=alpha, cmap=cmap, edgecolors='none')
    else:
        myscatter = ax.scatter(x_values, y_values, z_values, s=pointsize, alpha=alpha, cmap=cmap, edgecolors='none')

    if should_logx:
        plt.xscale('log', basex=log_base)
    if should_logy:
        plt.yscale('log', basey=log_base)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if zlabel:
        ax.set_zlabel(ylabel)
    if title:
        ax.set_title(title)
    if show_colorbar:
        plt.colorbar(myscatter)
    return figure


def add_legend_to_chart(ax, legend_on_chart=False, labels=None, rotate_labels=False):
    """
    add a legend to the chart. By default, overlapping the chart
    :param ax:
    :param legend_off_chart: move the legend entirely off.
    :return:
    """
    logger.debug("add_legend_to_chart")
    if legend_on_chart:
        logger.debug("Printing legend on chart")
        if labels:
            legend = ax.legend(labels, bbox_to_anchor=(1.1, 1.05), shadow=True, borderaxespad=0.)
        else:
            legend = ax.legend(bbox_to_anchor=(1.1, 1.05), shadow=True, borderaxespad=0.)
    else:
        logger.debug("Printing legend off chart")
        # shrink x axis by 40%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.25, box.width * 0.5, box.height * 0.5])
        if labels:
            legend = ax.legend(labels, loc='center left', bbox_to_anchor=(1, 0.5))
        else:
            legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
    frame = legend.get_frame()
    frame.set_facecolor('0.90')

    if rotate_labels:
        locs, ticklabels = plt.xticks()
        logger.debug("Rotating tick labels 90 degrees.")
        plt.setp(ticklabels, rotation=90)

    for label in legend.get_lines():
        label.set_linewidth(1.5)  # the legend line width


def big_hist_line(big_hist_data, xlabel='', ylabel=None, title='', plot_proportion=False):
    """

    :param big_hist_data:
    :param xlabel:
    :param ylabel:
    :param title:
    :param plot_proportion:
    :return:
    """
    return big_hist_multiline([big_hist_data], xlabel=xlabel, ylabel=ylabel, title=title,
                              plot_proportion=plot_proportion)


def big_hist_multiline(big_hist_datas, xlabel='', ylabel=None, title='', labels=None, plot_proportion=False):
    """
    Plot a multi-line plot from BigHistData data objects.
    Requires min_shown_value and max_shown_value are same for all the inputs
    :param big_hist_datas:
    :param xlabel:
    :param ylabel:
    :param title:
    :param labels:
    :param plot_proportion
    :return:
    """
    assert(len(set([big_hist_data.min_shown_value for big_hist_data in big_hist_datas])) == 1)
    assert(len(set([big_hist_data.max_shown_value for big_hist_data in big_hist_datas])) == 1)
    xvalues = big_hist_datas[0].generate_xvals()
    xvalueses = []
    yvalueses = []
    if not ylabel:
        if plot_proportion:
            ylabel = 'proportion'
        else:
            ylabel = 'count'
    for big_hist_data in big_hist_datas:
        xvalueses.append(xvalues)
        yvalues = big_hist_data.countdata
        if plot_proportion:
            yvalues = big_hist_data.generate_proportion_yvals()
        yvalueses.append(yvalues)
    return multiline(xvalueses, yvalueses, title=title, xlabel=xlabel, ylabel=ylabel, labels=labels)


def big_hist_multihist_proportion(big_hist_datas, title='', labels=None,
                                  bins=None):
    """
    Plot a multi-line plot from BigHistData data objects.
    Requires min_shown_value and max_shown_value are same for all the inputs
    :param big_hist_datas:
    :param xlabel:
    :param ylabel:
    :param title:
    :param labels:
    :param plot_proportion
    :return:
    """
    assert(len(set([big_hist_data.min_shown_value for big_hist_data in big_hist_datas])) == 1)
    assert(len(set([big_hist_data.max_shown_value for big_hist_data in big_hist_datas])) == 1)
    values_forhist = big_hist_datas[0].generate_xvals()
    xvalueses = []
    for big_hist_data in big_hist_datas:
        yvals = big_hist_data.generate_proportion_yvals()
        xvalues = []
        for i in xrange(0, len(yvals)):
            n_thishist_thisval = int(1000 * yvals[i])
            xvalues.extend([values_forhist[i]] * n_thishist_thisval)
        xvalueses.append(xvalues)
    bins = DEFAULT_HIST_BINS
    if len(values_forhist) < bins:
        bins = len(values_forhist) + 5
    return multihist_skinnybar(xvalueses, title=title, bins=bins, legend_labels=labels, normed=True)


class BigHistData:
    """
    A data structure for storing count data for big datasets, to be plotted as a line or multiline plot
    """
    def __init__(self, max_shown_value, min_shown_value=0):
        assert(max_shown_value > min_shown_value)
        self.min_shown_value = min_shown_value
        self.max_shown_value = max_shown_value
        self.countdata = [0] * (max_shown_value - min_shown_value + 1)
        self.min_real_value = float('inf')
        self.max_real_value = float('-inf')

    def add_value(self, value):
        cropped_value = max(self.min_shown_value, min(value, self.max_shown_value))
        self.countdata[cropped_value - self.min_shown_value] += 1
        self.min_real_value = min(value, self.min_real_value)
        self.max_real_value = max(value, self.max_real_value)

    def generate_xvals(self):
        return xrange(self.min_shown_value, self.max_shown_value+1)

    def generate_proportion_yvals(self):
        value_sum = sum(self.countdata)
        proportion_yvals = self.countdata
        if value_sum > 0:
            proportion_yvals = [float(y) / value_sum for y in self.countdata]
        return proportion_yvals

    def print_min_max_real_values(self):
        print("min=%f, max=%f" % (self.min_real_value, self.max_real_value))

    def __str__(self):
        result = "BigHistData:\nmin=%f, max=%f\n" % (self.min_real_value, self.max_real_value)
        for i in xrange(0, len(self.countdata)):
            if i > 0:
                result = result + ","
            result = result + str(self.countdata[i])
        result = result + '\n'
        return result
