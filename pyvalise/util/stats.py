#!/usr/bin/env python

"""Statistics. """

import logging
from scipy import stats as scistat
import numpy

__author__ = "Damon May"
__copyright__ = "Copyright (c) Damon May"
__license__ = ""
__version__ = ""

log = logging.getLogger(__name__)


def ttestp(vals1, vals2=None):
    """return the 2-sided p-value from a t-test with either 1 or 2 samples"""
    if vals2:
        return scistat.ttest_ind(vals1, vals2)[1]
    else:
        return scistat.ttest_1samp(vals1)[1]


def wilcoxp(vals1, vals2):
    """return the 2-sided p-value from a wilcox test.
    Todo: does this work?"""
    return scistat.ranksums(vals1, vals2)[1]


def get_scipy_distribution(data):
    """create a scipy.stats distribution from a univariate sample.
    data: a numpy.array"""
    kernel = scistat.gaussian_kde(data)

    class rv(scistat.rv_continuous):
        def _rvs(self, *x, **y):
            # don't ask me why it's using self._size 
            # nor why I have to cast to int
            return kernel.resample(int(self._size))

        def _cdf(self, x, *args):
            return kernel.integrate_box_1d(-numpy.Inf, x)

        def _pdf(self, x, *args):
            return kernel.evaluate(x)

    return rv(name='kdedist')


def list_to_array(valslist):
    valsarray = numpy.array(valslist, dtype=numpy.float)
    return valsarray

