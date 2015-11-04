#!/usr/bin/env python

"""
Cover methods for building charts with mayavi.
These don't play nicely with matplotlib. Maybe I can write them to the same pdf?
"""

from mayavi import mlab
import traits
import logging
from scipy import stats
import numpy as np


def show_density_plot(x_values, y_values, z_values, xlabel=None, ylabel=None, zlabel=None):
    """
    Cribbed from stackoverflow:
    http://stackoverflow.com/questions/25286811/how-to-plot-a-3d-density-map-in-python-with-matplotlib
    :param x_values:
    :param y_values:
    :param z_values:
    :param out_image_file: extension determines behavior. If .png, writes an image file. If .pdf, appends to pdf?
    :return:
    """
    figure = mlab.figure('DensityPlot')
    figure.scene.disable_render = True

    # calculate density
    xyz = np.vstack([x_values, y_values, z_values])
    kde = stats.gaussian_kde(xyz)

    # Evaluate kde on a grid
    xmin = min(x_values)
    ymin = min(y_values)
    zmin = min(z_values)
    xmax = max(x_values)
    ymax = max(y_values)
    zmax = max(z_values)
    xi, yi, zi = np.mgrid[xmin:xmax:30j, ymin:ymax:30j, zmin:zmax:30j]
    coords = np.vstack([item.ravel() for item in [xi, yi, zi]])
    density = kde(coords).reshape(xi.shape)

    # Plot scatter with mayavi
    grid = mlab.pipeline.scalar_field(xi, yi, zi, density)
    dens_min = density.min()
    dens_max = density.max()
    mlab.pipeline.volume(grid, vmin=dens_min, vmax=dens_max)# + .5*(dens_max-dens_min))

    mlab.axes()
    if xlabel:
        mlab.xlabel(xlabel)
    if ylabel:
        mlab.ylabel(ylabel)
    if zlabel:
        mlab.zlabel(zlabel)
    mlab.show()
