import os

import matplotlib.pyplot as plt
import cartopy
from cartopy.feature.nightshade import Nightshade

import earth_data
import constants
import geometry
from timer import timeit


def plot_clear(dtime, show=True):
    """
    Plot an orthographic projection of clear area

    dtime = datetime instance

    show = boolean, whether or not to display a plot

    Either shows a plot and returns nothing, or returns the axis on which
    the clear area is plotted, depending on the value of show
    """
    print('plotting clear area')

    # general earth setup
    ax = plt.axes(projection=constants.ortho)
    ax.add_feature(cartopy.feature.OCEAN, zorder=0)
    ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
    ax.add_feature(cartopy.feature.BORDERS, zorder=0, edgecolor='black')
    ax.add_feature(Nightshade(dtime, alpha=0.2))

    # add clear polys
    clear_polys, _ = earth_data.get_clear(dtime)
    print("adding", len(clear_polys), "clear polygons to plot")
    ax.add_geometries(clear_polys, crs=constants.lonlat, facecolor='g')

    # show the plot if desired
    if show:
        plt.show()
    else:
        # return the axis so we can plot more things overtop if we want
        return ax


@timeit
def plot_points(points, title, dtime, show=True):
    """
    Plots a set of observation points.

    points = [(lon, lat), (...), ...], a list of tuples with coordinates
    title = string, we'll name the file title.png

    dtime = datetime instance

    show = boolean, do you want to display the plot or save it as a png?
    """
    # plot clear_polys first
    ax = plot_clear(dtime, show=False)

    # then plot observations on the same axis
    polys = [geometry.obs_poly(lon, lat) for lon, lat in points]
    ax.add_geometries(polys, crs=constants.lonlat, edgecolor='r',
                      facecolor='', alpha=0.8)

    if show:
        plt.show()
    else:
        fname = os.path.join(constants.png_dir, "{}.png".format(title))
        plt.savefig(fname, dpi=300)
        plt.cla()
        plt.clf()
        print('saved image:', fname)


def plot_best_fitnesses(best_fitnesses):
    """Makes a plot of the best fitnesses seen during each iteration
    of the genetic algorithm. Currently just a basic matplotlib plot."""
    plt.plot(best_fitnesses)
    plt.xlabel("Iteration")
    plt.ylabel("Fitness")
    plt.show()
