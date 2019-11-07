import matplotlib.pyplot as plt
from os.path import join, realpath, dirname
import cartopy
from cartopy.feature.nightshade import Nightshade

import earth_data
import constants as c
import geometry as g
from timer import timeit


def plot_clear(dtime, clear_polys=None, show=True):
    """
    Show an orthographic projection of clear area

    dtime = datetime instance
    """
    print('plotting clear area')

    # general earth setup
    ax = plt.axes(projection=c.ortho)
    ax.add_feature(cartopy.feature.OCEAN, zorder=0)
    ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
    ax.add_feature(cartopy.feature.BORDERS, zorder=0, edgecolor='black')
    ax.add_feature(Nightshade(dtime, alpha=0.2))

    # add clear polys
    if clear_polys is None:
        clear_polys = earth_data.get_clear_polys(dtime)
    print("adding", len(clear_polys), "clear polygons to plot")
    ax.add_geometries(clear_polys, crs=c.lonlat, facecolor='g')

    # show the plot if desired
    if show:
        plt.show()

    # return the axis so we can plot more things overtop if we want
    return ax


@timeit
def plot_member(member, generation, clear_polys, dtime, show=True):
    """plots one set of FOVs"""
    # plot clear_polys first
    ax = plot_clear(dtime, clear_polys=clear_polys, show=False)

    # then plot member on the same axis
    member_polys = [g.obs_poly(lon, lat) for lon, lat in member]
    ax.add_geometries(
        member_polys, crs=c.lonlat, edgecolor='r', facecolor='', alpha=0.8)

    if show:
        plt.show()
    else:
        fname = join(c.png_dir, "gen_"+str(generation)+".png")
        plt.savefig(fname, dpi=300)
        plt.cla()
        plt.clf()
        print('saved image:', fname)


def plot_best_fitnesses(best_fitnesses):
    plt.plot(best_fitnesses)
    plt.xlabel("Iteration")
    plt.ylabel("Fitness")
    plt.show()
