"""
A module to be used as a template when trying new optimization methods.
"""
from datetime import datetime
import numpy as np

import earth_data
import geometry

def clear_poly(dtime):
    """
    Get a Polygon which represents all clear area on Earth.
    """
    clear_polys, big_clear_poly = earth_data.get_clear(dtime)
    return big_clear_poly


def cost(points, big_clear_poly):
    """
    Calculate the cost of a set of points,
    at a time where the clear area on Earth is represented by big_clear_poly
    """
    return -1 * geometry.fast_coverage(points, big_clear_poly)

def random_points():
    """
    Make some random set of points.
    """
    N = 100
    lons = 360 * np.random.rand(N) - 180
    lats = 180 * np.random.rand(N) - 90
    points = [(lon, lat) for lon, lat in zip(list(lons), list(lats))]
    # doesn't hurt to make sure all the points you pick are visible
    points = [(lon, lat) for lon, lat in points if geometry.visible(lon, lat)]
    if len(points) == 0:
        raise ValueError("no points visible")

    # points must be a 2d array
    points_arr = np.empty(len(points), dtype=tuple)
    for idx, _ in np.ndenumerate(points_arr):
        points_arr[idx] = points[idx[0]]
    points_arr = np.array([points_arr])
    return points_arr

def optimize():
    """
    Perform some optimization operation.
    """
    # pick a time
    dtime = datetime(2015, 5, 1)

    # pick some points
    points = random_points()

    # get cloud info
    big_clear_poly = clear_poly(dtime)

    # calculate cost
    cost_value = cost(points, big_clear_poly)

    # do something to the points...
    # optimize somehow...

    # return the optimized points
    return points

if __name__ == "__main__":
    optimize()
