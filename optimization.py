"""a module tobe used as a template when making new optimization methods"""
from datetime import datetime
import numpy as np

import earth_data
import geometry

def clear_poly(dtime):
    clear_polys, big_clear_poly = earth_data.get_clear(dtime)
    return big_clear_poly


def cost(points, big_clear_poly):
    return geometry.fov_coverage(points, big_clear_poly)

def random_points():
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
    # pick a time
    dtime = datetime(2015, 5, 1)

    # pick some points
    points = random_points()

    # get cloud info
    big_clear_poly = clear_poly(dtime)

    # calculate cost
    try:
        cost_value = cost(points, big_clear_poly)
    except Exception as e:
        print("Looks like this crashed, sorry...",
              "If you're feeling ambitious, look inside geometry.py",
              "and let me know if you figure out why that happened!",
              "But for now, just try running a different set of points!")

    # do something to the points...
    # optimize somehow...

    # return the optimized points
    return points

if __name__ == "__main__":
    optimize()




