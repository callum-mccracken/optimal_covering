"""a module tobe used as a template when making new optimization methods"""
from datetime import datetime
import numpy as np

import earth_data
import geometry

def clear_poly(dtime):
    clear_polys, big_clear_poly = earth_data.get_clear(dtime)
    return big_clear_poly


def cost(twoD_arr_of_points, big_clear_poly):
    return geometry.fov_coverage(twoD_arr_of_points, big_clear_poly)

def random_points():
    lons = 360 * np.random.rand() - 180
    lats = 180 * np.random.rand() - 90
    points = [(lon, lat) for lon, lat in zip(lons, lats)]
    # doesn't hurt to make sure all the points you pick are visible
    points = [point for point in points if geometry.visible(*point)]
    if len(points) == 0:
        raise ValueError("no points visible")
    # 2d array for cost function evaluation purposes
    return np.array([points])

def optimize():
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




