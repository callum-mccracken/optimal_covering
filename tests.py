"""contains functions for testing each function in each module"""
from constants import dtime
from datetime import datetime


def test_earth():
    print("working on: test_earth")
    import earth_data
    # check that land data works -- we should only need to run this once,
    # but best to make sure it works every time
    earth_data.get_land_mask()
    # UBC's coordinates, and coords on the opposite side of the world
    lon_lats = [[-123.2460, 49.2606], [70.12, -49.20]]
    now = datetime.now()
    for lon, lat in lon_lats:
        result = "\b" if earth_data.is_light(lon, lat, now) else "not"
        input("Manual check: Is it "+ result +" light right now at lon: " + \
              str(lon) + " lat: " + str(lat) + "? Hit enter when done: ")
    
    # download clouds for one month
    earth_data.download_merra_data(dtime.year, dtime.month)
    # get clear polygons for a specific time
    polys = earth_data.get_clear_polys(dtime)
    # not sure how to actually test this, I assume if it runs it's fine
    if polys is not None:
        print("passed test_earth")
    else:
        print("failed test_earth")


def test_utils():
    # if constants import, it's probably fine, not much in that file
    import constants
    # do all the geometry tests
    import geometry
    geometry.test_add_lon_lat()
    geometry.test_perpendicular_vector()
    geometry.test_set_magnitude()
    geometry.test_ortho_conversion()
    geometry.test_arcs_angles()
    geometry.test_xyz_conversion()
    geometry.test_ortho_angles()
    geometry.test_visible()
    geometry.test_polygon_area()
    geometry.test_obs_poly()
    geometry.test_fov_coverage()
    # test timer
    import time
    from timer import timeit
    @timeit
    def test_function():
        time.sleep(1)
    test_function()
    input("If the value of the timer is about 1 second, we passed this test!"
          "\nHit Enter to continue: ")
    print("passed geometry tests!")


def test_plotting():
    import plot_funcs
    from genetics import random_population
    from earth_data import get_clear_polys
    clear_polys = get_clear_polys(dtime)
    population = random_population(clear_polys)
    # this will also plot the clear polys, no need to test those separately
    plot_funcs.plot_member(population[0], 1, clear_polys, dtime)


def test_optimization():
    print("importing...")
    import genetics
    import plot_funcs
    best_points = genetics.do_genetics()
    print("Best spots to point:")
    print(best_points)

if __name__ == "__main__":
    #test_earth()
    #test_utils()
    #test_plotting()
    test_optimization()