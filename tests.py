"""contains functions for testing each function in each module"""
from optcov.utils.constants import dtime

# TODO: finish implimenting this


def test_land():
    from optcov import earth_data
    # we should only need to run this once,
    # but best to make sure it works every time
    earth_data.generate_landmask()
    
    # longitude, latitude of a hunk of desert -- this is definitely land
    lon = -119.16
    lat = 40.45
    print_arr = [lon, ",", lat, "is"]
    if not earth_data.is_land(lon, lat):
        print_arr.append("not")
    print_arr.append("land")
    print(*print_arr)

    # same test for ocean -- definitely not land
    lon = -160.38
    lat = 34.08
    print_arr = [lon, ",", lat, "is"]
    if not earth_data.is_land(lon, lat):
        print_arr.append("not")
    print_arr.append("land")
    print(*print_arr)


def test_cloud():
    from optcov import earth_data
    # download clouds for one month
    earth_data.download_merra_data(dtime.hour, dtime.minute)
    # get clear polygons for a specific time
    polys = earth_data.get_clear_polys(dtime)


def test_light():
    from datetime import datetime
    from optcov import earth_data
    # UBC's coordinates, opposite side of the world
    lon_lats = [[-123.2460, 49.2606], [70.12, -49.20]]
    dtime = datetime.now()
    for lon, lat in lon_lats:
        result = "" if earth_data.is_light(lon, lat, dtime) else "not"
        input("Manual check:",
              "Is it ", result, "light right now at lon:",
              lon, "lat:", lat, "? Hit enter when done: ")


def test_earth_data():
    test_land()
    test_cloud()
    test_light()
    # (note that we don't bother testing nightshade.py
    # since we pretty much copied that one from cartopy)


def test_constants():
    # if it imports, it's probably fine, nothing much happens in that file
    from optcov.utils import constants


def test_geometry():
    from optcov.utils import geometry
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


def test_timer():
    import time
    from optcov.utils.timer import timeit
    @timeit
    def test_function():
        time.sleep(1)
    test_function()
    print("if the value of the timer is about 1 second, we passed this test!")


def test_utils():
    test_constants()
    test_geometry()
    test_timer()


def test_plotting():
    from optcov import plot_funcs
    from optcov.genetics import random_population
    from optcov.earth_data import get_clear_polys
    clear_polys = get_clear_polys(dtime)
    population = random_population(clear_polys)
    # this will also plot the clear polys, no need to test those separately
    plot_funcs.plot_member(population[0], 1, clear_polys, dtime)


def test_optimization():
    from optcov import genetics
    genetics.perform_genetics(dtime)


if __name__ == "__main__":
    # test_earth_data()
    # test_utils()
    # test_plotting()
    test_optimization()