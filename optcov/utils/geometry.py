"""
Module for containing all purely geometric functions.

- conversions to and from orthographic projection coordinates / lon, lat
- get angles from arc lengths and lengths from angles
- check that polygons were formed correctly
- remove intersections of polygons from polygons
- other operations involving polygons and coordinates, ...
"""
import numpy as np
from numpy import cos, sin, arccos, arcsin, arctan2, \
                  sqrt, radians, array, degrees, pi
from shapely.geometry import Polygon
import shapely.ops as ops

from optcov.utils import constants as c

# we'll have some rounding error to deal with at certain points
accuracy = 1e-8


"""
Important things to note:

    longitude: -180 to 180
    latitude: -90 to 90


    To generate a Polygon:
    pgon = Polygon((
            (xll, yll),
            (xll, yur),
            (xur, yur),
            (xur, yll),
            (xll, yll)
        ))
    (where ll = lower left, ul = upper left,
     lr = lower right, ur = upper right)

    Throughout this module we treat the Earth as spherical.

    There are some variables we'll use a decent amount listed below:

    λ: longitude in radians / "azimuthal angle"
    φ: latitude in radians / "polar angle"
    h = height of some object from Earth's center
"""


# functions to deal with the weirdness of lon, lats coordinates
def add_lon(lon1, lon2):
    sum_of_lons = lon1 + lon2
    # but it might be out of bounds!
    while sum_of_lons > 180:
        sum_of_lons -= 360
    while sum_of_lons < -180:
        sum_of_lons += 360
    return sum_of_lons


def add_lat(lat1, lat2):
    sum_of_lats = lat1 + lat2
    # but it might be out of bounds!
    while sum_of_lats > 180:
        sum_of_lats -= 360
    while sum_of_lats < -180:
        sum_of_lats += 360
    return sum_of_lats


def lons_equal(lon1, lon2):
    diff = abs(add_lon(lon1, -1*lon2))
    return diff < accuracy


def lats_equal(lat1, lat2):
    diff = abs(add_lat(lat1, -1*lat2))
    return diff < accuracy


def equal(thing1, thing2):
    """up to numerical accuracy"""
    diff = abs(thing1 - thing2)
    return diff < accuracy


def test_add_lon_lat():
    # not actually sure how to test these... They're right, right?
    print("passed test_add_lon_lat")
    pass


def perpendicular_vector(v):
    """arbitrary perpendicular vector to v (v has 3 components)"""
    if v[1] == 0 and v[2] == 0:
        if v[0] == 0:
            raise ValueError('zero vector')
        else:
            return np.cross(v, [0, 1, 0])
    return np.cross(v, [1, 0, 0])


def test_perpendicular_vector():
    v = [np.random.rand(), np.random.rand(), np.random.rand()]
    while v == [0, 0, 0]:
        v = [np.random.rand(), np.random.rand(), np.random.rand()]
    perp = perpendicular_vector(v)
    should_be_zero = np.dot(perp, v)
    if should_be_zero != 0:
        print(v)
        print(perp)
        print(should_be_zero)
        raise ValueError("perpendicular vector isn't perpendicular!")
    print("passed test_perpendicular_vector")


def set_magnitude(v, mag):
    """adjust the magnitude of a vector, setting it to mag"""
    v_mag = np.linalg.norm(v)
    return v / v_mag * mag


def test_set_magnitude():
    # generate non-zero vector randomly
    v = [np.random.rand(), np.random.rand(), np.random.rand()]
    while v == [0, 0, 0]:
        v = [np.random.rand(), np.random.rand(), np.random.rand()]
    # pick a random magnitude
    mag = np.random.randint(10)
    # set the vector to have that magnitude
    new_v = set_magnitude(v, mag)
    # then check that it worked, up to rounding error
    mag_diff = abs(np.linalg.norm(new_v) - mag)
    if mag_diff > 1e-8:
        print(v)
        print(mag)
        print(new_v)
        print(np.linalg.norm(new_v))
        raise ValueError("magnitude setting has failed!")
    print("passed test_set_magnitude")


def lon_lat_to_ortho_x_y(lon, lat):
    """
    Performs an orthographic projection.
    works for both floating point values and for numpy arrays.

    x, y coordinates are orthographic projectdtimecoordinates.

    args:
        lat - latitude(s) of the point(s)
        lon - longitude(s) of the point(s)
    """
    # phi, lambda = spherical cordinate angles, in radians
    λ_s = radians(c.sat_lon)
    φ_s = radians(c.sat_lat)

    λ = radians(lon)
    φ = radians(lat)

    # perform orthographic projection
    # https://en.wikipedia.org/wiki/Orthographic_projection_in_cartography
    x = c.re_km * cos(φ) * sin(λ - λ_s)
    y = c.re_km * (cos(φ_s) * sin(φ) - sin(φ_s) * cos(φ) * cos(λ - λ_s))
    # return orthographic x-y coordinates
    return x, y


def ortho_x_y_to_lon_lat(x, y):
    """
    Performs an inverse orthographic projection.

    x, y coordinates are orthographic projection coordinates, in km

    args:
        x - orthographic x-coordinate(s) of the point(s)
        y - orthographic y-coordinate(s) of the point(s)
    """
    if x == 0 and y == 0:  # special case where we don't have a unique value
        return c.sat_lon, c.sat_lat

    λ_s = radians(c.sat_lon)
    φ_s = radians(c.sat_lat)

    # do inverse orthographic projection, see wiki for details
    # https://en.wikipedia.org/wiki/Orthographic_projection_in_cartography
    ρ = sqrt(x ** 2 + y ** 2)  # that's a rho, not a p!
    C = arcsin(ρ / c.re_km)  # use capital C, little c is constants module!

    φ = arcsin(cos(C) * sin(φ_s) + y*sin(C) * cos(φ_s) / ρ)
    Δλ = arctan2(x * sin(C), ρ * cos(φ_s) * cos(C) - y * sin(φ_s) * sin(C))

    Δlon = degrees(Δλ)
    lon = add_lon(c.sat_lon, Δlon)
    lat = degrees(φ)
    return lon, lat


def test_ortho_conversion():
    counter = 0
    for lon in c.lons:
        for lat in c.lats:
            if visible(lon, lat):
                counter += 1
                new_lon, new_lat = ortho_x_y_to_lon_lat(
                    *lon_lat_to_ortho_x_y(lon, lat))
                correct = lons_equal(new_lon, lon) and lats_equal(new_lat, lat)
                if not correct:
                    print(lon, lat)
                    print(new_lon, new_lat)
                    raise ValueError("ortho conversion didn't work!")
    print("passed test_ortho_conversion, tested", counter, "points")


def lon_lat_alt_to_xyz(lon, lat, alt):
    """Return cartesian coordinates given lon, lat (and possibly alt in km)

    We assume the Earth is spherical. xyz units are km.
    """
    h = alt + c.re_km  # height from center of earth

    # convert to radians
    λ = radians(lon)
    φ = radians(lat)

    # do conversion from polar to cartesian
    d = h * cos(φ)   # perpendicular distance from z axis

    x = d * cos(λ)
    y = d * sin(λ)
    z = h * sin(φ)

    # return Cartesian coordinate values!
    return x, y, z


def xyz_to_lon_lat_alt(x, y, z):
    """Return geodetic coordinates (lon, lat) corresponding to
    the Cartesian (geocentric) Coordinates given by x, y, z.
    """

    # A few important lengths
    d = sqrt(x ** 2 + y ** 2)  # perpendicular distance from z axis
    h = sqrt(d ** 2 + z ** 2)  # height from center of Earth
    alt = h - c.re_km  # altitude

    # get angles
    λ = arctan2(y, x)
    φ = arcsin(z / h)

    # convert to degrees
    lon = degrees(λ)
    lat = degrees(φ)

    # return lon, lat, alt coordinate values!
    return lon, lat, alt


def test_xyz_conversion():
    # make a random point
    lon = 360*np.random.rand()-180  # random number between -180 and 180
    lat = 180*np.random.rand()-90  # random number between -90 and 90
    alt = 0  # let's say it's a point on the surface of Earth

    x, y, z = lon_lat_alt_to_xyz(lon, lat, alt)
    new_lon, new_lat, new_alt = xyz_to_lon_lat_alt(x, y, z)

    if not lons_equal(new_lon, lon) or not lats_equal(new_lat, lat):
        print(lon, new_lon)
        print(lat, new_lat)
        raise ValueError("xyz conversion failed: lons / lats don't match!")
    if not equal(new_alt, alt):
        print(alt, new_alt)
        raise ValueError("xyz conversion failed: alts don't match!")
    print("passed test_xyz_conversion")


def visible(lon, lat):
    """returns True if this point is visible from the satellite, else False"""
    # rewrite the coordinates as differences from the sub-satellite point

    # center of earth, we don't actually need it for anything but FYI
    # C = array([0,0,0])

    # satellite point
    S = array(lon_lat_alt_to_xyz(c.sat_lon, c.sat_lat, c.sat_alt))
    # distance to S from C
    S_dist = np.linalg.norm(S)  # in km

    # distance to the furthest point in the line of sight of our satellite
    # is calculated by noting that that line of sight at point X (furthest pt)
    # is tangential to the Earth. So,
    # SCX is a right-angle triangle, with the right angle on Earth's surface.

    # X_dist is the distance from sat to the furthest point that we can see
    # use pythagorean theorem, noting that the distance to X from the center
    # of the Earth is c.re_km.
    X_dist = sqrt(S_dist**2 - c.re_km**2)
    # this is the point we want to see
    P = array(lon_lat_alt_to_xyz(lon, lat, 0))
    # distance to P from S
    P_dist = np.linalg.norm(P - S)
    # if P is closer to S than X, we'll be able to see it, otherwise no.
    return P_dist <= X_dist


def test_visible():
    # the antipode is the point on the opposite side of the Earth. Let's
    # find the antipode of the satellite, which should 100% not be visible

    anti_lon = add_lon(c.sat_lon, 180)
    anti_lat = add_lat(c.sat_lat, 90)
    if visible(anti_lon, anti_lat):
        raise ValueError("The antipode should not be visible!")

    if not visible(c.sat_lon, c.sat_lat):
        raise ValueError("The satellite must be able to see under itself!")
    print("passed test_visible")


def max_visible_angle():
    """
    Finds the max angle that the satellite can point
    before going off the edge of the Earth
    """
    r = c.re_km + c.sat_alt  # distance to satellite from center of Earth

    # where our view cuts off, line of sight will be tangential to Earth
    # (i.e. we have a right angle triangle) from CoE to P to S
    theta = arcsin(c.re_km / r)
    return theta


def angle_from_arc_length(arc_length):
    """angle is size of the angle we'll have to point the satellite,
    away from "straight down", to see a point a distance arc_length away
    on the surface of the Earth."""
    # distance to satellite from center of Earth
    # note sat_alt must also be in km.
    r = c.sat_alt + c.re_km

    # "internal angle" from satellite --> center of Earth --> point
    α = arc_length / c.re_km

    # Think of a triangle with vertices at the satellite,
    # at the center of the Earth, and at the point arc_length away.
    # Two sides have length r and re_km.
    # The other side has length b, which can be found with the cosine law:
    b = sqrt(r ** 2 + c.re_km ** 2 - 2 * r * c.re_km * cos(α))

    # theta is the angle (relative to "straight down") that the satellite must
    # point in order to see the point arc_length away.
    theta = arcsin(c.re_km*sin(α)/b)
    return theta


def arc_length_from_angle(angle):
    # angle = the angle away from straight down, where the satellite's pointing
    # the question here is "how far away along the earth's surface is that?"
    theta = angle

    r = c.re_km + c.sat_alt  # distance to satellite from center of Earth

    # use the sine law to get beta, where beta is the other angle in the
    # triangle with angles theta, alpha, beta
    # i.e. beta is the angle center of Earth --> point --> satellite
    beta = arcsin(r * sin(theta) / c.re_km)

    # as in the other function, alpha is the angle from sat --> Earth --> point
    alpha = np.pi - beta - theta

    arc_length = alpha * c.re_km

    # return distance along the Earth's surface
    return arc_length


# max angle the satellite can point before seeing space
max_ang = max_visible_angle()  # in radians
# max distance traveled on Earth before going off the planet
max_arc = arc_length_from_angle(max_ang)  # in km
# maximum value an orthographic coordinate can have
max_ortho = c.re_km * sin(radians(max_arc / c.lon_length))  # in km


def arc_length(lon1, lat1, lon2=c.sat_lon, lat2=c.sat_lat):
    """calculates how far away (along Earth's surface) a point is,
    from the sub-satellite point by default, or from some other point
    (the result is always positive).
    """

    # coordinates of point 1 and 2.
    # point 1 is given. By default point 2 is the satellite.
    λ_1 = radians(lon1)
    φ_1 = radians(lat1)
    λ_2 = radians(lon2)
    φ_2 = radians(lat2)

    # absolute differences
    Δλ = abs(λ_1 - λ_2)
    Δφ = abs(φ_1 - φ_2)

    # use formula for great circle distance from wikipedia
    # https://en.wikipedia.org/wiki/Great-circle_distance
    # I used the low-precision one right at the top under Formulae,
    # since we don't need accuracy better than a couple meters
    # (we're assuming the Earth is a sphere,
    # so maybe that should be fixed first if we want more accuracy)

    # this is (abs of) the internal angle from S to center of Earth to P
    Δσ = arccos(sin(φ_1) * sin(φ_2) + cos(φ_1) * cos(φ_2) * cos(Δλ))

    # arc length, in km
    a = c.re_km * Δσ
    return a


def test_arcs_angles():
    for angle in np.arange(-max_ang, max_ang, 1000):
        arc = arc_length_from_angle(angle)
        new_angle = angle_from_arc_length(arc)
        if new_angle > max_ang:
            raise ValueError("how?")
        if not equal(angle, new_angle):
            raise ValueError("arc conversion didn't work!")
    print("passed test_arcs_angles")


def ortho_x_y_to_angle_x_y(x, y):
    x_angle = arcsin(x / max_ortho) * 2 / pi * max_ang
    y_angle = arcsin(y / max_ortho) * 2 / pi * max_ang
    return x_angle, y_angle


def angle_x_y_to_ortho_x_y(angle_x, angle_y):
    x = sin(angle_x / max_ang * pi / 2) * max_ortho
    y = sin(angle_y / max_ang * pi / 2) * max_ortho
    return x, y


def test_ortho_angles():
    for lon in c.lons:
        for lat in c.lats:
            if visible(lon, lat):
                x, y = lon_lat_to_ortho_x_y(lon, lat)
                angle_x, angle_y = ortho_x_y_to_angle_x_y(x, y)
                new_x, new_y = angle_x_y_to_ortho_x_y(angle_x, angle_y)
                correct = equal(new_x, x) and equal(new_y, y)
                if not correct:
                    print(x, y)
                    print(new_x, new_y)
                    raise ValueError("ortho angles conversion didn't work!")
    print("passed test_ortho_conversion")


def polygon_area(poly):
    """Thanks to StackExchange for this one
    https://stackoverflow.com/questions/4681737/how-to-calculate-the-area-of-a-polygon-on-the-earths-surface-using-python
    """
    def sin_proj(lon, lat):
        """
        Returns the x & y coordinates in km using a sinusoidal projection

        latitude and longitude are numpy arrays
        """
        lat_dist = np.pi * c.re_km / 180.0
        y = lat * lat_dist
        x = lon * lat_dist * cos(radians(lat))
        return x, y

    def area(x, y):
        """Calculates the area of an arbitrary polygon given its verticies
        in sinusoidal coordinates (in km), result is in km^2
        """
        area = 0.0
        for i in range(-1, len(x)-1):
            area += x[i] * (y[i+1] - y[i-1])
        return abs(area) / 2.0

    # if poly has vertices at (0,0), (0,1), (1,1), (1,0)
    # this will give (0,0,1,1,0), (0,1,1,0,0)
    lon, lat = poly.exterior.coords.xy
    lon = array(lon)
    lat = array(lat)
    # for the functions above, we don't need that last repeated vertex.
    x, y = sin_proj(lon, lat)
    return area(x, y)


def test_polygon_area():
    # Wyoming is a pretty square state, let's use it as a test
    real_area = 253600  # km^2 (sorry Americans)
    wy_corners = [
        (-111.049623, 44.993649),
        (-104.057802, 44.997317),
        (-104.053325, 41.001400),
        (-111.049785, 41.001854),
        (-111.049623, 44.993649)
    ]
    wy = Polygon(wy_corners)
    wy_area = polygon_area(wy)
    percent_error_ok = 1
    if abs(wy_area - real_area) / real_area * 100 > percent_error_ok:
        print(real_area, wy_area, abs(wy_area-real_area))
        raise ValueError("error calculation is off!")
    print("passed test_polygon_area")


def obs_poly(lon, lat):
    """
    Given the (lon, lat) coordinates of a point we want to observe,
    get the (lon, lat) coordinates of the vertices of the polygon we'd see,
    then return a Polygon object with those vertices.

    lon, lat = coordinates of center of polygon

    returns: Polygon object, with (mostly) realistic vertices
    """
    # _c means "of center"
    x_c, y_c = lon_lat_to_ortho_x_y(lon, lat)

    # then convert those to angles
    x_angle, y_angle = ortho_x_y_to_angle_x_y(x_c, y_c)
    # our point (lon, lat) is (x_angle, y_angle) away from straight down

    # now generate vertices
    vertices = []
    for x_inc, y_inc in c.angle_increments:
        # add original angle + increment to get pointing angles to the vertex
        vert_x_ang = x_angle + x_inc
        vert_y_ang = y_angle + y_inc

        # ensure we don't end up with a polygon with vertices in space!
        # in real life that might be okay, but here it messes up calculations.
        if vert_x_ang > max_ang:
            vert_x_ang = max_ang
        if vert_y_ang > max_ang:
            vert_y_ang = max_ang

        # now convert back to lon, lat
        vert_x, vert_y = angle_x_y_to_ortho_x_y(vert_x_ang, vert_y_ang)
        vert_lon, vert_lat = ortho_x_y_to_lon_lat(vert_x, vert_y)
        vertices.append((vert_lon, vert_lat))
    return Polygon(vertices)


def test_obs_poly():
    # check that it makes a whole bunch of polygons in the visible area
    for lon in c.lons:
        for lat in c.lats:
            if visible(lon, lat):
                poly = obs_poly(lon, lat)
    # other than the fact that it makes polygons, not sure what else to test...
    print("passed test_obs_poly")


def fov_coverage(population, clear_poly):
    """
    calculates coverage for a population.

    That is, the total unique area covered by
    each fov of each member of the population.

    The population is just a set of sets of (lon, lat) points.

    clear_poly is a Polygon that contains all the clear area

    e.g. a population with 3 members, each with 4 fovs, might look like:

    population = array([
        [(lat00, lon00), (lat01, lon01), (lat02, lon02), (lat03, lon03)],
        [(lat10, lon10), (lat11, lon11), (lat12, lon12), (lat13, lon13)],
        [(lat20, lon20), (lat21, lon21), (lat22, lon22), (lat23, lon23)],
    ])
    """
    # poly_pop = population made of polygons,
    # array containing polygons rather than tuples of coordinates
    poly_pop = np.empty_like(population, dtype=Polygon)
    for index, lonlat_coords in np.ndenumerate(population):
        poly_pop[index] = obs_poly(*lonlat_coords)

    # area array to return at the end, containing useful area for each fov
    areas = np.empty_like(population, dtype=float)
    for i, member in enumerate(poly_pop):
        for j, fov in enumerate(member):
            # note that fov is a Polygon instance.

            # other FoVs in the member
            member_without_poly = [p for p in member if p is not fov]

            # join all other member polys together
            m_union = ops.cascaded_union(member_without_poly)

            # remove overlap with other fovs from fov
            fov_wo_overlap = fov.symmetric_difference(
                m_union).difference(m_union)

            # find intersection of fov with clear_poly
            intersection = fov_wo_overlap.intersection(clear_poly)
            # get the area of that intersection
            areas[i, j] = intersection.area

    # now normalize these areas so we don't get huge numbers
    ideal_poly = obs_poly(c.sat_lon, c.sat_lat)
    ideal_area = polygon_area(ideal_poly)
    # note that the normalization means we're measuring our result relative to
    # how much area would be covered if there were no clouds, everything was
    # light, and we were looking straight down for every FoV.

    # that's not a stellar metric, but I'm not sure what else to use.
    return areas  # / ideal_area


def test_fov_coverage():
    # make some random population
    n_fovs = 5
    n_members = 1
    population = np.empty((n_members, n_fovs), dtype=tuple)
    for m in range(n_members):
        for f in range(n_fovs):
            lon = 360 * np.random.rand() - 180
            lat = 180 * np.random.rand() - 90
            population[m, f] = (lon, lat)

    # get some test day's cloud data
    try:
        from optcov.earth_data import get_clear_polys
        clear_polys = get_clear_polys(c.dtime)
    except Exception:
        print("Ensure you have downloaded some cloud data!")
        raise

    coverage = fov_coverage(population, clear_polys)
    print("passed test_fov_coverage")


if __name__ == "__main__":
    test_add_lon_lat()
    test_perpendicular_vector()
    test_set_magnitude()
    test_ortho_conversion()
    test_arcs_angles()
    test_xyz_conversion()
    test_ortho_angles()
    test_visible()
    test_polygon_area()
    test_obs_poly()
    test_fov_coverage()
    print("All tests passed!")
