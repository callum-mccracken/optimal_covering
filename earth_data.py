import os
from os.path import join, exists, expanduser
import datetime as dt
from calendar import monthrange

import numpy as np
import netCDF4
import shapely.geometry
from shapely.ops import cascaded_union
from cartopy.feature.nightshade import Nightshade
import fiona
import cartopy.io.shapereader as shpreader
from shapely.prepared import prep

import constants as c
import geometry as g
from timer import timeit


def get_land_mask():
    """
    Creates a landmask for the entire world,
    i.e. a matrix with entries designating places as land or not land.
    """
    if not exists(c.land_mask_path):
        print("Creating landmask")
        # get land object for determining whether or not something is land
        geoms = fiona.open(shpreader.natural_earth(
            resolution='50m', category='physical', name='land'))
        land_geom = shapely.geometry.MultiPolygon(
            [shapely.geometry.shape(geom['geometry']) for geom in geoms])
        land = prep(land_geom)
        # make array of 1 for land 0 for non-land
        landmask = np.zeros([c.lats.size, c.lons.size])
        for i, lat in enumerate(c.lats):
            for j, lon in enumerate(c.lons):
                landmask[i, j] = land.contains(
                    shapely.geometry.Point(lon, lat))
        np.save(c.land_mask_path, landmask)
    return np.load(c.land_mask_path)


def get_clouds(dtime, interp_time=None):
    """
    Return cloud values at a given time.

    args:
        dtime - datetime: a time associated with cloud data
        interp_time - datetime: a time to use to interpolate cloud data
    """
    # get the filename for the date in question
    fname = c.cloud_path_format.format(
        Y=dtime.strftime("%Y"),
        m=dtime.strftime("%m"),
        strf8=dtime.strftime("%Y%m%d"))

    # check if it exists
    if not exists(fname):
        print("Cloud data does not exist... Downloading that month.")
        download_merra_data(dtime.year, dtime.month)

    print("Reading cloud data from MERRA file", fname)

    # we'll save cloud data in a Dataset object
    merra = netCDF4.Dataset(fname, 'r')

    # save MERRA datasets as attributes
    time = merra.variables["time"]
    cloud = merra.variables["CLOUD"]
    cloud.set_auto_maskandscale(False)
    cloud = cloud[...]

    # Get the times corresponding to each time step in the dataset
    t0 = dt.datetime(dtime.year, dtime.month, dtime.day) + dt.timedelta(
        seconds=int(time.begin_time))
    times = [t0 + dt.timedelta(
        seconds=int(time.time_increment)*(i+1))
        for i in range(int(time.shape[0]))]

    # try to avoid invalid values
    cloud[np.isnan(cloud)] = 1
    cloud[cloud > 1] = 1
    cloud[cloud < 0] = 0
    # take the avg value along the vertical axis
    cloud = np.mean(cloud, axis=1)

    # append the first column of lons to the end,
    # so that locations in the range (179.6875 E, 180 E]
    # can get assigned the proper cloud fraction values
    first_cld = np.array([
        np.vstack(cloud[k, :, 0]) for k in range(cloud.shape[0])])

    cloud = np.append(cloud, first_cld, axis=2)

    # now interpolate cloud values
    if interp_time is None:
        interp_time = times[0]
    # return the first timestep's array for low values
    if interp_time <= times[0]:
        return cloud[0, :, :]

    # return the last timestep's array for high values
    elif interp_time > times[-1]:
        return cloud[-1, :, :]

    else:
        # find the time steps that bound the new interp_time
        for k in range(len(times)-1):
            if times[k] < interp_time <= times[k+1]:
                lower = k
                upper = k+1

        # create the weights for the intepolation
        tlower = 0
        tupper = times[upper] - times[lower]
        t = interp_time - times[lower]

        wlower = (t - tupper)/(tlower-tupper)
        wupper = (t - tlower)/(tupper-tlower)

        # interpolate cloud values to the new time
        cld = wlower*cloud[lower, :, :] \
            + wupper*cloud[upper, :, :]

        # set any invalid values to 1 = most cloudy
        cld[np.isnan(cld)] = 1

        return cld


def clear_vertices(lon, lat):
    """calculate vertices of a cloud point, using set lon / lat resolution"""
    lo_lon = g.add_lon(lon, -c.lon_res / 2)
    lo_lat = g.add_lat(lat, -c.lat_res / 2)
    hi_lon = g.add_lon(lon,  c.lon_res / 2)
    hi_lat = g.add_lat(lat,  c.lat_res / 2)
    return [(lo_lon, lo_lat),
            (hi_lon, lo_lat),
            (hi_lon, hi_lat),
            (lo_lon, hi_lat),
            (lo_lon, lo_lat)]


def get_cloud_mask(dtime):
    """returns a 2d array of cloud fraction values"""
    cloud_mask_path = dtime.strftime(c.cloud_mask_format)
    if not exists(cloud_mask_path):
        cloud_fraction = get_clouds(dtime)
        # for now let's say we only care about whether or not it was cloudy,
        # so cast this as int, cloudy or not.
        cloud_mask = (cloud_fraction > c.cloud_thresh).astype(int)
        np.save(cloud_mask_path, cloud_mask)
    return np.load(cloud_mask_path)


def is_light(lons, lats, dtime, night=None):
    """test if lists / arrays of lons, lats are dark"""
    if night is None:
        night = Nightshade(dtime)
    dark_poly = list(night.geometries())[0]
    mesh_lons, mesh_lats = np.meshgrid(lons, lats)
    # check if each pair of points are inside the dark polygon
    result = np.empty_like(mesh_lons, dtype=int)
    for index, lon in np.ndenumerate(mesh_lons):
        lat = mesh_lats[index]
        # convert coordinate systems to the dark_poly's coords
        x, y = night.crs.transform_point(lon, lat, c.lonlat)
        light = not dark_poly.contains(shapely.geometry.Point(x, y))
        result[index] = light
    return result


@timeit
def get_light_mask(dtime):
    """get 2D array of points with 1 = light 0 = dark"""
    light_mask_file = dtime.strftime(c.light_mask_format)
    if not exists(light_mask_file):
        print("Creating light_mask, this might take 30 seconds")
        # this is the polygon that is the "shady" area
        night = Nightshade(dtime)
        # get all points not within the dark_poly
        light_mask = is_light(c.lons, c.lats, dtime, night)
        # save result to file
        np.save(light_mask_file, light_mask)
    return np.load(light_mask_file)


@timeit
def get_clear_polys(dtime):
    # 1 for clear, 0 for cloudy
    cloud_mask = get_cloud_mask(dtime)

    # get land_mask, 1 for land, 0 for water
    land_mask = get_land_mask()

    # check if we'd even be able to observe there (i.e. is each point lit up?)
    # 1 for light, 0 for dark
    light_mask = get_light_mask(dtime)

    # now multiply to get
    # 1 = good, observable region, 0 = bad, not observable
    useful_mask = cloud_mask * land_mask * light_mask

    # make an array to more conveniently store coordinates
    # access coords as [lon, lat] too
    lonlat_coords = np.empty(useful_mask.shape, dtype=tuple)
    for i, lat in enumerate(c.lats):
        for j, lon in enumerate(c.lons):
            lonlat_coords[i, j] = (lon, lat)

    # now create an array of Polygon objects based on the coords and mask
    ones = np.where(useful_mask == 1)
    useful_coords = lonlat_coords[ones]
    polygons = []

    for lon, lat in useful_coords:
        # these polys are nice to deal with since they're arranged in a grid
        # rather than created with a satellite's pointing geometry
        lonlat_vertices = clear_vertices(lon, lat)
        # ensure all vertices are actually visible
        vis = [g.visible(lon, lat) for lon, lat in lonlat_vertices]
        if all(vis):
            polygons.append(shapely.geometry.Polygon(lonlat_vertices))

    # cast as array for quick calculations later
    return np.array(polygons)


def download_merra_data(year, month):
    '''
    Notes:
    This function was created using these instructions from NASA:
    https://disc.gsfc.nasa.gov/data-access# note that
    MERRA data is only published for up to about 2 months ago/

    Inputs:
    year (string): the year you want to download
    month (string): the month you want to download
    No outputs, just saves the MERRA data to wherever you tell it to go.
    '''
    # save cwd so we can get back here at the end
    cwd = os.getcwd()

    # go to home directory
    os.chdir(expanduser("~"))

    # record credentials if they don't exist
    if not exists(".netrc"):
        user = "CallumMcCracken"  # input("Enter earthdata.nasa.gov username: ")
        password = "ThisIsMyV3ryStrongPassword"  # input("Enter password: ")
        os.system("touch .netrc")
        os.system(
            'echo "machine urs.earthdata.nasa.gov login {} '.format(user) +
            'password {}" >> .netrc'.format(password))
        os.system('chmod 777 .netrc')  # so future people can edit
        os.system('touch .urs_cookies')

    print("Downloading MERRA Data for", year, month)
    print("this will take approximately forever")

    # go to save directory
    if not os.path.exists(c.merra_dir):
        os.mkdir(c.merra_dir)
    os.chdir(c.merra_dir)
    # path format for files we want to download
    wget_format = (
        "https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/MERRA2/"
        "M2T3NVASM.5.12.4/{year}/{month:02}/MERRA2_400.tavg3_3d_asm_Nv.{year}"
        "{month:02}{day:02}.nc4.nc?CLOUD[0:7][0:71][0:360][0:575],lat[0:360],"
        "time[0:7],lon[0:575],lev")

    # download each file
    for day in range(1, monthrange(year, month)[1]+1):
        print("Working on day", day)
        os.system(
            "wget "
            "--load-cookies "+join(expanduser("~"), ".urs_cookies")+" "
            "--save-cookies "+join(expanduser("~"), ".urs_cookies")+" "
            "--keep-session-cookies "
            "--cut-dirs 3 "
            "-r -c -nH -nd -np "
            "-A MERRA2_400.tavg3_3d_asm_Nv.{}{:02}{:02}* ".format(
                year, month, day) +
            "--content-disposition " +
            wget_format.format(year=year, month=month, day=day))
        print("Done day", day)

    # rename / move downloaded files
    for fname in os.listdir():
        # for some reason all the files download with extension ".nc4.nc"
        # so rename all .nc4.nc files as just .nc4, then move them
        if fname.endswith("nc4.nc") and fname.contains("{}{:02}".format(year, month)):
            new_filename = fname[:-3]
            new_dir = join("{}".format(year), "{:02}".format(month))
            if not exists(str(year)):
                os.mkdir(str(year))
            if not exists(new_dir):
                os.mkdir(new_dir)
            os.system("mv "+fname+" "+join(new_dir, new_filename))
        # remove wget logs
        if "wget-log" in fname:
            os.remove(fname)

    # at the end return us to the cwd in case other code needs to be run
    os.chdir(cwd)


def get_clear(dtime):
    """really gets (clear_polys, big_clear_poly)

    clear_polys is a list of polygons which define where the clear area is

    big_clear_poly is the union of those

    dtime = datetime instance
    """
    clear_polys = get_clear_polys(dtime)
    big_clear_poly = cascaded_union(clear_polys)
    return clear_polys, big_clear_poly
