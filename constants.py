"""
Contains a bunch of constant parameters.

Satellite position, radius of earth, paths, etc.
"""
import numpy as np
from datetime import datetime
import cartopy.crs as ccrs  # this darn thing takes like 10 seconds to import
from os.path import realpath, join, dirname, exists
from os import mkdir

# the time will be constant over a run, I assume
dtime = datetime(2015, 4, 4, hour=22)

"""SATELLITE PARAMETERS"""

sat_lat = 90  # latitude of satellite, in degrees
sat_lon = -90  # longitude of satellite, in degrees
sat_alt = 10e3  # altitude of satellite, in km

# reference altitude, in km, used for calculating pixel distortion
ref_alt = 10e3

# size of one pixel when looking down and sat_alt == ref_alt, in radians
pixel_size = np.radians(1.0)

# number of pixels in each FoV, in x and y  (i.e. l/r, u/d)directions
xpixels = 1
ypixels = 1

# angular length of a FoV in each direction (when looking down)
x_increment = xpixels * pixel_size * sat_alt / ref_alt
y_increment = ypixels * pixel_size * sat_alt / ref_alt

# these are the adjustments we'll make with the satellite around a point
angle_increments = [(-x_increment/2, -y_increment/2),
                    (x_increment/2, -y_increment/2),
                    (x_increment/2, y_increment/2),
                    (-x_increment/2, y_increment/2),
                    (-x_increment/2, -y_increment/2)]

"""CLOUD PARAMETERS"""

# if MERRA data has cloud fraction greater than this, we say it's cloudy
cloud_thresh = 0.1


# lat and lon resolutions of cloud data and landmask data -- they must be equal
lon_res = 0.625
lat_res = 0.5

"""EARTH PARAMETERS"""

# radius of Earth, using WGS84 value
re_km = 6378.137

# lons, lats: clear area centers, the +1e-8 here is just so we get the endpoint
lons = np.arange(-180., 180. + 1e-8, lon_res)
lats = np.arange(-90, 90 + 1e-8, lat_res)

# arc length for a change of 1 degree:
lat_length = 2 * np.pi * re_km / 180  # 180 degrees of lat total
lon_length = 2 * np.pi * re_km / 360  # 360 degrees of lon total

# different projections to use at various points
ortho = ccrs.Orthographic(
    central_longitude=sat_lon,
    central_latitude=sat_lat)
lonlat = ortho.as_geodetic()

"""PATHS"""


# directory for cloud_mask, land_mask, and light_mask files
earth_dir = realpath(join(dirname(__file__), "earth_data"))

# directory for MERRA files
merra_dir = join(earth_dir, 'MERRA')


# path to save images
png_dir = realpath(join(dirname(__file__), 'images_and_output'))

# ensure everything we need now exists, we'll make output dir later
for d in [merra_dir, earth_dir]:
    if not exists(d):
        mkdir(d)

# formats of files we'll use for reading / writing
# note: strf8 means string format 8 characters long, like 20150312
# and %Y, %m, ... are datetime string format flags

cloud_path_format = join(
    merra_dir, '{Y}', '{m}', 'MERRA2_400.tavg3_3d_asm_Nv.{strf8}.nc4')

cloud_mask_format = join(earth_dir, 'cloud_mask_%Y_%m_%d_%H.npy')

land_mask_path = join(earth_dir, "land_mask_0.5_lat_0.625_lon.npy")

light_mask_format = join(earth_dir, "light_mask_%Y_%m_%d_%H.npy")
