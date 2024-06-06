# install.packages("sf")
# install.packages("terra")
# install.packages("spData")
# install.packages("spDataLarge", repos = "https://nowosad.r-universe.dev")



library(sf)            # classes and functions for vector data
#> Linking to GEOS  for geometry operations,
#> GDAL for reading and writing geographic data files
#> PROJ for representing and transforming projected coordinate reference systems
#> sf_use_s2() is TRUE - interface to Google’s spherical geometry library

library(terra)         # classes and functions for raster data
library(spData)        # load geographic data
library(spDataLarge)   # load larger geographic data

vignette(package = "sf") # see which vignettes are available
vignette("sf1")          # an introduction to the package


class(world)
names(world)
head(world)
str(world$geom)


plot(world)
# by default plots all attributes of df
summary(world$lifeExp)

world_mini= world[1:2,1:3]
world_mini

world_dfr = st_read(system.file("shapes/world.shp", package = "spData"))
# sf read functiin, 
world_tbl = read_sf(system.file("shapes/world.shp", package = "spData"))
# tidyverse read function
class(world_dfr)
#> [1] "sf"         "data.frame"
class(world_tbl)
#> [1] "sf"         "tbl_df"     "tbl"        "data.frame"

plot(world[3:6])
plot(world["pop"])

world_asia = world[world$continent =="Asia",]
asia = st_union(world_asia)

# add layers
plot(world["pop"],reset=F)
plot(asia,add=T,col="red")

plot(world["continent"], reset = FALSE)
cex = sqrt(world$pop) / 10000
world_cents = st_centroid(world, of_largest = TRUE)
plot(st_geometry(world_cents), add = TRUE, cex = cex)

# expand bounding box
india = world[world$name_long == "India", ]
plot(st_geometry(india), expandBB = c(0, 0.2, 0.1, 1), col = "gray", lwd = 3)
plot(st_geometry(world_asia), add = TRUE)

# creating sf object
lnd_point = st_point(c(0.1, 51.5))                 # sfg object (geometry)
lnd_geom = st_sfc(lnd_point, crs = "EPSG:4326")    # sfc object
lnd_attrib = data.frame(                           # data.frame object
  name = "London",
  temperature = 25,
  date = as.Date("2023-06-21")
)
lnd_sf = st_sf(lnd_attrib, geometry = lnd_geom)    # sf object

lnd_sf

# SFG - simple feature GEOMETRY
# sfg POINT

point1 = st_point(c(5, 2))
point2 = st_point(c(1, 3))

st_point(c(5,2))       # 2-dim XY
st_point(c(5,2,3))     # 3-dim XYZ
st_point(c(5,2,3),dim="XYM") # 3-dim XYM
st_point(c(5,2,3,6))        # 4-dim XYZM 

# SFC - simple feautre COLUMNS
# SFC = list of SFG + coordinate reference system
points_sfc = st_sfc(point1, point2)
points_sfc
#> Geometry set for 2 features 
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 1 ymin: 2 xmax: 5 ymax: 3
#> CRS:           NA
#> POINT (5 2)
#> POINT (1 3)

st_geometry_type(points_sfc)
# SFC usually consist of SFGs of one type - and so they inherit geometry type from SFG

st_crs(points_sfc)
#> Coordinate Reference System: NA
# all SFGs in SFC must have the same CRS

points_sfc_wgs = st_sfc(point1, point2, crs = "EPSG:4326")
st_crs(points_sfc_wgs)

# S2 - spherical geometry
sf_use_s2()

# sf_use_s2(FALSE) --> planar/projected geometry

india_buffer_with_s2 = st_buffer(india, 1) # 1 meter
sf_use_s2(FALSE)
#> Spherical geometry (s2) switched off
india_buffer_without_s2 = st_buffer(india, 1) # 1 degree
#> Warning in st_buffer.sfc(st_geometry(x), dist, nQuadSegs, endCapStyle =
#> endCapStyle, : st_buffer does not correctly buffer longitude/latitude data
#> dist is assumed to be in decimal degrees (arc_degrees).

plot(india_buffer_with_s2[,"iso_a2"])
plot(india_buffer_without_s2[,"iso_a2"])

sf_use_s2(TRUE)

# RASTER 
# Raster - continuus grid of cells (pixels)
# Header = CRS + extent of the matrix + origin of the grid
# matrix - one pixel (cell) can hold only single information
# The value can be either continuous or categorical 

# terra

raster_filepath = system.file("raster/srtm.tif", package = "spDataLarge")
my_rast = rast(raster_filepath)
class(my_rast)

my_rast

plot(my_rast)

single_raster_file = system.file("raster/srtm.tif", package = "spDataLarge")
single_rast = rast(raster_filepath)

new_raster = rast(nrows = 6, ncols = 6, 
                  xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5,
                  vals = 1:36)

plot(new_raster)

multi_raster_file = system.file("raster/landsat.tif", package = "spDataLarge")
multi_rast = rast(multi_raster_file)
multi_rast  #dimensons - nrow, ncol, NLYR
nlyr(multi_rast)
plot(multi_rast)

multi_rast1 = multi_rast[["landsat1"]]
multi_rast2 = multi_rast$landsat_2
multi_rast3 = subset(multi_rast, 3)
multi_rast4 = subset(multi_rast, "landsat_4")

multi_rast34 = c(multi_rast3, multi_rast4)
plot(multi_rast34)


## Coordinate Reference System = CRS
# geographic CRS
surface of the Earth is represented eirher by spherical or ellipsoidal surface.
Spherical - sphere with one radius
Elipsoidal - elipsoid with equatorial radius and polar radius. Equatorial > Polar
location is identified by  longitude and latitude (angular distance from point 0)
angular unit of measurment

The datum - which elipsoid to use , and what is the precise relationship 
between the coordinates and location on the earth surface

Geocentric datum (eg WGS84) - center of elipsoid is located in the Earth's' center of gravity
Local datum - elipsoidal surface is shifted to align with the surface of particular location (but misaligned with most of other locations on the Earth surface)

# Projected CRS
map projections, that are converting 3-dimensional geographic CRS models into X and Y values
Based on Cartesiazn coordinates, implicitly flat surface
linear unit of measurement (ie km)
The projection can be performed in many differnt ways, but always with deformation to one of these: area, direction, distance, shape
Projections are often named after the property they preserve ie equidistant projection preserves distance



# Units
luxembourg = world[world$name_long == "Luxembourg", ]
st_area(luxembourg) 
units::set_units(st_area(luxembourg), km^2)

# Exercises

summary(world$continent)
library(tidyverse)
world |> count(continent)

nigeria = world[world$name_long=="Nigeria",]
plot(world["pop"],reset=F)
plot(nigeria["pop"],add=T,col="red")

poland= world[world$name_long=="Poland",]
plot(poland)
plot(st_geometry(poland),expandBB=c(0,0.2,0.4,0.6),col="grey",lwd=3)
plot(st_geometry(world[world$continent=="Europe",]),add=T)

my_raster = rast(nrows = 10, ncols = 10, 
                  xmin = -2, xmax = 2, ymin = -2, ymax = 2,
                  vals = sample(1:10,100,replace = T))
plot(my_raster)

raster_filepath = system.file("raster/nlcd.tif", package = "spDataLarge")
my_rast2 = rast(raster_filepath)

my_rast2
crs(my_rast2)

plot(my_rast2)

