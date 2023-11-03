library(sf)
library(terra)
library(spData)
library(spDataLarge)
library(dplyr)

st_crs("EPSG:4326")
vector_filepath <- system.file("shapes/world.gpkg", package = "spData")
new_vector <- read_sf(vector_filepath)
st_crs(new_vector)
plot(world)
st_crs(new_vector)$IsGeographic
st_crs(new_vector)$units_gdal
st_crs(new_vector)$srid
st_crs(new_vector)$proj4string
new_vector <- st_set_crs(new_vector, "EPSG:4326")

raster_filepath <- system.file("raster/srtm.tif", package = "spDataLarge")
my_rast <- rast(raster_filepath)
cat(crs(my_rast))
crs(my_rast)
crs(my_rast) = "EPSG:26912"

london <- data.frame(lon = -0.1, lat = 51.5) |> st_as_sf(coords = c("lon", "lat"))
st_is_longlat(london)
st_crs(london)
london
plot(london)
london_geo <- st_set_crs(london, "EPSG:4326")
st_is_longlat(london_geo)
sf::sf_use_s2(F)
sf::sf_use_s2(T)
class(london)

london_buff_no_crs <- st_buffer(london, dist = 1)
plot(london_buff_no_crs)
london_buff_s2 <- st_buffer(london_geo, dist = 100000) #default value of max_cells = 1000
plot(london_buff_s2)
london_buff_s2_100_cells <- st_buffer(london_geo, dist = 100000, max_cells = 100)
plot(london_buff_s2_100_cells)
sf_use_s2(F)
london_buff_lonlat <- st_buffer(london_geo, dist = 1)
sf::sf_use_s2(T)
plot(london_buff_lonlat)

install.packages("geosphere")
library(geosphere)
geosphere::distGeo(c(-0.1,51.5), c(-1.1,51.5)) # intermeridian distance at London
london_proj <- data.frame(x = 530000, y = 180000) |> st_as_sf(coords = c("x", "y"), crs = "EPSG:27700")
st_crs(london_proj)
st_is_longlat(london_proj)

sf_use_s2(T)
london_buff_projected <- st_buffer(london_proj, 100000)
plot(london_buff_projected)

library(leaflet)
st_distance(london_geo, london_proj)

# For finding the UTM of any location on Earth

lonlat2UTM = function(lonlat){utm = (floor(lonlat[1] + 180 / 6) %% 60) + 1
if (lonlat[2]>0) {
  utm = 32600
} else {
  utm +32700
}
}

lonlat2UTM(c(174.4, -36.9))
lonlat2UTM(st_coordinates(london))

london_2 <- st_transform(london_geo, "EPSG:27700")
st_distance(london_2, london_proj)
st_crs(cycle_hire_osm)
crs_lnd <- st_crs(london_geo)
class(crs_lnd)
names(crs_lnd)
crs_lnd$Name
cat(crs_lnd$wkt)
crs_lnd$epsg
crs_lnd$proj4string

cycle_hire_osm_projected <- st_transform(cycle_hire_osm, "EPSG:27700")
st_crs(cycle_hire_osm_projected)
crs_lnd_new <- st_crs("EPSG:27700")
crs_lnd_new$Name
crs_lnd_new$proj4string
crs_lnd_new$epsg

cat_raster <- rast(system.file("raster/nlcd.tif", package = "spDataLarge"))
crs(cat_raster)
cat(crs(cat_raster))
unique(cat_raster)

cat_raster_wgs84 <- project(cat_raster, "EPSG:4326", method = "near")
plot(cat_raster_wgs84)
unique(cat_raster_wgs84)
nrow(cat_raster)
nrow(cat_raster_wgs84)

con_raster <- rast(system.file("raster/srtm.tif", package = "spDataLarge"))
crs(con_raster)
cat(crs(con_raster))
con_raster_ea <- project(con_raster, "EPSG:32612", method = "bilinear")
plot(con_raster_ea)
crs(con_raster_ea)

zion <- read_sf(system.file("vector/zion.gpkg", package = "spDataLarge"))
zion_centr <- st_centroid(zion)
zion_centr_wgs84 <- st_transform(zion_centr, "EPSG:4326")
st_as_text(st_geometry(zion_centr_wgs84))
cat(st_crs(zion)$wkt)

my_wkt = 'PROJCS["Custom_AEQD",
 GEOGCS["GCS_WGS_1984",
  DATUM["WGS_1984",
   SPHEROID["WGS_1984",6378137.0,298.257223563]],
  PRIMEM["Greenwich",0.0],
  UNIT["Degree",0.0174532925199433]],
 PROJECTION["Azimuthal_Equidistant"],
 PARAMETER["Central_Meridian",-113.0263],
 PARAMETER["Latitude_Of_Origin",37.29818],
 UNIT["Meter",1.0]]'
zion_aeqd <- st_transform(zion, my_wkt)

world_mollweide <- st_transform(world, crs = "+proj=moll")
plot(world_mollweide)
world_wintri <- st_transform(world, crs = "+proj=wintri")
plot(world_wintri)
world_laea2 <- st_transform(world, crs = "+proj=laea +x_0=0 +y_0=0 +lon_0=-74 +lat_0=40")
plot(world_laea2)
