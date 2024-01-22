library(sf)
library(terra)
library(dplyr)
library(spData)

download.file(url = "https://hs.pangaea.de/Maps/PeRL/PeRL_permafrost_landscapes.zip",
              destfile = "PeRL_permafrost_landscapes.zip", 
              mode = "wb")
unzip("PeRL_permafrost_landscapes.zip")
canada_perma_land = read_sf("PeRL_permafrost_landscapes/canada_perma_land.shp")
plot(canada_perma_land)
install.packages("FedData")
library(FedData)
install.packages("geodata")
library(geodata)
install.packages("osmdata")
library(osmdata)
install.packages("osmextract")
library(osmextract)
install.packages("rnaturalearth")

library(rnaturalearth)
india <- ne_countries(country = "India")
class(india)
plot(india)
geodata::gadm("India", level = 0, path = tempdir())
india_sf <- st_as_sf(india)

install.packages("rnoaa")
library(rnoaa)
install.packages("giscoR")
install.packages("GSODR")
library(GSODR)

library(geodata)
worldclim_prep <- worldclim_global("prec", res = 10, path = tempdir())

library(osmdata)
parks = opq(bbox = "leeds uk") |> add_osm_feature(key = "leisure", value = "park") |> osmdata_sf()

world2 <- spData::world
world3 <- read_sf(system.file("shapes/world.gpkg", package = "spData"))

install.packages("tidygeocoder")
library(tidygeocoder)
geo_df <- data.frame(address = "54 Frith St, London W1D 4SJ, UK")
geo_df <- geocode(geo_df, address, method = "osm")
geo_df
geo_sf <- st_as_sf(geo_df, coords = c("long", "lat"), crs = "EPSG:4326")

install.packages("httr")
library(httr)
base_url <- "http://www.fao.org"
endpoint <- "/figis/geoserver/wfs"
q = list(request = "GetCapabilities")
res = GET(url = modify_url(base_url, path = endpoint), query = q)
res$url
browseURL(res$url)
txt = content(res, "text")
xml = xml2::read_xml(txt)
xml

qf = list(request = "GetFeature", typeName = "area:FAO_AREAS")
file = tempfile(fileext = ".gml")
GET(url = base_url, path = endpoint, query = qf, write_disk(file))
fao_areas = read_sf(file)

install.packages("ows4R")
install.packages("datasets")

install.packages("tmap")
tmap_mode("view") # sets mapping mode to interactive
tm_shape("sf_object_name") +
tm_polygons("column_names_which_you_want_to_map")
install.packages("spdep")
install.packages("rgeoda")

sf_drivers <- st_drivers()
head(sf_drivers, n = 3)
summary(sf_drivers[-c(1:2)])

f <- system.file("shapes/world.gpkg", package = "spData")
world <- read_sf(f, quiet = T)
tanzania <- read_sf(f, query = 'SELECT * FROM world WHERE name_long = "Tanzania"')
tanzania_buf <- st_buffer(tanzania, 50000)
tanzania_buf_geom <- st_geometry(tanzania_buf)
class(tanzania_buf_geom)
tanzania_buf_wkt <- st_as_text(tanzania_buf_geom)
tanzania_neigh <- read_sf(f, wkt_filter = tanzania_buf_wkt)
plot(tanzania_buf_geom)
plot(tanzania_neigh)

cycle_hire_txt <- system.file("misc/cycle_hire_xy.csv", package = "spData")
cycle_hire_xy <- read_sf(cycle_hire_txt, options = c("X_POSSIBLE_NAMES=X", "Y_POSSIBLE_NAMES=Y"))

world_txt = system.file("misc/world_wkt.csv", package = "spData")
world_wkt = read_sf(world_txt, options = "GEOM_POSSIBLE_NAMES=WKT")
# the same as
world_wkt2 = st_read(world_txt, options = "GEOM_POSSIBLE_NAMES=WKT", 
                     quiet = TRUE, stringsAsFactors = FALSE, as_tibble = TRUE)
class(world_wkt2)

u = "https://developers.google.com/kml/documentation/KML_Samples.kml"
download.file(u, "KML_Samples.kml")
st_layers("KML_Samples.kml")

kml = read_sf("KML_Samples.kml", layer = "Placemarks")
raster_filepath <- system.file("raster/srtm.tif", package = "spDataLarge")
single_layer <- rast(raster_filepath)
multilayer_filepath <- system.file("raster/landsat.tif", package = "spDataLarge") 
multi_layer <- rast(multilayer_filepath)

myurl = "/vsicurl/https://zenodo.org/record/5774954/files/clm_snow.prob_esacci.dec_p.90_500m_s0..0cm_2000..2012_v2.0.tif"
snow = rast(myurl)
snow

rey <- data.frame(lon = -21.94, lat = 64.15)
snow_rey <- extract(snow, rey)
snow_rey

write_sf(obj = world, dsn = "world.gpkg")
write_sf(obj= world, dsn = "world_many_layers.gpkg", append = T)
st_write(obj = world, dsn = "world2.gpkg")

write_sf(cycle_hire_xy, "cycle_hire_xy.csv", layer_options = "GEOMETRY=AS_XY")
write_sf(world_wkt, "world_wkt.csv", layer_options = "GEOMETRY=AS_WKT")

writeRaster(single_layer, filename = "my_raster.tif", datatype = "INT2U")
writeRaster(x = single_layer, filename = "my_raster.tif",
            gdal = c("COMPRESS=NONE"), overwrite = TRUE)
writeRaster(x = single_layer, filename = "my_raster.tif",
            filetype = "COG", overwrite = TRUE)

png(filename = "lifeExp.png", width = 500, height = 350)
plot(world["lifeExp"])
dev.off()

library(tmap)
tmap_obj <- tm_shape(world) + tm_polygons(col = "lifeExp")
tmap_save(tmap_obj, filename = "lifeExp_tmap.png")

install.packages("mapview")
library(mapview)
mapview_obj <- mapview(world, zcol = "lifeExp", legend = T)
webshot::install_phantomjs()
mapshot(mapview_obj, file = "my_interactive_map.html")


