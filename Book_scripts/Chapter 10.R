library(sf)
library(terra)

install.packages("qgisprocess")
library(qgisprocess)
qgis_enable_plugins()
install.packages("remotes")
remotes::install_github("r-spatial/qgisprocess")

install.packages("Rsagacmd")
library(Rsagacmd)

install.packages("rgrass")
library(rgrass)

install.packages("rstac")
library(rstac)

install.packages("gdalcubes")
library(gdalcubes)

install.packages("mapedit")
library(mapedit)

input <- sf::read_sf(system.file("shape/nc.shp", package = "sf"))

result <- qgis_run_algorithm(
  "native:buffer",
  INPUT = input,
  DISTANCE = 1,
  DISSOLVE = T
)
result

output_sf <- sf::st_as_sf(result)
plot(sf::st_geometry(output_sf))

qgis_show_help("native:buffer")

qgis_plugins()
qgis_enable_plugins("grassprovider", quiet = T)
qgis_providers()

data("incongruent", "aggregating_zones", package = "spData")
incongr_wgs <- st_transform(incongruent, "EPSG:4326")
aggzone_wgs <- st_transform(aggregating_zones, "EPSG:4326")
qgis_algorithms()
qgis_search_algorithms("union")
qgis_show_help("native:union")

alg = "native:union"
union_arguments <- qgis_get_argument_specs(alg)
union_arguments
union_arguments$name
getwd()

union = qgis_run_algorithm(alg, INPUT = incongr_wgs, OVERLAY = aggzone_wgs)
union_sf <- st_as_sf(union)
plot(union_sf)
plot(st_union(incongr_wgs, aggzone_wgs))

# To clean silver polygons

qgis_search_algorithms("clean")
qgis_show_help("grass7:v.clean")
library(tidyr)
library(dplyr)
qgis_get_argument_specs("grass7:v.clean") |> select(name, description) |> dplyr::slice_head(n = 4)
clean <- qgis_run_algorithm("grass7:v.clean", input = union_sf, tool = "rmarea", threshold = 25000)
clean_sf <- st_as_sf(clean)
plot(clean_sf)

library(qgisprocess)
library(terra)
dem = system.file("raster/dem.tif", package = "spDataLarge")
qgis_search_algorithms("wetness") |> dplyr::select(provider_title, algorithm) |> head(2)

qgis_show_help("sagang:sagawetnessindex")
