library(sf)
library(terra)
library(dplyr)

# This chapter covers raster cropping and masking using vector objects;
# extracting raster values using different types of vector objects; 
# raster-vector conversion. 

srtm <- rast(system.file("raster/srtm.tif", package = "spDataLarge"))
plot(srtm)
zion <- read_sf(system.file("vector/zion.gpkg", package = "spDataLarge"))
plot(zion$geom)
zion_2 <- st_transform(zion, crs(srtm))
plot(zion_2$geom)
srtm_cropped <- crop(srtm, zion_2)
plot(srtm_cropped)
srtm_masked <- mask(srtm, zion_2)
srtm_final <- mask(srtm_cropped, zion_2)
all.equal(srtm_final, srtm_masked)
srtm_inverse_masked <- mask(srtm, zion_2, inverse = T)

data("zion_points", package = "spDataLarge")
elevation = terra::extract(srtm, zion_points)
zion_points <- cbind(zion_points, elevation)
plot(elevation)

zion_transect <- cbind(c(-113.2,-112.9), c(37.45,37.2)) |>
st_linestring() |> st_sfc(crs = crs(srtm)) |> st_sf(geometry = _)
zion_transect$geometry
zion_transect$id <- 1:nrow(zion_transect)
zion_transect <- st_segmentize(zion_transect, dfMaxLength = 250)
zion_transect <- st_cast(zion_transect, "POINT")
zion_transect <- zion_transect |> group_by(id) |> mutate(dist = st_distance(geometry)[,1])
zion_elev <- terra::extract(srtm, zion_transect)
zion_transect <- cbind(zion_transect, zion_elev)
plot(zion_transect$dist)

zion_srtm_values <- terra::extract(x = srtm, y = zion_2)
zion_srtm_values
group_by(zion_srtm_values, ID) |> summarize(across(srtm, list(min = min, max = max, mean = mean)))

nlcd <- rast(system.file("raster/nlcd.tif", package = "spDataLarge"))
class(nlcd)
zion2 <- st_transform(zion, st_crs(nlcd))
zion_nlcd <- terra::extract(nlcd, zion2)
zion_nlcd |> group_by(ID, levels) |> count()
plot(nlcd)
plot(zion2$geom)

install.packages("exactextractr")
library(exactextractr)

cycle_hire_osm = spData::cycle_hire_osm
cycle_hire_osm_projected <- st_transform(cycle_hire_osm, "EPSG:27700")
raster_template <- rast(ext(cycle_hire_osm_projected), resolution = 1000, crs = st_crs(cycle_hire_osm_projected)$wkt)

ch_raster1 <- rasterize(cycle_hire_osm_projected, raster_template) 
ch_raster2 <- rasterize(cycle_hire_osm_projected, raster_template, fun = "length")
ch_raster2
plot(ch_raster1)
plot(ch_raster2)
ch_raster3 <- rasterize(cycle_hire_osm_projected, raster_template, field = "capacity", fun = sum, na.rm = T)
plot(ch_raster3)

library(spData)
california = dplyr::filter(us_states, NAME == "California")
california_borders <- st_cast(california, "MULTILINESTRING")
raster_template2 = rast(ext(california), resolution = 0.5, crs = st_crs(california)$wkt)
plot(rasterize(california_borders, raster_template2, touches = T))
plot(rasterize(california, raster_template2, touches = T))

# Usually vectorization in R refers to bypassing for loops and such by the command 1:10/2

elev <- rast(system.file("raster/elev.tif", package = "spData"))
elev_point <- as.points(elev) |> st_as_sf()
plot(elev)
plot(elev_point)

dem = rast(system.file("raster/dem.tif", package = "spDataLarge"))
cl <- as.contour(dem) |> st_as_sf()
plot(dem, axes = F)
plot(cl, add = T)

grain = rast(system.file("raster/grain.tif", package = "spData"))
grain_poly <- as.polygons(grain, dissolve = T) |> st_as_sf()
plot(grain_poly)
plot(smoothr::smooth(grain_poly))

# Exercises:
library(sf)
library(terra)
library(spData)
zion_points_path <- system.file("vector/zion_points.gpkg", package = "spDataLarge")
zion_points <- read_sf(zion_points_path)
srtm = rast(system.file("raster/srtm.tif", package = "spDataLarge"))
ch = st_combine(zion_points) |> st_convex_hull() |> st_as_sf()
plot(ch)
plot(zion_points)

# E1

zion_points = st_transform(zion_points, crs(srtm))
srtm_cropped <- crop(srtm, zion_points)
plot(srtm_cropped)

ch = st_transform(ch, crs(srtm))
srtm_cropped_2 <- crop(srtm, ch)
plot(srtm_cropped_2)

srtm_mask <- mask(srtm, zion_points, inverse=F)
plot(srtm_mask)
srtm_mask_2 <- mask(srtm, ch)
plot(srtm_mask_2)

# E2

srtm_extracted <- terra::extract(srtm, zion_points)
srtm_extracted
zion_points_buffer <- st_buffer(zion_points, dist = 90)
srtm_extracted_2 <- terra::extract(srtm, zion_points_buffer, fun = mean, exact = F)
srtm_extracted_2

# DOUBT
library(exactextractr)
srtm_extracted_3 <- exact_extract(srtm,zion_points)
srtm_extracted_3
