install.packages("spData")
install.packages("spDataLarge", repos = "https://nowosad.r-universe.dev")
install.packages("remotes")
library(remotes)
remotes::install_github("r-tmap/tmap@v4")
remotes::install_github("geocompx/geocompkg", dependencies = TRUE)
library(sf)
sf :: sf_use_s2(T)
vignette(package = "sf")
vignette("sf1")
library(spData)
class(world)
names(world)
View(world)
plot(world)
summary(world["lifeExp"])
mini_world <- world[1:2, 1:3]
mini_world
mini_world$geom
mini_world_2 = world[1:2, 1:4]
mini_world_2
install.packages("styler")
install.packages("tibble")
world_df <- st_read(system.file("shapes/world.shp", package = "spData"))
class(world_df)
world_tbl <- read_sf(system.file("shapes/world.shp", package = "spData"))
class(world_tbl)
library(sp)
world_sp <- as(world,"Spatial")
world_sf <- st_as_sf(world_sp)
plot(world["pop"])
world_asia <- world[world$continent=="Asia",]
world_asia
asia <- st_union(world_asia)
plot(world["pop"], reset = F)
plot(asia, add=T, col="red")
# see ?graphics::plot and ?par
plot(world["continent"], reset = F)
cex = sqrt(world$pop)/10000
world_cent <- st_centroid(world, of_largest = T)
plot(st_geometry(world_cent), add = T, cex =cex)
india <- world[world$name_long == "India",]
plot(st_geometry(india),expandBB=c(0,0.2,0.1,1) ,lwd = 3, col = "grey")
plot(st_geometry(world_asia), add = T)

#Exercise
plot()

# To create an sf object
london_coords_pts <- st_point(c(0.1,51.5))
london_geometry <- st_sfc(london_coords_pts, crs = "EPSG:4326")
london_nongeographical_attributes <- data.frame(name = "London", temperature=25, date= as.Date("2017-06-21"))
london_sf <- st_sf(london_nongeographical_attributes, geometry = london_geometry)
london_sf

#creating a multipoint and linestring matrix
multipoint_matrix <- rbind(c(5,2), c(1,3), c(3,4), c(3,2))
st_multipoint(multipoint_matrix)
linestring_matrix <- rbind(c(1,5), c(4,4), c(4,1), c(2,2), c(3,2))
plot(st_linestring(linestring_matrix))

polygon_list = list(rbind(c(1,5), c(2,2), c(4,1), c(4,4), c(1,5)))
plot(st_polygon(polygon_list))

## POLYGON with a hole
polygon_border = rbind(c(1, 5), c(2, 2), c(4, 1), c(4, 4), c(1, 5))
polygon_hole = rbind(c(2, 4), c(3, 4), c(3, 3), c(2, 3), c(2, 4))
polygon_with_hole_list = list(polygon_border, polygon_hole)
plot(st_polygon(polygon_with_hole_list))
#> POLYGON ((1 5, 2 2, 4 1, 4 4, 1 5), (2 4, 3 4, 3 3, 2 3, 2 4))

multilinestring_1 <- list(rbind(c(1,1), c(2,2), c(3,3)), rbind(c(4,4), c(5,5), c(6,6)))
multilinestring1 = st_multilinestring(multilinestring_1)
multilinestring_2 <- list(rbind(c(7,7), c(8,8), c(9,9)), rbind(c(10,10), c(11,11), c(12,12)))
multilinestring2 = st_multilinestring(multilinestring_2)
multilinestring_sfc <- st_sfc(multilinestring1, multilinestring2)
st_geometry_type(multilinestring_sfc)
st_crs(multilinestring_sfc)

m = matrix(1:8, ncol = 2)
sfheaders::sfg_linestring(obj=m)
v = c(1,2)
check <- sfheaders::sfg_linestring(obj=v)
st_geometry_type(check)
df <- data.frame(x=1:4, y=4:1)
sfheaders::sf_polygon(obj=df)
df

sf_use_s2()
plot(india_buffer_with_s2 <- st_buffer(india,1))
sf_use_s2(F)
plot(india_buffer_with_s2 <- st_buffer(india,1))

# Creating a Spatraster object
library(terra)
raster_filepath <- system.file("raster/srtm.tif", package = "spDataLarge")
my_rast <- rast(raster_filepath)
class(my_rast)
my_rast
inMemory(my_rast)
help("terra-package")
plot(my_rast)
install.packages("tmap", "rasterVis")
new_raster <- rast(nrows = 6, ncols = 6, xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5, vals = 1:36)
multi_raster_file <- system.file("raster/landsat.tif", package = "spDataLarge")
multi_rast <- rast(multi_raster_file)
multi_rast
nlyr(multi_rast)
multi_rast3 <- subset(multi_rast, 3)
multi_rast4 <- subset(multi_rast, "landsat_4")
multi_rast34 <- c(multi_rast3, multi_rast4)
install.packages("PROJ")
library(PROJ)
sf_proj_info(type = "proj")
attributes(st_area(india))
install.packages("units")
library(units)
units::set_units(st_area(india), km^2)
