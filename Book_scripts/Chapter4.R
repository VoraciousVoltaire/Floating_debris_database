install.packages("fuzzyjoin")
library(spData)
library(terra)
library(sf)
library(dplyr)
canterbury <- nz |> filter(Name == "Canterbury")
plot(canterbury_height <- nz_height[canterbury,])
View(nz)
View(nz_height)
plot(nz)
plot(nz_height)
nz_height[canterbury,,op = st_disjoint]
nz_height[canterbury,2,op = st_disjoint]
weird_name <- st_intersects(x = nz_height, y = canterbury)
class(weird_name)
weird_name
weird_logical <- lengths(weird_name) > 0 
canterbury_2 <- nz_height[weird_logical,]
plot(canterbury_2)
canterbury_3 <- nz_height |> st_filter(y= canterbury, .predicate = st_intersects)

# Creating an sfc object
polygon_matrix <- cbind(
  x= c(0,0,1,1,0),
  y= c(0,1,1,0.5,0)
)
polygon_sfc <- st_sfc(st_polygon(list(polygon_matrix)))
polygon_sfc

line_sfc <- st_sfc(st_linestring(cbind(x=c(0.4,1), y=c(0.2,0.5))))
point_df <- data.frame(x=c(0.2,0.7,0.4), y=c(0.1,0.2,0.8))
point_sf <- st_as_sf(point_df, coords=c("x", "y"))
st_intersects(point_sf, polygon_sfc)
st_intersects(point_sf, polygon_sfc, sparse = F)
st_within(point_sf, polygon_sfc)
st_touches(point_sf, polygon_sfc)
st_disjoint(point_sf, polygon_sfc, sparse=F)[,1]
st_distance(point_sf, polygon_sfc)
st_is_within_distance(point_sf, polygon_sfc, dist = 0.135, sparse = F)[,1]

# Overwhelming DE-9IM 
xy2sfc <- function(x,y) st_sfc(st_polygon(list(cbind(x,y))))
x <- xy2sfc(x=c(0,0,1,1,0), y=c(0,1,1,0.5,0))
y <- xy2sfc(x=c(0.7,0.7,0.9,0.7), y=c(0.8,0.5,0.5,0.8))
st_relate(x,y)
st_queen <- function(x,y) st_relate(x,y,pattern="F***T****")
st_rook <- function(x,y) st_relate(x,y, pattern="F***1****")
grid = st_make_grid(x,n=3)
grid_sf = st_sf(grid)
grid_sf$queens <- lengths(st_queen(grid,grid[5])) > 0
plot(grid,col=grid_sf$queens)
grid_sf$rooks <- lengths(st_rook(grid, grid[5])) > 0
plot(grid, col=grid_sf$rooks)

# spatial joining
set.seed(2018)
(bb = st_bbox(world))
random_df <- data.frame(x=runif(n=10,min=bb[1],max=bb[3]),
                        y=runif(n=10,min=bb[2],max=bb[4]))
random_points <- random_df |> st_as_sf(coords=c("x","y"), crs="EPSG:4326")
class(random_points)
world_random <- world[random_points,]
nrow(world_random)
random_joined <- st_join(random_points, world["name_long"])
plot(world_random)
View(world_random)
plot(random_joined)
plot(st_geometry(cycle_hire), col="blue")
plot(st_geometry(cycle_hire_osm), col="red", add=T, pch=3)     
any(st_touches(cycle_hire,cycle_hire_osm, sparse=F))
sel = st_is_within_distance(cycle_hire,cycle_hire_osm,dist=units::set_units(20,"m"))
summary(lengths(sel)>0)
z<-st_join(cycle_hire, cycle_hire_osm, st_is_within_distance, dist=units::set_units(20,"m"))
nrow(cycle_hire)
nrow(z)
z <- z |> group_by(id) |> summarize(capacity=mean(capacity))
nrow(z)==nrow(cycle_hire)
plot(cycle_hire_osm["capacity"])
plot(z["capacity"])

nz_agg <- aggregate(x=nz_height, by=nz,FUN=mean)
st_geometry(nz_agg)==st_geometry(nz)
plot(nz_agg)
nz_agg_2 <- st_join(x=nz,y=nz_height) |> group_by(Name) |> summarize(elevation=mean(elevation,na.rm=T))
plot(nz_agg_2)

iv <- incongruent["value"]
agg_aw <- st_interpolate_aw(iv, aggregating_zones, extensive = T)
agg_aw$value
nz_highest <- nz_height |> slice_max(n=1, order_by=elevation)
canterbury_centroid <- st_centroid(canterbury)
st_distance(nz_highest, canterbury_centroid)
View(nz_height)
co = filter(nz, grepl("Canter|Otag", Name))
st_distance(nz_height[1:3,], co)
plot(st_geometry(co)[2])
plot(st_geometry(nz_height)[2:3], add=T)

elev=rast(system.file("raster/elev.tif", package="spData"))
grain= rast(system.file("raster/grain.tif", package="spData"))

id = cellFromXY(elev, xy = matrix(c(0.1,0.1), ncol=2))
id
terra::extract(elev, matrix(c(0.1,0.1), ncol=2))

clip = rast(xmin = 0.9,xmax = 1.8,ymin= -0.45, ymax= 0.45, resolution = 0.3, vals = rep(1,9))
elev[clip]
terra::extract(elev, ext(clip))
plot(clip)
plot(elev)
plot(elev[clip])
elev[1:2, drop = F]

rmask = elev
values(rmask) = sample(c(NA,T), 36, replace = T)
elev[rmask, drop=F]
mask(elev, rmask)
elev[elev<20] = NA

plot(elev + elev)
plot(elev^2)
plot(log(elev))
plot(elev>5)

rcl = matrix(c(0,12,1,12,24,2,24,36,3), ncol=3, byrow=T) # reclassification matrix
rcl
recl = classify(elev, rcl = rcl)
recl

multi_raster <- rast(system.file("raster/landsat.tif", package = "spDataLarge"))
View(multi_raster)
ndvi_fun <- function(nir,red){(nir-red)/(nir+red)}
ndvi_rast <- lapp(multi_raster[[c(4,3)]], fun = ndvi_fun)
plot(ndvi_rast)
plot(multi_raster)

r_focal <- focal(elev, w = matrix(1, nrow=3,ncol=3), fun = min)
plot(r_focal)

z = zonal(elev, grain, fun = "mean")
z

install.packages("geodata")
library(geodata)
aut <- geodata::elevation_30s(country="AUT", path = tempdir())
ch <- geodata::elevation_30s(country="CHE", path = tempdir())
aut_ch <- merge(aut,ch)
plot(aut_ch)
plot(aut)
plot(ch)


