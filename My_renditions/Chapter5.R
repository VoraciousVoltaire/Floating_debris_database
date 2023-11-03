library(spDataLarge)
library(sf)

# SIMPLIFICATION

seine_simp <- st_simplify(seine, dTolerance = 2000)
plot(seine_simp)
plot(seine)
object.size(seine)
object.size(seine_simp)

us_states_2163 <- st_transform(us_states, "EPSG:2163")
us_states_simp1 <- st_simplify(us_states_2163, dTolerance = 100000)
plot(us_states_2163)
plot(us_states_simp1)

install.packages("rmapshaper")
library(rmapshaper)
us_states_simp2 <- rmapshaper::ms_simplify(us_states_2163, keep = 0.01, keep_shapes = T)
plot(us_states_simp2)

install.packages("smoothr")
library(smoothr)
us_states_simp3 <- smoothr::smooth(us_states_2163, method = "ksmooth", smoothness = 6)

nz_centroid <- st_centroid(nz)
seine_centroid <- st_centroid(seine)
nz_pos <- st_point_on_surface(nz)
seine_pos <- st_point_on_surface(seine)

plot(nz_pos)
plot(seine)
plot(seine_pos)
seine_buff_5km <- st_buffer(seine, dist = 5000)
seine_buff_50km <- st_buffer(seine, dist = 50000)
plot(seine_buff_5km)
plot(seine_buff_50km)
nz_sfc <- st_geometry(nz)
class(nz_sfc)
nz_shift <- nz_sfc + c(0,100000)
plot(nz_shift)
plot(nz_sfc)
nz_centroid_sfc <- st_centroid(nz_sfc)
nz_scale <- (nz_sfc - nz_centroid_sfc)*0.5 + nz_centroid_sfc

rotation <- function(a){r = a*pi/180
matrix(c(cos(r), sin(r), -sin(r), cos(r)), nrow = 2, ncol = 2)}
nz_rotate = (nz_sfc - nz_centroid_sfc) * rotation(30) + nz_centroid_sfc
plot(nz_rotate)
nz_scale_sf <- st_set_geometry(nz, nz_scale)
plot(nz_scale_sf)

b = st_sfc(st_point(c(0,1)), st_point(c(1,1)))
b = st_buffer(b, dist=1)
plot(b, border = "grey")
text(x=c(-0.5,1.5), y=1, labels = c("x","y"), cex=3)

x = b[1]
y = b[2]
x_and_y <- st_intersection(x,y)
plot(b,border="grey")
plot(x_and_y, col = "lightgrey", border = "grey", add = T)

bb = st_bbox(st_union(x,y))
box = st_as_sfc(bb)
class(bb)
class(box)
set.seed(2017)
p = st_sample(x = box, size = 10)
p_xy1 = p[x_and_y]
plot(box, border = "grey", lty = 3)
plot(x, add=T, border="grey")
plot(y, add=T, border="grey")
plot(p, add=T)
plot(p_xy1, cex=3, col="red", add=T)
text(x=c(-0.5,1.5), y=1, labels=c("x","y"), cex=2)
p_xy2 <- st_intersection(x_and_y, p)
sel_p_xy <- st_intersects(p,x,sparse=F)&st_intersects(p,y,sparse=F)
p_xy3 <- p[sel_p_xy]
p_xy1 == p_xy2 
p_xy2 == p_xy3

regions <- aggregate(x=us_states[,"total_pop_15"], by=list(us_states$REGION), FUN=sum, na.rm=T)
regions_2 <- us_states |> group_by(REGION) |> summarize(pop = sum(total_pop_15, na.rm = T))
plot(regions)
plot(regions_2)
us_west <- us_states[us_states$REGION=="West",]
plot(us_west$geometry)
us_west_union <- st_union(us_west)
plot(us_west_union)
texas = us_states[us_states$NAME=="Texas",]
plot(texas$geometry)
texas_union=st_union(us_west_union, texas)
plot(texas_union)

# Type trasnformations

multipoint <- st_multipoint(matrix(c(1,3,5,1,3,1), ncol = 2))
plot(multipoint)
linestring <- st_cast(multipoint, "LINESTRING")
plot(linestring, add = T)
polyg <- st_cast(multipoint, "POLYGON")
plot(polyg, add = T)
multipoint_2 <- st_cast(linestring, "MULTIPOINT")
multipoint_3 <- st_cast(polyg, "MULTIPOINT")
all.equal(multipoint, multipoint_2)
multipoint == multipoint_2

multiline_string_list = list(matrix(c(1,4,5,3), nrow = 2),
                             matrix(c(4,4,4,1), nrow = 2),
                             matrix(c(2,4,2,2), nrow = 2))
multiline_string <- st_multilinestring(multiline_string_list)
multiline_string
plot(multiline_string)
multiline_string_sf <- st_sf(geom = st_sfc(multilinestring))
multiline_string_sf
class(st_sfc(multilinestring))
linestring_sf2 <- st_cast(multiline_string_sf, "LINESTRING")
linestring_sf2
linestring_sf2$name = c("Riddle Rd", "Marshall Ave", "Foulke St")
linestring_sf2$length = st_length(linestring_sf2)
linestring_sf2

elev <- rast(system.file("raster/elev.tif", package = "spData"))
clip = rast(xmin = 0.9, xmax = 1.8, ymin = -0.45, ymax = 0.45, resolution = 0.3, vals = rep(1,9) )
plot(elev[clip, drop = F])
elev
clip

elev_2 <- extend(elev, c(1,2))
plot(elev_2)
plot(elev)
elev_3 <- elev + elev_2
elev_4<- extend(elev, elev_2)
origin(elev_4)
origin(elev_4) = c(0.25,0.25)

dem = rast(system.file("raster/dem.tif", package = "spDataLarge"))
dem_agg <- aggregate(dem, fact = 5, fun = mean)
plot(dem)
plot(dem_agg)
dem_disagg <- disagg(dem_agg, fact = 5, method = "bilinear")
identical(dem, dem_disagg)
plot(dem_disagg)
compareGeom(dem,dem_disagg)
extend(dem, dem_disagg)
compareGeom(dem, dem_disagg)
all.equal(dem, dem_disagg)

target_rast <- rast(xmin = 794650, xmax = 798250, ymin = 8931750, ymax = 8935350, 
                    resolution = 300, crs = "EPSG:32717")
dem_resampl <- resample(dem, y = target_rast, method = "bilinear")
plot(dem_resampl)
plot(target_rast)
sf::gdal_utils()

# Exercises

# E1

nz_simplify <- st_simplify(nz, dTolerance = 100000)
plot(nz_simplify$geom)
library(rmapshaper)
nz_simplify_2 <- ms_simplify(nz, keep_shapes = T, keep = 0.2)
plot(nz_simplify_2$geom)
# st_simplify() works on the Douglas-Peuker algorithm which removes topological geometry resulting in holey and overlapping aerial units; 
# ms_simplify() works on the Visvalingam algorithm.

# E2

