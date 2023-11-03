library(sf)
library(terra)
library(dplyr)
library(spData)
library(spDataLarge)
remotes::install_github("r-tmap/tmap@v4")
library(tmap)
install.packages("shiny")
library(shiny)
library(leaflet)
library(ggplot2)
nz_elev <- rast(system.file("raster/nz_elev.tif", package = "spDataLarge"))
install.packages("cartogram")
library(cartogram)

tm_shape(nz) + tm_fill()
tm_shape(nz) + tm_borders()
tm_shape(nz) + tm_borders() + tm_fill()    
map_nz <- tm_shape(nz) + tm_polygons()
class(map_nz)
map_nz
tm_shape(nz) + tm_symbols()

map_nz1 <- map_nz + tm_shape(nz_elev) + tm_raster(alpha = 0.9)
map_nz1

nz_water <- st_union(nz) |> st_buffer(22200) |> st_cast(to = "LINESTRING")
map_nz2 <- map_nz1 + tm_shape(nz_water) + tm_lines()
map_nz2
map_nz3 <- map_nz2 + tm_shape(nz_height) + tm_symbols()
map_nz3
tmap_arrange(map_nz1, map_nz2, map_nz3)

ma1 <- tm_shape(nz) + tm_polygons(col = "red")
ma2 <- tm_shape(nz) + tm_polygons(col = "violet", alpha = 0.2)
ma3 <- tm_shape(nz) + tm_polygons(col = "blue")
ma4 <- tm_shape(nz) + tm_polygons(lwd = 3)
ma5 <- tm_shape(nz) + tm_polygons(lty = 2)
ma6 <- tm_shape(nz) + tm_polygons(border.col = "maroon", alpha = 0.9, col = "orange", lwd = 3, lty = 2)
tmap_arrange(ma1,ma2,ma3,ma4,ma5,ma6)

plot(st_geometry(nz), col = nz$Land_area)
tm_shape(nz) + tm_fill(col = "Land_area")
ma7 <- tm_shape(nz) + tm_polygons(col = "Median_income")

# Doesn't work because of outdated package in book
ma8 <- tm_shape(nz) + tm_polygons(col = "Median_income", fill.scale = tm_scale_bar(breaks = c(0,30000,40000,50000)))  
ma8
ma9 <- tm_shape(nz) + tm_polygons(col = "Median_income", fill.scale = tm_scale_bar(bg.color = "Green"))
tmap_arrange(ma7, ma8, ma9)
?tm_scale_bar

install.packages("classInt")
install.packages("tmaptools")
library(tmaptools)          

# Help youtube video

install.packages("osmdata")
library(osmdata)
installed.packages("osmplotr")
library(osmplotr)
# spread function in tidyr takes the dataset as the first argument, key = column_name as the second argument which specifies which column do you want to spread out such that the values inside those respective columns become sub-columns themselves and the third argument is value = attribute_which_is_to_be_stored_inside_the_sub-column_matrix.
# view mode in the tmap_mode is for interactive maps. tmap_mode = "plot" is for static maps.
# Doesn't work for some resaon becakse the video is old.
roads <- extract_osm_objects(key = "highway", bbox = c(xmin = -96.36947, 
                                                       ymin = 30.56764,
                                                       xmax = -96.34981,
                                                       ymax = 30.60611), sf = T)
# tm_bubbles(size = 0.25, col=c("column_name_1","column_name_2) +
# tm_compass(type = "arrow", text.size = 0.7, position = c(0.8,0.2)) +
# tm_layput(legend.position = c(0.6,0.5), legend.bg.color = "white", legend.frame = "black")
# tm_facets(by = "column_name", free.cords = F) free.cords asks if you wish the resolution of the zoomed in panels (resolution of panels) to be different for different coordinates. 

# Doesn't work:
legend_title <- expression("Area")
map_nz <- tm_shape(nz) + tm_polygons(col = "Land_area", 
                                      fill.legend = tm_legend(title = legend_title))
map_nz

map_nz + tm_layout(scale = 3)
map_nz + tm_layout(bg.color = "lightblue")
map_nz + tm_layout(frame = FALSE)

urb_1970_2030 = urban_agglomerations |> 
  filter(year %in% c(1970, 1990, 2010, 2030))

tm_shape(world) +
  tm_polygons() +
  tm_shape(urb_1970_2030) +
  tm_symbols(col = "black", size = "population_millions") +
  tm_facets(by = "year", nrow = 2)
nz_region = st_bbox(c(xmin = 1340000, xmax = 1450000,
                      ymin = 5130000, ymax = 5210000),
                    crs = st_crs(nz_height)) |> 
  st_as_sfc()
nz_height_map = tm_shape(nz_elev, bbox = nz_region) + tm_polygons() +
  tm_shape(nz_height) + tm_symbols(shape = 2, col = "red", size = 1) +
  tm_scale_bar(position = c("left", "bottom"))
nz_height_map
nz_map = tm_shape(nz) + tm_polygons() +
  tm_shape(nz_height) + tm_symbols(shape = 2, col = "red", size = 0.1) + 
  tm_shape(nz_region) + tm_borders(lwd = 3) +
  tm_layout(bg.color = "lightblue")

us_states_map = tm_shape(us_states, projection = "EPSG:2163") + tm_polygons() + 
  tm_layout(frame = FALSE)
us_states_map

hawaii_map = tm_shape(hawaii) +
  tm_polygons() + 
  tm_layout(frame = FALSE, bg.color = NA, 
            title.position = c("LEFT", "BOTTOM"))
alaska_map = tm_shape(alaska) +
  tm_polygons() +
  tm_layout(frame = FALSE, bg.color = NA)
us_states_map
print(hawaii_map, vp = grid::viewport(0.35, 0.1, width = 0.2, height = 0.1))
print(alaska_map, vp = grid::viewport(0.15, 0.15, width = 0.3, height = 0.3))

# Animated maps

install.packages("gganimate")
urb_anim = tm_shape(world) + tm_polygons() + 
  tm_shape(urban_agglomerations) + tm_symbols(size = "population_millions") +
  tm_facets(by = "year", nrow = 1, ncol = 1, free.coords = FALSE)
urb_anim
install.packages("gifski")
tmap_animation(urb_anim, filename = "urb_anim.gif", delay = 25)

install.packages("mapview")
install.packages("mapdeck")
installed.packages("leaflet.extras")

tmap_mode("view")
map_nz

map_nz + tm_basemap(server = "OpenTopoMap")

world_coffee <- left_join(world, coffee_data, by = "name_long")
facets = c("coffee_production_2016", "coffee_production_2017")
tm_shape(world_coffee) + tm_polygons(facets) + tm_facets(nrow = 1, sync = T)

tmap_mode("plot")
mapview::mapview(nz)
library(mapview)
oberfranken = subset(franconia, district == "Oberfranken")
oberfranken
trails |>
  st_transform(st_crs(oberfranken)) |>
  st_intersection(oberfranken) |>
  st_collection_extract("LINESTRING") |>
  mapview(color = "red", lwd = 3, layer.name = "trails") +
  mapview(franconia, zcol = "district") + breweries

install.packages("googleway")
library(googleway)
library(mapdeck)
set_token(Sys.getenv("MAPBOX"))
crash_data = read.csv("https://git.io/geocompr-mapdeck") # LINK DOESN'T WORK
crash_data = na.omit(crash_data)
ms = mapdeck_style("dark")
mapdeck(style = ms, pitch = 45, location = c(0, 52), zoom = 4) |>
  add_grid(data = crash_data, lat = "lat", lon = "lng", cell_size = 1000,
           elevation_scale = 50, colour_range = hcl.colors(6, "plasma"))

pal = colorNumeric("RdYlBu", domain = cycle_hire$nbikes)
leaflet(data = cycle_hire) |> 
  addProviderTiles(providers$CartoDB.Positron) |>
  addCircles(col = ~pal(nbikes), opacity = 0.9) |> 
  addPolygons(data = lnd, fill = FALSE) |> 
  addLegend(pal = pal, values = ~nbikes) |> 
  setView(lng = -0.1, 51.5, zoom = 12) |> 
  addMiniMap()

# leaflet youtube video
leaflet() %>% addTiles() %>% setView(lng = -81.2001, lat = 28.6024, zoom = 1)
providers
# esri natgeoworld map and cartodb are good default ones to work with 
library(dplyr)
library(tidyr)
leaflet() %>% addProviderTiles(provider = "OpenStreetMap.France") |> setView(lat=48.150623780019785, lng=-1.6208811102497036, zoom= 15)

leaflet() %>% addProviderTiles(provider = "Stamen") |> setView(lat=48.150623780019785, lng=-1.6, zoom = 15)

leaflet() %>% addProviderTiles(provider = "HikeBike") |> setView(lat=48.150623780019785, lng=-1.6, zoom = 15)
leaflet (
  options = leafletOptions(
    minZoom = 12,
    maxZoom = 18,
    dragging = T
  )
) %>%
  addProviderTiles(
    "Esri.WorldImagery"
  ) %>% setView(
    lng = -81.2001,
    lat = 26.604,
    zoom = 15
  ) %>%
  setMaxBounds(
    lng1 = -81.2001 + 0.05,
    lng2 = -81.2001 - 0.05,
    lat1 = 26.604 + 0.05,
    lat2 = 26.604 - 0.05
  )

v_longitude <- c(-81.1926, -81.200489, -81.2006)
v_latitude <- c(28.608, 28.601, 28.6005)
v_label <- c("Stadium","Student Union","Statistics")
leaflet() |> addProviderTiles("Esri.WorldImagery") |> setView(lng = -81.2001, lat = 26.604, zoom = 15) |>
  addMarkers(lng = v_longitude, lat = v_latitude, label = v_label)
leaflet() |> addProviderTiles("Esri.WorldImagery") |> setView(lng = -81.2001, lat = 26.604, zoom = 15) |>
  addMarkers(lng = v_longitude, lat = v_latitude, popup = v_label)

# Use clearBounds() and clearMarkers() to clear existing boundaries and markers respectively.

data(list = "breweries91", package = "leaflet")
breweries91 <- data.frame(breweries91, sp::coordinates(obj = breweries91))
colnames(breweries91)[colnames(breweries91) == "latitude"] <- "lat"
colnames(breweries91)[colnames(breweries91) == "longitude"] <- "lng"

colorNumeric_founded <- colorNumeric(palette = c("red", "blue"), domain = range(breweries91$founded, na.rm = T))
colorBin_founded <- colorBin(palette = c("red", "blue"), domain = breweries91$founded, bins = 2)
colorQuantile_founded <- colorQuantile(palette = c("red", "blue"), domain = breweries91$founded)
v_villages <- names(sort(x = table(breweries91$village), decreasing = T))
colorFactor_village <- colorFactor(
  palette = c("#440154FF", "#443A83FF", "#31688EFF", "#21908CFF", "#35B779FF", "#8FD744FF", rep("#FDE725FF",18)), 
  domain = breweries91$founded, 
  levels = v_villages)

leaflet(data = breweries91) |> addProviderTiles(provider = "Stamen") |>
  addCircles(lng = ~lng, lat = ~lat, label = ~brewery, radius = 19, color = "red")
leaflet(data = breweries91) |> addProviderTiles(provider = "Stamen") |>
  addMarkers(lng = ~lng, lat = ~lat, label = ~brewery)
leaflet(data = breweries91) |> addProviderTiles(provider = "Stamen") |>
  addLabelOnlyMarkers(lng = ~lng, lat = ~lat, label = ~brewery)
leaflet(data = breweries91) |> addProviderTiles(provider = "Stamen") |>
  addCircleMarkers(lng = ~lng, lat = ~lat, label = paste(breweries91$brewery, breweries91$founded),
                   color = ~colorNumeric_founded(breweries91$founded))
leaflet(data = breweries91) |> addProviderTiles(provider = "Stamen") |>
  addCircleMarkers(lng = ~lng, lat = ~lat, label = ~paste(brewery,founded),
                   color = ~colorFactor_village(breweries91$village),
                   opacity = 1, fillOpacity = 1) |> addLegend(pal = colorFactor_village, values = v_villages, title = "Village", position = "topleft", opacity = 1)

install.packages("leaflet.extras")
library(leaflet.extras)

leaflet() |> addTiles() |> setView(lng = -81.2001, lat = 28.6024, zoom = 15) |> addSearchOSM()

median_lng <- median(x = breweries91$lng, na.rm = T)
median_lat <- median(x = breweries91$lat, na.rm = T)
breweries91$Quadrant <- NA
breweries91$Quadrant[breweries91$lng <= median_lng & breweries91$lat <= median_lat] <- "SW"
breweries91$Quadrant[breweries91$lng > median_lng & breweries91$lat <= median_lat] <- "SE"
breweries91$Quadrant[breweries91$lng <= median_lng & breweries91$lat > median_lat] <- "NW"
breweries91$Quadrant[breweries91$lng > median_lng & breweries91$lat > median_lat] <- "NE"
colorFactor_Quadrant <- colorFactor(palette = c("red", "green", "blue", "black"), levels = c("SW", "SE", "NW", "NE"))

leaflet() |> addProviderTiles(providers$Esri.WorldStreetMap) |>
  addCircleMarkers(
    data = breweries91[breweries91$Quadrant == "SW",],
    lat = ~lat,
    lng = ~lng,
    col = ~colorFactor_Quadrant(Quadrant),
    group = "SW"
  ) |>
  addCircleMarkers(
    data = breweries91[breweries91$Quadrant == "SE",],
    lat = ~lat,
    lng = ~lng,
    col = ~colorFactor_Quadrant(Quadrant),
    group = "SE" ) |>
  addCircleMarkers(
    data = breweries91[breweries91$Quadrant == "NW",],
    lat = ~lat,
    lng = ~lng,
    col = ~colorFactor_Quadrant(Quadrant),
    group = "NW" ) |>
  addCircleMarkers(
    data = breweries91[breweries91$Quadrant == "NE",],
    lat = ~lat,
    lng = ~lng,
    col = ~colorFactor_Quadrant(Quadrant),
    group = "NE" ) |>
  addLayersControl(overlayGroups = c("SW", "SE", "NW", "NE"))

leaflet() |> setView(lng = median_lng, lat = median_lat, zoom = 10) |>
  addTiles(group = "OSM") |>
  addProviderTiles("CartoDB", group = "Carto") |>
  addProviderTiles("Esri", group = "Esri") |>
  addLayersControl(baseGroups = c("OSM", "Carto", "Esri"))

# use clusterOptions() = markerClusterOptions() inside addCircleMarkers to club clusters for proximal data points

install.packages("shiny")
library(shiny)
library(leaflet)
library(spData)
ui = fluidPage(
  sliderInput(inputId = "life", "Life expectancy", 49, 84, value = 80),
  leafletOutput(outputId = "map")
)
server = function(input, output){
  output$map <- renderLeaflet(
    {
      leaflet() %>% 
        addPolygons(data = world[world$lifeExp < input$life, ])
    }
  )
}
shinyApp(ui, server)

setwd("/Users/ameydanole/Desktop/r")
runApp("app.r")
library(shiny)    # for shiny apps
library(leaflet)  # renderLeaflet function
library(spData)   # loads the world dataset 
ui_2 = fluidPage(
  titlePanel = "Mother sir",
  sliderInput(inputId = "Age", "Places to visit at the age of", 5, 83, value = 43),
  leafletOutput(outputId = "map")
)
server_2 = function(input, output) {
  output$map = renderLeaflet({
    leaflet() |> 
      addPolygons(data = world[world$lifeExp > input$Age, ])})
}
shinyApp(ui_2, server_2)

g <- st_graticule(nz, lon = c(170,175), lat = c(-45, -40, -35))
plot(nz_elev, graticule = g, axes = T, col = "blue")
terra::plot(nz_elev/1000, add = T, axes = F)
plot(st_geometry(nz), add = T)

library(ggplot2)
g1 <- ggplot() + geom_sf(data = nz, aes(fill = Median_income)) + 
  geom_sf(data = nz_height) + scale_x_continuous(breaks = c(170,175))
g1

install.packages("plotly")
library(plotly)
g2 <- plotly::ggplotly(g1)
g2
