# Compare behavioural state estimates among models---- 
# credits: Joshua Cullen

# Packages
library(tidyverse)
library(lubridate)
library(amt)
library(sf)
library(rnaturalearth)
library(plotly)
library(ggplot2)

# Prepare data for amt package requirements
trial <- relevant_new_data_2 |> filter(colony == "Bjørnøya")
dat.track <- make_track(trial, lon, lat, timestamp, crs = 4326, all_cols = T)

# Calculating isopleth levels of MCPs: at what percent you want to estimate the MCP: doesn't include all points
dat.mcp <- hr_mcp(dat.track, levels = c(0.5,0.95,1)) # 50% is the usual mark for denoting core area of use; 
# 95% isopleth is used to estimate total area of use to exlude any extreme outliers 
plot(dat.mcp, col = c('red', 'green', 'blue'))
norway <- ne_countries(scale = 50, country = "Norway", returnclass = 'sf')
coupled_countries <- ne_countries(scale = 50, country= c("Norway", "Sweden", "Finland", "Russia"), returnclass = 'sf')

# ggplot with couple spatial vector : for population
ggplot() +
  geom_sf(data = coupled_countries) +
  geom_point(data = trial, aes(lon, lat), alpha = 0.05, size = 1) +
  geom_sf(data = dat.mcp$mcp, aes(colour = factor(level)), fill = 'transparent', size = 0.75) +
  scale_colour_viridis_d(direction = -1) +
  theme_bw() +
  coord_sf(xlim = c(-66, 72), ylim = c(50, 83))

# Calculating MCPs per ID (at 95% level)
# Need to turn dataframe into a list to map hr_mcp function and then recombine

dat.id.mcp <- dat.track %>% 
  split(.$individ_id) %>% # each element of the list is an individual id
  map(hr_mcp, levels = 0.95) %>% # what this basically says is that map hr_mcp function for each element that is split into unique IDs
  map(pluck, 1) %>% # just plucking or extracing the spatial vector layer
  do.call(rbind, .) #combining all them polygons together

dat.id.mcp <- dat.id.mcp %>% 
  mutate(id = rownames(.), .before = level) # adding a column for id by pasting all rownames in that column and placing it before
# the 'level' column


# Using ggplotly so that we can also zoom in

ggplotly(
  ggplot() +
    geom_sf(data = coupled_countries) +
    geom_point(data = trial, aes(lon, lat, colour = factor(individ_id)), alpha = 0.05, size = 1) +
    geom_sf(data = dat.id.mcp, aes(colour = id), fill = 'transparent', size = 0.75) +
    scale_colour_viridis_d() +
    theme_bw() +
    coord_sf(xlim = c(-66, 72), ylim = c(50, 83))
)

# Export data for ez use
setwd('/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/My_renditions/')
save(dat.id.mcp, file = "MCP_fits.RData")

# Kernel density estimate vid
trast <- make_trast(dat.track) # why can't I set a desired resolution here?
trast

# Calculating KDE using href
dat.kde.ref <- hr_kde(dat.track, trast = trast, h = hr_kde_ref(dat.track), levels = c(0.5, 0.95))
plot(dat.kde.ref, col = c("red", "blue"))
kde.href.contours <- hr_isopleths(dat.kde.ref)

# ggplot with couple
ggplot() +
  # geom_tile(data = plastics_spdf_df, aes(x = x, y = y, fill = Value), alpha = 0.8) + # because produces a very distorted graph
  geom_sf(data = coupled_countries) +
  geom_point(data = trial, aes(lon,lat), alpha = 0.25, size = 1, color = "chartreuse") +
  scale_fill_viridis_c() +
  theme_bw() +
  coord_sf(xlim = c(-66, 72), ylim = c(50, 83))

# Make isopleth contours at requested levels
h_ref_ggplot <- ggplot() +
  geom_sf(data = coupled_countries) +
  geom_path(data = trial, aes(lon,lat,group = individ_id), alpha = 0.25, size = 0.3) +
  geom_sf(data = kde.href.contours, aes(colour = factor(level)), fill = 'transparent', size = 0.75) +
  scale_fill_viridis_c() +
  coord_sf(xlim = c(-66, 72), ylim = c(50, 83))
h_ref_ggplot

# Plug-in method:
h_pi <- hr_kde_pi(dat.track)
dat.kde.pi <- hr_kde(dat.track, trast = trast, h = h_pi, levels = c(0.5,0.95))
dat.kde.pi
plot(dat.kde.pi, col = c("red", "blue"))
kde.hpi.contours <- hr_isopleths(dat.kde.pi)

# same as before
h_pi_ggplot <- ggplot() +
  geom_sf(data = coupled_countries) +
  geom_path(data = trial, aes(lon,lat,group = individ_id), alpha = 0.25, size = 0.3) +
  geom_sf(data = kde.hpi.contours, aes(colour = factor(level)), fill = 'transparent', size = 0.75) +
  scale_fill_viridis_c() +
  coord_sf(xlim = c(-66, 72), ylim = c(50, 83))

# comparing both methods: very close
ggarrange(h_ref_ggplot, h_pi_ggplot, ncol = 1, nrow = 2)

# Map kde per individual id----
# href
dat.id.kde.href <- dat.track %>%
  split(.$individ_id) %>%
  map(~hr_kde(.x,
              trast = make_trast(.x),
              h = hr_kde_ref(.x),
              levels = c(0.5, 0.95))) %>%
  # map(hr_isopleths) 
  do.call(rbind, .)
 
dat.id.kde.href <- as.data.frame(dat.id.kde.href) %>%
  mutate(id = rownames(.), .before = estimator)
ggplot() +
  geom_sf(data = coupled_countries) +
  geom_path(data = trial, aes(lon, lat, group = individ_id), alpha = 0.25, size = 0.3) +
 # geom_sf(data = dat.id.kde.href, aes(colour = factor(level)), fill = 'transparent', size = 0.75) +
  coord_sf(xlim = c(-66, 72), ylim = c(50, 83)) +
  facet_wrap(~individ_id)

# h_pi
dat.id.kde.hpi <- dat.track %>%
  split(.$individ_id) %>%
  map(~hr_kde(.x,
              trast = make_trast(.x),
              h = hr_kde_pi(.x),
              levels = c(0.5, 0.95))) %>%
  map(hr_isopleths) %>% # this worked!
  do.call(rbind, .)
dat.id.kde.hpi <- dat.id.kde.hpi %>%
  mutate(individ_id = rownames(.), .before = level)
ggplot() +
  geom_sf(data = coupled_countries) +
  geom_path(data = trial, aes(lon, lat, group = individ_id), alpha = 0.25, size = 0.3) +
  geom_sf(data = dat.id.kde.hpi, aes(color = factor(level)), fill = 'transparent', size = 0.75) +
  coord_sf(xlim = c(-66, 72), ylim = c(50, 83)) +
  facet_wrap(~individ_id)
# the best graph so far 
# To do: overlay the location of Bjønøya atop this

# Dynamic Brownian Bridge Models----

install.packages("move")
library(move) # used for running dynamic BBMs

# Creating a column for location errors
trial <- trial %>% mutate(location_error = ifelse(
  loc_type == "IRMA", 180, 200
), .before = colony)
trial_sf <- st_as_sf(trial, coords = c("lon", "lat"), crs = 4326)
trial_sf_merc <- st_transform(trial_sf, CRS('+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs'))
trial_sf_merc_split <- trial_sf_merc %>% mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2])
trial_merc_df <- data.frame(trial_sf_merc_split)[,-6]
trial_2 <- trial_merc_df %>% split(.$individ_id)


# Creating 'move' objects that dBBM can work on 
dat.list <- vector("list", length(trial_2))
contours <- vector("list", length(trial_2))



# Creating a for loop for all individuals
# Vector memory exhausted 
for(i in length(dat.list)){
  print(paste("ID", trial_2[[i]]$individ_id[1]))
  dat.move <- move(x = trial_2[[i]]$x, y = trial_2[[i]]$y, time = trial_2[[i]]$timestamp, data = trial_2[[i]],
                   proj = CRS('+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs'),
                   animal = trial_2[[i]]$individ_id)
  x.ext <- diff(dat.move@bbox[1,])
  rast.ext <- ifelse(x.ext < 250, 3, 0.3)
  dat.list[[i]] <- brownian.bridge.dyn(object = dat.move, raster = 0.25, location.error = "location_error", 
                                       margin = 9, window.size = 29, ext = rast.ext)
  res <- raster2contour(dat.list[[i]], levels = c(0.5, 0.95))
  contours[[i]] <- st_as_sf(res)
  
  if(st_geometry_type(contours[[i]])[1] == 'LINESTRING'){
    contours[[i]] <- st_cast(contours[[i]], "POLYGON")
  } else{
    contours[[i]] <- st_cast(contours[[i]], "MULTIPOLYGON") %>%
      st_make_valid()
  }
}

sessionInfo()

