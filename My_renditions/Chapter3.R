library(sf)
library(terra)
library(dplyr)
methods(class="sf")

st_sf(data.frame(n<-world$name_long), g=world$geom)
class(world)
dim(world)
world_df <- st_drop_geometry(world)
class(world_df)
ncol(world_df)
i_small <- world$area_km2 < 10000
summary(i_small)
small_countries <- world[i_small,]
small_countries
names(world)
filter(world, area_km2<10000)
world_new <- world |> filter(continent=="Asia") |> select(name_long,continent) |> slice(1:5)
world_new
world_aggregate <- aggregate(pop~continent, FUN = sum, data = world, na.rm = T)
class(world_aggregate)
world_aggregate
world_aggregate_2 <- aggregate(world["pop"], list(world$continent), FUN = sum, na.rm = T)
world_aggregate_2
world_aggregate_3 <- world |> group_by(continent) |> summarize(population = sum(pop, na.rm = T))
world_aggregate_3
world_aggregate_4 <- world |> group_by(continent) |> summarize(population = sum(pop, na.rm = T), area = sum(area_km2), N=n())
world_aggregate_4
world_aggregate_5 <- world |> st_drop_geometry() |> group_by(continent) |> summarize(Pop = sum(pop, na.rm = T), N = n(), Area = sum(area_km2)) |> mutate(density = Pop/Area) |> slice_max(Pop, n=3) |> arrange(desc(N))
world_aggregate_5

#---- Joining

world_coffee <- left_join(world,coffee_data)
names(world_coffee)
plot(world_coffee$coffee_production_2016)
plot(world_coffee["coffee_production_2016"])
View(world_coffee)

coffee_renamed <- rename(coffee_data, nm = name_long) 
left_join(world, coffee_renamed, join_by (name_long == nm))

world_coffee_inner <- inner_join(world, coffee_data)
nrow(world_coffee_inner)
setdiff(coffee_data$name_long, world$name_long)

drc <- stringr::str_subset(world$name_long,"Dem*.+Congo")
drc
coffee_data$name_long[grepl("Congo,", coffee_data$name_long)] = drc
world_coffee_match <- inner_join(world,coffee_data)
nrow(world_coffee_match)
coffee_world <- inner_join(coffee_data,world)
class(coffee_world)
coffee_world_new <- st_as_sf(coffee_world)
class(coffee_world_new)

new_world <- world
new_world$popdensity <- new_world$pop / new_world$area_km2
View(new_world)
new_world_2 <- world |> mutate(pop_dens = pop / area_km2)
View(new_world_2)
new_world_3 <- world |> transmute(pop_dens = pop/area_km2)
View(new_world_3)

world_unite <- world |> tidyr::unite("conreg", continent:region_un, remove = T) 
View(world_unite)
world_separate <- world_unite |> tidyr::separate(conreg, c("continent","region_un"), sep=":")

world |> rename(name = name_long)
new_names <- c("i", "n", "c", "r", "s", "t", "a", "p", "l", "gP", "geom")
world_new_names <- world |> setNames(new_names)
View(world_new_names)

# ---- Raster objects

elev <- rast(nrow = 6, ncol = 6, xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5, vals = 1:36)
elev

grain_order <- c("clay", "silt", "sand")
grain_char <- sample(grain_order, 36, replace = T)
grain_fact <- factor(grain_char, levels = grain_order)
grain <- rast(nrows = 6, ncols = 6, xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5, vals = grain_fact )
cats(grain)

levels(grain) <- data.frame(value= c(0,1,2), wetness = c("wet", "moist", "dry"))
levels(grain)
cats(grain)
plot(grain)
coltb <- data.frame(value=1:3, col=rainbow(3))
coltb
has.colors(grain)
coltab(grain) <- coltb
plot(grain)

# Raster sub-setting
elev[1,1]
elev[1]
two_layers <- c(grain,elev)
class(two_layers)
two_layers[1]
elev[1,1] = 0
elev[]
two_layers[1] = cbind(c(1), c(4))
two_layers[]

global(elev,sd)
summary(elev)
summary(c(elev,grain))
hist(elev)

# Exercises
library(sf)
library(dplyr)
library(terra)
library(spData)
data("us_states")
data("us_states_df")
us_states_name <- us_states |> select(NAME)
us_states_name
class(us_states_name)
us_states_name_2 <- us_states[,c("NAME", "geometry")]
us_states_name_2
View(us_states)
method_1 <- us_states |> select("total_pop_10", "total_pop_15")
method_2 <- us_states |> select(contains("pop"))
method_3 <- us_states |> select(matches("pop"))

midwest_world <- us_states |> group_by(REGION="Midwest") |> filter(as.numeric(AREA) > 250000 & as.numeric(total_pop_15) > 5000000)
View(midwest_world)
sum_world <- us_states |> group_by(REGION) |> summarize(N = n())
View(sum_world)
sum_world_2 <- us_states |> group_by(REGION) |> 
  summarize(min_2015 = min(total_pop_15), max_2015 = max(total_pop_15), 
            sum_pop_2015 = sum(total_pop_15))
View(sum_world_2)
class(left_join(us_states, us_states_df, join_by(NAME == state)))
us_states_rename <- rename(us_states, state = NAME)
anti_join(us_states_rename, us_states_df) #### Doubt E.8

E9 <- us_states_rename |> 
  group_by(state) |>
  summarize(sum1 = sum(total_pop_10)/AREA, sum2 = sum(total_pop_15)/AREA) |>
  mutate(new = ((sum2-sum1)/sum1) * 100)
plot(E9$new)
