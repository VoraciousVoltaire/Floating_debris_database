# Figuring out the deal with Alkefjellet_10

# Loading data----
datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/input_data/"
mylocs <- readRDS(paste0(datadir, "SEATRACK_FUGLA_20220307_v2.3_FA.rds"))
summary_info <- readRDS(paste0(datadir,"summaryTable.rds"))

# Creating an sf object----

indiv_merged_df <- merge(mylocs, summary_info, by = "ring") 
nbs_mylocs <- indiv_merged_df %>% 
  filter(!grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp)) %>%
  dplyr::mutate(year = year(timestamp))
names(nbs_mylocs)[names(nbs_mylocs) == "ring"] <- "individ_id"
nbs_mylocs_sf <- nbs_mylocs %>% 
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) 

# Assorting by months----

all_data <- nbs_mylocs
df <- all_data[!is.na(all_data$timestamp), ]
df <- df %>% dplyr::mutate(month = month(timestamp))
months <- sort(unique(df$month))

dir_kernels <- dir_kernels <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/outputs/Alkefjellet_10/"
proj_wgs84 <- sp::CRS(sp::proj4string(land))
land <- as(world, "Spatial")

  sub <- df %>% filter(colony == "Alkefjellet")
  
    tracks_wgs <- sub[sub$month == 10,]
    
    if(nrow(tracks_wgs) > 4){
      
      if(min(tracks_wgs$lon) <= -179 ){ lon_min <- -180
      } else {lon_min <- floor(min(tracks_wgs$lon))-1 }
      
      if(max(tracks_wgs$lon) >= 179){ lon_max <- 180
      } else { lon_max <- ceiling(max(tracks_wgs$lon))+1 }
      
      if(min(tracks_wgs$lat) <= -89 ){ lat_min <- -90 
      } else { lat_min <- floor(min(tracks_wgs$lat))-1 }
      
      if(max(tracks_wgs$lat) >= 89){ lat_max <- 90
      } else { lat_max <- ceiling(max(tracks_wgs$lat))+1 }
      
      so.grid <- expand.grid( LON = seq(lon_min, lon_max, by=1), 
                              LAT = seq(lat_min, lat_max, by=1))
      
      sp::coordinates(so.grid) <- ~LON+LAT
      sp::proj4string(so.grid) <- sp::proj4string(land)
      
      mean_loc <- geosphere::geomean(cbind(tracks_wgs$lon,tracks_wgs$lat))
      DgProj <- sp::CRS(paste0("+proj=laea +lon_0=",mean_loc[1],
                               " +lat_0=",mean_loc[2])) 
      
      so.grid.proj <- sp::spTransform(so.grid, CRS = DgProj)
      coords <- so.grid.proj@coords
      
      c <- min(coords[,1])-1000000   ## to check my min lon
      d <- max(coords[,1])+1000000   ## to check my max lon
      
      e <- min(coords[,2])-1000000   ## to check my min lat
      f <- max(coords[,2])+1000000   ## to check my max lat
      
      a <- seq(c, d, by=10000)
      b <- seq(e, f, by=10000)
      null.grid <- expand.grid(x = a,y = b)
      sp::coordinates(null.grid) <- ~x+y
      sp::gridded(null.grid) <- TRUE
      
      sp::coordinates(tracks_wgs) <- ~lon+lat
      sp::proj4string(tracks_wgs) <- sp::proj4string(land)
      tracks <- sp::spTransform(tracks_wgs, CRS = DgProj)
      
      tracks$month <- factor(tracks@data$month)
      
      kudl <- adehabitatHR::kernelUD(tracks[,"month"], 
                                     grid = null.grid, h = 200000)  ## smoothing factor equals 200 km for GLS data
      
      vud <- adehabitatHR::getvolumeUD(kudl)
      
      fud <- vud[[1]]
      hr95 <- as.data.frame(fud)[,1]
      hr95 <- as.data.frame(fud)[,1]
      hr95 <- data.frame(hr95)
      coordinates(hr95) <- coordinates(fud)
      sp::gridded(hr95) <- TRUE
      
      kde_spixdf <- adehabitatHR::estUDm2spixdf(kudl)
      kern95 <- kde_spixdf
      
      stk_100 <- raster::stack(kern95)
      stk_95 <- raster::stack(hr95)
      
      sum_all_100 <- stk_100[[1]]
      sum_all_95 <- stk_95[[1]]
      
      sum_all_raw <- sum_all_100*sum_all_95
      
      rast <- sum_all_raw/sum(raster::getValues(sum_all_raw))
      rast[rast == 0] <- NA
      
      x.matrix <- is.na(as.matrix(rast))
      colNotNA <- which(colSums(x.matrix) != nrow(rast))
      rowNotNA <- which(rowSums(x.matrix) != ncol(rast))
      
      croppedExtent <- raster::extent(rast, 
                                      r1 = rowNotNA[1]-2, 
                                      r2 = rowNotNA[length(rowNotNA)]+2,
                                      c1 = colNotNA[1]-2, 
                                      c2 = colNotNA[length(colNotNA)]+2)
      
      cropped <- raster::crop(rast, croppedExtent)
      cropped[is.na(cropped)] <- 0
      
      
      mask_proj <- sp::spTransform(land, DgProj)   ## changing projection
      mask_proj_pol <- as(mask_proj, "SpatialPolygons")   ## converting SpatialPolygonsDataFrame to SpatialPolygons
      plot(cropped)
      plot(mask_proj_pol, add = T)
      
      # The error lies here; can't figure out why 
      rast_mask_na <- raster::mask(cropped, mask = mask_proj_pol, inverse = T) 
      plot(rast_mask_na)
      
      rast_mask <- rast_mask_na
      rast_mask[is.na(rast_mask)] <- 0
      rast_mask_sum1 <- rast_mask/sum(raster::getValues(rast_mask))
      ?raster::mask
      #rast_mask[rast_mask == 0] <- NA
      rast_mask_final <- raster::mask(rast_mask_sum1, mask_proj_pol, inverse = T)
      rast_mask_final
      rast_mask_final2 <- rast_mask_final 
      
      mask_wgs84 <- raster::projectRaster(rast_mask_final2, crs = proj_wgs84, over = F)
       
      raster::writeRaster(mask_wgs84, filename = paste0(dir_kernels,"Alkefjellet_10.tif"), # you could add month_number here to get KDE for each month; I still don't get what they aggregated those months
                          format = "GTiff", overwrite = TRUE)
      
      png(filename = paste0(dir_kernels,"Alkefjellet_10.png"))
      plot(mask_wgs84, main = paste0("",i))
      plot(land, add = T, col = "#66000000")
      dev.off()
      
    }
  
    
    
    # Raster multiplication---- (next step after the issue is resolved)
    
    dir_demClasses <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/outputs/Alkefjellet_10/"
    dir_1by1 <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/outputs/Alkefjellet_10/" # edit later
    setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/")
    plastics <- raster::raster("00_PlasticsRaster.tif")
    setwd(paste0(dir_demClasses))
    
    r <- raster()
    r
    
    yelblus <- c(brewer.pal(n = 9, name = "YlGnBu"),"#00172e")
    cols <- colorRampPalette(yelblus)(255)
    colsviri <- cols[20:255]
    
    colsinf <- rev(inferno(200))
    
    plastics2 <- plastics
    plastics2[is.na(plastics2)] <- 0 
    p_sum1    <- plastics2/sum(raster::getValues(plastics2))
    p_sum1[is.na(plastics)] <- NA
    
    RES <- res(plastics) 
    R <- 6371007.2 
    lat <- raster::yFromRow(plastics, 1:nrow(plastics)) # latitude of the centroid of each cell (in degrees, need to be converted in radians)
    area <- (sin(pi/180*(lat + RES[2]/2)) - sin(pi/180*(lat - RES[2]/2))) * (RES[1] * pi/180) * R^2
    r_area <- raster::setValues(plastics, rep(area, each=ncol(plastics))) # gives the area of each grid cell in meters 
    plot(r_area, col=colsviri)
    
    files <- list.files(dir_demClasses, full.names = TRUE, pattern=".*\\.tif$"); files
    
    nas <- as.data.frame(files)
    nas$name <- NA
    nas$vals <- NA
    nas$nas <- NA
    nas$files <- NULL
    
    # colony_exposure_score_vector <- vector()
    # exposure_scores_df <- data.frame(matrix(NA, nrow = 12))
    # cbind.fill <- function(...){
    #   nm <- list(...) 
    #   nm <- lapply(nm, as.matrix)
    #   n <- max(sapply(nm, nrow)) 
    #   do.call(cbind, lapply(nm, function (x) 
    #     rbind(x, matrix(, n-nrow(x), ncol(x))))) 
    
    
    dat <- data.frame()
    result <- c()
    
    for (i in 1:length(files)){
      
      a <- raster(files[i])
      name <- a@data@names[1]
      nas$name[i] <- name
      a[is.na(a)] <- 0 
      b <- sum(raster::getValues(a)) 
      # b: why is this defined 
      a_proj <- raster::projectRaster(a, plastics, method = "bilinear")
      a_proj2 <- a_proj * r_area / 100000000 # rescaling the values in each cell
      a_proj2[is.na(a_proj2)] <- 0 
      c <- sum(raster::getValues(a_proj2)) 
      # c: why is this defined
      na_a_proj2 <- a_proj2
      na_a_proj2[na_a_proj2 > 0] <- 1
      na_a_proj2[is.na(na_a_proj2)] <- 0
      nas$vals[i] <- sum(na_a_proj2@data@values)
      
      na_p_sum1 <- p_sum1
      na_p_sum1[is.na(na_p_sum1)] <- 1
      na_p_sum1[na_p_sum1 != 1] <- 0
      na_over <- na_a_proj2 * na_p_sum1
      nas$nas[i] <- sum(na_over@data@values) 
      
      a_proj2[is.na(plastics)] <- NA
      
      raster_name_1 <- gsub(dir_demClasses, "", files[i])
      print(raster_name_1) 
      
      raster_name_2 <- paste0(dir_1by1, raster_name_1)
      
      raster::writeRaster(a_proj2, filename=raster_name_2, format="GTiff", overwrite=TRUE)
      
      over <- a_proj2 * p_sum1
      
      over_score <- over
      summary(over_score@data@values)
      over_score[is.na(over_score)] <- 0
      exposure_score <- round(sum(raster::getValues(over_score))*1000000,4)
      
      png(paste0(dir_1by1,"/maps/",name,".png"), width=1399,height=455)
      par(mfrow=c(1,2))
      plot(a_proj2,main=paste0(name," distribution"),col=colsviri,legend=F)
      plot(over,main=paste0("Exposure score = ", exposure_score),
           col=colsinf,legend=F)
      
      dev.off()
      
      if(sum(na_over@data@values)> 0){
        png(paste0(dir_1by1,"/na_maps/",name,".png"), width=1399,height=455)
        par(mfrow=c(1,2))
        plot(a_proj2,main=name,col=colsviri,legend=F)
        plot(na_over,main=paste0("Exposure Nas = ",sum(na_over@data@values)),
             col=colsinf,legend=F)
        dev.off()
      }
      
      name_split <- strsplit(name,"_")[[1]]
      month <- name_split[length(name_split)]
      population <- paste(name_split[1])
      
      check <- cbind(b,c)
      result <- cbind(population, month, exposure_score)
      dat <- rbind(dat, as.data.frame(result))
      print(i)
    }
    
    write.csv(dat, "exposure_scores_by_month.csv",
              row.names = F)  
    head(nas)
    nas$percent_na <- nas$nas/nas$vals*100
    
    #2 don't exist
    nas_no_na <- subset(nas,vals > 0)
    mean(nas_no_na$percent_na)
