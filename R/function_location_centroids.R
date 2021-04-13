function_centroids <- function(){
  # Import population grid file
  pop <- sf::st_read("https://gisco-services.ec.europa.eu/grid/grid_2km_surf.gpkg") 
  # Set to data table format
  pop_test <- as.data.table(pop)[, .(X_LLC, Y_LLC, TOT_P_2011)]
  # Initialise raster
  z_init <- rasterFromXYZ(pop_test, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
  
  # Import Shapefile map of French departements
  temp <- tempfile()
  download.file("http://osm13.openstreetmap.fr/~cquest/openfla/export/departements-20180101-shp.zip",
                temp)
  
  ## Extract the shapefile into a temporary folder
  unzip(temp, exdir = tempdir())
  unlink(temp)
  ## Read the shapefile 
  map <- sf::st_read(paste0(tempdir(), "/departements-20180101.shp"),
                     quiet = T)
  map$nuts3 <- as.character(map$nuts3)
  # Fix inconsistencies between the map and nuts3 data
  map[as.character(map$code_insee) == "13", "nuts3"] <- "FR824"
  map[as.character(map$code_insee) == "73", "nuts3"] <- "FR717"
  map[as.character(map$code_insee) == "69D", "nuts3"] <- "FR716"
  map[as.character(map$code_insee) == "74", "nuts3"] <- "FR718"
  map$nuts3 <- as.factor(map$nuts3)
  #merge units as desired by group variable
  st_union_by = function(geo, group) {
    y2 <- list()
    #loop over by groups and merge units
    for (i in unique(group)) {
      #which units
      z = geo[group == i]
      #merge
      y = Reduce(st_union, z)
      y2[[i]] = y
    }
    st_sfc(y2)
  }
  map$group <- 1
  # Group by nuts3 region
  map[!duplicated(map$nuts3),]$group <- seq_along(unique(map$nuts3))
  map[duplicated(map$nuts3),]$group <- 
    map[!duplicated(map$nuts3) & !is.na(map$nuts3) &
          map$nuts3 == map[duplicated(map$nuts3),]$nuts3,]$group
  # Remove duplicates
  map <- map[!duplicated(map$nuts3),]
  map1 <- map
  map1$geometry <- st_union_by(map$geometry, map$group)
  map[map$group == 64, ]$geometry <- map1[map1$group == 64, ]$geometry
  # Convert coordinate and rasterize
  map <- st_transform(map, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
  z <- rasterize(map, z_init)
  ## Compute weighted x and y coordinates within each rasterized region
  xx <- zonal(init(z_init, v="x")*z_init, z) / zonal(z_init,z)
  yy <- zonal(init(z_init, v="y")*z_init, z) / zonal(z_init,z)
  ## Combine results in a matrix
  res <- cbind(xx[,2],yy[,2]) %>% as.data.table
  res[, nuts3 := z@data@attributes[[1]]$nuts3]
  colnames(res) <- c("long", "lat", "nuts3")
  # Remove na values
  res <- res[!is.na(res$long),]
  # Convert to an sf object
  sf_res <- st_as_sf(res, coords = c(1,2), crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
  # Convert coordinates
  res <- st_transform(sf_res, crs = "+proj=longlat +datum=WGS84 +no_defs")
  # Generate and save the final coordinates of the centroids
  final_res <- st_coordinates(res) %>% as.data.table
  final_res[, nuts3 := res$nuts3]
  saveRDS(final_res, file = "Data/location_centroids.RDS")
  
}