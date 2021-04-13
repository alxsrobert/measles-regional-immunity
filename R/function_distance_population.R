importation_distance <- function(regs, neighbours = F){
  if(neighbours == F){
    ## Set distance data set (each element = distance row-column)
    # Import population centroid file, generated with function_centroid.R
    loc_centroid <- readRDS(file = "Data/location_centroids.RDS")
    # Keep only regions of interest
    loc_centroid <- loc_centroid[is.element(nuts3, regs),]
    loc_centroid <- loc_centroid[!duplicated(rownames(loc_centroid)),]
    # Create data table with ditance between every region
    dt_centroid <- as.data.table(matrix(nrow = length(regs)**2, ncol = 7))
    colnames(dt_centroid) <-  
      c("NUTS3_1", "NUTS3_2", "long1", "lat1", "long2", "lat2", "distance")
    # Import longitude and lattitude of each region
    dt_centroid[, NUTS3_1 := rep(loc_centroid$nuts3, nrow(loc_centroid))]
    dt_centroid[, long1 := rep(loc_centroid$X, nrow(loc_centroid))]
    dt_centroid[, lat1 := rep(loc_centroid$Y, nrow(loc_centroid))]
    dt_centroid[, NUTS3_2 := rep(loc_centroid$nuts3, each = nrow(loc_centroid))]
    dt_centroid[, long2 := rep(loc_centroid$X, each = nrow(loc_centroid))]
    dt_centroid[, lat2 := rep(loc_centroid$Y, each = nrow(loc_centroid))]
    # Generate distance between population centroids
    dt_centroid[, distance := distGeo(p1 = dt_centroid[, .(long1, lat1)],
                                      p2 = dt_centroid[, .(long2, lat2)])/1000]
    # Create the distance matrix
    distance_matrix <- matrix(dt_centroid$distance,
                              ncol = length(unique(dt_centroid$NUTS3_1)))
    colnames(distance_matrix) <- rownames(distance_matrix) <- 
      dt_centroid$NUTS3_2 %>% unique
  } else{
    # Import shapefile of the French NUTS3 region
    temp <- tempfile()
    download.file("http://osm13.openstreetmap.fr/~cquest/openfla/export/departements-20180101-shp.zip",
                  temp)
    
    ## Extract the shapefile into a temporary folder
    unzip(temp, exdir = tempdir())
    ## Read the shapefile and create a map
    map <- sf::st_read(paste0(tempdir(), "/departements-20180101.shp"), 
                       quiet = T)
    map$nuts3 <- as.character(map$nuts3)
    # Fix inconsistencies between the map and nuts3 data
    map[as.character(map$code_insee) == "13", "nuts3"] <- "FR824"
    map[as.character(map$code_insee) == "73", "nuts3"] <- "FR717"
    map[as.character(map$code_insee) == "69D", "nuts3"] <- "FR716"
    map[as.character(map$code_insee) == "74", "nuts3"] <- "FR718"
    map$nuts3 <- as.factor(map$nuts3)
    ## Compute the number of degrees between the regions
    list_deg1 <- st_intersects(map, map)
    names(list_deg1) <- as.character(map$nuts3)
    list_deg1 <- lapply(list_deg1, function(X) return(as.character(map$nuts3)[X]))
    # Initialise the empty distance matrix
    distance_matrix <- matrix(0, ncol = length(regs), nrow = length(regs))
    colnames(distance_matrix) <- rownames(distance_matrix) <- regs
    # for each region, set the distance to its neighbours to 1 
    for(i in colnames(distance_matrix)) {
      distance_matrix[list_deg1[[i]], i] <- 
        distance_matrix[i, list_deg1[[i]]] <- 1
    }
    # Number of degrees between a region and itself is 0
    diag(distance_matrix) <- 0
  }
  
  return(distance_matrix)
}

importation_pop_area <- function(corres, regs, all_dates){
  ## Load age stratified population
  age <- read.csv("Data/age_structure.csv", header = T) %>% 
    as.data.table
  # Initialise age data table containing the region, region id, year, 
  # total number of inhabitant, and surface
  dt_age <- data.table(region = age$region, id = age$number, year = age$year, 
                         total = age$tot, area = age$area)
  # Change the region entry to the Geocode used in all datasets
  dt_age[, region := corres[as.character(dt_age$id),GeoCode]]
  # Select only the regions of interest
  dt_age <- dt_age[is.element(region, regs), ] 
  # Select only the dates of interest
  all_dates <- all_dates[lubridate::year(all_dates) >= min(dt_age$year)]
  # Create the population matrix at each date
  pop_mat <- sapply(as.Date(all_dates), function(X){
    return(dt_age[year == lubridate::year(X), total])
  }) %>% t
  # Create the surface matrix at each date
  area_mat <- sapply(as.Date(all_dates), function(X){
    return(dt_age[year == lubridate::year(X), area])
  }) %>% t
  colnames(pop_mat) <- colnames(area_mat) <- regs
  rownames(pop_mat) <- rownames(area_mat) <- as.character(all_dates)
  
  return(list(pop = pop_mat, area = area_mat))
}
