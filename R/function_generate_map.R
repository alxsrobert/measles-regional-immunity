generate_map <- function(){
  temp <- tempfile()
  download.file("http://osm13.openstreetmap.fr/~cquest/openfla/export/departements-20180101-shp.zip",
                temp)
  
  ## Extract the shapefile into a temporary folder
  unzip(temp, exdir = tempdir())
  ## Read the shapefile and create one map for each model
  map <- sf::st_read(paste0(tempdir(), "/departements-20180101.shp"),
                      quiet = TRUE)
  map$nuts3 <- as.character(map$nuts3)
  map[as.character(map$code_insee) == "13", "nuts3"] <- "FR824"
  map[as.character(map$code_insee) == "73", "nuts3"] <- "FR717"
  map[as.character(map$code_insee) == "69D", "nuts3"] <- "FR716"
  map[as.character(map$code_insee) == "74", "nuts3"] <- "FR718"
  map$nuts3 <- as.factor(map$nuts3)
  
  ## Delete the temporary file
  unlink(temp)
  return(map)
}
