admin0 <- readRDS('M:/Eradication/rds/admin0.RDS')
admin0 <- admin0[which(admin0$ISO=='SEN'), ]
sf::st_as_sf(admin0)
sf::st_crs(admin0) # EPSG:4326

props <- rep(NA, 3)
species <- c('arabiensis', 'funestus', 'gambiae_ss')

for (i in seq_along(species)) {
  s <- species[[i]]
  dat <- terra::rast(paste0('data/vectors/2010_Anopheles_', s, '.tif'))
  dat <- terra::project(dat, "EPSG:4326", method = "bilinear")
  dat <- raster::raster(dat)
  
  props[[i]] <- raster::extract(
    dat,
    admin0,
    method='simple',
    fun=mean,
    na.rm=T,
    weights=F
  )
}

df <- data.frame(
  species = species,
  props = props,
  normalised_props = props / sum(props)
)
