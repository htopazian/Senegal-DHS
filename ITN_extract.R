# ----------------
# ITN extract
# ----------------

library(sf)
library(tidyverse)
library(raster)
library(terra)


# ITN -----------------
# data from malaria atlas project
# read in admin shapefile from malaria shared folder
admin0 <- readRDS("M:/Eradication/rds/admin0.RDS")
admin0 <- admin0[which(admin0$ISO=='SEN'), ]       # extract Senegal

admin1 <- readRDS("M:/Eradication/rds/admin1.RDS")
admin1 <- admin1[which(admin1$ISO=='SEN'), ]       # extract Senegal

# read in raster files of ITN coverage from malaria atlas project
files <- list.files(path = "./data/ITN_2000_use_mean/", pattern = "*.tif", full.names = TRUE)
dat_list <- lapply(files, function (x) terra::rast(x))

crs(admin0); crs(dat_list[[1]]) # admin 0 is EPSG 4326, ITN is EPSG 9122
# reproject the raster files to EPSG 4326
dat_list <- lapply(dat_list, function (x) terra::project(x, "EPSG:4326", method = "bilinear"))

dat_list <- lapply(files, function (x) raster::raster(x))

# extract mean ITN coverage
rastextract <- function(x){
  admin0 %>%
  mutate(year = files[x],
         year = gsub("./data/ITN_2000_use_mean/ITN_", "", year),
         year = gsub("_use_mean.tif", "", year),
         ITN = raster::extract(dat_list[[x]], admin0, method = "simple", fun = mean, na.rm=T, weights=F))
}

ITNextract <- map_dfr(c(1:length(dat_list)), rastextract)

ITNextract <- ITNextract %>% as_tibble() %>% dplyr::select(year, ITN)


# print coefficients
write.csv(ITNextract, './results/ITN.csv', row.names=F)
