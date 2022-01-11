# ----------------
# ITN
# ----------------

library(sf)
library(tidyverse)


# DHS DOWNLOAD -----------------


# read in admin shapefile from malaria shared folder
admin0 <- readRDS("M:/Eradication/rds/admin0.RDS")
admin0 <- admin0[which(admin0$ISO=='SEN'), ]       # extract Senegal

admin1 <- readRDS("M:/Eradication/rds/admin1.RDS")
admin1 <- admin1[which(admin1$ISO=='SEN'), ]       # extract Senegal