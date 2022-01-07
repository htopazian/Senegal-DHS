# ----------------
# RAINFALL
# ----------------

library(umbrella) # remotes::install_github('mrc-ide/umbrella@gee', force=T)
library(rgee) # remotes::install_github('r-spatial/rgee', force=T)
library(reticulate)
library(purrr)
library(sf)


# RAINFALL ---------------------------------------------------------------------

reticulate::use_python('C:/Users/htopazia/Anaconda3/python.exe', required=T)

# sign into google earth account
ee_Initialize(user='htopazian@gmail.com', drive=T) 

# specify time-range
dhs <- readRDS('./data/dhs_full.rds')
table(dhs$hv007) # surveys are from 1992 to 2019

# when run from 1992 to 2019 it times out. Take yearly estimates
start_date <- lapply(seq(1992, 2019, 1), function(x){paste0(x, '-01-01')})
end_date <- lapply(seq(1992, 2019, 1), function(x){paste0(x, '-12-31')})
dates <- data.frame(cbind(start_date, end_date))


# read in admin shapefile from malaria shared folder
admin0 <- readRDS("M:/Eradication/rds/admin0.RDS")
admin0 <- admin0[which(admin0$ISO=='SEN'), ]       # extract Senegal

admin1 <- readRDS("M:/Eradication/rds/admin1.RDS")
admin1 <- admin1[which(admin1$ISO=='SEN'), ]       # extract Senegal


# extract rainfall function
extract_rain <- function(start_date, end_date){
  # extract
  daily_rain <- pull_daily_rainfall(sf = admin0, 
                                    start_date = start_date, end_date = end_date)
  # process
  process(gee_output = daily_rain, ISO, Country)
}

# extract
rainfall <- map2_dfr(start_date, end_date, extract_rain)

# pull out fit profile
out <- map_dfr(seq(1, nrow(rainfall), 1), function(x){do.call(data.frame, rainfall[x,]$profile)})
out$year <- ceiling(1:nrow(out)/365)

# pull out rainfall data
out2 <- map_dfr(seq(1, nrow(rainfall), 1), function(x){do.call(data.frame, rainfall[x,]$rainfall_data)})
out2$year <- as.numeric(substr(out2$date, 1, 4))
out2$t <- sequence(rle(out2$year)$lengths)

# plot
library(ggplot2)
ggplot() + 
  geom_line(data = out2, aes(x=t, y=rainfall), col = 'red', alpha = 0.3) +
  geom_line(data = out, aes(x=t*365, y=y, group=year), col='blue', alpha = 0.3) + 
  labs(x='days', y='rainfall (mm)', 
       title='daily rainfall (red) and fit (blue), by year, 1992-2019') + 
  theme_classic()

ggsave('./results/rainfall.pdf', width=15, height=10, units="cm")


# get coefficients
coeff <- map_dfr(seq(1, nrow(rainfall), 1), function(x){do.call(data.frame, rainfall[x,]$coefficients)})
coeff$year <- seq(1992,2019,1)

# print coefficients
write.csv(coeff, './results/rain_coeff.csv', row.names=F)


# now repeat for each admin1 unit
# pull just data for clusters within a single admin unit
# dhs_mini <- st_join(admin1[which(admin1$admin1=="Kedougou"),], dhs, join = st_contains)

