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
raw_rainfall <- do.call('rbind', rainfall$rainfall_data)
rainfall_fit <- fourier_model(
  raw_rainfall$rainfall,
  seq(nrow(raw_rainfall)) %% 365 / 365
)

t <- seq(nrow(raw_rainfall)) / 365
profile <- umbrella:::fourier_predict(rainfall_fit, t)

out <- data.frame(
  x=c(t, t),
  y=c(raw_rainfall$rainfall, profile$y),
  label=rep(c('rainfall', 'fit'), each=length(t))
)

# plot
library(ggplot2)
ggplot() + 
  geom_line(data = out, aes(x=x, y=y, color=label, group=label, alpha=.3)) +
  labs(x='t', y='rainfall (mm)', 
       title='daily rainfall and fit between 1992-2019') + 
  theme_classic()

ggsave('./results/rainfall.pdf', width=15, height=10, units="cm")

# print coefficients
write.csv(rainfall_fit$coefficients, './results/rain_coeff.csv', row.names=F)


# now repeat for each admin1 unit
# pull just data for clusters within a single admin unit
# dhs_mini <- st_join(admin1[which(admin1$admin1=="Kedougou"),], dhs, join = st_contains)

