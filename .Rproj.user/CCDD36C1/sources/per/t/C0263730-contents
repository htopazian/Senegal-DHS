# load libraries
library(rdhs)
library(sf)
library(data.table)
library(umbrella) # remotes::install_github('mrc-ide/umbrella@gee', force=T)
library(rgee) # remotes::install_github('r-spatial/rgee', force=T)
library(reticulate)


# DHS --------------------------------------------------------------------------
# log-in
rdhs::set_rdhs_config(
  email = "h.topazian@imperial.ac.uk",
  project = 'Calibration of an Individual-based model of Malaria transmission',
  config_path = "rdhs.json",
  global = FALSE)

# load available surveys
survs <- rdhs::dhs_surveys(countryIds = 'SN')

# select those that meet this criteria
datasets <- rdhs::dhs_datasets(
  surveyIds = survs$SurveyId,
  fileFormat = "FL", # flat flie
  fileType = "PR")   # individual recode

# download datasets
downloads <- rdhs::get_datasets(datasets$FileName)

vars <- c(
  'hhid',  # Case identification number
  'hv001', # Cluster number
  'hv006', # Month
  'hv007', # Year
  'hv024', # Region
  'hv042', # Household selected for hemoglobin
  'hv103', # Slept under ITN last night
  'hc1',   # Child's age in months
  'hml32', # Final result of malaria from blood smear test
  'hml35', # Result of malaria rapid test
  'sh214r', # Lab malaria test result for MIS surveys
  'hv005', # Household sample weight
  'hv023'  # stratification weight
)

questions <- rdhs::search_variables(
  datasets$FileName,
  variables = vars
)


# extract variables of interest and add geographic information
extract <- rdhs::extract_dhs(questions, add_geo = TRUE)

dhs <- data.table::rbindlist(extract[c(1:5,7,9:10,12:13)], fill=T, use.names=T) 
# 6, 8, and 11 have problems with hv023 var - it is not in an int+label format
extract[6]$SNPR61FL$SurveyId[1] # 2010
extract[8]$SNPR70FL$SurveyId[1] # 2014
extract[11]$SNPR7ZFL$SurveyId[1] # 2017

dhs$hv023 <- as.integer(dhs$hv023) # remove label

dhs2 <- data.table::rbindlist(extract[c(6,8,11)], fill=T, use.names=T) # add in

dhs <- merge(dhs, dhs2, deparse.level=0, all = T) # get full dataset

# re-code hml32 smear results
# 0 = negative, 1 = positive, 6 = undetermined, 7 = sample not found in database
table(dhs$hml32, useNA = 'always')
dhs$hml32 <- ifelse(dhs$hml32==7, NA, dhs$hml32)


# plot
# read in admin shapefile from malaria shared folder
admin0 <- readRDS("M:/Eradication/rds/admin0.RDS")
admin0 <- admin0[which(admin0$ISO=='SEN'), ]       # extract Senegal

admin1 <- readRDS("M:/Eradication/rds/admin1.RDS")
admin1 <- admin1[which(admin1$ISO=='SEN'), ]       # extract Senegal

dhsmap <- st_as_sf(dhs[which(dhs$LATNUM!=0),], # removing missing coordinates
                 coords = c("LONGNUM", "LATNUM"), 
                 crs = 4326)

# look at clusters - blue points indicate survey clusters in the 2019 survey
# plot(admin0['ISO'], lwd = 2, col = 'khaki', main = NA)
plot(st_geometry(admin0), lwd = 2, col = 'khaki', main = NA)
plot(st_geometry(admin1[which(admin1$admin1=="Kedougou"),]), lwd = 2, col = 'lightgreen', add = T, density=20)
plot(st_geometry(dhsmap[which(dhsmap$SurveyId=="SN2019DHS"), ]), pch = 3, col = 'cornflowerblue', add = T)


# calculate weighted prevalence
library(survey)

# specify complex survey design: id = cluster number, strata = stratification weight, weights = household sample weight
# sample weights need to be divided by 1000000 before using https://dhsprogram.com/data/Guide-to-DHS-Statistics/Analyzing_DHS_Data.htm
dhs$hv005 <- dhs$hv005 / 1000000 

SID <- unique(dhs$SurveyId)

# find prevalence by year
dhsprev <- function(SID){
  
  # choose one survey
  dhs_mini <- dhs[which((dhs$SurveyId==SID)),]
  
  # create complex survey design
  DHSdesign <- svydesign(id=~hv001, strata=~hv023, weights=~hv005, data=dhs_mini, nest=T) 

  # descriptive statistics
  # without survey.lonely.psu options, function fails b/c of single clusters
  options(survey.lonely.psu="adjust")

  # calculated weighted totals - # who have a positive test and a negative test
  m <- svytotal(~interaction(as.factor(hv006), as.factor(hml32)), DHSdesign, na.rm=T)
  m <- as.data.frame(m)
  m <- data.frame(rep(dhs_mini$hv007[1],nrow(m)), rep(sort(unique(dhs_mini$hv006)),2), c(rep(0,nrow(m)/2),rep(1,nrow(m)/2)), m$total, m$SE)
  names(m) <- c('year', 'month', 'ptest', 'total', 'SE')
  
  return(m)
}

# some survey IDs are missing values
dhs <- dhs[which(!is.na(dhs$hv023)), ] # remove missing hv023 strata values
SID <- unique(dhs$SurveyId)[c(3:4,7:9)]

m <- map_dfr(SID, dhsprev)

# all survey IDs do not work - why? 
# 2008 MIS - MIS survey, use other variable
# 1997 DHS - hml32 all NA
# 2019 DHS - hml32 all NA
# 2018 DHS - hml32 all NA
# 1993 DHS - hml32 all NA
# 2005 DHS - hml32 all NA
# 2006 MIS`- no data for either test variable

dhs_mini <- dhs[which((dhs$SurveyId=='SN2006MIS')),]
table(dhs_mini$hml32, useNA = 'always')

misprev <- function(SID){
 
  dhs_mini <- dhs[which((dhs$SurveyId==SID)),]
  
  DHSdesign <- svydesign(id=~hv001, strata=~hv023, weights=~hv005, data=dhs_mini, nest=T) 
  
  options(survey.lonely.psu="adjust")
  
  m <- svytotal(~interaction(as.factor(hv006), as.factor(sh214r)), DHSdesign, na.rm=T)
  m <- as.data.frame(m)
  m <- data.frame(rep(dhs_mini$hv007[1],nrow(m)), rep(sort(unique(dhs_mini$hv006)),2), c(rep(1, nrow(m)/2), rep(0, nrow(m)/2)), m$total, m$SE) # order reversed for MIS
  names(m) <- c('year', 'month', 'ptest', 'total', 'SE')
  
  return(m)
}

SID <- unique(dhs$SurveyId)[c(1)]
m2 <- map_dfr(SID, misprev)

m3 <- merge(m, m2, all=T)

m3 <- reshape(data=m3,
              idvar=c("year", "month"),
              v.names = c('total','SE'),
              timevar = "ptest",
              direction="wide")

m3$prev <- m3$total.1 / (m3$total.0 + m3$total.1)


# print weighted prevalence estimates
write.csv(m3, 'C:/Users/htopazia/Desktop/prevmonth.csv', row.names=F)

# timing
t <- table(dhs[which(dhs$SurveyId==SID)]$SurveyId, dhs$hv006)
table(dhs[which(dhs$hv024==1)]$hv007, dhs[which(dhs$hv024==1)]$hv006)
write.csv(t, 'C:/Users/htopazia/Desktop/timing.csv', row.names=F)


# RAINFALL ---------------------------------------------------------------------

reticulate::use_python('C:/Users/htopazia/Anaconda3/python.exe', required=T)

# sign into google earth account
ee_Initialize(user='htopazian@gmail.com', drive=T) 

# specify time-range
table(dhs$hv007) # surveys are from 1992 to 2019
start_date <- "2019-01-01" # when started from 1992 it times out
end_date <- "2019-12-31"

start_date <- lapply(seq(1992, 2019, 1), function(x){paste0(x, '-01-01')})
end_date <- lapply(seq(1992, 2019, 1), function(x){paste0(x, '-12-31')})
dates <- data.frame(cbind(start_date, end_date))

library(purrr)
extract_rain <- function(start_date, end_date){
  # extract
  daily_rain <- pull_daily_rainfall(sf = admin0, 
                                    start_date = start_date, end_date = end_date)
  # process
  process(gee_output = daily_rain, ISO, Country)
}

rainfall <- map2_dfr(start_date, end_date, extract_rain)

out <- map_dfr(seq(1, nrow(rainfall), 1), function(x){do.call(data.frame, rainfall[x,]$profile)})
out$year <- ceiling(1:nrow(out)/365)
out2 <- map_dfr(seq(1, nrow(rainfall), 1), function(x){do.call(data.frame, rainfall[x,]$rainfall_data)})
out2$year <- as.numeric(substr(out2$date, 1, 4))
out2$t <- sequence(rle(out2$year)$lengths)

# plot
library(ggplot2)
ggplot() + 
  geom_line(data = out2, aes(x=t, y=rainfall), col = 'red', alpha = 0.3) +
  geom_line(data = out, aes(x=t*365, y=y, group=year), col='blue', alpha = 0.3) + 
  labs(x='days', y='rainfall (mm)', title='daily rainfall (red) and fit (blue), by year, 1992-2019') + 
  theme_classic()

# get coefficients
coeff <- map_dfr(seq(1, nrow(rainfall), 1), function(x){do.call(data.frame, rainfall[x,]$coefficients)})
coeff$year <- seq(1992,2019,1)

# print coefficients
write.csv(coeff, 'C:/Users/htopazia/Desktop/coeff.csv', row.names=F)


# now repeat for each admin1 unit
# pull just data for clusters within a single admin unit
# dhs_mini <- st_join(admin1[which(admin1$admin1=="Kedougou"),], dhs, join = st_contains)


# ITN --------------------------------------------------------------------------



# PLOT IT ALL ------------------------------------------------------------------





