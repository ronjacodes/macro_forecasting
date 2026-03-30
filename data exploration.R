renv::init()

library(jsonlite)

rawdata <- fromJSON("swiss_nowcast_metadata/swiss_nowcast_data.json", 
                    flatten = TRUE)
head(rawdata)
names(rawdata)
summary(rawdata)