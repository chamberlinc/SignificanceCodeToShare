#####
### THIS CODE WILL CREATE ONE WHITENOISE SERIES OF LENGTH 1000. THE OUTPUT OF THIS CODE IS ALREADY PROVIDED AS "whitenoise_compiled.RData"
#####

#####
### Requires Input!
#####
wd <- "<PATH/TO/FOLDER/HERE/SignificanceCodeToShare>"
# wd <- "/Users/cac118/OneDrive - Duke University/FFTvHHT/code/SignificanceCodeToShare"

#####
### Setup
#####
library(ggplot2)
library(hht)
set.seed(64)
setwd(wd)

#####
### Create white noise dataseries
#####

### This code creates a white noise series 1000 observations in length. I have not done a Monte Carlo of this for time reasons.
### I've also attached timestamps to it to make the data/code easily transferable to my real datasets

whitenoise <- rnorm(1000, 0, 1) # White noise is a normal distribution right?
DateTime_EST <- seq.POSIXt(strptime("2010-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST"), strptime("2011-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST"), "hour")[1:1000]
### the following vectors express the time into the series. This is helpful to quantify periods and frequencies in an understandable way (e.g. to find peaks that occur at frequencies of 365 times per year)
seconds <- as.numeric(DateTime_EST) - as.numeric(head(DateTime_EST, 1)) # how many seconds from the beginning of the series
minutes <- seconds / 60
hours <- minutes / 60
days <- hours / 24
weeks <- days / 7
years <- days / 365

timeseries <- data.frame(DateTime_EST, seconds, minutes, hours, days, weeks, years, 
                         whitenoise)
#####
### quick view all data series
#####
  title <- "whitenoise"

  ggplot(timeseries,
         aes_string(x = "DateTime_EST",
                    y = title)) +
    geom_line() +
    labs(y = "Signal", x = "Date") +
    theme_classic(base_size = 12) +
    ggtitle(title)
  
  spectrum(timeseries$whitenoise) #shows the FFT spectrum


#####
### Make EEMDs
#####

### Takes ~ 10 minutes
  ### I saved the output so you can skip this step
   trials <- 1000 #This is going to take a while!
  nimf <- 10
  trials.dir <- paste0(title, "_eemd")
  noise.amp <- 0.1
  
  eemd.result <- EEMD(timeseries[,i], timeseries$years, trials = trials, noise.amp = noise.amp, nimf = nimf, trials.dir = trials.dir)
  eemd.result <- EEMDCompile(trials.dir, trials, nimf)
  resift.rule <- "max.var"
  resift.result <- EEMDResift(eemd.result, resift.rule)
  named.result <- c(title, resift.result)

  save(named.result, paste0(wd, "/whitenoise_compiled.RData"))
