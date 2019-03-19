#####
### THIS CODE WILL CREATE ONE WHITENOISE SERIES OF LENGTH 1000. THE OUTPUT OF THIS CODE IS ALREADY PROVIDED AS "whitenoise_compiled.RData"
#####

#####
### Requires Input!
#####

#####
### Setup
#####
library(ggplot2)
library(hht)
set.seed(64)

#####
### Create white noise dataseries
#####

### This code creates a white noise series 1000 observations in length. I have not done a Monte Carlo of this for time reasons.
### I've also attached timestamps to it to make the data/code easily transferable to my real datasets

whitenoise <- rnorm(1000, 0, 1) #MB: Yes, white noise is a normal distribution with mean of 0, the noisy data has to be normalized (sd = 1)
Sys.setenv(TZ='America/Montreal') #Define system time zone
DateTime_EST <- seq.POSIXt(strptime("2010-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST"), strptime("2011-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST"), "hour")[1:1000]
### the following vectors express the time into the series. This is helpful to quantify periods and frequencies in an understandable way (e.g. to find peaks that occur at frequencies of 365 times per year)
seconds <- as.numeric(DateTime_EST) - as.numeric(DateTime_EST[1]) # how many seconds from the beginning of the series
minutes <- seconds / 60
hours <- minutes / 60
days <- hours / 24
weeks <- days / 7
years <- days / 365

timeseries <- data.frame(DateTime_EST, seconds, minutes, hours, days, weeks, years, 
                         whitenoise)
#MB: Not sure why line 27 to 34 is necessary if what you just computed is not in you final output (named.result)
#I thought the vector for EEMD had to be a time serie object (using ts()), for example use code below to have 
#a series that start Jan 1st 2010 and has an hourly time step
#whitenoise <- ts(whitenoise, start = c(as.numeric(strftime("2010-01-01", "%j")),1), frequency = 24)
#plot(whitenoise, main = 'White noise', xlab = 'Day of the year')

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
  ### I saved the output so you can skip this step, MB: I did not run the code below (but did read it)
  trials <- 1000 #MB: Du et al. 2019 and Bekka and Berrouche 2013 suggest methods to find optimal no of trials (and noise amplitude)
  nimf <- 10
  trials.dir <- paste0(title, "_eemd")
  noise.amp <- 0.1 #MB:I think recommandation from Wu and Huang 2004 is 0.2 
  
  eemd.result <- EEMD(timeseries[,8], timeseries$years, trials = trials, noise.amp = noise.amp, nimf = nimf, trials.dir = trials.dir) #I changed timeseries[,i], to timeseries[,8]  
  resift.rule <- "max.var"
  #MB: Before running EEMDResift I run : EEMDCompile("trials_directory", no of trial, no of IMF)
  resift.result <- EEMDResift(eemd.result, resift.rule)
  named.result <- c(title, resift.result)

  save(named.result)
