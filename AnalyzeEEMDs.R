#####
### THIS CODE WILL ANALYZE THE SIGNFICANCE OF EACH IMF FOLLOWING WU & HUANG 2004
#####

#####
### Requires Input!
#####
# wd <- "<PATH/TO/FOLDER/HERE/SignificanceCodeToShare>"
wd <- "/Users/cac118/OneDrive - Duke University/FFTvHHT/code/SignificanceCodeToShare"

dataset <- "Ichetucknee_compiled.RData" #options are "whitenoise_compiled.RData", "Ichetucknee_compiled.RData", or "Potomac_compiled.RData"

#####
### Setup
#####
library(gridExtra)
library(zoo)
library(foreach)
library(dplyr)
library(xts)
source(paste0(wd, "/helpfulfunctions.R"))

if(dataset == "whitenoise_compiled.RData") {
  set.seed(64) ### See "MakeWhiteNoise.R"
  whitenoise <- rnorm(1000, 0, 1) 
  DateTime_EST <- seq.POSIXt(strptime("2010-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST"), strptime("2011-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST"), "hour")[1:1000]
  seconds <- as.numeric(DateTime_EST) - as.numeric(head(DateTime_EST, 1)) # how many seconds from the beginning of the series
  minutes <- seconds / 60
  hours <- minutes / 60
  days <- hours / 24
  weeks <- days / 7
  years <- days / 365
  
  timeseries <- data.frame(DateTime_EST, seconds, minutes, hours, days, weeks, years, 
                           whitenoise)
  } else if(dataset == "Potomac_compiled.RData") {
  load(file = paste0(wd, "/01646500_Nitrate.RData"))
  attr(dat$DateTime_UTC, "tzone") <- "EST"; names(dat)[1] <- "DateTime_EST"
  
  timeseries <- dat %>%
    filter(DateTime_EST > "2015-01-01 00:00:00" & DateTime_EST < "2016-01-01 00:00:00") %>%
    mutate(seconds = as.numeric(DateTime_EST) - as.numeric(head(DateTime_EST, 1)),
           minutes = seconds / 60,
           hours = minutes / 60,
           days = hours / 24,
           weeks = days / 7,
           years = days / 365) %>%
    filter(!is.na(DateTime_EST)) %>%
    mutate(Potomac = na.approx(Nitrate_mgL)) %>%
    select(DateTime_EST, seconds, minutes, hours, days, weeks, years, Potomac)
  } else if(dataset == "Ichetucknee_compiled.RData") {
    timeseries <- read.csv(paste0(wd, "/Ich_short.csv")) %>%
      mutate(DateTime_EST = as.POSIXct(strptime(Date, format = "%m/%d/%Y %H:%M", tz = "EST")), 
             seconds = as.numeric(DateTime_EST) - as.numeric(head(DateTime_EST, 1)),
             minutes = seconds / 60,
             hours = minutes / 60,
             days = hours / 24,
             weeks = days / 7,
             years = days / 365) %>%
      filter(!is.na(DateTime_EST)) %>%
      mutate(Ichetucknee = na.approx(NO3.mg.L)) %>%
      select(DateTime_EST, seconds, minutes, hours, days, weeks, years, Ichetucknee)
  }

#####
### View Results
#####
load(paste0(wd, "/", dataset))

my_PlotIMFs(named.result[2:15], main = named.result[1]) # I just changed this slightly to add the title

investigate.IMF(6, named.result, timeseries$DateTime_EST) # Look closer at individual IMFs


#####
### Find mean periods of IMFs
#####

### "periods" is roughly equivalent to Table 1 in Wu & Huang
periods <- as.data.frame(foreach(j = 1:ncol(named.result$imf), .combine = "rbind", .packages = c("pracma", "dplyr")) %do% {
  x <- timeseries$DateTime_EST
  y <- named.result$imf[,j]
  
  peaks.places <- as.data.frame(findpeaks(y, zero = '+', nups = 1))  %>%
    mutate(Height = V1, 
           time_max = x[V2], 
           time_begin = x[V3], 
           time_end = x[V4]) %>%
    select(Height, time_max)
  
  times <- peaks.places$time_max
  
  periods <- difftime(times, c(head(times, 1), times[-length(times)]), units = "secs")
  mean.period_hr <- mean(periods) / (60 * 60)
  var.period_hr <- sd(periods) / (60 * 60)
  median.period_hr <- median(periods) / (60 * 60)
  range.period_hr <- (max(periods) - min(periods)) / (60 * 60)
  npeaks <- length(times)
  
  c(j, mean.period_hr, var.period_hr, median.period_hr, range.period_hr, npeaks)
}); names(periods) <- c("IMF", "MeanPeriod_hrs", "VariancePeriod_hrs", "MedianPeriod_hrs", "RangePeriod_hrs", "NumberPeaks")
View(periods)

### This code produces histograms of the periods of each IMF
period.hists <- foreach(j = 1:ncol(named.result$imf), .packages = c("pracma", "ggplot2")) %do% {
  x <- timeseries$DateTime_EST
  y <- named.result$imf[,j]
  
  peaks.places <- as.data.frame(findpeaks(y, zero = '+', nups = 1))  %>%
    mutate(Height = V1, 
           time_max = x[V2], 
           time_begin = x[V3], 
           time_end = x[V4]) %>%
    select(Height, time_max)
  
  times <- peaks.places$time_max
  
  periods <- as.data.frame(difftime(times, c(head(times, 1), times[-length(times)]), units = "secs")); names(periods) <- "Period"
  
  ggplot(periods, 
         aes(Period / (60*60))) + 
    geom_histogram(color = "black", fill = "grey80") +
    labs(y = "Count", x = "Period (hrs)") +
    theme_classic(base_size = 12) +
    ggtitle(paste(named.result[1], "IMF", j))
  
}
grid.arrange(grobs = period.hists)

#####
### Significance Test
#####

### This code tries to reproduce Figure 5 in Wu & Huang 2004.
peaks <- as.data.frame(foreach(j = 1:ncol(named.result$imf), .packages = c("dplyr", "pracma", "zoo"), .combine = "rbind") %do% {
  ts <- named.result$imf[,j] #
  spectra <- spectrum(ts, plot = F, detrend = T) #Does a fourier transform of the IMF
  ### Next three are the results of spectra that we need
  period <- (1 / spectra$freq)   # period
  fdensity <- spectra$spec #spectral density
  freq <- spectra$freq # frequency
  
  dlnperiod <- c(diff(rev(log(period))), 0) #This is the dlnT that is used in equation 2.4. It is in reverse to have positive numbers.
  integral <- sum(dlnperiod * rev(fdensity))# This is equation 2.4. The spectral density ("Fourier spectrum") is reversed to line up with the lnT that is also reversed.
  
  dfreq <- c(diff(freq), 0) # This is domega that is in equation 2.5
  energy <- sum(ts^2)/length(ts) # this is equation 2.2. C_n(j) I believe is the raw IMF, and N is the number of frequencies in the fourier spectrum
  # energy <- sum(fdensity * dfreq) / length(fdensity) # This could be an alternative way of getting energy based on the most lefthand part of equation 2.5. S_omega,n is the fourier specrum, and d_omega is the difference in frequency

  meanperiod <- integral / (sum(rev(fdensity/period) * dlnperiod)) #This is equation 2.6 to get the mean period of the IMF

  c(j, meanperiod, energy)

}); names(peaks) <- c("IMF", "Period", "Energy")
# grid.arrange(grobs = peaks)

peaks$Energy <- peaks$Energy * (1 / sum(peaks$Energy)) # This normalizes the energy so that it always equals one. I think this is kosher, it definitely does not work well without it, but would love feedback.

### fun.1 and fun.2 are the upper and lower boundaries given by equation 3.9 w/ 99th percentiles. I believe this is valid so long as you are interested in periods roughly equivalent to the mean period. So probably only for IMF 3-6 perhaps?
### for these functions I think that N = the number of IMFs, but maybe it equals the number of frequencies in the fourier transform?
fun.1 <- function(x) {
  k <- 0.675
  -x + k * sqrt(2 / ncol(named.result$imf) * exp(x / 2))
}

fun.2 <- function(x) {
  k <- 0.675
  -x - k * sqrt(2 / ncol(named.result$imf) * exp(x / 2))
}

ggplot(peaks,
       aes(x = log(Period),
           y =log(Energy),
           color = as.factor(IMF))) +
  geom_point()  +
  stat_function(fun = fun.1) +
  stat_function(fun = fun.2)

