
###############################################
# METHOD 1 SCRIPT FOR RUNNING ON A SINGLE SITE 
###############################################

#run installPackages.R script first

library("vegan")
library("rioja")
library("plyr")
library("analogue")

# Sort the directory

mainDir <- "~/EcoRe3"
# subDir <- "ForRes"
# dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
# setwd(file.path(mainDir, subDir))
setwd(mainDir)
source("ecore3_functions.R")


###############################################
# Script
###############################################

dataFile <- "TAS1110_POLLEN_NEW.csv"
dataDir <- "data"
fileToLoad <- file.path(mainDir, dataDir, dataFile)

poll <- read.csv2(fileToLoad, header = TRUE)
head(poll)
ages <- poll$AGE
poll <- poll[,-1]
siteObj <- list(core = poll, ages = ages, dataset = "TAS1110")
rm(poll, ages)

# Plot the data
pollDat <- data.frame(age = siteObj$ages, chooseTaxa(siteObj$core, n.occ = 5, max.abun = 5))
Stratiplot(age ~ ., data = pollDat, 	type = c("h", "l", "g"), sort = "wa")

# Calculate the BC distances
polBC <- calcBC(site= siteObj )

# Make the BC distance plot
plotBC(polBC)

# Choose how many perturbations you want and find out where they are in the record
distEvents <- extractDist.perc(calcBC.obj = polBC, perc = 0.95, BCselect = 2)

# find out which events are in consecutive samples and only select non-consecutive ones
dEvents <- c(diff(distEvents$sample), 1000) # the 1000 is just an arbitrary number added to make the calculation work so the last sample can be included
sglEvnts <- distEvents[dEvents > 2,]
# write.csv(sglEvnts, file = "sglEvnts.csv")

# Plot the full PCA and the full dataset to check the disturbance points make sense
plotCore(site = siteObj, evntObj = sglEvnts)
plotCore(site = siteObj, evntObj = sglEvnts, strat= FALSE)

# Then you can inspect turnover using PCA within a given window and also define manually the time points for the return time/recovery rate (e.g. by inspection of dist1 or counting the number of samples)
resInteractive(site = siteObj, window = 500, dist = sglEvnts$ages[1], method = "BC")

# The above interactive code is essentially a wrapper function for the following
# dist1 <- localOrd(site = siteObj, dist = sglEvnts$ages[1])
# localOrd.plot(dist1)
# dist1.re3 <- re3(dist1, bottom = 1597, end = 1544, site = siteObj, method = "BC")


################# Then needs to be repeated for all other disturbance events in the 'sglEvnts' object

dist1 <- resInteractive(site = siteObj, window = 1200, dist = sglEvnts$ages[1], method = "BC")
1597
1544
dev.off()
dist2 <- resInteractive(site = siteObj, window = 1200, dist = sglEvnts$ages[2], method = "BC")
5658
5658
dev.off()
dist3 <- resInteractive(site = siteObj, window = 1200, dist = sglEvnts$ages[3], method = "BC")
7501
6881
dev.off()
dist4 <- resInteractive(site = siteObj, window = 1200, dist = sglEvnts$ages[4], method = "BC")
7964
7909
dev.off()
dist5 <- resInteractive(site = siteObj, window = 1200, dist = sglEvnts$ages[5], method = "BC")
11536
11536
dev.off()
# Combine them all into a complete table
hdg <- rbind(dist1, dist2, dist3, dist4, dist5)
write.csv(hdg, file = "hdg.csv")
# Do a quick cluster analysis on the pollen data to identify which pollen zones each of the disturbance events belong to
hdg <- doClust(site = siteObj, hdg )

# Make the plot
plot.hdg(hdg = hdg)

rm(list= ls())


