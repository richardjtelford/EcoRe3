mainDir <- "~/EcoRe3"
setwd(mainDir)

################################################################
# Taking it forward to improve the analysis?
################################################################

# Also started to write some breakpoint models that can do this. Load
source("ecore3_functions.R") # and then you try to apply this on some data and fit the model using segmented.

par(mfrow = c(2,2))
# There are number of different shapes a disturbance model can take
# For some simulated data, you can model a classic disturbance model
simDat <- makeDummyTS(x = 1:30, brk = data.frame(kD = 5, k1 = 10, k2 = 20), coef = data.frame(d0 = 0, d1 = -2, d2 = 3, d3 = -1))
# plotSimDat(simDat)
fit <- fit.bkPt(x = simDat$x, y = simDat$z, bkpt = list(x = simDat$brk))
plot.bkPt(fit)

# Or try a ramped response
simDat <- makeDummyTS(x = 1:30, brk = data.frame(kD = 10, k1 = 20), coef = data.frame(d0 = 0, d1 = -2, d2 = 3))
# plotSimDat(simDat)
fit <- fit.bkPt(x = simDat$x, y = simDat$z, bkpt = list(x = simDat$brk))
plot.bkPt(fit)

# Or a shallower response
simDat <- makeDummyTS(x = 1:30, brk = data.frame(kD = 5, k1 = 25), coef = data.frame(d0 = 0, d1 = -2, d2 = 3))
# plotSimDat(simDat)
fit <- fit.bkPt(x = simDat$x, y = simDat$z, bkpt = list(x = simDat$brk))
plot.bkPt(fit)

# Or maybe there is just a trend following the disturbance
simDat <- makeDummyTS(x = 1:30, brk = data.frame(kD = 20), coef = data.frame(d0 = 0, d1 = -0.5))
# plotSimDat(simDat)
fit <- fit.bkPt(x = simDat$x, y = simDat$z, bkpt = list(x = 10))
plot.bkPt(fit)


######################################################################
### The problem is that these models are difficult to fit on palaeodata. 
######################################################################

###############################################
# The code from earlier to get the data ready
###############################################

dataFile <- "TAS1110_POLLEN_NEW.csv"
poll <- read.csv2(dataFile, header = TRUE)
head(poll)
ages <- poll$AGE
poll <- poll[,-1]
siteObj <- list(core = poll, ages = ages, dataset = "TAS1110")

# Calculate the BC distances
polBC <- calcBC(site= siteObj )

# Make the BC distance plot
plotBC(polBC)

# Choose how many perturbations you want and find out where they are in the record
distEvents <- extractDist.perc(calcBC.obj = polBC, perc = 0.95, BCselect = 2)

# find out which events are in consecutive samples and only select non-consecutive ones
dEvents <- c(diff(distEvents$sample), 1000) # the 1000 is just an arbitrary number added to make the calculation work so the last sample can be included
sglEvnts <- distEvents[dEvents > 2,]

# Plot the full PCA and the full dataset to check the disturbance points make sense
plotCore(site = siteObj, evntObj = sglEvnts)

###############################################
# Now use the segmented package instead
###############################################

dist1 <- localOrd(site = siteObj, dist = sglEvnts$ages[1])
localOrd.plot(dist1)
fitDist1 <- fitDist(	x = dist1$PC$ages,
						y = dist1$PC$PC1,
						bk = c(dist1$dist, 1597, 1544))
						
# But it (sometimes) works better here, for example
dist3 <- localOrd(site = siteObj, dist = sglEvnts$ages[3], window= 1200)
localOrd.plot(dist3)
fitDist3 <- fitDist(	x = dist3$PC$ages,
						y = dist3$PC$PC1,
						bk = c(dist3$dist, 7501, 6881))

dist2 <- localOrd(site = siteObj, dist = sglEvnts$ages[2], window= 500)
localOrd.plot(dist2)
fitDist2 <- fitDist(	x = dist2$PC$ages,
						y = dist2$PC$PC1,
						bk = c(dist2$dist, 5658))

######## So with high resolution data, it may be possible to use breakpoint models to automaticaly detect and define disturbances. What other data types (palaeo/ non-palaeo) could we apply these methods to? This a potential discussion point.

# Development of Bayesian versions of this I think would be preferential (e.g. set the priors on the breakpoint nodes, might allow for a bit more rigidity in the models that you would like to fit)





