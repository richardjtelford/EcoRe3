install.packages("rioja")
install.packages("vegan")
install.packages("plyr")
install.packages("analogue")

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


### In some records we have indepedent disturbance data, which makes it easier to look at resistance in more detail

dataFile <- "TAS1110_CHAR_NEW.csv"
dataDir <- "data"
fileToLoad <- file.path(mainDir, dataDir, dataFile)
char <- read.csv2(fileToLoad, header = TRUE)
head(char)
colnames(char) <- c("age", "macro", "ageOs", "micro", "X")
char<- char[!is.na(char$age), ]
char <- list(ages = char$age, macro = char$macro)

dataFile <- "TAS1110_POLLEN_NEW.csv"
dataDir <- "data"
fileToLoad <- file.path(mainDir, dataDir, dataFile)
poll <- read.csv2(fileToLoad, header = TRUE)
head(poll)
ages <- poll$AGE
poll <- poll[,-1]
siteObj <- list(core = poll, ages = ages, dataset = "TAS1110")
rm(poll, ages)
dev.off()
# Step 1: Identify fire events against background variability (NB, realise there are other methods available here, e.g. Higuera's CHAR analysis, but using a simple percentiile approach for now for simplicity)

distChar <- charEvnts(char)$dist
makeCharPlot(char = char, dist.df = distChar)

# Then we can use the same methods as previously to investigate the transitions
d1 <- resInteractive(site, dist = distChar$start[1], window= 500, method = "BC", interactive = TRUE, bottom = NA, end = NA)
4334
4264
d2 <- resInteractive(site, dist = distChar$start[2], window= 500, method = "BC", interactive = TRUE, bottom = NA, end = NA)
4970
4615
d3 <- resInteractive(site, dist = distChar$start[3], window= 500, method = "BC", interactive = TRUE, bottom = NA, end = NA)
5181
5110
d4<- resInteractive(site, dist = distChar$start[4], window= 500, method = "BC", interactive = TRUE, bottom = NA, end = NA)
5388
5319
d5 <- resInteractive(site, dist = distChar$start[5], window= 500, method = "BC", interactive = TRUE, bottom = NA, end = NA)
5658
5658
d6 <- resInteractive(site, dist = distChar$start[6], window= 500, method = "BC", interactive = TRUE, bottom = NA, end = NA)
5658
5658
# Unsure what to do with this one
d7 <- resInteractive(site, dist = distChar$start[7], window= 500, method = "BC", interactive = TRUE, bottom = NA, end = NA)
6492
6492
d8 <- resInteractive(site, dist = distChar$start[8], window= 500, method = "BC", interactive = TRUE, bottom = NA, end = NA)
6789
6693
# Becomes interesting when you scrutinise this event
d9 <- resInteractive(site, dist = distChar$start[9], window= 500, method = "BC", interactive = TRUE, bottom = NA, end = NA)
7964
7909

# Need to fix this error
d10 <- resInteractive(site, dist = distChar$start[10], window= 500, method = "BC", interactive = TRUE, bottom = NA, end = NA)


### Then combine them all into one dataframe
hdg <- rbind(d1, d2, d3, d4, d5, d6, d7, d8, d9)
hdg <- doClust(site = site, hdg )

distChar <- distChar[-10,] # because I haven't fixed this error
hdgCh <- cbind(hdg, distChar)


par(mfrow = c(3,2))
plot(hdgCh$maxMag, hdgCh$resistance, ylab = "Change in state", xlab = "Charcoal.Z", bg = hdgCh$clust, pch = 21, type = "b")
plot(hdgCh$resistance, hdgCh$return, xlab = "Change in state", ylab = "Return Time", bg = hdgCh$clust, pch = 21, cex = sqrt(hdgCh$maxMag), type = "b")
plot(hdgCh$resistance, hdgCh$recovery, xlab = "Change in state", ylab = "Recovery Time", bg = hdgCh$clust, pch = 21, type = "b", cex = sqrt(hdgCh$maxMag))
plot(hdgCh$resistance, hdgCh$r.resilience, xlab = "Change in state", ylab = "r.Resilience", bg = hdgCh$clust, pch = 21, type = "b", cex = sqrt(hdgCh$maxMag))
	plot(hdgCh$recovery, hdgCh$resilience, xlab = "Recovery", ylab = "Resilience", bg = hdgCh$clust, pch = 21, type = "b", cex = sqrt(hdgCh$maxMag))
	plot(hdgCh$recovery, hdgCh$r.resilience, xlab = "Recovery", ylab = "r.Resilience", bg = hdgCh$clust, pch = 21, type = "b", cex = sqrt(hdgCh$maxMag))

