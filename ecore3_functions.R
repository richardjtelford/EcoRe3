

###############################################
## Calculate Bray Curtis Distance
###############################################

# test edit

calcBC <- function(site) {
	
	core <- site$core
	ages <- site$ages
	dataset <- site$dataset
	
	# Do the B-C metric
	turn <- as.matrix(vegdist(core, method="bray"))
	BC1 <- diag(turn[-1, -ncol(turn)]) # compared to next sample
	BC2 <- diag(turn[-(1:2), -ncol(turn)]) # compared to second to next sample

	return(list(ages= ages, BC1 = BC1, BC2= BC2, dataset = dataset))
}

###############################################
## Make BC inspection plot
###############################################

plotBC<- function(BCobject, print.pdf= FALSE)	{
	
	# Make the decision plot
	with(BCobject, {
		if(print.pdf == TRUE) pdf(paste("decisionPlot_", dataset, ".pdf", sep ="" ))
		par(mfrow = c(3,2), mar = c(3,3,1,1), mgp = c(1.5, .5, 0))
		plot(ages[-1], BC1, type = "h")
		plot(ages[-c(1,2)], BC2, type = "h")
		hist(BC1)
		hist(BC2)
		qqnorm(BC1)
		qqline(BC1)
		qqnorm(BC2)
		qqline(BC2)
		par(mfrow = c(1,1))
		if(print.pdf == TRUE) dev.off()
		
	})
}

###############################################
## Extract the high BC events
###############################################

extractDist <- function(calcBC.obj, nDist, BCselect){
	
if(BCselect == 1) {
	want <- sort(calcBC.obj$BC1, decreasing = TRUE)[1:nDist]
	wantTime <- match(want, calcBC.obj$BC1)
	} else {
		want <- sort(calcBC.obj$BC2, decreasing = TRUE)[1:nDist]
	wantTime <- match(want, calcBC.obj$BC2)
	}
	agesCrop <- calcBC.obj$ages[-BCselect]
	distEvents <- 	agesCrop[wantTime]
	return(distEvents)
}

###############################################
## Extract the high BC events using the nth percentiles
###############################################

extractDist.perc <- function(calcBC.obj, perc = 0.9, BCselect= 2){
	
	if(BCselect == 1) BC.hist <- calcBC.obj$BC1 else BC.hist <- calcBC.obj$BC2
	
	names(BC.hist) <- calcBC.obj$ages[- c(1: BCselect)]
	thresh <- quantile(BC.hist, probs = perc)
	turnover <- BC.hist[BC.hist > thresh]
	distEvents <- data.frame(ages = as.numeric(names(turnover)),
								BC = turnover,
								sample = match(as.numeric(names(turnover)), calcBC.obj$ages ))
								
								
	return(distEvents)
}


###############################################
# Do a PCA of the whole dataset and plot it
###############################################

plotCore <- function(site, evntObj, strat= TRUE){
	
	require("analogue")
	require("vegan")
	
	core <- site$core
	ages <- site$ages
	pca <- rda(sqrt(core))
	PC <- data.frame(scores(pca)$sites )
	if(strat== TRUE){
		pollDat <- data.frame(age = ages, chooseTaxa(core, n.occ = 5, max.abun = 5))
		Stratiplot(age ~ ., data = pollDat, 	type = c("h", "l", "g"), sort = "wa")
	} else {
		m <- rbind(c(1,1,3), c(2,2,3))
		layout(m)
		plot(ages, PC$PC1, pch= 21, bg = "darkgrey", type = "b")
		points(c(min(ages), evntObj$ages, max(ages)), rep(max(PC$PC1), nrow(evntObj)+2), pch = 17, col = c("black", rep("red", nrow(evntObj)+2), "black"), type = "b")
		plot(ages, PC$PC2, pch= 21, bg = "darkgrey", type = "b")
		points(c(min(ages), evntObj$ages, max(ages)), rep(max(PC$PC2), nrow(evntObj)+2), pch = 17, col = c("black", rep("red", nrow(evntObj)+2), "black"), type = "b")
		plot(pca)
	}
}


###############################################
# Estimate turnover using PCA within a given window
###############################################

localOrd <- function( site, dist, window = 500){
	require("analogue")
	core <- site$core
	ages <- site$ages
	want <- ages > dist-window  & ages < dist + window
	coreWant <- sqrt(core[want,]) # Note the sqrt transformation
	ageWant <- ages[want]
	
	pollDat <- data.frame(ageWant = ageWant, chooseTaxa(coreWant^2, n.occ = 5, max.abun = 5))
	Stratiplot(ageWant ~ ., data = pollDat,
							type = c("h", "l", "g"), sort = "wa")

	#quartz() # NB won't work on windows, use window()/windows()?
	pca <- rda(coreWant)
	PC <- data.frame(scores(pca)$sites)
	bstick <- screeplot(pca, bstick = TRUE)
	
	resultPC <- data.frame(ages = ageWant, PC)
	listPC <- list(PC = resultPC, dist = dist, window = window)
	return(listPC)
}



###############################################
# Plot local Ordination
###############################################

localOrd.plot <- function(localOrd.obj){
	
	PC <- localOrd.obj$PC
	dist <- localOrd.obj$dist
	
	# par(mfrow = c(2, 1))
	plot(PC$ages, PC$PC1, pch= 21, bg = "grey", col = "black", type= "b")
	abline(v = dist, col = "red", lty = 2)
	# plot(PC$ages, PC$PC2, pch= 21, bg = "grey", col = "black")
	# abline(v = dist, col = "red", lty = 2)
	par(mfrow = c(1, 1))	
}


###############################################
# Calculate resistance, recovery, resilience with manual input of timings
###############################################
	
re3 <- function(localOrd.obj, bottom, end_, site, method = c("BC", "PCA")){
	
	PC <- localOrd.obj$PC
	PC1 <- PC$PC1
	#start <- localOrd.obj$dist
	# This bit needed for matching up pollen and charcoal data
	age <- site$ages
	diffAge <- abs(age- localOrd.obj$dist)
	start <- 	age[diffAge == min(diffAge)]	
	
	returnTime <- start - end_
	recoveryTime <- bottom - end_

	if(method == "BC"){
		core <- site$core
		age <- site$ages
		ageResist <- match(c(start, bottom), age)
		ageResil <- match(c(start, end_), age)
		ageRecov <- match(c(bottom, end_), age)
		
		resistance <- vegdist(sqrt(core[ageResist,]), method="bray")[1]
		recovery <- vegdist(sqrt(core[ageRecov,]), method="bray")[1]
		resilience <-  vegdist(sqrt(core[ageResil,]), method="bray")[1]
		r.resilience <-  resilience / resistance
		 
	} else if(method == "PCA") {
		
	  PreD <- abs((PC1[PC$ages == start] - PC1[PC$ages == bottom]))
	  D <- 0
	  PostD <-  abs(PC1[PC$ages == end] - PC1[PC$ages == bottom])
	  resistance <-  D / PreD
		recovery <- PostD / D
		resilience <- PostD / PreD
	  r.resilience <- ((PostD-D)/(PreD- D)) * (1 -(D/PreD))
		}
		
	re3.obj <- data.frame(start = start, bottom = bottom, end = end, 
	                      resistance = resistance, resilience = resilience, r.resilience = r.resilience, recovery = recovery, 
	                      return = returnTime, recoveryTime = recoveryTime)
	return(re3.obj)
}


###############################################
# Function to interactively calculate resistance and recovery axes and enter values manually
###############################################


resInteractive <- function(site, dist, window, method = BC, interactive = TRUE, bottom = NA, end = NA){
	
	distResult <- localOrd(site= site, dist = dist, window)
	localOrd.plot(distResult)
	if(interactive == TRUE){
	print(distResult)
	bottom <- readline("What is the age of the disturbance minimum?   ")
	end <- readline("What is the age of the recovery date?     ")
	bottom <- as.numeric(unlist(strsplit(bottom, ",")))
	end <- as.numeric(unlist(strsplit(end, ",")))
	}
	re3Result <- re3(distResult, bottom = bottom, end = end, site = site, method = method)
	return(re3Result)

} 




###############################################
# Cluster analysis to identify the vegetation zones and append to hdg table
###############################################

doClust <- function (site, hdg){
	diss <- vegdist(sqrt(site$core/100), method = "bray")
	clust <- chclust(diss)
	sigZone <- bstick(clust,10)
	sigZone$diff <- sigZone$dispersion - sigZone$bstick
	nG <- 1
	for(i in 1:length(sigZone$diff)) ifelse(sigZone$diff[i] > 0 , nG <- sigZone$nGroups[i], break)
	
	site$clustZone <- cutree(clust, k = nG)
	want <- match(hdg$start, site$ages)
	hdg$clust <- site$clustZone[want]
	return(hdg)
}


###############################################
# Hodgson Re3 plots
###############################################

plot.hdg <- function(hdg){
	par(mfrow = c(3,2))
	plot(hdg$resistance, hdg$recoveryTime, xlab = "Change in state", ylab = "Recovery Time", bg = hdg$clust, pch = 21, type = "b", cex = 3)
	plot(hdg$resistance, hdg$return, xlab = "Change in state", ylab = "Return Time", bg = hdg$clust, pch = 21, type = "b", cex = 3)
	plot(hdg$resistance, hdg$r.resilience, xlab = "Change in state", ylab = "r.Resilience", bg = hdg$clust, pch = 21, type = "b", cex = 3)
	plot(hdg$recoveryTime, hdg$resilience, xlab = "Recovery", ylab = "Resilience", bg = hdg$clust, pch = 21, type = "b", cex = 3)
	plot(hdg$recoveryTime, hdg$r.resilience, xlab = "Recovery", ylab = "r.Resilience", bg = hdg$clust, pch = 21, type = "b", cex = 3)
}


###############################################
# Try fitting a disturbance model (uses the breakpoints close to it)
###############################################

fitDist <- function(x, y, bk){
	
	xN <- (x - max(x))*-1
	bkN <- (bk - max(x))*-1
	fit <- fit.bkPt(x = xN, y = PC1, bkpt = list(x = bkN))
	plot.bkPt(fit)
	return(fit)
}

################################################
# Fit a breakpoint model with given number of breakpoints
################################################

fit.bkPt <- function(x, y, bkpt){
	require(segmented)
	fit.lm <- lm(y ~ x )
	fit.seg <- segmented(fit.lm, seg.Z = ~ x, psi = bkpt)
	fit.lm <- update(fit.lm, .~.-x)
	fit.seg1 <- update(fit.seg)
	pred <- predict(fit.seg1)
	BIC.fit <- BIC(fit.seg1)
	results <- list(x = x, y = y, pred.y = pred, BIC = BIC.fit, model = fit.seg1)
	return(results)
	}



###############################################
# Plot the breakpoint model
###############################################

plot.bkPt <- function(fit.bkPt.obj){
	x <- fit.bkPt.obj$x
	y <- fit.bkPt.obj$y
	predX <- data.frame(x = seq(min(x), max(x), 0.1))
		
	plot(x, y, pch = 20, col = "darkgrey")	
	points(x, y)
	lines(predX$x, predict(fit.bkPt.obj$model, newdata = predX), col = "blue")
	points.segmented(fit.bkPt.obj$model, col = "red", pch = 13)
	
	#	abline(v = c(simDat$brk), col = "red", lty = 2)	
	
}


################################################
# Simulate some dummy data with breakpoints
################################################

makeDummyTS <- function(x, brk, coef){
	
	with(cbind(brk, coef), {
		nbrk <- length(brk)
		if(nbrk == 3) {		
		 	y <- ifelse(x-kD<= 0, d0,
				ifelse(x-kD > 0 & x-k1 <= 0, d0 + d1* x,
					ifelse(x-k1 > 0 & x-k2 <= 0, d0 + d1* x + d2*(x-k1),
						d0 + d1* x + d2*(x-k1) + d3*(x-k2) ))) + rnorm(length(x), 0, 3)
		} else if(nbrk == 2){
				y <- ifelse(x-kD<= 0, d0,
				ifelse(x-kD > 0 & x-k1 <= 0, d0 + d1* x,
					 d0 + d1* x + d2*(x-k1))) + rnorm(length(x), 0, 3)
		} else		
				y <- ifelse(x-kD<= 0, d0, d0 + d1* x)  + rnorm(length(x), 0, 3)	
	results <- list(x = x, y = y, brk = brk, coef = coef, z = scale(y))
	return(results)
	})
		
} 

################################################
# Plot the dummy data with breakpoints
################################################

plotSimDat <- function(simDat){
	x <- simDat$x
	y <- simDat$z
	plot(x, y, pch = 20, col = "darkgrey")	
	points(x, y)
	abline(v = c(simDat$brk), col = "red", lty = 2)
	
}


##############################################
# Char events
##############################################

charEvnts <- function(char, perc = 0.9){
  require("plyr")
  ages <- char$ages
  macro <- char$macro
  
  
  # find the 90th percentile for the charcoal data
  hist(macro, breaks = length(macro)/10)
  cutOff <- quantile(macro, c(perc))
  events <- which(macro > cutOff)
  
  # Plot them
  plot(ages, macro, type = "h")
  points(ages[events], macro[events], pch = "x", col= "red", cex = 0.75)
  
  # find out which events are in consecutive samples
  dEvents <- c(diff(events), 1000)
  sglEvnts <- events[dEvents > 3]
  points(ages[sglEvnts], macro[sglEvnts], pch = 20, col= "blue")
  points(ages[sglEvnts], macro[sglEvnts], pch = 1)
  
  #sglEvents gives you the start time the of the events
  
  #Get the duration of the event	
  nEvents <- rep(1, length(events))
  N <- 1
  for(i in 1:(length(dEvents) - 1)){
    if(dEvents[i] < 3) {
      nEvents[i+1] <- N
    } else{
      N <- N + 1
      nEvents[i+1] <- N
    }
  }
  points(ages[events], macro[events], pch = 20, col= (nEvents+1), cex = 0.75)
  
  # Then make an event table
  evnt.df <- data.frame(	sample = events,
                         age = ages[events],
                         peak = macro[events],
                         peakZ = scale(macro)[events],
                         nEvnt = as.factor(nEvents)
  )
  
  #get the maximum value and other values for an event. Creates your "disturbance table"
  
  dist.df <- ddply(evnt.df, .(nEvnt), summarize, maxMag = max(peakZ), meanMag = mean(peakZ), start = max(age), end = min(age), duration = max(age)- min(age), startSmp = max(sample))
  
  charEvntObj <- list(dist = dist.df, evnt = evnt.df)
  return(charEvntObj)
}

################################################
# function to make the standard charcoal plot
################################################

makeCharPlot <- function(char, dist.df, xLim= NULL){
  
  plot(char$age, char$macro, type = "h", xlim = xLim)
  for(i in 1:nrow(dist.df)){	
    st <- dist.df$start[i]
    end <- dist.df$end[i] 
    polygon(x = c(st, end, end, st), y = c(0,0, 5,5), density = NA, col = rgb(0,0,255, max = 255, alpha = 50))
    abline(v = dist.df$start[i], lty = 3)	
  }
}







# # 
# fit.bkPt <- function(x, y, bkpt){
	# nBkpt <- length(bkpt$x)
		# if(nBkpt == 1) { res <- fit.bkPt1(x, y, bkpt)
	# } 	else if(nBkpt == 2) {
			# res <- fit.bkPt2(x, y, bkpt)
	# } else  res<- fit.bkPt3(x, y, bkpt)
	# return(res)
# }


# fit.bkPt1 <- function(x, y, bkpt){
	# require(segmented)
	# fit.lm <- lm(y ~ x )
	# fit.seg <- segmented(fit.lm, seg.Z = ~ x, psi = bkpt)
	# fit.lm <- update(fit.lm, .~.-x)
	# fit.seg1 <- update(fit.seg)
	# pred <- predict(fit.seg1)
	# BIC.fit <- BIC(fit.seg1)
	# results <- list(x = x, y = y, pred.y = pred, BIC = BIC.fit, model = fit.seg1)
	# return(results)
	# }



# fit.bkPt2 <- function(x, y, bkpt){
	# require(segmented)
	# fit.lm <- lm(y ~ x )
	# fit.seg <- segmented(fit.lm, seg.Z = ~ x, psi = bkpt)
	# fit.lm <- update(fit.lm, .~.-x)
	# fit.seg1 <- update(fit.seg)
	# pred <- predict(fit.seg1)
	# BIC.fit <- BIC(fit.seg1)
	# results <- list(x = x, y = y, pred.y = pred, BIC = BIC.fit, model = fit.seg1)
	# return(results)
	# }



# fit.bkPt3 <- function(x, y, bkpt){
	# require(segmented)
	# fit.lm <- lm(y ~ x )
	# fit.seg <- segmented(fit.lm, seg.Z = ~ x, psi = bkpt)
	# fit.lm <- update(fit.lm, .~.-x)
	# fit.seg1 <- update(fit.seg)
	# pred <- predict(fit.seg1)
	# BIC.fit <- BIC(fit.seg1)
	# results <- list(x = x, y = y, pred.y = pred, BIC = BIC.fit, model = fit.seg1)
	# return(results)
	# }


	# plot(ageWant, DC$DCA1, pch= 21, bg = "grey", col = "black")
	# abline(v = dist, col = "red", lty = 2)
	# plot(ageWant, DC$DCA2, pch= 21, bg = "grey", col = "black")
	# abline(v = dist, col = "red", lty = 2)
	# plot(ageWant, DC$DCA3, pch= 21, bg = "grey", col = "black")
	# abline(v = dist, col = "red", lty = 2)
	# plot(ageWant, DC$DCA4, pch= 21, bg = "grey", col = "black")
	# abline(v = dist, col = "red", lty = 2)




