#-----------------------------------------------------------------------------#
# Mini meta on spread rates in exp-evo studies: analysis
# Date created: Sept 26 2018
# When using lnCVR codes, please cite: 
# Nakagawa et al. (2015) Methods in Ecology and Evolution, 6, 143-152.
#-----------------------------------------------------------------------------#


# Load packages and data:
#-----------------------------------------------------------------------------#
setwd("C:/Users/tm9/Dropbox/Manuscripts/EcoEvo range expansions/metaanalysis")

# Remove any old objects from the R environment
rm(list=ls())

# Load packages
library(metafor)
library(dplyr)

# Load functions
source("spread_functions.R")

# Load data
spdata <- read.csv("spread_data.csv", stringsAsFactors = FALSE, strip.white = T)

# Calculate speed and effect sizes:
#-----------------------------------------------------------------------------#

# Estimate speed and SD from distance, time, and SD
spdata$evo_speed      <- calcspeed(spdata$evo_mean_dist, spdata$gens)
spdata$evo_sd_speed   <- calcspeed(spdata$evo_sd_dist, spdata$gens)
spdata$noevo_speed    <- calcspeed(spdata$noevo_mean_dist, spdata$gens)
spdata$noevo_sd_speed <- calcspeed(spdata$noevo_sd_dist, spdata$gens)

# Calculate effect sizes: lnRR (mean contrast of speed between evolution vs. no evolution treatments) and its sampling variance
spdata$lnRR <- calclnRR(spdata$evo_speed, spdata$noevo_speed)
spdata$samplnRR <- calcsamplnRR(spdata$evo_sd_speed, spdata$noevo_sd_speed, spdata$evo_speed, spdata$noevo_speed, spdata$evo_N, spdata$noevo_N)

# Calculate effet sizes: lnCVR (variance contrast of speed between evolution vs. no evolution treatments) and its sampling variance; see Nakagawa et al. ref above 
spdata$lnCVR <- calclnCVR(spdata$evo_sd_speed, spdata$noevo_sd_speed, spdata$evo_speed, spdata$noevo_speed, spdata$evo_N, spdata$noevo_N)
spdata$samplnCVR <- calcsamplnCVR(spdata$evo_sd_speed, spdata$noevo_sd_speed, spdata$evo_speed, spdata$noevo_speed, spdata$evo_N, spdata$noevo_N)


# Model 1: all control-evolved data (but without protist core-edge)
#-----------------------------------------------------------------------------#

# Delete protist data
nonprot.data <- spdata[-c(7:8),] # remaining = 10 data points

# Effect size ID
nonprot.data$effectsize_ID <- as.factor(seq(1, length(nonprot.data$lnRR), 1))

# MLMA run for overall lnRR (mean) and lnCVR (variance)
mlma.lnRR <- rma.mv(yi = lnRR, V = samplnRR, data = nonprot.data, random = list(~1|ï..paper, ~1|effectsize_ID))
mlma.lnCVR <- rma.mv(yi = lnCVR, V = samplnCVR, data = nonprot.data, random = list(~1|ï..paper, ~1|effectsize_ID))

# Model summaries
summary(mlma.lnRR)
summary(mlma.lnCVR)

# Calculate Jensen's inequalities for each model to estimate exponential outcomes
Jen.mlma.lnRR <- ((sum(mlma.lnRR$sigma2))*0.5)
Jen.mlma.lnCVR <- ((sum(mlma.lnCVR$sigma2))*0.5)

# Estimate exponential outcomes for MLMA models
exp(mlma.lnRR$ci.lb+Jen.mlma.lnRR)
exp(mlma.lnRR$b+Jen.mlma.lnRR)
exp(mlma.lnRR$ci.ub+Jen.mlma.lnRR) 
exp(mlma.lnRR$se+Jen.mlma.lnRR) 

exp(mlma.lnCVR$ci.lb+Jen.mlma.lnCVR)
exp(mlma.lnCVR$b+Jen.mlma.lnCVR)
exp(mlma.lnCVR$ci.ub+Jen.mlma.lnCVR) 
exp(mlma.lnCVR$se+Jen.mlma.lnCVR) 

# save model outputs
saveRDS(mlma.lnRR, file = "./mod_outputs/mlma.lnRR.rds")
saveRDS(mlma.lnCVR, file = "./mod_outputs/mlma.lnCVR.rds")


# Model 2: subset to only homogeneous, favourable data
#-----------------------------------------------------------------------------#

# Delete protist data
cont.data <- nonprot.data[-c(2,4:6),] # remaining = 10 data points

# Effect size ID
cont.data$effectsize_ID <- as.factor(seq(1, length(cont.data$lnRR), 1))

# MLMA run for overall lnRR (mean) and lnCVR (variance)
mlma.lnRR.2 <- rma.mv(yi = lnRR, V = samplnRR, data = cont.data, random = list(~1|effectsize_ID))
mlma.lnCVR.2 <- rma.mv(yi = lnCVR, V = samplnCVR, data = cont.data, random = list(~1|effectsize_ID))

# Model summaries
summary(mlma.lnRR.2)
summary(mlma.lnCVR.2)

# Calculate Jensen's inequalities for each model to estimate exponential outcomes
Jen.mlma.lnRR.2 <- ((sum(mlma.lnRR.2$sigma2))*0.5)
Jen.mlma.lnCVR.2 <- ((sum(mlma.lnCVR.2$sigma2))*0.5)

# Estimate exponential outcomes for MLMA models
exp(mlma.lnRR.2$ci.lb+Jen.mlma.lnRR.2)
exp(mlma.lnRR.2$b+Jen.mlma.lnRR.2)
exp(mlma.lnRR.2$ci.ub+Jen.mlma.lnRR.2) 
exp(mlma.lnRR$se+Jen.mlma.lnRR) 

exp(mlma.lnCVR.2$ci.lb+Jen.mlma.lnCVR.2)
exp(mlma.lnCVR.2$b+Jen.mlma.lnCVR.2)
exp(mlma.lnCVR.2$ci.ub+Jen.mlma.lnCVR.2) 
exp(mlma.lnCVR.2$se+Jen.mlma.lnCVR.2) 

# save model outputs
saveRDS(mlma.lnRR.2, file = "./mod_outputs/mlma.lnRR.2.rds")
saveRDS(mlma.lnCVR.2, file = "./mod_outputs/mlma.lnCVR.2.rds")

# Model 3: model 1 but with landscape as covariate
#-----------------------------------------------------------------------------#

# MLMA run for overall lnRR (mean) and lnCVR (variance)
mlma.lnRR.3 <- rma.mv(yi = lnRR, mods = ~landscape, V = samplnRR, data = nonprot.data, random = list(~1|ï..paper, ~1|effectsize_ID))
mlma.lnCVR.3 <- rma.mv(yi = lnCVR, mods = ~landscape, V = samplnCVR, data = nonprot.data, random = list(~1|ï..paper, ~1|effectsize_ID))

# Model summaries
summary(mlma.lnRR.3)
summary(mlma.lnCVR.3)

# Calculate Jensen's inequalities for each model to estimate exponential outcomes
Jen.mlma.lnRR.3 <- ((sum(mlma.lnRR.3$sigma2))*0.5)
Jen.mlma.lnCVR.3 <- ((sum(mlma.lnCVR.3$sigma2))*0.5)

# Estimate exponential outcomes for MLMA models
exp(mlma.lnRR.3$ci.lb[1]+Jen.mlma.lnRR.3)
exp(mlma.lnRR.3$b[1]+Jen.mlma.lnRR.3)
exp(mlma.lnRR.3$ci.ub[1]+Jen.mlma.lnRR.3) 
exp(mlma.lnRR.3$se[1]+Jen.mlma.lnRR.3) 

exp(mlma.lnCVR.3$ci.lb[1]+Jen.mlma.lnCVR.3)
exp(mlma.lnCVR.3$b[1]+Jen.mlma.lnCVR.3)
exp(mlma.lnCVR.3$ci.ub[1]+Jen.mlma.lnCVR.3) 
exp(mlma.lnCVR.3$se[1]+Jen.mlma.lnCVR.3) 

# save model outputs
saveRDS(mlma.lnRR.3, file = "./mod_outputs/mlma.lnRR.3.rds")
saveRDS(mlma.lnCVR.3, file = "./mod_outputs/mlma.lnCVR.3.rds")


# Figure
#-----------------------------------------------------------------------------#

# arrange data by reverse alphabetical order of author
nonprot.data <- arrange(nonprot.data, desc(nonprot.data$ï..paper))

# set plot area
par(mfrow = c(1, 2), mar=c(5,3,5,3)) 

# (a) plot area for lnRR
plot(NULL, NULL, 
     xlim = c(-0.2, 1.4),
     ylim = c(0, 10),
     xaxt = "n",
     yaxt = "n",
     axes = F,
     xlab = "",
     ylab = "",
     yaxs="i")
axis(1, seq(-0.2, 1.4, 0.4), cex.axis = 1.5) # define x-axis
abline(v = 0) # zero y axis

# polygon
xvect <- c(mlma.lnRR$ci.lb, mlma.lnRR$ci.ub, mlma.lnRR$ci.ub, mlma.lnRR$ci.lb) # x vect coords for polygon
yvect <- c(0, 0, 10, 10) # y vect coords for polygon
polygon(xvect, yvect,
        col =  adjustcolor("lightcyan", alpha.f = 0.6), border = NA) # define polygon

# plot points lnRR
ypositions <- seq(2, 9, 1) # y position points for study estimates 
points(nonprot.data$lnRR, 
       ypositions,
       pch = 16, col = "dark gray", cex = 1.5) # plot mean study estimates
arrows(y0 = ypositions,
       x0 = nonprot.data$lnRR,
       x1 = nonprot.data$lnRR+nonprot.data$samplnRR,
       y1 = ypositions, 
       code = 0, col = "dark gray") # upper var
arrows(y0 = ypositions,
       x0 = nonprot.data$lnRR,
       x1 = nonprot.data$lnRR-nonprot.data$samplnRR,
       y1 = ypositions, 
       code = 0, col = "dark gray") # lower var
points(mlma.lnRR$b, 1,
       pch = 18, col = "dodgerblue2", cex = 2) # plot overall model estimate
arrows(y0 = 1,
       x0 = mlma.lnRR$b,
       x1 = mlma.lnRR$ci.lb,
       y1 = 1, 
       code = 0, col = "dodgerblue2") # LCI
arrows(y0 = 1,
       x0 = mlma.lnRR$b,
       x1 = mlma.lnRR$ci.ub,
       y1 = 1, 
       code = 0, col = "dodgerblue2") # UCI

# (b) plot area lnCVR
plot(NULL, NULL, 
     xlim = c(-1.0, 1.0),
     ylim = c(0, 10),
     xaxt = "n",
     yaxt = "n",
     axes = F,
     xlab = "",
     ylab = "",
     yaxs="i")
axis(1, seq(-1.0, 1.0, 0.4), cex.axis = 1.5) # define x-axis
abline(v = 0) # define y-axis

# polygon
xvect <- c(mlma.lnCVR$ci.lb, mlma.lnCVR$ci.ub, mlma.lnCVR$ci.ub, mlma.lnCVR$ci.lb) # x vect coords for polygon
yvect <- c(0, 0, 10, 10) # y vect coords for polygon
polygon(xvect, yvect,
        col =  adjustcolor("lightcyan", alpha.f = 0.6), border = NA) # define polygon

# plot points lnCVR
ypositions <- seq(2, 9, 1)
points(nonprot.data$lnCVR, 
       ypositions,
       pch = 16, col = "dark gray", cex = 1.5)
arrows(y0 = ypositions,
       x0 = nonprot.data$lnCVR,
       x1 = nonprot.data$lnCVR+nonprot.data$samplnCVR,
       y1 = ypositions, 
       code = 0, col = "dark gray")
arrows(y0 = ypositions,
       x0 = nonprot.data$lnCVR,
       x1 = nonprot.data$lnCVR-nonprot.data$samplnCVR,
       y1 = ypositions, 
       code = 0, col = "dark gray")
points(mlma.lnCVR$b, 1,
       pch = 18, col = "dodgerblue2", cex = 2) # plot overall model estimate
arrows(y0 = 1,
       x0 = mlma.lnCVR$b,
       x1 = mlma.lnCVR$ci.lb,
       y1 = 1, 
       code = 0, col = "dodgerblue2") # LCI
arrows(y0 = 1,
       x0 = mlma.lnCVR$b,
       x1 = mlma.lnCVR$ci.ub,
       y1 = 1, 
       code = 0, col = "dodgerblue2") # UCI

## Tom's figure idea
ben_cont <- subset(nonprot.data,landscape=="benign_continuous")
RRmodel <- mlma.lnRR.3
CVRmodel <- mlma.lnCVR.3

win.graph()
par(mfrow=c(2,1),mar=c(1,4,1,1))
plot(NULL, NULL, 
     ylim = c(-0.2, 1.4),
     xlim = c(0.5, 4.5),
     xaxt = "n",
     yaxt = "n",
     axes = F,
     ylab = "Effect on mean speed",cex.lab=1,
     xlab = "")
axis(2, seq(-0.2, 1.4, 0.2), cex.axis = 0.6) # define x-axis
abline(v=c(1.5,2.5,3.5,4.5),col="lightgrey")
abline(h = 0,lty=2) # zero y axis
yvect <- c(RRmodel$ci.lb[1], 
           RRmodel$ci.ub[1], 
           RRmodel$ci.ub[1], 
           RRmodel$ci.lb[1]) # x vect coords for polygon
xvect <- c(0, 0, 5, 5) # y vect coords for polygon
polygon(xvect, yvect,
        col =  adjustcolor("dodgerblue", alpha.f = 0.2), border = NA) # define polygon
abline(h = RRmodel$b["intrcpt",],col="dodgerblue",lwd=2)
arrows(x0 = c(0.8,2:4),
       y0 = ben_cont$lnRR-ben_cont$samplnRR,
       y1 = ben_cont$lnRR+ben_cont$samplnRR,
       x1 = c(0.8,2:4), 
       code = 0)
points(c(0.8,2:4),ben_cont$lnRR,pch=21,bg="dodgerblue")
## Arabidopsis landscape treatments
arrows(x0 = c(1,1.1,1.2),
       y0 = nonprot.data$lnRR[2:4]-nonprot.data$samplnRR[2:4],
       y1 = nonprot.data$lnRR[2:4]+nonprot.data$samplnRR[2:4],
       x1 = c(1,1.1,1.2), 
       code = 0)
points(c(1,1.1,1.2),nonprot.data$lnRR[2:4],pch=22:24,bg="gray")
## Tribolium harsh
arrows(x0 = 2.1,
       y0 = nonprot.data$lnRR[7]-nonprot.data$samplnRR[7],
       y1 = nonprot.data$lnRR[7]+nonprot.data$samplnRR[7],
       x1 = 2.1, 
       code = 0)
points(2.1,nonprot.data$lnRR[7],bg="gray",pch=21)

#####################################################
plot(NULL, NULL, 
     ylim = c(-1, 1),
     xlim = c(0.5, 4.5),
     xaxt = "n",
     yaxt = "n",
     axes = F,
     ylab = "Effect on speed CV",cex.lab=1,
     xlab = "")
axis(2, seq(-1, 1, 0.25), cex.axis = 0.6) # define x-axis
abline(v=c(1.5,2.5,3.5,4.5),col="lightgrey")
abline(h = 0,lty=2) # zero y axis
yvect <- c(CVRmodel$ci.lb[1], 
           CVRmodel$ci.ub[1], 
           CVRmodel$ci.ub[1], 
           CVRmodel$ci.lb[1]) # x vect coords for polygon
xvect <- c(0, 0, 5, 5) # y vect coords for polygon
polygon(xvect, yvect,
        col =  adjustcolor("dodgerblue", alpha.f = 0.2), border = NA) # define polygon
abline(h = CVRmodel$b["intrcpt",],col="dodgerblue",lwd=2)
arrows(x0 = c(0.8,2:4),
       y0 = ben_cont$lnCVR-ben_cont$samplnCVR,
       y1 = ben_cont$lnCVR+ben_cont$samplnCVR,
       x1 = c(0.8,2:4), 
       code = 0)
points(c(0.8,2:4),ben_cont$lnCVR,pch=21,bg="dodgerblue")
## Arabidopsis landscape treatments
arrows(x0 = c(1,1.1,1.2),
       y0 = nonprot.data$lnCVR[2:4]-nonprot.data$samplnCVR[2:4],
       y1 = nonprot.data$lnCVR[2:4]+nonprot.data$samplnCVR[2:4],
       x1 = c(1,1.1,1.2), 
       code = 0)
points(c(1,1.1,1.2),nonprot.data$lnCVR[2:4],pch=22:24,bg="gray")
## Tribolium harsh
arrows(x0 = 2.1,
       y0 = nonprot.data$lnCVR[7]-nonprot.data$samplnCVR[7],
       y1 = nonprot.data$lnCVR[7]+nonprot.data$samplnCVR[7],
       x1 = 2.1, 
       code = 0)
points(2.1,nonprot.data$lnCVR[7],bg="gray",pch=21)


############# flip axes
pointsize <-1.6
pdf(file="metaanalysis_fig",height = 4, width = 7)

win.graph()
par(mfrow=c(1,2),mar=c(5,1,2,1))
plot(NULL, NULL, 
     xlim = c(-0.2, 1.4),
     ylim = c(0.5, 4.5),
     yaxt = "n",
     xaxt = "n",
     axes = F,
     xlab = "Effect on mean speed",cex.lab=1.4,
     ylab = "")
axis(1, seq(-0.2, 1.4, 0.2), cex.axis = 1) # define x-axis
abline(h=c(1.5,2.5,3.5,4.5),col="lightgrey")
abline(v = 0,lty=2) # zero y axis
xvect <- c(RRmodel$ci.lb[1], 
           RRmodel$ci.ub[1], 
           RRmodel$ci.ub[1], 
           RRmodel$ci.lb[1]) # x vect coords for polygon
yvect <- c(0, 0, 5, 5) # y vect coords for polygon
polygon(xvect, yvect,
        col =  adjustcolor("dodgerblue", alpha.f = 0.2), border = NA) # define polygon
abline(v = RRmodel$b["intrcpt",],col="dodgerblue",lwd=2)
arrows(y0 = c(0.8,2:4),
       x0 = ben_cont$lnRR-ben_cont$samplnRR,
       x1 = ben_cont$lnRR+ben_cont$samplnRR,
       y1 = c(0.8,2:4), 
       code = 0)
points(ben_cont$lnRR,c(0.8,2:4),pch=21,bg="dodgerblue",cex=pointsize)
## Arabidopsis landscape treatments
arrows(y0 = c(1,1.1,1.2),
       x0 = nonprot.data$lnRR[2:4]-nonprot.data$samplnRR[2:4],
       x1 = nonprot.data$lnRR[2:4]+nonprot.data$samplnRR[2:4],
       y1 = c(1,1.1,1.2), 
       code = 0)
points(nonprot.data$lnRR[2:4],c(1,1.1,1.2),pch=22:24,bg="gray",cex=pointsize)
## Tribolium harsh
arrows(y0 = 2.1,
       x0 = nonprot.data$lnRR[7]-nonprot.data$samplnRR[7],
       x1 = nonprot.data$lnRR[7]+nonprot.data$samplnRR[7],
       y1 = 2.1, 
       code = 0)
points(nonprot.data$lnRR[7],2.1,bg="gray",pch=21,cex=pointsize)
title("A",font=3,adj=0)
legend("topright",title="Landscape",pt.cex=1.2,bg="white",cex=0.8,
       legend=c("Continuous, benign",
                "Continuous, harsh",
                "Fragmented (small gaps)",
                "Fragmented (medium gaps)",
                "Fragmented (large gaps)"),
       pch=c(21,21,22:24),
       pt.bg=c("dodgerblue","gray","gray","gray","gray"))

#####################################################
plot(NULL, NULL, 
     xlim = c(-1, 1),
     ylim = c(0.5, 4.5),
     xaxt = "n",
     yaxt = "n",
     axes = F,
     xlab = "Effect on speed CV",cex.lab=1.4,
     ylab = "")
axis(1, seq(-1, 1, 0.25), cex.axis = 1) # define x-axis
abline(h=c(1.5,2.5,3.5,4.5),col="lightgrey")
abline(v = 0,lty=2) # zero y axis
xvect <- c(CVRmodel$ci.lb[1], 
           CVRmodel$ci.ub[1], 
           CVRmodel$ci.ub[1], 
           CVRmodel$ci.lb[1]) # x vect coords for polygon
yvect <- c(0, 0, 5, 5) # y vect coords for polygon
polygon(xvect, yvect,
        col =  adjustcolor("dodgerblue", alpha.f = 0.2), border = NA) # define polygon
abline(v = CVRmodel$b["intrcpt",],col="dodgerblue",lwd=2)
arrows(y0 = c(0.8,2:4),
       x0 = ben_cont$lnCVR-ben_cont$samplnCVR,
       x1 = ben_cont$lnCVR+ben_cont$samplnCVR,
       y1 = c(0.8,2:4), 
       code = 0)
points(ben_cont$lnCVR,c(0.8,2:4),pch=21,bg="dodgerblue",cex=pointsize)
## Arabidopsis landscape treatments
arrows(y0 = c(1,1.1,1.2),
       x0 = nonprot.data$lnCVR[2:4]-nonprot.data$samplnCVR[2:4],
       x1 = nonprot.data$lnCVR[2:4]+nonprot.data$samplnCVR[2:4],
       y1 = c(1,1.1,1.2), 
       code = 0)
points(nonprot.data$lnCVR[2:4],c(1,1.1,1.2),pch=22:24,bg="gray",cex=pointsize)
## Tribolium harsh
arrows(y0 = 2.1,
       x0 = nonprot.data$lnCVR[7]-nonprot.data$samplnCVR[7],
       x1 = nonprot.data$lnCVR[7]+nonprot.data$samplnCVR[7],
       y1 = 2.1, 
       code = 0)
points(nonprot.data$lnCVR[7],2.1,bg="gray",pch=21,cex=pointsize)
title("B",font=3,adj=0)

dev.off()
