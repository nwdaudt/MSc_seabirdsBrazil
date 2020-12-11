###################################################################################
## Seabird assemblage at the mouth of Amazon River and its relationship with 
## environmental characteristics
##
## Daudt et al. (2019) Journal of Sea Research, 155, 101826
## 
## Supplementary material (S1) – R code
##
## Code by Nicholas W Daudt
## R version 3.4.2
####################################################################################

## Load required packages ####

# install.packages("devtools")
# devtools::install_github("gavinsimpson/ggvegan")
# install.packages("countreg", repos="http://R-Forge.R-project.org")

library(tidyverse)
library(matrixStats)
library(vegan)
library(ggplot2)
library(ggvegan)
library(countreg)

## Open data ####

allData <- read.csv("./allData.csv", sep = ";", dec = ",", header = TRUE)
allData$Date <- as.Date.character(allData$Date)

## Subset of the respective types of seabird-counts (in this case, 'ccTotal' 
## refers to 'continuous-count' (cc))
ccTotal <- subset(x = allData, subset = CensusDescription %in% c("Contínuo"))
        # 'Contínuo' means 'continuous-count'.

## In 'ccTotal' subset, seabird families were grouped to subsequent analysis
ccTotal$abundTot <- rowSums(ccTotal[,c(18:42)]) #abundTot = total abundance

ccTotal$Procellariidae <- 
  ccTotal$Calonectris.borealis + ccTotal$Puffinus.puffinus + ccTotal$Puffinus.sp.

ccTotal$Hydrobatidae <- 
  ccTotal$Oceanites.oceanicus + ccTotal$Oceanodroma.leucorhoa + ccTotal$Hydrobatidae.sp

ccTotal$Stercorariidae <- 
  ccTotal$Stercorarius.parasiticus + ccTotal$Stercorarius.sp.

ccTotal$Laridae <- ccTotal$Leucophaeus.atricilla

ccTotal$Sternidae <- ccTotal$Sterna.hirundo + ccTotal$Sterna.sp.

## For multivariate CCA, delete the counts that sum 'zero'
ccTotalMultivar <- ccTotal[-ccTotal$abundTot !=0,]

## Subset just birds
ccBirdsGroups <- subset(x = ccTotalMultivar, select = Procellariidae:Sternidae)
ccBirdsGroups$abundTotGr <- rowSums(ccBirdsGroups[, c(1:5)])
ccBirdsGroups <- ccBirdsGroups[-ccBirdsGroups$abundTotGr !=0,]

ccBirdsGroups <- subset(x = ccBirdsGroups, select = Procellariidae:Sternidae)

## Subset environmental variables
ccEnvMultivar <- subset(x = ccTotalMultivar, select = SSS:BAT)

## Summary ####

## Density from 'continuous-counts' (cc). 
## For 'snapshot-counts' (ci) the same procedure was made.
ccAreaTotal <- matrixStats::sum2(ccTotal$Area)
matrixStats::colSums2(as.matrix(ccBirds))/ccAreaTotal

## Frequency of occurrence (%FO) from 'continuous-counts' (cc). 
## For 'ship-followers', 'snapshot-counts', and 'point-counts' the same procedure.
zero.cc <- as.vector(print(matrixStats::colCounts(x = as.matrix(ccBirds), 
                                                  value = 0)))
cc.fo <- as.vector(100 - ((zero.cc/118) * 100))
cc.fo

## Exploratory analysis ####
## Need {lattice} package
## Code based on:
## Alain F. Zuur, et al. (2010) 
## A protocol for data exploration to avoid common statistical problems.
## Methods Ecol. Evo. doi: 10.1111/j.2041-210X.2009.00001.x

## CCA ####

ccaGroups <- vegan::cca(ccBirdsGroups ~ SST + DIST + CHL + BAT, 
                        data = ccEnvMultivar[-c(1),])

anova(ccaGroups, by = "axis", permutations = 99)
anova(ccaGroups, by = "term", permutations = 99)

ccaPlot <- ggvegan::autoplot(ccaGroups) + 
  ggplot2::theme_bw() + ggplot2::expand_limits(x = 1.5)

ccaPlot

## Multiple-regression linear models (GLM) ####

## Zero-part: Binomial GLM
Procella.bin <- glm(factor(ccTotal$Procellariidae>0) ~ SST + CHL + ((BAT)*-1), 
                    data = ccTotal, family = binomial)
summary(Procella.bin)
# and so on, for each group...

## Checking zero-part model fitting
plot(Procella.bin)
# and so on, for each group...

## Count-part: Zero-truncated GLM
Procella.zt.p <- 
  countreg::zerotrunc(ccTotal$Procellariidae ~ SST + CHL + ((BAT)*-1), 
                           data = ccTotal, subset = ccTotal$Procellariidae>0)
summary(Procella.zt.p)
# and so on, for each group...

## Checking count-part model fitting
countreg::qqrplot(Procella.zt.p)
countreg::rootogram(Procella.zt.p)
# and so on, for each group...