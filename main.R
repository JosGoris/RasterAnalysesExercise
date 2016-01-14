#Hakken en Zagen
#13 January 2016
#Mark ten Vregelaar, Jos Goris

#Start with empty environment
rm(list=ls())

#Import libraries
library(raster)
library(rasterVis)

ifolder <- './data/'
ofolder <- './output/'
dir.create(ifolder, showWarnings = FALSE)
dir.create(ofolder, showWarnings = FALSE)

# Download the data
dataURL <- "https://github.com/GeoScripting-WUR/AdvancedRasterAnalysis/archive/gh-pages.zip"
inputZip <- list.files(path=ifolder, pattern= '^.*\\.zip$')

if (length(inputZip) == 0){ ##only download when not alrady downloaded
  download.file(url = dataURL, destfile = 'data/landsatData.zip', method = 'wget')
} 

unzip('data/landsatData.zip', exdir=ifolder)

# load data
load("data/AdvancedRasterAnalysis-gh-pages/data/GewataB2.rda")
load("data/AdvancedRasterAnalysis-gh-pages/data/GewataB3.rda")
load("data/AdvancedRasterAnalysis-gh-pages/data/GewataB4.rda")
load("data/AdvancedRasterAnalysis-gh-pages/data/GewataB1.rda")
load("data/AdvancedRasterAnalysis-gh-pages/data/GewataB5.rda")
load("data/AdvancedRasterAnalysis-gh-pages/data/GewataB7.rda")
load("data/AdvancedRasterAnalysis-gh-pages/data/vcfGewata.rda")
load("data/AdvancedRasterAnalysis-gh-pages/data/trainingPoly.rda")


#Brick the layers
gewata <- brick(GewataB1, GewataB2, GewataB3, GewataB4, GewataB5, GewataB7)

#Remove all values above 100 for 
vcfGewata[vcfGewata > 100] <- NA

#Add VCF to brick
covs <- addLayer(gewata, vcfGewata)

#Pair the bands to the vcf
pairs(covs)

#Rename the columns
names(covs) <- c('band1', 'band2', 'band3', 'band4', 'band5', 'band7', 'vcf')

#Make covs a data frame and omit the NA's
valuetable <- getValues(covs) 
valuetable <- na.omit(valuetable)
valuetable <- as.data.frame(valuetable)

#Make the linear model and show the summary
regression <- lm(vcf ~ band1 + band2 + band3 + band4 + band5 + band7, data = valuetable)
summary(regression)

#Predict the VCF
predVCF <- predict(covs, model=regression, na.rm=TRUE)

#Remove the negative values
predVCF[predVCF<0] <- NA

#Plot the results
mycolorkey <- list(at=seq(0, 100, 5), labels=list(labels=seq(0, 100, 10), at=seq(0, 100, 10)), space = 'bottom')
levelplot(stack(vcfGewata, predVCF), col.regions = colorRampPalette(c('brown','yellow','darkgreen'))(255), colorkey = mycolorkey,
          names.attr = c('Original VCF','Predicted VCF'), main = list('Tree cover in Gewata',cex=1.8,vjust=-1),scales=list(draw=FALSE),xlab=list('Tree cover (%)',vjust=7))

#Calculate the difference between the tree cover and the predicted tree cover
difference <- (predVCF - vcfGewata)

#Calculate the root mean square error
RMSE <- sqrt((cellStats((difference)^2, stat = mean, na.rm = TRUE)))
names(RMSE) <- 'RMSE Total Raster'

#Add column 'Code' to the trainingPoly 
trainingPoly@data$Code <- as.numeric(trainingPoly@data$Class)

#Rasterize the trainingPoly
zones <- rasterize(trainingPoly, difference, field='Code')

#Mask the classes
differenceMasked <- mask(difference, zones)

#Compute the mean for the difference^2 for the zones
zones <- zonal((differenceMasked)^2, zones)

#Calculate the RSME for the zones
RMSEzones <- sqrt(zones)

#plot the result, remove the column zone, name the rows
rownames(RMSEzones) <- c('Croplands','Forest','Wetlands')
RMSEzones <- RMSEzones[,-1]
barplot(c(RMSEzones, RMSE), main = 'RMSE per zone', ylim = c(0, 12), col = c('lightgreen','darkgreen','lightblue','pink') )

