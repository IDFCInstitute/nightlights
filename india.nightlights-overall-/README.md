# VIIRS Night-time India Analysis

---

## India Analysis
VIIRS Nighttime data obtained from: https://ngdc.noaa.gov/eog/viirs/download_dnb_composites.html (2017 monthly data, June month, Tile3_75N060E)
This tile is the one  containing data that excludes any stray light.The file ending in avg_rade9 contains average radiances; this has been used here.

India Shapefile: http://www.gadm.org/country (.rds file type)
Here, we use the admin0 layer that has the country's shapefile.This shapefile is read into the wgs84 coordinate reference system, the standard for GPS. The same is done for the raster data from the tiff file.
By assigning the same spatial projection, we can work across shapefiles and raster files.
India Population Data: http://www.census2011.co.in/states.php (used as a .csv file)


```{r}

library(doParallel)
library(foreach)
library(raster)
library(sp)
library(rgdal)
library(ggmap)
library(plotly)
library(ggplot2)
library(devtools)
library(spatstat)
library(maptools)
getwd()
setwd("C:/Users/admin1/Desktop/imagery")

imagery = "[path]/imagery"

##Obtain a list of TIF files, load in the first file in list. Stored as 'imagery' folder
tifs = list.files(imagery,pattern = "\\.tif")
rast <- raster("image.tif")

wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
projection(rast) <- CRS(wgs84)
 
      
bound<-readRDS("IND_adm0.rds")

projection(bound) <- CRS(wgs84)


require(rgdal)
# Read SHAPEFILE.shp from the current working directory (".")
shape <- readOGR(dsn = "C:/Users/admin1/Desktop/imagery", layer = "IND_adm0")
projection(shape) <- CRS(wgs84)
s<-shapefile("IND_adm0")



pop<- read.csv("pca.csv")
pop$NAME <- as.character(pop$State) 

```
Mapping India's nighttime data by cropping the shapefile.
We first set the graph layout with no margins (mai), with 1 row and 1 column (mfrow), with a navy blue background (bg). A data frame is created for the coordinates which will be later called from source (google) in the looped code. The map we create focuses on comparing patterns across cities using the same color coding. Thus, one interval is generated based on a random sample of pixels (generated via set.seed) from across the country. A k-means clustering algorithm is used to find natural intervals within the radiance distribution. For each cluster of 10 pixels, we extract the maximum radiance. The country is then mapped with a navy to yellow color palette with intervals from the k-means clustering. To map the country, the extent specifies the spatial bounding box (the frame around the country) as +/-18 degree longitude and +/-18 degree latitude, from which all parts other than India are masked.

```{r}
country<- c("India,In")

par(mai=c(0,0,0,0),mfrow = c(1,1),bg='#001a4d', bty='n')
coords <- data.frame()
#Run clustering
set.seed(123) #set seed for reproducibility
sampled <- sample(rast, 20000) #sample 20,000 pixels
clusters <- 10 ##10 clusters
clust <- kmeans(sampled,clusters)$cluster
combined <- as.data.frame(cbind(sampled,clust))
brk <- sort(aggregate(combined[,1], list(combined[,2]), max)[,2])

##Loop through the country
for(i in 1:length(country)){
  
  temp_coord <- geocode(country[i], source = "google")
coords <- rbind(coords,temp_coord)
   
  e <- extent(temp_coord$lon - 18, temp_coord$lon + 18,
              temp_coord$lat - 18, temp_coord$lat + 18)
  rc <- crop(rast, e)    
  a<-mask(rc,s)
  
  #Plots
  plot(a, breaks=brk, col=colorRampPalette(c("#001a4d","#0066FF", "yellow"))(clusters), 
       legend=F,yaxt='n',xaxt='n',frame = T, asp=1.5)
  plot(bound, add=T, border="white")
  text(temp_coord$lon ,temp_coord$lat + 0.15,
       substr(country[i],1,regexpr(",",country[i])-1), 
       col="white", cex=1.25)
  
  rm(combined)
}
```
![indiaplot](https://user-images.githubusercontent.com/31407895/30101311-719c3140-9309-11e7-9853-cc27fad0b895.png)

To analyse the distribution of these radiances as seen on the graph, we now create a function that extracts the radiance values (rast) based on the shape of a geographic area from  shapefile (s), doing so for 1 to the i-th shape. The data frame hence created, contains Geoid/Locational name ( as Geoids for India are not available), the longitude, latitude and radiance values.The data is the processed first for India (country), and turned into a series of histograms using a combination of ggplot and plotly.

```{r}


masq <- function(s,rast,i){
  #Extract one polygon based on index value i
  polygon <- s[i,] #extract one polygon
  extent <- extent(polygon) #extract the polygon extent 
  
  #Raster extract
  outer <- crop(rast, extent) #extract raster by polygon extent
  inner <- mask(outer,polygon) #keeps values from raster extract that are within polygon
  
  #Convert cropped raster into a vector
  #Specify coordinates
  coords <- expand.grid(seq(extent@xmin,extent@xmax,(extent@xmax-extent@xmin)/(ncol(inner)-1)),
                        seq(extent@ymin,extent@ymax,(extent@ymax-extent@ymin)/(nrow(inner)-1)))
  #Convert raster into vector
  rastdata <- as.vector(inner)
  
  #package data in neat dataframe
  rastdata <- cbind(as.character(s@data$NAME_ENGLI[i]),coords, rastdata) 
  colnames(rastdata)<-c("GEOID","lon","lat","avg_rad") #note that 
  rastdata <- rastdata[!is.na(rastdata$avg_rad),] #keep non-NA values only
  
  return(rastdata)
}


skt<-c("India")

  radiances <- data.frame() 
 
for(i in length(skt)){
  
    print(skt)
    
    #Extract i polygon
      shp_temp <- shape[shape@data$NAME_ENGLI==skt,]
    
        loc = as.character(shp_temp@data$NAME_ISO)[1]
    
    #Extract the radiances, append to radiances placeholder
      rad <- masq(shp_temp,rast,1)$avg_rad 
      temp <- data.frame(loc = as.character(paste(loc,"(TNL = ",round(sum(rad),0),")",sep="")), avg_rad = rad) 
      radiances <- rbind(radiances,temp)
}
attach(radiances)
  #Use ggplot to create histogram.
  ggplot(radiances, aes(x=log(avg_rad))) +
    geom_histogram(position="identity",stat="bin", binwidth = 0.9, na.rm=FALSE, alpha=0.6) +
    facet_grid(. ~ loc)

#Remove all axes labels for style
    x <- list(
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
    )
    y <- list(
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
    ) 
    
#Initiate a plotly graph without axes
  ggplotly()  %>% layout(xaxis=x, yaxis=y)

```
|Average Radiance Summary Statistics|Radiance Value|
|-----------------------------------|--------------|
|Minimum  |	-0.562|
|Mean    	| 0.328|
|Median   |	0.45|
|Maximum	 | 19843.096|
|Standard Error| 0.001079989|

The summary statistics clearly indicate the variance of the data collected. The VIIRS data is inclusive and sensitive to the Day/Night Band (DNB) which is panchromatic, sensitive to visible and near-infrared (NIR) wavelengths ranging from daylight down to the low levels of radiation observed at night.  When the DNB radiances are calculated during the production of the SDR files, the clear-sky offsets are subtracted, which may result in negative radiance values being reported. These values are kept in the SDR files (i.e. not replaced with SOUB fill) even though negative radiances are physically impossible because some physical features visible in DNB imagery would be lost if the negative values were excluded by the user. (See NOAA Techincal Report, VIIRS Imagery EDR User's Guide, 2013). [[http://rammb.cira.colostate.edu/projects/npp/VIIRS_Imagery_EDR_Users_Guide.pdf](url)]

On further exploring the radiance dataset for India, we look at radiance distribution within different ranges through boxplots.

![boxplots](https://user-images.githubusercontent.com/31407895/31211452-c61ce066-a9b6-11e7-9b2e-c6ab62752c64.PNG)

The nature of such dispersed data is why we then plot the distribution of logarithimic values of radiances.

![india histo](https://user-images.githubusercontent.com/31407895/30099611-db695e82-9303-11e7-93c8-776280785b00.png)

We also see the radiance distribution as k-density distribution curve. The histogram and the curve both follow the same shape. VIIRS data is already an average monthly composite term. Thus, the k-denisty plots this reported average radiance of each pixel.

```{r}
#To get the k-density plot to see the distribution curve.
lograd<-log(avg_rad)
d <- density(lograd, na.rm=T)
plot(d, main="Kernel Density of Average Radiance")
polygon(d, col="red",border="blue")

```
![india kdensity](https://user-images.githubusercontent.com/31407895/30099626-ebe2f818-9303-11e7-87fc-a6345cd8fe88.png)

Extra:
When the extremely high radiance values from the dataset are mapped, we see that these occur due to technical issues ( in a pattern on the map) or on geographic places that could have high reflectance (such as waterbodies).
![image](https://user-images.githubusercontent.com/31407895/31211755-493acea2-a9b9-11e7-8673-a9bb9ef7378a.png)
![image](https://user-images.githubusercontent.com/31407895/31211765-5b91051c-a9b9-11e7-9ab4-5f83820f73b0.png)
![image](https://user-images.githubusercontent.com/31407895/31211773-750abd8a-a9b9-11e7-9076-63581b020a47.png)




 



