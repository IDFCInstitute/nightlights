# radiance-calibrated-ols

## India's Nightlights Analysis (DMSP-OLS)


DMSP Nighttime data obtained from: https://ngdc.noaa.gov/eog/dmsp/downloadV4composites.html (2013 data, F182013) These files contain annual copmosite data regarding Average Visible, Stable Lights, & Cloud Free Coverages. The products are 30 arc second grids, spanning -180 to 180 degrees longitude and -65 to 75 degrees latitude.

India Shapefile: http://www.gadm.org/country (.rds file type) Here, we use the admin0 layer that has the country's shapefile.This shapefile is read into the wgs84 coordinate reference system, the standard for GPS. The same is done for the raster data from the tiff file. By assigning the same spatial projection, we can work across shapefiles and raster files. India Population Data: http://www.census2011.co.in/states.php (used as a .csv file)

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
getwd()
setwd("C:/Users/admin1/Desktop/ols")

imagery = "[path]/ols"

##Obtain a list of TIF files, load in the first file in list. Stored as 'imagery' folder
#tifs = list.files(ols,pattern = "\\.tif")- didn't work
rast <- raster("F16_20100111-20110731_rad_v4.avg_vis.tif")
rast

wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
projection(rast) <- CRS(wgs84)
 
      
bound<-readRDS("IND_adm0.rds")

projection(bound) <- CRS(wgs84)
names(bound)

require(rgdal)
# Read SHAPEFILE.shp from the current working directory (".")
shape <- readOGR(dsn = "C:/Users/admin1/Desktop/imagery", layer = "IND_adm0")
projection(shape) <- CRS(wgs84)
s<-shapefile("IND_adm0")

pop<- read.csv("pca.csv")
pop$NAME <- as.character(pop$State) 

```
Mapping India's nighttime data by cropping the shapefile. We first set the graph layout with no margins (mai), with 1 row and 1 column (mfrow), with a navy blue background (bg). A data frame is created for the coordinates which will be later called from source (google) in the looped code. The map we create focuses on comparing patterns across cities using the same color coding. Thus, one interval is generated based on a random sample of pixels (generated via set.seed) from across the country. A k-means clustering algorithm is used to find natural intervals within the radiance distribution. For each cluster of 10 pixels, we extract the maximum radiance. The country is then mapped with a navy to yellow color palette with intervals from the k-means clustering. To map the country, the extent specifies the spatial bounding box (the frame around the country) as +/-18 degree longitude and +/-18 degree latitude, from which all parts other than India are masked.


```{r}
country<-c("India,In")

par(mai=c(0,0,0,0),mfrow = c(1,1),bg='#001a4d', bty='n')
coords <- data.frame()
#Run clustering
set.seed(123) #set seed for reproducibility
sampled <- sample(rast, 20000) #sample 20,000 pixels
clusters <- 15 ##10 clusters
clust <- kmeans(sampled,clusters)$cluster
combined <- as.data.frame(cbind(sampled,clust))
brk <- sort(aggregate(combined[,1], list(combined[,2]), max)[,2])

##Loop through the country
for(i in 1:length(country)){
  
  temp_coord <- geocode(country[i], source = "google")
coords <- rbind(coords,temp_coord)
   
  e <- extent(temp_coord$lon - 13, temp_coord$lon + 18,
              temp_coord$lat - 17, temp_coord$lat + 17)
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
![india rad calib](https://user-images.githubusercontent.com/31407895/31481289-3ef97a24-af40-11e7-974e-2aa4d7dfcc59.png)
To analyse the distribution of these radiances as seen on the graph, we now create a function that extracts the radiance values (rast) based on the shape of a geographic area from shapefile (s), doing so for 1 to the i-th shape. The data frame hence created, contains Geoid/Locational name ( as Geoids for India are not available), the longitude, latitude and radiance values.The data is the processed first for India (country), and turned into a series of histograms using a combination of ggplot and plotly. This is done for the raw radiance values as well as the log values of these radiances. As already noted by other researchers (e.g. (Henderson, et al., 2012)), the amount of dim lit pixels is unreasonably low probably due to the empirical operation of noise removal. (See also  Pestalozzi, 2012.)
[[http://www.worldatnight.ethz.ch/content/doc/Nicola_Pestalozzi_Master_Thesis.pdf](url)]

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
  
  summary(rastdata)
  
  #package data in neat dataframe
  rastdata <- cbind(as.character(s@data$NAME_ENGLI[i]),coords, rastdata) 
  colnames(rastdata)<-c("GEOID","lon","lat","stable_lights") #note that 
  rastdata <- rastdata[!is.na(rastdata$stable_lights),] #keep non-NA values only
  
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
      rad <- masq(shp_temp,rast,1)$stable_lights 
      temp <- data.frame(loc = as.character(paste(loc,"(TNL = ",round(sum(rad),0),")",sep="")), stable_lights = rad) 
      radiances <- rbind.data.frame(radiances,temp)
}
attach(radiances)  
  
  #stderr <- function(stable_lights) sqrt(var(stable_lights,na.rm=TRUE)/length(na.omit(stable_lights)))

#stderr(stable_lights)


#Use ggplot to create histogram.
  ggplot(radiances, aes(x=stable_lights)) +
    geom_histogram(position="identity",stat="bin", binwidth = 1, na.rm=FALSE, alpha=0.6) +
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

![india histo rad calib](https://user-images.githubusercontent.com/31407895/31481365-803cec78-af40-11e7-9ee2-c0ff361925db.png)

|Average Radiance Summary Statistics|Radiance Value|
|-----------------------------------|--------------|
|Minimum  |	0.000|
|Mean    	| 4.740|
|Median   |	3.865|
|Maximum	 |581.221|
|Standard Error| 0.00499166|

![boxplot india](https://user-images.githubusercontent.com/31407895/31481819-ff49325e-af42-11e7-8ccb-0b380b859a7a.png)
The k-density plots are also seen for the distribution of randiances.
```{r}

d <- density(stable_lights, na.rm=T)
plot(d, main="Kernel Density of Radiance")
polygon(d, col="red",border="blue")
```

![k density calib](https://user-images.githubusercontent.com/31407895/31481847-1da8c23c-af43-11e7-9e37-d3d07d0e42b5.png)

