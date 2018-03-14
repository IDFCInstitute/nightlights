# DMSP OLS Night-time State-wise Analysis

## India's State-wise Nightlights Analysis (DMSP-OLS)


DMSP Nighttime data obtained from: https://ngdc.noaa.gov/eog/dmsp/downloadV4composites.html (2013 data, F182013) These files contain annual copmosite data regarding Average Visible, Stable Lights, & Cloud Free Coverages. The products are 30 arc second grids, spanning -180 to 180 degrees longitude and -65 to 75 degrees latitude.


India Shapefile: http://www.gadm.org/country (.rds file type)

India Population Data: http://www.census2011.co.in/states.php (used as a .csv file)
```{r}
library(doParallel)
library(foreach)
library(reshape2)
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
rast <- raster("F182013.v4c_web.stable_lights.avg_vis.tif")
rast

wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
projection(rast) <- CRS(wgs84)
 
      
bound<-readRDS("IND_adm1.rds")

projection(bound) <- CRS(wgs84)
names(bound)

require(rgdal)
# Read SHAPEFILE.shp from the current working directory (".")
shape <- readOGR(dsn = "C:/Users/admin1/Desktop/ols", layer = "IND_adm1")
projection(shape) <- CRS(wgs84)
s<-shapefile("IND_adm1")

pop<- read.csv("pca.csv")
pop$NAME <- as.character(pop$State) 


```
Now, after selecting the states to be compared, we geocode the states and create their spatial bounding box. This is then extended via longitude and latitide to include a greater statespatial box and cropped from the raster tile. This cropped raster formes the sample size within which we first identify clusters of 15 and create a separate data frame which we can then plot.

```{r}
states<- c("Maharashtra,India","Bihar,India","Gujarat,India","West Bengal,India","Kerala,India","Madhya Pradesh,India")
par(mai=c(0,0,0,0),mfrow = c(3,3),bg='#001a4d', bty='n')

coords <- data.frame() ##place holder

for(i in 1:length(states)){
  ##Coords
  temp_coord <- geocode(states[i], source = "google")
  coords <- rbind(coords,temp_coord)
  
  e <- extent(temp_coord$lon - 4, temp_coord$lon + 4,
              temp_coord$lat - 2, temp_coord$lat + 2)
 rc<- crop(rast,e)   
  
  ##Rescale brackets
  sampled <- as.vector(rc)
  clusters <- 15
  clust <- kmeans(sampled,clusters)$cluster
  combined <- as.data.frame(cbind(sampled,clust))
  brk <- sort(aggregate(combined[,1], list(combined[,2]), max)[,2])
  
  #Plots
  plot(rc, breaks=brk, col=colorRampPalette(c("#001a4d","#0066FF", "yellow"))(clusters), 
       legend=F,yaxt='n',xaxt='n',frame = T, asp=1.5); plot(bound, add=T, border="white")
  text(temp_coord$lon ,temp_coord$lat + 0.15,
       substr(states[i],1,regexpr(",",states[i])-1), 
       col="white", cex=1.25)
  
  rm(combined)
}
```
![ols states1](https://user-images.githubusercontent.com/31407895/31212233-e397e108-a9bc-11e7-845f-268a8d7f0f91.PNG)

For inter-state comparisons, now the entire raster is the sample size from which clusters are identified. Now, these clusters are sorted by radiance. Post this, individual spatial bound boxes for each state are made and plotted.
```{r}
par(mai=c(0,0,0,0),mfrow = c(3,3),bg='#001a4d', bty='n')

#Run clustering
set.seed(53542) #set seed for reproducibility
sampled <- sample(rast, 20000) #sample 20,000 pixels
clusters <- 15 ##15 clusters
clust <- kmeans(sampled,clusters)$cluster
combined1 <- as.data.frame(cbind(sampled,clust))
brk <- sort(aggregate(combined1[,1], list(combined1[,2]), max)[,2])

##Loop through each city
for(i in 1:length(states)){
  
  temp_coord <- coords[i,] ##re-use the coordinates 
  e <- extent(temp_coord$lon - 4, temp_coord$lon + 4,
              temp_coord$lat - 2, temp_coord$lat + 2)
  rc <- crop(rast, e)    
  
  #Plots
  plot(rc, breaks=brk, col=colorRampPalette(c("#001a4d","#0066FF", "yellow"))(clusters), 
       legend=F,yaxt='n',xaxt='n',frame = F, asp=1.5)
  plot(bound, add=T, border="white")
  text(temp_coord$lon ,temp_coord$lat + 0.15,
       substr(states[i],1,regexpr(",",states[i])-1), 
       col="white", cex=1.25)
  
  rm(combined1)
}
```
![ols states2](https://user-images.githubusercontent.com/31407895/31212245-fd8cb3fe-a9bc-11e7-83c5-7bc746fc0707.PNG)

However, since the data is low resolution, significant differences in inter and intra state nighttime visibility cannot be made out.

To create a histogram for the 6 state comparison:
Average radiances are extracted from the raster and shape files for the states and plotted. However, theses plots do not take into account the area of each state and analysis based on the histograms should be made keeping this fact in mind. Here, unlike the VIIRS dataset, we use the radiance values as their range is not large to require logarithmic transformation. As already noted by other researchers (e.g. (Henderson, et al., 2012)), the amount of dim lit pixels is unreasonably low probably due to the empirical operation of noise removal. (See also  Pestalozzi, 2012.)
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
  
  #package data in neat dataframe
  rastdata <- cbind(as.character(s@data$NAME_1[i]),coords, rastdata) 
  colnames(rastdata)<-c("GEOID","lon","lat","stable_lights") #note that 
  rastdata <- rastdata[!is.na(rastdata$stable_lights),] #keep non-NA values only
  
  return(rastdata)
}
 
skt<-c(1:36)

  radiances <- data.frame() 
 
for(i in skt){
  
    print(skt)
    
    #Extract i polygon
      shp_temp <- shape[shape@data$ID_1==i,]
    
   #if(regexpr(" ",as.character(shp_temp@data$NAME_1)[1])[1]==-1){
        loc = as.character(shp_temp@data$NAME_1)[1]
      #} else{
        #loc = (as.character(shp_temp@data$NAME_1)[0])
      #}
    
    #Extract the radiances, append to radiances placeholder
      rad <- masq(shp_temp,rast,1)$stable_lights 
      temp <- data.frame(loc = as.character(paste(loc,"(TNL = ",round(sum(rad),0),")",sep="")), stable_lights = rad) 
      radiances <- rbind(radiances,temp)
      #print(radiances)
}

#Use ggplot to create histograms by States.
  p<-ggplot(radiances, aes(x=stable_lights)) +
    geom_histogram(position="identity", alpha=0.6)
    
    p + facet_wrap( ~ loc, scales="free_y")

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
|Minimum  |	0.000|
|Mean    	| 0.000|
|Median   |	5.067|
|Maximum	 | 7.000|
|Standard Error| 0.004247306|

![all states dmsp ols hist](https://user-images.githubusercontent.com/31407895/31531688-2ae9784c-b006-11e7-9835-2206b54a3d3e.png)

We now have a scatter plot of Total Nighttime Light and Population to see the underlying correlation between the two. For this, extracted GEOIDs and average radiance data frame extracted is merged with the population data frame and then plotted. Based on the plot, there exists a positive correlation between the two parameters. The Total Nighttime Light is the sum of all the average radiance values of each pixel of the state. The procedure is then used to plot Total Nighttime Light and Population Density and State GDP, respectively. For State GDP, 2013 figures are taken at Constant Prices(2004-05) from MOSPI.
  
```{r}

library(sp)
  registerDoParallel(cores=2)
    extract <- foreach(i=1:nrow(s@data),.packages='raster', .combine=rbind) %dopar% {
        rastdata <- masq(shape,rast,i)
       data.frame(GEOID = rastdata$GEOID[1],sum = sum(rastdata$stable_lights))
    }
   #extract$GEOID <- as.numeric(as.character(extract$GEOID))
    
   
   
  # D <- as.numeric(as.character(extract$GEOID))
    
  ##Join in data
    
  joined<-merge(extract, pop[,c("Density","NAME","Population")],by.x="GEOID",by.y="NAME")
   
    colnames(joined) <- c("State","TNL","Density","Population")
   joined<- read.csv("joined.csv")
   joined1<-na.omit(joined)
    attach(joined1)
    aTNL<-log(joined1$TNL)
    bPop<- log(joined1$Density)
    cState<-joined1$State
    m1<- lm(GDP2013~TNL, data=joined1)
    summary(m1)
    
    m2<-lm(bPop~aTNL, na.rm=T)
    summary(m2)
    
    
    joinedfd<-cbind(joined1,aTNL,bPop,cState)
    attach(joinedfd)
     t <- list(
  family = "sans serif",
  size = 9,
  color = toRGB("grey50"))
        
     
      p<-plot_ly(joinedfd, 
            x = aTNL,
            y = bPop, 
            text = ~cState,
             mode = "markers+text",    marker = list(size = 3),
            color = TNL,colors="PuOr") %>%
      add_markers()  %>%
        add_lines(y = fitted(loess(m2)),
            line = list(color = '#07A4B5'),
           # text="GDP=(4.672e+04)+(2.316e-01)*TNL",
            name = "Loess Smoother", showlegend = T) %>%
  add_text(textfont = t, textposition = "bottom right") %>%
                                                                                   #layout(text="GDP=(4.672e+04)+(2.316e-01)*TNL, adjusted R^s= 0.5478, p-value=3.264e-05")%>%
            layout(title="Total Nighttime Light vs. Population Density ",xaxis = list(title = 'log(TNL)'),
         yaxis = list(title = 'log(Population Density)'), showlegend = F)
    
p

```
![tnl gdp ols log](https://user-images.githubusercontent.com/31407895/31212998-92e5dcb0-a9c1-11e7-8ddb-1e632417f1f7.png)
![tnl pop ols log](https://user-images.githubusercontent.com/31407895/31213000-92f7aada-a9c1-11e7-874f-bec38e909fac.png)
![tnl density ols log](https://user-images.githubusercontent.com/31407895/31212999-92ef9872-a9c1-11e7-82d1-7f375dc49c1c.png)


|Dependant Variable|Intercept|Independant Variable(TNL)|
|------------------|---------|------------------|
|GDP 2013|7.811e+03| 2.632e-01 \*\*\*\|
|Population|1.193e+07| 4.669e+01\*\*\*|
Population Density|1.266e+03\*| -5.747e-04|

|Dependant Variable|Intercept|Independant Variable(log(TNL))|
|------------------|---------|------------------|
|log(GDP 2013)|  1.72772\* | 0.77411***|
|log(Population)|  5.6398**| 0.8842***|
|log(Population Density)| 4.92551 \*|   0.06886|

Signif. codes:  0 ‘\*\*\*’ 0.001 ‘\*\*’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



