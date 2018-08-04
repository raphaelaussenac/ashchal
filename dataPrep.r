# packages
library(rgdal)
library(ggplot2)
library(ggmap)

# Open shapefiles
frax <- readOGR("~/Téléchargements/Mailles_2012/Volume17_wgs84.shp")
regions <- readOGR("~/Téléchargements/regions/regions-20180101.shp")

# Convert Volume shapefile into dataframe
frax <- as.data.frame(frax@data)
colnames(frax)[colnames(frax)=="xl"] <- "X"
colnames(frax)[colnames(frax)=="yl"] <- "Y"

# plot
ggplot(data=frax[frax$V2 != 0, ])+
  geom_point(aes(X, Y, col=V2))+
  coord_fixed()


library(ggmap)

### centre de la carte
lat <- c(42, 51.5)
lon <- c(-5, 10)
### télécharger fond de carte chez google
map <- get_map(location = c(lon = mean(lon), lat = mean(lat)), zoom = 5,
               maptype = "toner-background", source = "google")
### integrer fond de carte google dans ggplot
bckgrd <- ggmap(map)+
  scale_x_continuous(limits = lon, expand = c(0, 0)) +
  scale_y_continuous(limits = lat, expand = c(0, 0))

### parametres carte
theme=theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.background = element_rect(fill = 'white'),
            legend.position = c(1,0),
            legend.justification=c(1,0),
            text = element_text(size=12),
            axis.text.x = element_text(size=10),
            legend.key = element_blank())

####################################################
##                        Map                     ##
####################################################
# all sites
bckgrd+
  theme_bw()+
  theme+
  xlab("longitude") + ylab("latitude")+
  ggtitle("All sites with at least one core of one of the species")+
  geom_point(data=frax[frax$V2 != 0, ], aes(X, Y, col=V2))



