rm(list=ls(all=TRUE))

######################################################
# Load packages
######################################################

library(rgdal)
library(ggplot2)
library(ggmap)

######################################################
# Load shapefile and convert into dataframe
######################################################

# load shapefile
frax <- readOGR("~/Dropbox/chalarose/ashchal/Volume17_wgs84.shp")
fraxR <- frax
# Convert Volume shapefile into dataframe
frax <- as.data.frame(frax@data)
colnames(frax)[colnames(frax)=="xl"] <- "X"
colnames(frax)[colnames(frax)=="yl"] <- "Y"

######################################################
# Map settings
######################################################

# Location
lat <- c(42, 51.5)
lon <- c(-5, 8.5)

# Download background
map <- get_map(location = c(lon = mean(lon), lat = mean(lat)), zoom = 5,
               maptype = "toner-background", source = "google")

# set background
bckgrd <- ggmap(map)+
  scale_x_continuous(limits = lon, expand = c(0, 0)) +
  scale_y_continuous(limits = lat, expand = c(0, 0))

# theme settings
theme=theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.background = element_rect(fill = 'white'),
            legend.position = c(0,0),
            legend.justification=c(-0.05,-0.05),
            text = element_text(size=12),
            axis.text.x = element_text(size=10),
            legend.key = element_blank())

####################################################
# Plot Map
####################################################

bckgrd+
  theme_bw()+
  theme+
  xlab("longitude") + ylab("latitude")+
  geom_point(data=frax[frax$V2 != 0, ], aes(X, Y, col=Volume), shape = 15, size = 3)+
  # scale_colour_gradient(limits = c(0,600), colours=c("lightgreen","green4","black"))+
  scale_colour_gradient2(low = "lightgreen", mid = "green4", high = "black", midpoint = 200)+
  # scale_color_manual(values=c("palegreen1", "seagreen1", "seagreen3", "seagreen4"))+
  geom_polygon(data = fortify(fraxR), aes(long, lat, group = group), color = "grey", fill = "grey", alpha = 0.1)

####################################################
# Save Map and volume data
####################################################

ggsave("~/Dropbox/chalarose/ashchal/volume.pdf")
vol <- frax
save(vol, file = '~/Dropbox/chalarose/ashchal/volume.rdata')


############################################################################################################
############################################################################################################
# Load temperautre txt
######################################################

temp <- read.table("~/Dropbox/chalarose/ashchal/MMTP_2007-2017_coord.txt", header = TRUE, sep = '\t')

# coordinates names
colnames(temp)[colnames(temp) == "x_wgs84"] <- "X"
colnames(temp)[colnames(temp) == "y_wgs84"] <- "Y"

bckgrd+
  theme_bw()+
  theme+
  xlab("longitude") + ylab("latitude")+
  geom_point(data=temp[temp$MMTP_2007 >0,], aes(X, Y, col=MMTP_2007), shape = 15, size = 3, alpha = 0.8)+
  scale_colour_gradient2(low = "blue", mid = "green", high = "red", midpoint = 20)+
  geom_polygon(data = fortify(fraxR), aes(long, lat, group = group), color = "grey", fill = "grey", alpha = 0.1)

####################################################
# Save Map and volume data
####################################################

ggsave("~/Dropbox/chalarose/ashchal/temp.pdf")
save(temp, file = '~/Dropbox/chalarose/ashchal/temp.rdata')

############################################################################################################
############################################################################################################
# Load chalara csv
######################################################


##### !!!!!!!!!!! coordonnées en WGS84 !!!!!!!!!!!!!!!!

chal <- read.csv("~/Dropbox/chalarose/ashchal/chalThierry.csv", header = TRUE, sep = ';')

plot(chal$xl2e, chal$yl2e)
for (i in 2008:2018){
  lines(chal[chal$annee == i, 'xl2e'], chal[chal$annee == i, 'yl2e'], col = i-2007, pch = 16, type = 'p')
  Sys.sleep(1)
}
