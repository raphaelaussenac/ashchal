rm(list=ls(all=TRUE))

######################################################
# Load packages
######################################################

library(rgdal)
library(ggplot2)
library(ggmap)

######################################################
# Load data
######################################################

# load Volume shapefile
frax <- readOGR("~/Dropbox/chalarose/ashchal/Volume17_wgs84.shp")
fraxR <- frax
# Convert Volume shapefile into dataframe
vol <- as.data.frame(frax@data)
colnames(vol)[colnames(vol)=="xl"] <- "X"
colnames(vol)[colnames(vol)=="yl"] <- "Y"
vol <- vol[, c("X", "Y", "Volume")]

# load temperature data
temp <- read.table("~/Dropbox/chalarose/ashchal/MMTP_2007-2017_coord.txt", header = TRUE, sep = '\t')
# coordinates names
colnames(temp)[colnames(temp) == "x_wgs84"] <- "X"
colnames(temp)[colnames(temp) == "y_wgs84"] <- "Y"
temp$x_l2e <- NULL
temp$y_l2e <- NULL

# load chalara data
chal <- read.csv("~/Dropbox/chalarose/ashchal/chalThierry.csv", header = TRUE, sep = ';')
# coordinates names
colnames(chal)[colnames(chal) == "xwgs84"] <- "X"
colnames(chal)[colnames(chal) == "ywgs84"] <- "Y"
chal <-  chal[, c("X", "Y", "annee")]
chal$annee <- as.numeric(chal$annee)

######################################################
# merge Temperature and Volume data
######################################################

# assign a unique number to each cell in the volume dataframe
vol$ID <- paste(vol$X, vol$Y, sep = "")

# assign a unique number to each cell in the temperature dataframe
temp$ID <- paste(temp$X, temp$Y, sep = "")

# merge dataframes
df <- merge(vol[, c("ID", "Volume")], temp, by = "ID")
df$annee <- 999

######################################################
# assign DSF data to cells in the "Temp-Volume" grid
######################################################

# assign a unique ID to each cell in the "DSF" grid
chal$ID <- paste(chal$X, chal$Y, sep = "")

# bind the 2 grids for distance calculation
combinedCoord <- rbind(df[, c("ID", "X", "Y")], chal[, c("ID", "X", "Y")])
row.names(combinedCoord) <- combinedCoord$ID

# calculate the distances among points (cells)
distMat <- as.data.frame(as.matrix(dist(x = combinedCoord[,c("X", "Y")], method = "euclidean", diag = FALSE, upper = FALSE)))

# lines = DSF grid
distMat$ID <- rownames(distMat)
distMat <- distMat[distMat$ID %in% chal$ID, ]
distMat$ID <- NULL

# columns = "Temp-Volume" grid
distMat <- as.data.frame(t(distMat))
distMat$ID <- rownames(distMat)
distMat <- distMat[distMat$ID %in% df$ID, ]
distMat$ID <- NULL
distMat <- as.data.frame(t(distMat))

# assign DSF data to cells in the "Temp-Volume" grid
for (i in chal$ID){
  # find the nearest "Temp-Volume" cell (from the DSF cell)
  a <- as.data.frame(t(distMat[i,]))
  a$ID <- rownames(a)
  a <- as.data.frame(a[order(a[,1]),])
  a <- a[a[,1]>0,]
  nearest <- a[1, "ID"]

  # year of contamination
  annee <- chal[chal$ID == i,"annee"]
  # assign to "Temp-Volume" cell
  df[df$ID == nearest, "annee"] <- annee
}

# keep only inland cells
df$annee <- as.numeric(df$annee)
df1 <- df[is.na(df$annee),] # keep only NA
df2 <- df[!is.na(df$annee),] # keep only != NA
df3 <- df2[df2$annee != 999,] # remove 999
df4 <- rbind(df1, df3) # reassemble
df <- df4

########
########
######## Attention, en procédant de la sorte on perd des
######## données DSF (car les cellules sont plus grandes)
########
########


######################################################
# Plot
######################################################

# Map settings
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

# Plot chalara
bckgrd+
  theme_bw()+
  theme+
  xlab("longitude") + ylab("latitude")+
  geom_point(data = df, aes(X, Y, col = annee), shape = 15, size = 3, alpha = 0.8)+
  scale_colour_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 2012)+
  geom_polygon(data = fortify(fraxR), aes(long, lat, group = group), color = "grey", fill = "grey", alpha = 0.1)

ggsave("~/Dropbox/chalarose/ashchal/chal.pdf")

# plot Volume
  bckgrd+
    theme_bw()+
    theme+
    xlab("longitude") + ylab("latitude")+
    geom_point(data=df[df$Volume > 0,], aes(X, Y, col=Volume), shape = 15, size = 3)+
    # scale_colour_gradient(limits = c(0,600), colours=c("lightgreen","green4","black"))+
    scale_colour_gradient2(low = "lightgreen", mid = "green4", high = "black", midpoint = 200)+
    # scale_color_manual(values=c("palegreen1", "seagreen1", "seagreen3", "seagreen4"))+
    geom_polygon(data = fortify(fraxR), aes(long, lat, group = group), color = "grey", fill = "grey", alpha = 0.1)

ggsave("~/Dropbox/chalarose/ashchal/volume.pdf")

# plot temperature
    bckgrd+
      theme_bw()+
      theme+
      xlab("longitude") + ylab("latitude")+
      geom_point(data=df, aes(X, Y, col=MMTP_2007), shape = 15, size = 3, alpha = 0.8)+
      scale_colour_gradient2(low = "blue", mid = "green", high = "red", midpoint = 20)+
      geom_polygon(data = fortify(fraxR), aes(long, lat, group = group), color = "grey", fill = "grey", alpha = 0.1)

ggsave("~/Dropbox/chalarose/ashchal/temp.pdf")

####################################################
# Save data
####################################################

save(df, file = '~/Dropbox/chalarose/ashchal/df.rdata')
