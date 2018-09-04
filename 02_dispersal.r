rm(list=ls(all=TRUE))

######################################################
# Load packages
######################################################

library(ggplot2)
library(ggmap)

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

######################################################
# import volume data
######################################################

load("~/Dropbox/chalarose/ashchal/df.rdata")
world <- df

# remove Corsica Island
world <- world[world$X < 8.5,]

# assign ID to each point
row.names(world) <- seq(10001, nrow(world)+10000, 1)
world$pointID <- row.names(world)

# add the infection status
world$infected <- 0

######################################################
# distance matrix
######################################################

distMat <- as.data.frame(as.matrix(dist(x = world[,c("X", "Y")], method = "euclidean", diag = FALSE, upper = FALSE)))

######################################################
# dispersal function
######################################################

# set the maximum dispersal distance
mdd <- 2

# set the lambda parameter of the dispersal function
# (to be optimised)
lambda <- 2.5

# set the bt0 and bt1 parameters of the temperature function
bt0 <- -18
bt1 <- 0.8

# set the bb0 and bb1 parameters of the Volume function
bb0 <- 3
bb1 <- -0.2

# create an annual loop
for (annee in 2008:2015){
  # list of infected cells
  if (annee == 2008){
    world[!is.na(world$annee) & world$annee == 2008, "infected"] <- 1
  }
  infectCells <- rownames(world[world$infected == 1,])
  # for each infected cell, create the list of neighboring cells
  # within the maximum dispersal distance
  for (i in infectCells){
    # create the list of neighboring cells
    neighCells <- t(distMat[i,distMat[i,] < mdd])
    # remove those cells that are already infected
    # for that we first merge the infection and distance information
    neighCells <- merge(world, neighCells, by = "row.names")
    # then we keep only healthy cells
    # if there are some left
    if (length(neighCells[neighCells$infected == 0, 1])>0){
      neighHealCells <- neighCells[neighCells$infected == 0, ]
      colnames(neighHealCells)[ncol(neighHealCells)] <- "dist"

      # assign each cell a probability of being contaminated using
      # the dispersal and temperature function:
      # first we must define the temperature we are going to use
      # (previous year --> annee - 1)
      colnames(neighHealCells)[substr(colnames(neighHealCells), 6, 9) == (annee - 1)] <- "temperature"
      # then we calculate the probability of being contaminated
      # exp(-lambda*dist) * (1/(1+exp(bt0+bt1*temperature))) * (1/(1+exp(bb0+bb1*volume)))
      neighHealCells$proba <- exp(-lambda*neighHealCells$dist) * (1/(1+exp(bb0+bb1*neighHealCells$Volume))) * (1/(1+exp(bt0+bt1*neighHealCells$temperature)))
      # draw a random value [0;1]
      neighHealCells$rdmvalue <- runif(nrow(neighHealCells), 0, 1)

      # if proba>rdmvalue then the cell is infected
      neighHealCells$test <- neighHealCells$proba - neighHealCells$rdmvalue
      # create a list of infected neighbors
      infectNeigh <- neighHealCells[neighHealCells$test > 0, "pointID"]
      # update the infected state of cells in the "world" data frame
      world[world$pointID %in% infectNeigh, "infected"] <- 1
    }
  }
  # plot
  print(bckgrd+
  theme_bw()+
  theme+
  geom_point(data=world[world$Volume != 0, ], aes(X, Y, col=Volume), shape = 15, size = 3)+
  scale_colour_gradient2(low = "lightgreen", mid = "green4", high = "black", midpoint = 200)+
  geom_point(data=world[world$infected == 1, ], aes(X, Y), col = "red", shape = 16, size = 2)+
  # coord_fixed()+
  ggtitle(annee))
  print(annee)
}

#
# ######################################################
# # set the dispersal function
# ######################################################
# lambda <- 2.5
# curve(exp(-lambda*x), from = 0, to = 2)
#
# ######################################################
# # set the temperature logistic function
# ######################################################
# bt0 <- -18
# bt1 <- 0.8
# curve(1/(1+exp(bt0+bt1*x)), from = 0, to = 40)
#
# ######################################################
# # set the volume exponential function
# ######################################################
# bb0 <- 3
# bb1 <- -0.2
# curve(1/(1+exp(bb0+bb1*x)), from = 0, to = 100)
