rm(list=ls(all=TRUE))

######################################################
# Load packages
######################################################

library(ggplot2)

######################################################
# create a new world!
######################################################

# central coordinates of cells
X <- Y <- seq(0,300, 8)
world <- expand.grid(X = X, Y = Y)
row.names(world) <- seq(10001, nrow(world)+10000, 1)
world$pointID <- row.names(world)

# assign random state to cells (1 = infected)
world$infected <- 0
world[world$X >= 150 & world$Y >= 250, "infected"] <- 1 # world[round(runif(1, 0, nrow(world)), digits = 0), "infected"] <- 1

# assign temperature to cells (min temp of the hottest 6 hours of the year)
world$temperature <- (world$Y-300) / -10

# assign biomass to cells
world$biomass <- world$X * world$Y / 1000

# plot chalarose
ggplot(data=world, aes(X, Y))+
geom_raster(aes(fill=infected), interpolate=F)+
scale_fill_gradientn(limits = c(0,1), colours=c("green4", "black"))

# plot temperature
ggplot(data=world, aes(X, Y))+
geom_raster(aes(fill=temperature), interpolate=F)+
scale_fill_gradientn(limits = c(0,30), colours=c("blue","green","red"))

# plot biomass
ggplot(data=world, aes(X, Y))+
geom_raster(aes(fill=biomass), interpolate=F)+
scale_fill_gradientn(limits = c(0,100), colours=c("lightgreen","green4","black"))

######################################################
# distance matrix
######################################################

distMat <- as.data.frame(as.matrix(dist(x = world[,c("X", "Y")], method = "euclidean", diag = FALSE, upper = FALSE)))

######################################################
# dispersal function
######################################################

# set the maximum dispersal distance
mdd <- 200

# set the lambda parameter of the dispersal function
# (to be optimised)
lambda <- 0.1

# set the bt0 and bt1 parameters of the temperature function
bt0 <- -35
bt1 <- 0.8

# set the bb0 and bb1 parameters of the biomass function
bb0 <- 3
bb1 <- -0.5

# create an annual loop
for (yr in 2000:2010){
  # list of infected cells
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
      # exp(-lambda*dist) * (1/(1+exp(bt0+bt1*temperature))) * (1/(1+exp(bb0+bb1*temperature)))
      neighHealCells$proba <- exp(-lambda*neighHealCells$dist) * (1/(1+exp(bt0+bt1*neighHealCells$temperature))) * (1/(1+exp(bb0+bb1*neighHealCells$biomass)))
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
  print(ggplot(data=world, aes(X, Y))+
  geom_raster(aes(fill=infected), interpolate=F)+
  scale_fill_gradientn(limits = c(0,1), colours=c("lightgreen","green4","black"))+
  ggtitle(yr))
  print(yr)
}










#
# ######################################################
# # set the dispersal function
# ######################################################
# lambda <- 0.05
# curve(exp(-lambda*x), from = 0, to = 100)
#
# ######################################################
# # set the temperature logistic function
# ######################################################
# bt0 <- -15
# bt1 <- 0.8
# curve(1/(1+exp(bt0+bt1*x)), from = 0, to = 40)
#
# ######################################################
# # set the biomass exponential function
# ######################################################
# bb0 <- 3
# bb1 <- -0.5
# curve(1/(1+exp(bb0+bb1*x)), from = 0, to = 40)
