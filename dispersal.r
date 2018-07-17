rm(list=ls(all=TRUE))

######################################################
# Load packages
######################################################

library(ggplot2)

######################################################
# create a new world!
######################################################

# central coordinates of cells
X <- Y <- seq(0,1000, 8)
world <- expand.grid(X = X, Y = Y)
row.names(world) <- seq(10001, nrow(world)+10000, 1)
world$pointID <- row.names(world)

# assign random state to cells (1 = infected)
world$infected <- 0
world[round(runif(3, 0, nrow(world)), digits = 0), "infected"] <- 1

# plot initial conditions
ggplot(data=world, aes(X, Y))+
geom_raster(aes(fill=infected), interpolate=F)+
scale_fill_gradientn(limits = c(0,1), colours=c("lightgreen","green4","black"))

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
      # the dispersal function: exp(-lambda*dist)
      neighHealCells$proba <- exp(-lambda*neighHealCells$dist)
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











######################################################
# set the dispersal function
######################################################
lambda <- 0.1
curve(exp(-lambda*x), from = 0, to = 100)

######################################################
# set the temperature function logistique
######################################################
b0 <- 15   # seuil de température a partir duquel on passe à un proba de 1
b1 <- -0.9  # doit être négatif pour que les temperatures élevées aient une proba associé de 0
curve(1/(1+exp(b0+b1*x)), from = 0, to = 40)
