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

# assign random state to cells (1 = infected)
world$infected <- round(runif(nrow(world), 0, 1), digits = 0)

# initial conditions
ggplot(data=world, aes(X, Y))+
geom_raster(aes(fill=infected), interpolate=F)+
scale_fill_gradientn(limits = c(0,1), colours=c("lightgreen","green4","black"))

######################################################
# set the dispersal function
######################################################
lambda <- 0.1
curve(exp(-lambda*x), from = 0, to = 100)
