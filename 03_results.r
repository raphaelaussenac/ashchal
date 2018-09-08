rm(list=ls(all=TRUE))

######################################################
# Load packages
######################################################

library(ggplot2)
library(ggmap)
library(reshape2)
library(gridExtra)

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
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=12),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.key = element_blank())

######################################################
# import and prepare data
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

# add a column to recieve simulated dates of infection
world$simulAnnee <- NA

######################################################
# distance matrix
######################################################

distMat <- as.data.frame(as.matrix(dist(x = world[,c("X", "Y")], method = "euclidean", diag = FALSE, upper = FALSE)))

######################################################
# parameters' initialisation from results
######################################################

load("~/Dropbox/chalarose/ashchal/paramBackup.rdata")
# set the maximum dispersal distance
mdd <- 2

# set the lambda parameter of the dispersal function
# (to be optimised)
lambda <- paramBackup[1, "lambda"]

# set the bt0 and bt1 parameters of the temperature function
bt0 <- paramBackup[1, "bt0"]
bt1 <- paramBackup[1, "bt1"]

# set the bb0 and bb1 parameters of the Volume function
bb0 <- paramBackup[1, "bb0"]
bb1 <- paramBackup[1, "bb1"]

param <- data.frame(mdd, lambda, bt0, bt1, bb0, bb1)
param$pt <- NA
param$iter <- NA

# save parameters of all runs
paramBackup <- data.frame("mdd" = NA, "lambda" = NA, "bt0" = NA, "bt1" = NA, "bb0" = NA, "bb1" = NA, "pt" = NA, "iter" = NA)

######################################################
# RUN !!!!
######################################################

iter <- 1
nb <- 1
# for (iter in 1:1){  # Number of iterations
  # for (nb in 1:1){  # number of individuals
    # initial state for each individual (individual = model)
    worldInd <- world
    # create an annual loop
    for (annee in 2008:2017){
      # list of infected cells
      if (annee == 2008){
        worldInd[!is.na(worldInd$annee) & worldInd$annee == 2008, "infected"] <- 1
      }
      if (annee == 2010){ # new contamination in the North
        worldInd[!is.na(worldInd$annee) & worldInd$annee == 2010 & worldInd$Y > 49.5, "infected"] <- 1
      }
      infectCells <- rownames(worldInd[worldInd$infected == 1,])
      # for each infected cell, create the list of neighboring cells
      # within the maximum dispersal distance
      for (i in infectCells){
        # create the list of neighboring cells
        neighCells <- t(distMat[i,distMat[i,] < mdd])
        # remove those cells that are already infected
        # for that we first merge the infection and distance information
        neighCells <- merge(worldInd, neighCells, by = "row.names")
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
          neighHealCells$proba <- exp(-param[nb, "lambda"]*neighHealCells$dist) * (1/(1+exp(param[nb, "bb0"]+param[nb, "bb1"]*neighHealCells$Volume))) * (1/(1+exp(param[nb, "bt0"]+param[nb, "bt1"]*neighHealCells$temperature)))
          # draw a random value [0;1]
          neighHealCells$rdmvalue <- runif(nrow(neighHealCells), 0, 1)

          # if proba>rdmvalue then the cell is infected
          neighHealCells$test <- neighHealCells$proba - neighHealCells$rdmvalue
          # create a list of infected neighbors
          infectNeigh <- neighHealCells[neighHealCells$test > 0, "pointID"]
          # update the infected state of cells in the "worldInd" data frame
          worldInd[worldInd$pointID %in% infectNeigh, "infected"] <- 1
          # assign the date of infection to newly infected cells
          worldInd[worldInd$pointID %in% infectNeigh, "simulAnnee"] <- annee + 1
        }
      }
      # plot
      print(bckgrd+
      # theme_bw()+
      theme+
      geom_point(data=worldInd, aes(X, Y, col=annee), shape = 15, size = 3)+
      scale_colour_gradient2(low = "orange", mid = "green", high = "blue", midpoint = 2012)+
      geom_point(data=worldInd[worldInd$infected == 1, ], aes(X, Y), col = "red", shape = 16, size = 1)+
      # coord_fixed()+
      ggtitle(annee))
    }

  ######################################################
  # calculate the performance of the run
  ######################################################

  # index to maximise: score
  # (minimising the sum of difference (obs - predicted dates)
  # this encourages the model not to disseminate)

  # first we select all points with a DSF date of infection
  tabOptim <- worldInd[!is.na(worldInd$annee), c("ID", "annee", "simulAnnee")]
  # difference between observed and predicted dates
  tabOptim$diff <- NA
  tabOptim[!is.na(tabOptim$simulAnnee), "diff"]<- sqrt((tabOptim[!is.na(tabOptim$simulAnnee), "annee"] - tabOptim[!is.na(tabOptim$simulAnnee), "simulAnnee"])^2)
  # assign score depending on the situation
  tabOptim$pt <- NA
  tabOptim[!is.na(tabOptim$diff) & tabOptim$diff == 0, "pt"] <- 10
  tabOptim[!is.na(tabOptim$diff) & tabOptim$diff == 1, "pt"] <- 9
  tabOptim[!is.na(tabOptim$diff) & tabOptim$diff == 2, "pt"] <- 8
  tabOptim[!is.na(tabOptim$diff) & tabOptim$diff == 3, "pt"] <- 7
  tabOptim[!is.na(tabOptim$diff) & tabOptim$diff == 4, "pt"] <- 6
  tabOptim[!is.na(tabOptim$diff) & tabOptim$diff == 5, "pt"] <- 5
  tabOptim[!is.na(tabOptim$diff) & tabOptim$diff == 6, "pt"] <- 4
  tabOptim[!is.na(tabOptim$diff) & tabOptim$diff == 7, "pt"] <- 3
  tabOptim[!is.na(tabOptim$diff) & tabOptim$diff == 8, "pt"] <- 2
  tabOptim[!is.na(tabOptim$diff) & tabOptim$diff == 9, "pt"] <- 1
  tabOptim[!is.na(tabOptim$diff) & tabOptim$diff > 9, "pt"] <- 0
  tabOptim[is.na(tabOptim$diff), "pt"] <- 0
  # sum of the scores
  param[nb, "pt"]  <- sum(tabOptim$pt)
  # indicate iteration
  param[nb, "iter"]  <- iter

######################################################
# Plot results
######################################################

# last map
map <- bckgrd+
theme_bw()+
theme+
geom_point(data=worldInd[!is.na(worldInd$annee),],aes(X,Y, col=annee), shape = 15, size = 3)+
scale_colour_gradient2(low = "orange", mid = "green", high = "blue", midpoint = 2012)+
geom_point(data=worldInd[worldInd$infected == 1, ], aes(X, Y), col = "red", shape = 16, size = 1)+
ggtitle(annee)

# observed - predicted lag
results <- as.data.frame(t(table(tabOptim$diff)))
results <- results[, 2:3]
colnames(results) <- c("lag", "freq")
results$percent <- results$freq * 100 / nrow(tabOptim)
results$cumulPercent <- cumsum(results$percent)

lag <- ggplot(data = results, aes(x = lag, y = freq))+
geom_bar(stat="identity")+
geom_text(aes(label = paste(round(cumulPercent, 2), "%")), col = "red", vjust = 1.5, size = 2)+
theme_bw()

# Convergence plot
load("~/Dropbox/chalarose/ashchal/paramBackup.rdata")

vecMax <- c()
for (i in sort(unique(paramBackup$iter))){
  a <- max(paramBackup[paramBackup$iter %in% c(1:i), 'pt'])
  vecMax <- c(vecMax, a)
}

vecQt <- c()
for (i in sort(unique(paramBackup$iter))){
  a <- quantile(paramBackup[paramBackup$iter %in% c(1:i), 'pt'], 0.7)
  vecQt <- c(vecQt, a)
}

vecMed <- c()
for (i in sort(unique(paramBackup$iter))){
  a <- a <- quantile(paramBackup[paramBackup$iter %in% c(1:i), 'pt'], 0.5)
  vecMed <- c(vecMed, a)
}

vec <- data.frame()
vec <- as.data.frame(t(rbind(vec, vecMax)))
vec <- cbind(vec, vecMed, vecQt)
vec$iter <- c(1:nrow(vec))
colnames(vec) <- c("max", "0.7qt", "0.5qt", "iter")
vec <- melt(vec, id.vars = c("iter"))

converg <- ggplot(data = vec, aes(iter, value, group = variable))+
geom_line(aes(col = variable))+
theme_bw()+
ylab("score")+
xlab("iteration")+
theme(legend.position="bottom", legend.title=element_blank(), legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

# Trace plot
longparam <- melt(paramBackup, id.vars = c("mdd", "pt", "iter"))
trace <- ggplot(data = longparam, aes(iter, value))+
xlab("iteration")+
# geom_point(shape = 16, size = 1)+
geom_smooth(col = "orange")+
facet_wrap(~ variable, ncol = 1, scales = "free", strip.position = "right")+
theme(axis.title.y = element_blank(), strip.background = element_rect(colour = "white", fill = "white"), legend.title=element_blank(), legend.position = "bottom", panel.spacing = unit(0, "lines"))

# set the dispersal function
# lambda <- 2.5
disp <- ggplot(data.frame(x=c(0, 2)), aes(x)) +
stat_function(fun=function(x) exp(-2.5*x), col = "grey")+
stat_function(fun=function(x) exp(-lambda*x), col = "orange")+
ggtitle("dispersal")+
theme(axis.title = element_blank())

# set the temperature logistic function
# bt0 <- -18
# bt1 <- 0.8
tempC <- ggplot(data.frame(x=c(0, 40)), aes(x)) +
stat_function(fun=function(x) 1/(1+exp(-18+0.8*x)), col = "grey")+
stat_function(fun=function(x) 1/(1+exp(bt0+bt1*x)), col = "orange")+
ggtitle("temperature")+
theme(axis.title = element_blank())

# set the volume exponential function
# bb0 <- 3
# bb1 <- -0.2
vol <- ggplot(data.frame(x=c(0, 100)), aes(x)) +
stat_function(fun=function(x) 1/(1+exp(3+-0.2*x)), col = "grey")+
stat_function(fun=function(x) 1/(1+exp(bb0+bb1*x)), col = "orange")+
ggtitle("biomass")+
theme(axis.title = element_blank())

# multiplot
lay <- rbind(c(1, 1, 3, 5),
             c(1, 1, 3, 6),
             c(4, 4, 2, 7))
grid.arrange(map, lag, trace, converg, disp, tempC, vol, layout_matrix = lay)



######################################################
# Exploration synthesis
######################################################
par(mfrow=c(2,2))
# dispersal
curve(exp(-2.5*x), col = "black", from = 0, to = 2, main = "dispersal function", ylab = "", xlab = "lat - long (degree?)")
for (i in 1:nrow(paramBackup)){
  curve(exp(-paramBackup[i, "lambda"]*x), add = TRUE, col = "grey")
}
curve(exp(-lambda*x), add = TRUE, col = "red")

# temperature
curve(1/(1+exp(-18+0.8*x)), col = "black", from = 0, to = 45, main = "temperature limitation function", ylab = "", xlab = "temperature (°C)")
for (i in 1:nrow(paramBackup)){
  curve(1/(1+exp(paramBackup[i, "bt0"]+paramBackup[i, "bt1"]*x)), add = TRUE, col = "grey")
}
curve(1/(1+exp(bt0+bt1*x)), add = TRUE, col = "red")

# biomass
curve(1/(1+exp(3+-0.2*x)), col = "black", from = 0, to = 100, main = "biomass limitation function", ylab= "", xlab = "biomass (m3/ha)")
for (i in 1:nrow(paramBackup)){
  curve(1/(1+exp(paramBackup[i,"bb0"]+paramBackup[i, "bb1"]*x)), add = TRUE, col = "grey")
}
curve(1/(1+exp(bb0+bb1*x)), add = TRUE, col = "red")
 par(mfrow=c(1,1))
