rm(list=ls(all=TRUE))

######################################################
# Load packages
######################################################

library(ggplot2)
library(ggmap)
library(doParallel)

######################################################
# Map settings
######################################################

# # Location
# lat <- c(42, 51.5)
# lon <- c(-5, 8.5)
#
# # Download background
# map <- get_map(location = c(lon = mean(lon), lat = mean(lat)), zoom = 5,
#                maptype = "toner-background", source = "google")
#
# # set background
# bckgrd <- ggmap(map)+
#   scale_x_continuous(limits = lon, expand = c(0, 0)) +
#   scale_y_continuous(limits = lat, expand = c(0, 0))
#
# # theme settings
# theme=theme(panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank(),
#             strip.background = element_blank(),
#             panel.background = element_rect(fill = 'white'),
#             legend.position = c(0,0),
#             legend.justification=c(-0.05,-0.05),
#             text = element_text(size=12),
#             axis.text.x = element_text(size=10),
#             legend.key = element_blank())

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
# parameters' initialisation
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

# for X individuals
nbInd <- 10
mdd <- rnorm(nbInd,mdd, mdd*0.1)
lambda <- rnorm(nbInd,lambda, lambda*0.1)
bt0 <- rnorm(nbInd,bt0, bt0*-0.1)
bt1 <- rnorm(nbInd,bt1, bt1*0.1)
bb0 <- rnorm(nbInd,bb0, bb0*0.1)
bb1 <- rnorm(nbInd,bb1, bb1*-0.1)

param <- data.frame(mdd, lambda, bt0, bt1, bb0, bb1)
param$sumDiff <- NA
param$iter <- NA

# save parameters of all runs
paramBackup <- data.frame("mdd" = NA, "lambda" = NA, "bt0" = NA, "bt1" = NA, "bb0" = NA, "bb1" = NA, "sumDiff" = NA)

######################################################
# RUN !!!!
######################################################


start_time <- Sys.time()

for (iter in 1:2){  # Number of iterations
  for (nb in 1:nbInd){  # number of individuals
    # create an annual loop
    for (annee in 2008:2017){
      # list of infected cells
      if (annee == 2008){
        world[!is.na(world$annee) & world$annee == 2008, "infected"] <- 1
      }
      if (annee == 2010){ # new contamination in the North
        world[!is.na(world$annee) & world$annee == 2010 & world$Y > 49.5, "infected"] <- 1
      }
      infectCells <- rownames(world[world$infected == 1,])
      # for each infected cell, create the list of neighboring cells
      # within the maximum dispersal distance
      for (i in infectCells){
        # create the list of neighboring cells
        neighCells <- t(distMat[i,distMat[i,] < param[nb, "mdd"]])
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
          neighHealCells$proba <- exp(-param[nb, "lambda"]*neighHealCells$dist) * (1/(1+exp(param[nb, "bb0"]+param[nb, "bb1"]*neighHealCells$Volume))) * (1/(1+exp(param[nb, "bt0"]+param[nb, "bt1"]*neighHealCells$temperature)))
          # draw a random value [0;1]
          neighHealCells$rdmvalue <- runif(nrow(neighHealCells), 0, 1)

          # if proba>rdmvalue then the cell is infected
          neighHealCells$test <- neighHealCells$proba - neighHealCells$rdmvalue
          # create a list of infected neighbors
          infectNeigh <- neighHealCells[neighHealCells$test > 0, "pointID"]
          # update the infected state of cells in the "world" data frame
          world[world$pointID %in% infectNeigh, "infected"] <- 1
          # assign the date of infection to newly infected cells
          world[world$pointID %in% infectNeigh, "simulAnnee"] <- annee + 1
        }
      }
      # # plot
      # print(bckgrd+
      # theme_bw()+
      # theme+
      # geom_point(data=world[world$Volume != 0, ], aes(X, Y, col=Volume), shape = 15, size = 3)+
      # scale_colour_gradient2(low = "lightgreen", mid = "green4", high = "black", midpoint = 200)+
      # geom_point(data=world[world$infected == 1, ], aes(X, Y), col = "red", shape = 16, size = 2)+
      # # coord_fixed()+
      # ggtitle(annee))
      print(annee)
      print(nb)
    }

  ######################################################
  # calculate the performance of the run
  ######################################################

  # index to minimise: sum of the differences between
  # DSF dates of infection and simulated dates

  # first we select all points with a DSF date of infection
  tabOptim <- world[!is.na(world$annee), c("ID", "annee", "simulAnnee")]
  # then we delete those points for which there are no simulated infection dates
  # (should be a tiny number of cells if the model runs on wide temporal window)
  tabOptim <- tabOptim[!is.na(tabOptim$simulAnnee),]
  tabOptim$diff <- sqrt((tabOptim$annee - tabOptim$simulAnnee)^2)
  # sum of the differences
  param[nb, "sumDiff"]  <- sum(tabOptim$diff)
  # indicate iteration
  param[nb, "iter"]  <- iter

  # save parameters and performance
  paramBackup <- rbind(paramBackup, param[nb, ])
  }

  ######################################################
  # Genetic operators (selection, crossover, mutation)
  ######################################################

  # selection
  param <- param[order(param$sumDiff),]
  selected <- param[1:4,]

  # crossover
  crossed <- param[5:7,]
  # create an empty table to recieve randomised values
  cross <- as.data.frame(matrix(ncol = 7, nrow = 3))
  colnames(cross) <- colnames(param)
  # for each parameter we randomly exchange values among individuals
  # first create a list of possible random values
  rdm <- c(1, 2, 3)
  for (np in 1:ncol(crossed-1)){  # each column (parameter)
    for (ni in 1:nrow(cross)){  # each line (individual)
      draw <- round(runif(1, 1, length(rdm)), 0)  # random number in rdm
      cross[ni, np] <- crossed[rdm[draw], np]  # write the randmly selected value in the cross table
      rdm <- rdm[-draw]  # delete that value in the remaining choices
    }
    rdm <- c(1, 2, 3)
  }
  crossed <- cross

  # mutation
  mutated <- param[8:10,]
  # add random value to positive parameters
  for (np in c("mdd", "lambda", "bt1","bb0")){
    mutated[, np] <- mutated[, np] + rnorm(3, 0, mutated[, np] * 0.3)
  }
  # add random value to negative parameters
  for (np in c("bt0", "bb1")){
    mutated[, np] <- mutated[, np] + rnorm(3, 0, mutated[, np] * -0.3)
  }

  # combine new set of parameters
  param <- rbind(selected, crossed, mutated)
  param$sumDiff <- NA

}
end_time <- Sys.time()
# time elapsed
end_time - start_time

# results
paramBackup <- paramBackup[-1,]
hist(paramBackup$sumDiff)


save(paramBackup, file = '~/Dropbox/chalarose/ashchal/paramBackup.rdata')

# ####################################################
# ## Parallel
# ####################################################
# # number of workers
# registerDoParallel(cores = 2)
# getDoParWorkers()
#
# start_time <- Sys.time()
# out <- foreach(i = 1:nrow(param) ,.combine='rbind') %dopar% {
#   individualScore()
# }
# end_time <- Sys.time()
#
# # time elapsed
# end_time - start_time
#
#
#
# #
# out <- out[!is.na(out$mdd),]
#
# out1 <- out[order(out$sumDiff), ]
#
# # save(out,file='./out.rdata')







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