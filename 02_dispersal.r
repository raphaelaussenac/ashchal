rm(list=ls(all=TRUE))

######################################################
# Load packages
######################################################

library(ggplot2)
library(ggmap)
# library(doParallel)

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
lambda <- 0.54

# set the bt0 and bt1 parameters of the temperature function
bt0 <- -25.06
bt1 <- 1.07

# set the bb0 and bb1 parameters of the Volume function
bb0 <- 2.11
bb1 <- -0.08

# for X individuals
nbInd <- 10
# first numbers are approximate optimums
lambda <- c(lambda, 0.067, rnorm(nbInd-1, lambda, sqrt((lambda * 0.5)^2)))
bt0 <- c(bt0, -25.06, rnorm(nbInd-2, bt0, sqrt((bt0 * 0.3)^2)))
bt1 <- c(bt1, 1.070, rnorm(nbInd-2, bt1, sqrt((bt1 * 0.3)^2)))
bb0 <- c(bb0, 2.7123, rnorm(nbInd-2, bb0, sqrt((bb0 * 0.3)^2)))
bb1 <- c(bb1, -0.06437, rnorm(nbInd-2, bb1, sqrt((bb1 * 0.3)^2)))

param <- data.frame(mdd, lambda, bt0, bt1, bb0, bb1)
param$pt <- NA
param$iter <- NA

# save parameters of all runs
paramBackup <- data.frame("mdd" = NA, "lambda" = NA, "bt0" = NA, "bt1" = NA, "bb0" = NA, "bb1" = NA, "pt" = NA, "iter" = NA)

######################################################
# RUN !!!!
######################################################

# number of iteration
nbIter <- 5

# progress bar
pb = txtProgressBar(min = 0, max = nbIter, initial = 0, style = 3)

start_time <- Sys.time()
for (iter in 1:nbIter){  # Number of iterations
  for (nb in 1:nbInd){  # number of individuals
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
      # # plot
      # print(bckgrd+
      # theme_bw()+
      # theme+
      # geom_point(data=worldInd, aes(X, Y, col=annee), shape = 15, size = 3)+
      # scale_colour_gradient2(low = "orange", mid = "green", high = "blue", midpoint = 2012)+
      # geom_point(data=worldInd[worldInd$infected == 1, ], aes(X, Y), col = "red", shape = 16, size = 1)+
      # # coord_fixed()+
      # ggtitle(annee))
      print(paste("iter = ", iter, "/" , nbIter , " ; annee: ", annee, " ; individu: ", nb, "/" , nbInd , sep = ""))
    }

  ######################################################
  # calculate the performance of the run
  ######################################################

  # index to maximise: score
  tabOptim <- worldInd[, c("ID", "annee", "simulAnnee")]
  # difference between observed and predicted dates
  tabOptim$diff <- NA
  tabOptim[!is.na(tabOptim$simulAnnee), "diff"]<- sqrt((tabOptim[!is.na(tabOptim$simulAnnee), "annee"] - tabOptim[!is.na(tabOptim$simulAnnee), "simulAnnee"])^2)
  # assign score depending on the lag between observed and predicted dates
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
  # model doesn't infect a cell that has no chalara on the field --> 0
  tabOptim[is.na(tabOptim$annee) & is.na(tabOptim$simulAnnee), "pt"] <- 0
  # doesn't infect a cell that has chalara --> 0
  tabOptim[!is.na(tabOptim$annee) & is.na(tabOptim$simulAnnee), "pt"] <- -10
  # infect a cell that has no chalara on the field --> -10
  tabOptim[is.na(tabOptim$annee) & !is.na(tabOptim$simulAnnee), "pt"] <- -10
  # sum of the scores
  param[nb, "pt"]  <- sum(tabOptim$pt)
  # indicate iteration
  param[nb, "iter"]  <- iter

  # save parameters and performance
  paramBackup <- rbind(paramBackup, param[nb, ])
  print(param)
  }

  ######################################################
  # Genetic operators (selection, crossover, mutation)
  ######################################################

  # selection
  # First define the probability (for each ind) of being
  # part of the reproducing population (proba is
  # is proportional to score)
  if (min(param$pt) < 0){
    param$proba <- (param$pt + sqrt(min(param$pt)^2)) / (max(param$pt) + sqrt(min(param$pt)^2))
  } else {
    param$proba <- param$pt / max(param$pt)
  }

  selected <- data.frame()
  # till we have nbInd individuals
  while (nrow(selected) < nbInd){
    # Then drawing with replacement
    # for that we first assign random numbers [0; 1]
    # to each individual
    param$rdmnb <- runif(nrow(param), 0, 1)
    # If the proba > rdmnb, the individual is selected
    # for reproduction
    param$select <- param$proba - param$rdmnb
    selected <- rbind(selected, param[param$select > 0, ])
  }
  selected <- selected[1:nbInd,1:8]

  # crossover
  # crossover probability = 0.8 (applies to allels)
  # divide the reproducing population into 2 sets of individuals
  male <- selected[1:(nbInd/2), c("lambda", "bt0","bt1", "bb0", "bb1")]
  female <- selected[((nbInd/2)+1):nbInd, c("lambda", "bt0","bt1", "bb0", "bb1")]

  crossed <- data.frame()
  for (i in 1:nrow(male)){  # for each couple
    # create couples
    couple <- rbind(male[i,], female[i,])
    couple <- as.data.frame(t(couple))
    # assign random number [0, 1] to each parameter
    rdmnb <- runif(5, 0, 1)
    couple <- cbind(couple, rdmnb)
    couple$crossProba <- 0.8
    couple$diff <- couple$crossProba - couple$rdmnb
    for (j in 1:nrow(couple)){  # for each allel
      # crossover if rdmnb < crossProba
      if (couple[j, "diff"] > 0){
        a <- couple[j, 1]
        couple[j, 1] <- couple[j, 2]
        couple[j, 2] <- a
      }
    }
    couple <- as.data.frame(t(couple))
    couple <- couple[1:2,]
    crossed <- rbind(crossed, couple)
  }

  # mutation
  # mutation probability = 0.2 (applies to allels)
  for (i in colnames(crossed)){
    for (j in 1:nrow(crossed)){
      rdmnb <- runif(1, 0, 1)
      if (rdmnb < 0.2){
        crossed[j, i] <- crossed[j, i] + rnorm(1, 0, sqrt((crossed[j, i] * 0.3)^2))
      }
    }
  }
  mutated <- crossed

  # elitist strategy suspended becaus best fitted most often already selected
  #
  # # elitist strategy
  # # repace one of the new individuals by the best individual of
  # # the previous generation
  # rdmnb <- round(runif(1, 1, ncol(mutated)), 0)
  # elitInd <- param[param$pt == max(param$pt), c("lambda", "bt0","bt1", "bb0", "bb1")]
  # mutated[rdmnb, ] <- elitInd

  # format new parameters
  mutated$pt <- NA
  mutated$iter <- NA
  a <- runif(nrow(mutated), mdd, mdd)
  mutated <- cbind(a, mutated)
  colnames(mutated)[1] <- "mdd"
  param <- mutated

  ######################################################
  # Follow progression
  ######################################################

  # progress bar
  # Sys.sleep(0.01)
  setTxtProgressBar(pb,iter)
  middle_time <- Sys.time()
  print(middle_time - start_time)

  # Algorithme convergence plot
  vecMax <- c()
  for (i in sort(unique(paramBackup$iter))){
    a <- max(paramBackup[paramBackup$iter %in% c(1:i), 'pt'], na.rm = T)
    vecMax <- c(vecMax, a)
  }

  vecQt <- c()
  for (i in sort(unique(paramBackup$iter))){
    a <- quantile(paramBackup[paramBackup$iter %in% c(1:i), 'pt'], 0.7, na.rm = T)
    vecQt <- c(vecQt, a)
  }

  vecMed <- c()
  for (i in sort(unique(paramBackup$iter))){
    a <- quantile(paramBackup[paramBackup$iter %in% c(1:i), 'pt'], 0.5, na.rm = T)
    vecMed <- c(vecMed, a)
  }

  plot(vecMax, type = "l", col = "red", ylim = c(min(vecMed, vecMax, vecQt), max(vecMed, vecMax, vecQt)))
  lines(vecMed, type = "l", col = "orange")
  lines(vecQt, type = "l", col = "grey")

}
end_time <- Sys.time()
# time elapsed
end_time - start_time

# results
paramBackup <- paramBackup[-1,]
paramBackup <- paramBackup[order(-paramBackup$pt),]

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
# out1 <- out[order(out$pt), ]
#
# # save(out,file='./out.rdata')
