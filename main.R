# This is the main file for multiple files
# starting with GA_***.R
#
# develop a genetic algorithm for optimized grid adaptation
# and Analog Ensemble technique
#
library(sp)
library(GA)
library(maps)
library(gstat)
library(ncdf4)
library(rgdal)
library(fields)
library(raster)
library(deldir)
library(spatstat)
library(maptools)
library(stringr)
library(RColorBrewer)
library(velox)
library(RAnEn)

# Run this if you haven't installed RAnEn
#
# install.packages('RAnEn_2.1.1.tar.gz', type = 'source', repo = NULL)

# whether to penalize distance in GA evolution
distance.penalty <- F

# visualize data
visualize <- T

# save every stage GA population
options(save.population.stage = F)

# pre-computed AnEn
precomputed.AnEn <- 'AnEn.nc'

# specify outout for GA figures
# options(GA_output_form = 'jpeg')
options(GA_output_form = NULL)

# develop unstructured grid on test day index and flt
analog.day.index <- 1
flt.index <- 1

# resample brewer colors
brewer.cols <- colorRampPalette(
  brewer.pal(11, 'Spectral')[11:1])(100)

source('functions.R')

parallel <- T


########################################################
#                     prepare data                     #
########################################################
source('GA_real-data.R')

vec.obs <- as.vector(observations[
  1, , test.ID.start + analog.day.index - 1, flt.index])
vec.anen <- as.vector(apply(analogs[
  , analog.day.index, flt.index, , 1, drop = F], 1, mean))
vec.fcst <- as.vector(forecasts[
  forecast.var, , test.ID.start + analog.day.index - 1, flt.index])

# create raster objects
x.min <- range(xs)[1]
x.max <- range(xs)[2]
y.min <- range(ys)[1]
y.max <- range(ys)[2]

rast.base <- raster(
  nrows = ny, ncols = nx,
  xmn = x.min, xmx = x.max,
  ymn = y.min, ymx = y.max)
rast.obs <- raster(
  nrows = ny, ncols = nx,
  xmn = x.min, xmx = x.max,
  ymn = y.min, ymx = y.max,
  vals = vec.obs)
rast.anen <- raster(
  nrows = ny, ncols = nx,
  xmn = x.min, xmx = x.max, 
  ymn = y.min, ymx = y.max,
  vals = vec.anen)
rast.fcst <- raster(
  nrows = ny, ncols = nx,
  xmn = x.min, xmx = x.max, 
  ymn = y.min, ymx = y.max,
  vals = vec.fcst)

# flip rasters upside down
rast.obs <- flip(rast.obs, 2)
rast.anen <- flip(rast.anen, 2)
rast.fcst <- flip(rast.fcst, 2)


########################################################
#                    algorithm design                  #
########################################################
#------------------ genetic algorithm -----------------#
#
source('GA_core.R')

suggestions <- matrix(NA, nrow = 0, ncol = num.selected*2)
for (i in 1:pop.size) {
  indv <- as.vector(rbind(
    runif(num.selected, min = x.min , max = x.max),
    runif(num.selected, min = y.min , max = y.max)))
  suggestions <- rbind(suggestions, indv)
}
save(suggestions, file = 'suggestions.RData')

# GA with real value optimization
GA <- ga.new(
  type = 'real-valued',
  
  fitness = fitness,
  interpolation = interpolation, n_nei = n_nei,
  rast.base = rast.base, rast.obs = rast.obs,
  rast.anen = rast.anen, profile = F,
  distance.penalty = distance.penalty,
  
  monitor = plot.map,
  selection = func_selection,
  mutation = func_mutation,
  
  popSize = pop.size, maxiter = num.iter, seed = 1,
  pmutation = pmutation, pcrossover = pcrossover,
  lower = rep(c(x.min, y.min), num.selected),
  upper = rep(c(x.max, y.max), num.selected),
  keepBest = T, parallel = parallel, suggestions = suggestions)

ga.fitness <- GA@fitness
ga.pts <- matrix(GA@solution, ncol = 2, byrow = T)


########################################################################
#                         Pure random                                  #
########################################################################
# Prevent algorithm output
options(GA_output_form = NULL)

library(parallel)
library(foreach)
library(doParallel)

cl <- makeCluster(detectCores())
registerDoParallel(cl)

Population <- suggestions
lower = rep(c(x.min, y.min), num.selected)
upper = rep(c(x.max, y.max), num.selected)


# Calculate fitness
Fitness <- foreach(i. = 1:nrow(Population), .combine = "c") %dopar% 
{
  library(sp);library(GA);library(maps);library(RAnEn)
  library(gstat);library(ncdf4);library(rgdal)
  library(fields);library(raster);library(deldir)
  library(spatstat);library(maptools)
  
  fitness(Population[i., ],
          interpolation = interpolation, n_nei = n_nei,
          rast.base = rast.base, rast.obs = rast.obs,
          rast.anen = rast.anen, profile = F,
          distance.penalty = distance.penalty)
}

for (iter in 1:num.iter) {
  print(paste("Iteration", iter))

  # Random selection
  sel <- sample(1:length(Fitness), size = pop.size, replace = TRUE)
  Population = Population[sel, , drop = FALSE]
  Fitness = Fitness[sel]
  
  # Save data
  if (getOption("save.population.stage")) {
    file = paste('random-data-iteration_',
                 str_pad(iter, 3, pad = '0'), '.RData', sep = '')
    save(Population, Fitness, file = file)
  }
  
  # Random mutation
  num.mutation <- floor(pop.size*pmutation)
  s <- sample(1:pop.size, size = num.mutation, replace = T)
  mutate <- parent <- Population[s, ]
  dampeningFactor <- 1.1 - iter/num.iter
  value.max <- apply(parent, 1, max, rm.na = T) - 
    apply(parent, 1, min, rm.na = T)
  
  # each gene should have a possibility to be mutated
  gene.prob <- .5
  selected <- runif(ncol(parent)) < gene.prob
  
  # mutate genes
  for (row in 1:nrow(mutate)) {
    mutate[row, selected] <- rnorm(
      length(which(selected)), mean = parent[row, selected],
      sd = value.max/2/ncol(mutate)*dampeningFactor)
    
    selected <- mutate[row, ] < lower
    mutate[row, selected] <- lower[selected]
    selected <- mutate[row, ] > upper
    mutate[row, selected] <- upper[selected]
  }
  
  mutate.Fitness <- foreach(i. = 1:nrow(mutate), .combine = "c") %dopar% 
  {
    library(sp);library(GA);library(maps);library(RAnEn)
    library(gstat);library(ncdf4);library(rgdal)
    library(fields);library(raster);library(deldir)
    library(spatstat);library(maptools)
    
    fitness(mutate[i., ],
            interpolation = interpolation, n_nei = n_nei,
            rast.base = rast.base, rast.obs = rast.obs,
            rast.anen = rast.anen, profile = F,
            distance.penalty = distance.penalty)
  }
  
  Population <- rbind(Population, mutate)
  Fitness <- c(Fitness, mutate.Fitness)
}

stopCluster(cl)


###########################################################################
#                            Random w/ Elitism                            #
###########################################################################
library(parallel)
library(foreach)
library(doParallel)

cl <- makeCluster(detectCores())
registerDoParallel(cl)

Population <- suggestions
lower = rep(c(x.min, y.min), num.selected)
upper = rep(c(x.max, y.max), num.selected)


# Calculate fitness
Fitness <- foreach(i. = 1:nrow(Population), .combine = "c") %dopar% 
{
  library(sp);library(GA);library(maps);library(RAnEn)
  library(gstat);library(ncdf4);library(rgdal)
  library(fields);library(raster);library(deldir)
  library(spatstat);library(maptools)
  
  fitness(Population[i., ],
          interpolation = interpolation, n_nei = n_nei,
          rast.base = rast.base, rast.obs = rast.obs,
          rast.anen = rast.anen, profile = F,
          distance.penalty = distance.penalty)
}

for (iter in 1:num.iter) {
  print(paste("Iteration", iter))
  
  #--------------------------------- Elitism ------------------------------#
  # Random selection
  best <- order(Fitness, decreasing = T)[1]
  best.fitess <- Fitness[best]
  best.indv <- Population[best, ]
  
  sel <- sample(1:length(Fitness), size = pop.size-1, replace = TRUE)
  Population = Population[sel, , drop = FALSE]
  Fitness = Fitness[sel]
  
  # Add the best from last iteration
  Population <- rbind(Population, best.indv)
  Fitness <- c(Fitness, best.fitess)
  #------------------------------------------------------------------------#
  
  # Save data
  if (getOption("save.population.stage")) {
    file = paste('elitism-data-iteration_',
                 str_pad(iter, 3, pad = '0'), '.RData', sep = '')
    save(Population, Fitness, file = file)
  }
  
  # Random mutation
  num.mutation <- floor(pop.size*pmutation)
  s <- sample(1:pop.size, size = num.mutation, replace = T)
  mutate <- parent <- Population[s, ]
  dampeningFactor <- 1.1 - iter/num.iter
  value.max <- apply(parent, 1, max, rm.na = T) - 
    apply(parent, 1, min, rm.na = T)
  
  # each gene should have a possibility to be mutated
  gene.prob <- .5
  selected <- runif(ncol(parent)) < gene.prob
  
  # mutate genes
  for (row in 1:nrow(mutate)) {
    mutate[row, selected] <- rnorm(
      length(which(selected)), mean = parent[row, selected],
      sd = value.max/2/ncol(mutate)*dampeningFactor)
    
    selected <- mutate[row, ] < lower
    mutate[row, selected] <- lower[selected]
    selected <- mutate[row, ] > upper
    mutate[row, selected] <- upper[selected]
  }
  
  mutate.Fitness <- foreach(i. = 1:nrow(mutate), .combine = "c") %dopar% 
  {
    library(sp);library(GA);library(maps);library(RAnEn)
    library(gstat);library(ncdf4);library(rgdal)
    library(fields);library(raster);library(deldir)
    library(spatstat);library(maptools)
    
    fitness(mutate[i., ],
            interpolation = interpolation, n_nei = n_nei,
            rast.base = rast.base, rast.obs = rast.obs,
            rast.anen = rast.anen, profile = F,
            distance.penalty = distance.penalty)
  }
  
  Population <- rbind(Population, mutate)
  Fitness <- c(Fitness, mutate.Fitness)
}

stopCluster(cl)


# Visualization #
source('different-DOUG.R')
