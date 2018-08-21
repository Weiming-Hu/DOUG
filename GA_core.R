# GA core functions for experiments
#
# sourced by GA_develop-unstructured-grid.R
#
library(GA)

source('functions.R')

plot.map <- function(object, digits = getOption("digits"), ...) {
  require(raster)
  
  fitness <- object@fitness
  
  best.index <- order(object@fitness, decreasing = T)[1]
  worst.index <- order(object@fitness)[1]
  
  best.coords <- matrix(object@population[best.index, ],
                        ncol = 2, byrow = T)
  
  worst.coords <- matrix(object@population[worst.index, ],
                         ncol = 2, byrow = T)
  
  xlab <- paste('Best  fitness =', round(fitness[best.index], 3), ';',
                'Worst fitness =', round(fitness[worst.index], 3))
  
  flag <- getOption('GA_output_form')
  if (!is.null(flag)) {
    require(stringr)
    if (flag == 'pdf') {
      pdf(paste('GA-iteration-plot_',
                str_pad(object@iter, 3, pad = '0'), '.pdf', sep = ''),
          compress = F, height = 10, width = 20)
    } else {
      jpeg(paste('GA-iteration-plot_',
                 str_pad(object@iter, 3, pad = '0'), '.jpeg', sep = ''),
           height = 10, width = 20,
           units = 'in', res = 200)
    }
  }
  
  par(mfrow = c(1, 2), mar =  c(5, 4, 4, 5) + 0.1)
  
  # plot(rast.obs, xlab = xlab,
  #      main = paste('Best and Worst of Iteration',
  #                   object@iter, '(Base map: observation)'),
  #      col = brewer.cols)
  plot(rast.anen, xlab = xlab,
       col = brewer.cols)
  
  points(best.coords, pch = 16, cex = 0.8, col = 'black')
  #points(worst.coords, pch = 1, cex = 1, col = 'grey', lwd = 1.5)
  legend('top', col = 'black', pch = 16,
         legend = 'best', horiz = T)
  
  plot(object, ...)
  
  if (!is.null(flag)) {
    dev.off()
  }
  
}

fitness <- function(
  indv, interpolation, n_nei,
  rast.base, rast.obs, rast.anen,
  profile = F,
  distance.penalty = F) {
  # defines the fitness for an individual
  #
  # the input indv is a vector with xs and ys of the 
  # selected points. For example, for a vector with 
  # length of 6, the values are:
  #             (x1, y1, x2, y2, x3, y3)
  #
  
  if (profile) {
    time.start <- Sys.time()
  }
  
  coords <- matrix(data = indv, ncol = 2, byrow = T)
  coords <- cbind(coords, colFromX(rast.base, coords[, 1]))
  coords <- cbind(coords, rowFromY(rast.base, coords[, 2]))
  colnames(coords) <- c('x', 'y', 'col', 'row')
  
  # add four corners
  four.corners <- matrix(c(
    1, 1, nrow(rast.base), 1, 1, ncol(rast.base),
    nrow(rast.base), ncol(rast.base)),
    byrow = T, ncol = 2)
  four.corners <- cbind(
    four.corners, xFromCol(rast.base, four.corners[, 2]))
  four.corners <- cbind(
    four.corners, yFromRow(rast.base, four.corners[, 1]))
  four.corners <- four.corners[, c(3, 4, 1, 2)]
  colnames(four.corners) <- c('x','y','col','row')
  coords <- rbind(coords, four.corners)
  
  
  if (distance.penalty) {
    distances <- nndist(coords[, c(1, 2)])
    distances.mean <- mean(distances)
  } else {
    distances.mean <- 0
  }
  
  # use vector indexing with multiple values
  # to speed up raster value indexing
  #
  # the conversion from row/col to 1-D index is
  #     index = (row - 1) * ncol + col
  #
  # selected.vals <- values(rast.anen)[
  #   ((coords[, 4] - 1) * ncol(rast.anen) + coords[, 3])]
  
  if (profile) {
    reformat.end <- Sys.time()
  }
  
  if (interpolation == 'idw') {
    require(spatstat)
    require(gstat)
    selected.vals <- values(rast.anen)[
      ((coords[, 4] - 1) * ncol(rast.anen) + coords[, 3])]
    
    df <- data.frame(x = coords[, 'x'],
                     y = coords[, 'y'],
                     z = selected.vals)
    mod <- gstat(id = 'z', formula = z~1, locations = ~x+y,
                 data = df, nmax = 7)
    rast.int <- interpolate(rast.base, mod, debug.level = 0)
    
  } else if (interpolation == 'nni') {
    require(RAnEnExtra)
    selected.vals <- values(rast.anen)[
      ((coords[, 4] - 1) * ncol(rast.anen) + coords[, 3])]
    
    rast.int <- RAnEnExtra::nni(
      coords[, 1], coords[, 2],
      selected.vals, rast.base,
      n = n_nei, searchtype = 'priority') 
    
  } else if (interpolation == 'voronoi') {
    require(rgeos)
    require(deldir)
    
    offset <- res(rast.base) * .1
    extent <- extent(rast.base)
    extent[1] <- extent[1] + offset[1]
    extent[2] <- extent[2] - offset[1]
    extent[3] <- extent[3] + offset[2]
    extent[4] <- extent[4] - offset[2]
    extent <- as(extent, "SpatialPolygons")
    
    polys <- voronoipolygons(coords)
    polys <- gIntersection(polys, extent, byid = T)
    
    rast.int <- meshing(rast.origin = rast.anen,
                        mesh.polys = polys,
                        return.type = 'raster')
    
  } else if (interpolation == 'delaunay') {
    require(spatstat)
    require(maptools)
    
    coords <- data.frame(coords)
    W <- owin(xrange = range(coords[, 1]),
              yrange = range(coords[, 2]))
    polys <- as(delaunay(as.ppp(coords, W=W)),
                "SpatialPolygons")
    rast.int <- meshing(rast.origin = rast.anen,
                        mesh.polys = polys,
                        return.type = 'raster')
  } else {
    stop(paste('Invalid interpolation method:',
               interpolation))
  }
  
  if (profile) {
    int.end <- Sys.time()
  }
  
  rmse <- sqrt(mean(values( rast.obs - rast.int )^2,
                    na.rm = T))
  
  if (profile) {
    compute.end <- Sys.time()
    
    reformat.time <- reformat.end - time.start
    int.time <- int.end - reformat.end
    compute.time <- compute.end - int.end
    function.time <- compute.end - time.start
    
    cat(paste("function:\t\t ( 100.00 % )", 
              format(function.time)), '\n')
    cat(paste("reformat:\t\t (",
              round(as.numeric(reformat.time) /
                      as.numeric(function.time) *
                      100, 2), "% )",
              format(reformat.time)), '\n')
    cat(paste("interpolation:\t\t (",
              round(as.numeric(int.time) /
                      as.numeric(function.time) *
                      100, 2), "% )",
              format(int.time)), '\n')
    cat(paste("error calculation:\t (",
              round(as.numeric(compute.time) /
                      as.numeric(function.time) *
                      100, 2), "% )",
              format(compute.time)), '\n')
  }
  
  return(distances.mean - rmse)
}

mutate.guido <- function (object, parent, ...) 
{
  # use tourname selection to determine which individual to mutate
  s <- sample(1:object@popSize, size = 3)
  parent <- which.max(object@fitness[s])
  mutate <- parent <- as.vector(object@population[parent, ])
  
  n <- length(parent)
  
  # Define dampening factor
  dampeningFactor <- 1.1 - object@iter/object@maxiter
  
  # define max step
  value.max <- (max(parent) - min(parent))  
  
  # each gene should have a possibility to be mutated
  gene.prob <- .5
  
  # genes that would get mutated
  selected <- runif(n) < gene.prob
  
  # mutate genes
  mutate[selected] <- rnorm(length(which(selected)), mean = parent[selected], sd = value.max/2/n*dampeningFactor)
  
  # make sure points are within boundary
  selected <- mutate < object@lower
  mutate[selected] <- object@lower[selected]
  selected <- mutate > object@upper
  mutate[selected] <- object@upper[selected]
  
  return(mutate)
}

ga.new <- function (type = c("binary", "real-valued", "permutation"), 
                    fitness, ..., lower, upper, nBits, population = gaControl(type)$population, 
                    selection = gaControl(type)$selection, crossover = gaControl(type)$crossover, 
                    mutation = gaControl(type)$mutation, popSize = 50, pcrossover = 0.8, 
                    pmutation = 0.1, elitism = base::max(1, round(popSize * 
                                                                    0.05)), updatePop = FALSE, postFitness = NULL, maxiter = 100, 
                    run = maxiter, maxFitness = Inf, names = NULL, suggestions = NULL, 
                    optim = FALSE, optimArgs = list(method = "L-BFGS-B", poptim = 0.05, 
                                                    pressel = 0.5, control = list(fnscale = -1, maxit = 100)), 
                    keepBest = FALSE, parallel = FALSE, monitor = if (interactive()) gaMonitor else FALSE, 
                    seed = NULL) 
{
  call <- match.call()
  type <- match.arg(type, choices = eval(formals(ga)$type))
  if (!is.function(population)) 
    population <- get(population)
  if (!is.function(selection)) 
    selection <- get(selection)
  if (!is.function(crossover)) 
    crossover <- get(crossover)
  if (!is.function(mutation)) 
    mutation <- get(mutation)
  if (missing(fitness)) {
    stop("A fitness function must be provided")
  }
  if (!is.function(fitness)) {
    stop("A fitness function must be provided")
  }
  if (popSize < 10) {
    warning("The population size is less than 10.")
  }
  if (maxiter < 1) {
    stop("The maximum number of iterations must be at least 1.")
  }
  if (elitism > popSize) {
    stop("The elitism cannot be larger that population size.")
  }
  if (pcrossover < 0 | pcrossover > 1) {
    stop("Probability of crossover must be between 0 and 1.")
  }
  if (is.numeric(pmutation)) {
    if (pmutation < 0 | pmutation > 1) {
      stop("If numeric probability of mutation must be between 0 and 1.")
    }
    else if (!is.function(population)) {
      stop("pmutation must be a numeric value in (0,1) or a function.")
    }
  }
  callArgs <- list(...)
  if (any("min" %in% names(callArgs))) {
    lower <- callArgs$min
    callArgs$min <- NULL
    warning("'min' arg is deprecated. Use 'lower' instead.")
  }
  if (any("max" %in% names(callArgs))) {
    upper <- callArgs$max
    callArgs$max <- NULL
    warning("'max' arg is deprecated. Use 'upper' instead.")
  }
  if (missing(lower) & missing(upper) & missing(nBits)) {
    stop("A lower and upper range of values (for 'real-valued' or 'permutation' GA) or nBits (for 'binary' GA) must be provided!")
  }
  switch(type, binary = {
    nBits <- as.vector(nBits)[1]
    lower <- upper <- NA
    nvars <- nBits
    if (is.null(names)) names <- paste0("x", 1:nvars)
  }, `real-valued` = {
    lnames <- names(lower)
    unames <- names(upper)
    lower <- as.vector(lower)
    upper <- as.vector(upper)
    nBits <- NA
    if (length(lower) != length(upper)) stop("lower and upper must be vector of the same length!")
    nvars <- length(upper)
    if (is.null(names) & !is.null(lnames)) names <- lnames
    if (is.null(names) & !is.null(unames)) names <- unames
    if (is.null(names)) names <- paste0("x", 1:nvars)
  }, permutation = {
    lower <- as.vector(lower)[1]
    upper <- as.vector(upper)[1]
    nBits <- NA
    nvars <- length(seq.int(lower, upper))
    if (is.null(names)) names <- paste0("x", 1:nvars)
  })
  if (is.null(suggestions)) {
    suggestions <- matrix(nrow = 0, ncol = nvars)
  }
  else {
    if (is.vector(suggestions)) {
      if (nvars > 1) 
        suggestions <- matrix(suggestions, nrow = 1)
      else suggestions <- matrix(suggestions, ncol = 1)
    }
    else {
      suggestions <- as.matrix(suggestions)
    }
    if (nvars != ncol(suggestions)) 
      stop("Provided suggestions (ncol) matrix do not match number of variables of the problem!")
  }
  if (is.logical(monitor)) {
    if (monitor) 
      monitor <- gaMonitor
  }
  if (is.null(monitor)) 
    monitor <- FALSE
  if (optim) {
    optimArgs.default <- eval(formals(ga)$optimArgs)
    optimArgs.default$control[names(optimArgs$control)] <- optimArgs$control
    optimArgs$control <- NULL
    optimArgs.default[names(optimArgs)] <- optimArgs
    optimArgs <- optimArgs.default
    rm(optimArgs.default)
    if (any(optimArgs$method == c("L-BFGS-B", "Brent"))) {
      optimArgs$lower <- lower
      optimArgs$upper <- upper
    }
    else {
      optimArgs$lower <- -Inf
      optimArgs$upper <- Inf
    }
    optimArgs$poptim <- min(max(0, optimArgs$poptim), 1)
    optimArgs$pressel <- min(max(0, optimArgs$pressel), 
                             1)
    optimArgs$control$maxit <- as.integer(optimArgs$control$maxit)
    if (is.null(optimArgs$control$fnscale)) 
      optimArgs$control$fnscale <- -1
    if (optimArgs$control$fnscale > 0) 
      optimArgs$control$fnscale <- -1 * optimArgs$control$fnscale
  }
  if (is.logical(parallel)) {
    if (parallel) {
      parallel <- startParallel(parallel)
      stopCluster <- TRUE
    }
    else {
      parallel <- stopCluster <- FALSE
    }
  }
  else {
    stopCluster <- if (inherits(parallel, "cluster")) 
      FALSE
    else TRUE
    parallel <- startParallel(parallel)
  }
  on.exit(if (parallel & stopCluster) stopParallel(attr(parallel, 
                                                        "cluster")))
  `%DO%` <- if (parallel && requireNamespace("doRNG", quietly = TRUE)) 
    doRNG::`%dorng%`
  else if (parallel) 
    `%dopar%`
  else `%do%`
  if (!is.null(seed)) 
    set.seed(seed)
  i. <- NULL
  fitnessSummary <- matrix(as.double(NA), nrow = maxiter, 
                           ncol = 6)
  colnames(fitnessSummary) <- names(gaSummary(rnorm(10)))
  bestSol <- if (keepBest) 
    vector(mode = "list", length = maxiter)
  else list()
  Fitness <- rep(NA, popSize)
  object <- new("ga", call = call, type = type, lower = lower, 
                upper = upper, nBits = nBits, names = if (is.null(names)) 
                  character()
                else names, popSize = popSize, iter = 0, run = 1, maxiter = maxiter, 
                suggestions = suggestions, population = matrix(), elitism = elitism, 
                pcrossover = pcrossover, pmutation = if (is.numeric(pmutation)) 
                  pmutation
                else NA, fitness = Fitness, summary = fitnessSummary, 
                bestSol = bestSol)
  if (maxiter == 0) 
    return(object)
  Pop <- matrix(as.double(NA), nrow = popSize, ncol = nvars)
  ng <- min(nrow(suggestions), popSize)
  if (ng > 0) {
    Pop[1:ng, ] <- suggestions
  }
  if (popSize > ng) {
    Pop[(ng + 1):popSize, ] <- population(object)[1:(popSize - 
                                                       ng), ]
  }
  object@population <- Pop
  for (iter in seq_len(maxiter)) {
    object@iter <- iter
    if (!parallel) {
      for (i in 1:length(Fitness)) 
        if (is.na(Fitness[i])) {
          fit <- do.call(fitness, c(list(Pop[i, ]), callArgs))
          if (updatePop) 
            Pop[i, ] <- attributes(fit)[[1]]
          Fitness[i] <- fit
        }
    }
    else {
      Fitness <- foreach(i. = 1:length(Fitness), .combine = "c") %DO% 
      {
        if (is.na(Fitness[i.])) 
          do.call(fitness, c(list(Pop[i., ]), callArgs))
        else Fitness[i.]
      }
    }
    object@population <- Pop
    object@fitness <- Fitness
    fitnessSummary[iter, ] <- gaSummary(object@fitness)
    object@summary <- fitnessSummary
    if (is.function(monitor)) {
      monitor(object)
    }
    if (optim & (type == "real-valued")) 
      if (optimArgs$poptim > runif(1)) {
        i <- sample(1:popSize, size = 1, prob = optimProbsel(Fitness, 
                                                             q = optimArgs$pressel))
        opt <- try(suppressWarnings(do.call(stats::optim, 
                                            c(list(fn = fitness, par = Pop[i, ], method = optimArgs$method, 
                                                   lower = optimArgs$lower, upper = optimArgs$upper, 
                                                   control = optimArgs$control), callArgs))), 
                   silent = TRUE)
        if (is.function(monitor)) {
          if (!inherits(opt, "try-error")) 
            cat("\b\b | Local search =", format(opt$value, 
                                                digits = getOption("digits")))
          else cat(" |", opt[1])
        }
        if (!inherits(opt, "try-error")) {
          Pop[i, ] <- opt$par
          Fitness[i] <- opt$value
        }
        object@population <- Pop
        object@fitness <- Fitness
        fitnessSummary[iter, ] <- gaSummary(object@fitness)
        object@summary <- fitnessSummary
      }
    if (keepBest) {
      object@bestSol[[iter]] <- unique(Pop[Fitness == 
                                             max(Fitness, na.rm = TRUE), , drop = FALSE])
    }
    if (is.function(postFitness)) {
      object <- do.call(postFitness, c(object, callArgs))
      Fitness <- object@fitness
      Pop <- object@population
    }
    if (iter > 1) 
      object@run <- garun(fitnessSummary[seq(iter), 1])
    if (object@run >= run) 
      break
    if (max(Fitness, na.rm = TRUE) >= maxFitness) 
      break
    if (object@iter == maxiter) 
      break
    ord <- order(Fitness, decreasing = TRUE)
    PopSorted <- Pop[ord, , drop = FALSE]
    FitnessSorted <- Fitness[ord]
    
    if (is.function(selection)) {
      sel <- selection(object)
      Pop <- sel$population
      Fitness <- sel$fitness
    }
    else {
      sel <- sample(1:popSize, size = popSize, replace = TRUE)
      Pop <- object@population[sel, ]
      Fitness <- object@fitness[sel]
    }
    object@population <- Pop
    object@fitness <- Fitness
    
    # save the selected population
    save.population.stage <- getOption("save.population.stage")
    
    if (is.logical(save.population.stage)) {
      require(stringr)
      if (save.population.stage) {
        if (identical(selection, random.weiming)) {
          file = paste('random-data-iteration_',
                       str_pad(object@iter, 3, pad = '0'), '.RData', sep = '')
        } else if (identical(selection, tour.weiming)) {
          file = paste('GA-data-iteration_',
                       str_pad(object@iter, 3, pad = '0'), '.RData', sep = '')
        } else {
          file = paste('unknown-selection-data-iteration_',
                       str_pad(object@iter, 3, pad = '0'), '.RData', sep = '')
        }
        save(object, file = file)
        if (identical(selection, random.weiming)) {
          object@fitness <- rep(0, length(object@fitness))
        }
      }
    }
    
    # Crossover
    if (is.function(crossover) & pcrossover > 0) {
      nmating <- floor(popSize/2)
      mating <- matrix(sample(1:(2 * nmating),
                              size = (2 * nmating)), ncol = 2)
      for (i in seq_len(nmating)) {
        
        # Dampning pcrossover
        if ((pcrossover*(1- object@iter/object@maxiter)) > runif(1)) {
          parents <- mating[i, ]
          Crossover <- crossover(object, parents)
          
          # Pop[parents, ] <- Crossover$children
          # Fitness[parents] <- Crossover$fitness
          
          Pop <- rbind(Pop, Crossover$children)
          Fitness <- c(Fitness, NA)
        }
      }
      object@population <- Pop
      object@fitness <- Fitness
    }
    pm <- if (is.function(pmutation)) 
      pmutation(object)
    else pmutation
    if (is.function(mutation) & pm > 0) {
      for (i in seq_len(popSize)) {
        if (pm > runif(1)) {
          Mutation <- mutation(object, i)
          
          Pop <- rbind(Pop, Mutation)
          Fitness <- c(Fitness, NA)
          
          # Pop[i, ] <- Mutation
          # Fitness[i] <- NA
        }
      }
      object@population <- Pop
      object@fitness <- Fitness
    }
    if (elitism > 0) {
      ord <- order(object@fitness, na.last = TRUE)
      u <- which(!duplicated(PopSorted, margin = 1))
      Pop[ord[1:elitism], ] <- PopSorted[u[1:elitism], 
                                         ]
      Fitness[ord[1:elitism]] <- FitnessSorted[u[1:elitism]]
      object@population <- Pop
      object@fitness <- Fitness
    }
    # if (is.function(monitor)) {
    #   cat("\n")
    #   flush.console()
    # }
  }
  if (optim & (type == "real-valued")) {
    optimArgs$control$maxit <- rev(optimArgs$control$maxit)[1]
    i <- which.max(object@fitness)
    opt <- try(suppressWarnings(do.call(stats::optim, c(list(fn = fitness, 
                                                             par = object@population[i, ], method = optimArgs$method, 
                                                             lower = optimArgs$lower, upper = optimArgs$upper, 
                                                             control = optimArgs$control), callArgs))), silent = TRUE)
    if (is.function(monitor)) {
      if (!inherits(opt, "try-error")) 
        cat(" | Final local search =", format(opt$value, 
                                              digits = getOption("digits")))
      else cat(" |", opt[1])
    }
    if (!inherits(opt, "try-error")) {
      object@population[i, ] <- opt$par
      object@fitness[i] <- opt$value
    }
  }
  if (is.function(monitor)) {
    cat("\n")
    flush.console()
  }
  object@summary <- na.exclude(object@summary)
  attr(object@summary, "na.action") <- NULL
  object@fitnessValue <- max(object@fitness, na.rm = TRUE)
  valueAt <- which(object@fitness == object@fitnessValue)
  solution <- object@population[valueAt, , drop = FALSE]
  if (nrow(solution) > 1) {
    eps <- gaControl("eps")
    solution <- unique(round(solution/eps) * eps, margin = 1)
  }
  colnames(solution) <- parNames(object)
  object@solution <- solution
  if (keepBest) 
    object@bestSol <- object@bestSol[!sapply(object@bestSol, 
                                             is.null)]
  return(object)
}

setGeneric(name = "parNames", 
           def = function(object, ...) { standardGeneric("parNames") }
)

setMethod("parNames", "ga",
          function(object, ...)
          { 
            names <- object@names
            nvars <- ncol(object@population)
            if(length(names) == 0)
            { names <- paste("x", 1:nvars, sep = "") }
            return(names)
          })

tour.weiming <- function (object, k = 3, ...) 
{
  sel <- rep(NA, object@popSize)
  for (i in 1:object@popSize) {
    s <- sample(1:object@popSize, size = k)
    sel[i] <- s[which.max(object@fitness[s])]
  }
  out <- list(population = object@population[sel, , drop = FALSE], 
              fitness = object@fitness[sel])
  return(out)
}

random.weiming <- function(object) {
  sel <- sample(1:length(object@fitness), size = object@popSize, replace = TRUE)
  out <- list(population = object@population[sel, , drop = FALSE], 
              fitness = object@fitness[sel])
  return(out)
}


#######################################
#          Parameter Setting          #
#######################################
n_nei <- 3
parallel <- T
num.iter <- 100
pop.size <- 100
pmutation <- .5
pcrossover <- 0
num.selected <- 200
interpolation <- 'voronoi'

func_mutation <- mutate.guido
func_selection <- tour.weiming
