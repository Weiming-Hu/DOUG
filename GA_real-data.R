# use the real data
# and compute AnEn
#
# sourced by GA_develop-unstructured-grid.R
#
# read the actual NAM data
require(ncdf4)
require(RAnEn)

print("Reading data from NetCDF files ...")

forecasts.file <- 'forecasts.nc'
observations.file <- 'analysis.nc'

forecasts.nc <- nc_open(forecasts.file)
observations.nc <- nc_open(observations.file)

xs <- ncvar_get(forecasts.nc, 'xs')
ys <- ncvar_get(forecasts.nc, 'ys')

var.names <- ncvar_get(forecasts.nc, 'Names')

nx <- ncatt_get(forecasts.nc, 'Data', 'nx')$value
ny <- ncatt_get(forecasts.nc, 'Data', 'ny')$value

forecasts <- ncvar_get(
  forecasts.nc, 'Data')

forecasts[1, , , ] <- forecasts[1, , , ] - 273.15

observations <- ncvar_get(
  observations.nc, 'Data', collapse_degen = F) - 273.15

nc_close(forecasts.nc)
nc_close(observations.nc)

test.ID.start <- 337 # this is 01/01/2016
test.ID.end <- 400

members.size <- 10
circulars <- 4
rolling <- -2
verbose <- 2
quick <- T
cores <- 8

search.extension <- F
output.search.stations <- F
output.metric <- F
extend.obs <- F
num.nearest <- 8

# this is the variable in the forecasts that would be
# used to compare with AnEn and the analysis
#
forecast.var <- 1

if ("precomputed.AnEn" %in% ls()) {
  if (!file.exists(precomputed.AnEn)) {
    command <- paste(
      '/Users/wuh20/github/CAnalogsV2/install/bin/canalogs -N',
      '--forecast-nc', forecasts.file,
      '--observation-nc', observations.file,
      '--test-ID-start', test.ID.start-1,
      '--test-ID-end', test.ID.end-1,
      '--members-size', members.size,
      '--rolling', rolling,
      '--cores', cores,
      '-o', precomputed.AnEn
    )
    if (quick) command <- paste(command, '--quick')
    system(command)
  }
  
  print(paste('read from precomputed AnEn NetCDF',
              precomputed.AnEn))
  AnEn.nc <- nc_open(precomputed.AnEn)
  analogs <- ncvar_get(AnEn.nc, 'Data', collapse_degen = F) - 273.15
  stations <- ncvar_get(AnEn.nc, 'MembersStation', collapse_degen = F)
  days <- ncvar_get(AnEn.nc, 'MembersTime', collapse_degen = F)
  nc_close(AnEn.nc)
  
  combine <- array(dim = c(dim(analogs), 3))
  combine[, , , , 1] <- analogs
  combine[, , , , 2] <- stations + 1
  combine[, , , , 3] <- days + 1
  analogs <- combine
  rm(stations, days, combine)
  
} else {
  AnEn <- compute_analogs(
    forecasts = forecasts, observations = observations,
    # stations_ID = stations.ID,
    # observation_ID = observation.ID,
    test_ID_start = test.ID.start, test_ID_end = test.ID.end,
    # train_ID_start = train.ID.start, train_ID_end = train.ID.end,
    members_size = members.size,
    circulars = circulars,
    # weights = weights, 
    rolling = rolling,
    output_search_stations = output.search.stations,
    output_metric = output.metric,
    search_extension = search.extension,
    quick = quick, cores = cores, verbose = verbose,
    #distance = distance,
    num_nearest = num.nearest,
    xs = as.numeric(xs), ys = as.numeric(ys),
    observation_from_extended_station = extend.obs)
  
  analogs <- AnEn$analogs
}

bias.corrected.anen <- 'bias_corrected_anen.RData'
if (file.exists(bias.corrected.anen)) {
  print('Read bias correced AnEn ...')
  load(bias.corrected.anen)
} else {
  print("Bias correction")
  analogs <- bias_correction(
    analogs = analogs, observations = observations, forecasts = forecasts,
    forecast.ID = 1, test.ID.start = test.ID.start, group.func = mean,
    na.rm = T, keep.bias = F)$bias.corrected.analogs
  save(analogs, file = bias.corrected.anen)
}