# plot different types of meshing grid 
# This is the main file for multiple files
# starting with GA_***.R
#
# develop a genetic algorithm for optimized grid adaptation
# and Analog Ensemble technique
#
library(sp)
library(GA)
library(maps)
library(RAnEn)
library(gstat)
library(ncdf4)
library(rgdal)
library(fields)
library(raster)
library(deldir)
library(spatstat)
library(maptools)
library(RColorBrewer)

# visualize data
visualize <- T

# output figures
output.figure <- F

# develop unstructured grid on test day index and flt
analog.day.index <- 1
flt.index <- 1

# resample brewer colors
brewer.cols <- colorRampPalette(
  brewer.pal(11, 'Spectral')[11:1])(100)


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


# interpolate AnEn on GA grid
# ditto for the GA solution
ga.pts <- matrix(GA@solution, ncol=2, byrow=T)
ga.pts <- cbind(ga.pts, colFromX(
  rast.base, ga.pts[, 1]))
ga.pts <- cbind(ga.pts, rowFromY(
  rast.base, ga.pts[, 2]))
colnames(ga.pts) <- c('x', 'y', 'col', 'row')

spd.voronoi <- interpolate.mesh(
  coords = ga.pts, rast.from = rast.anen, method = 'voronoi',
  return.type = 'SpatialPolygonsDataFrame')


# define the points on the four edges
four.corners <- matrix(c(
  1, 1, nrow(rast.base), 1, 1, ncol(rast.base),
  nrow(rast.base), ncol(rast.base)),
  byrow = T, ncol = 2)
four.corners <- cbind(
  four.corners, xFromCol(rast.base, four.corners[, 2]))
four.corners <- cbind(
  four.corners, yFromRow(rast.base, four.corners[, 1]))
four.corners <- four.corners[, c(3, 4, 1, 2)]
ga.pts.w.corners <- rbind(ga.pts, four.corners)

spd.delaunay <- interpolate.mesh(
  coords = ga.pts.w.corners, rast.from = rast.anen, method = 'delaunay',
  return.type = 'SpatialPolygonsDataFrame')

if (visualize) {
  
  if (output.figure) {
    jpeg('different-DOUG.jpeg',
         width = 10, height = 5,
         res = 300, units = 'in')
  }
  
  cex.sub <- 1.5
  
  par(mfrow = c(1, 3), mar = c(3, 1, 1, 1) + .1,
      oma = c(6, 2, 0, 1))
  zlim <- range(c(values(rast.obs),
                  spd.voronoi@data$z,
                  spd.delaunay@data$z),
                na.rm = T)
  plot(rast.obs, col = brewer.cols, zlim = zlim, legend = F)
  title(sub = '(a)', line = 2, cex.sub = cex.sub)
  
  plot(rast.obs, col = 'white', zlim = zlim, legend = F)
  plot(spd.voronoi, add = T, col = map.to.colors(
    brewer.cols, spd.voronoi@data$z))
  title(sub = '(b)', line = 2, cex.sub = cex.sub)
  
  plot(rast.obs, col = 'white', zlim = zlim,
       legend = F)
  plot(spd.delaunay, add = T, col = map.to.colors(
    brewer.cols, spd.delaunay@data$z))
  title(sub = '(c)', line = 2, cex.sub = cex.sub)
  
  par(mfrow=c(1,1), new=FALSE, oma=c(0,0,0,0))
  plot(rast.obs, legend.only = T, legend.width = .5,
       legend.shrink = .9,
       horizontal = T, zlim = zlim, col = brewer.cols)
  mtext('°C', side = 1, adj = 1, at = c(-72, 38))
  
  if (output.figure) dev.off()
}
