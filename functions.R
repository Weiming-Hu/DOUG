voronoipolygons <- function(crds) {
  # this function builds Voronoi polygons and
  # return SpatialPolygons
  #
  # referece 
  # http://carsonfarmer.com/2009/09/voronoi-polygons-with-r/
  #
  require(deldir)
  require(sp)
  
  z = deldir(crds[, 'x'], crds[, 'y'])
  w = tile.list(z)
  polys = vector(mode='list', length=length(w))
  
  for (i in seq(along=polys)) {
    pcrds = cbind(w[[i]]$x, w[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1,])
    polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(i))
  }
  
  SP = SpatialPolygons(polys)
  return (SP)
}

meshing <- function(rast.origin, mesh.polys,
                    return.type = 'SpatialPolygonsDataFrame',
                    use.velox = T) {
  # this function generates raster with the mesh specified in
  # mesh.polys. Raster values inside the same mesh will be reassigned
  # a new value that is the mean of the vertex value
  #
  z <- as.numeric(lapply(mesh.polys@polygons, function (poly) {
    vertice <- poly@Polygons[[1]]@coords[-1, ]
    return(mean(extract(rast.origin, vertice),
                na.rm = T))
  }))
  
  if (return.type == 'SpatialPolygonsDataFrame') {
    mesh.spd <- SpatialPolygonsDataFrame(
      mesh.polys, data = data.frame(z = z),
      match.ID = F)
    return(mesh.spd)
  } else if (return.type == 'raster') {
    if (use.velox) {
      require(velox)
      velox.rast <- velox(rast.origin)
      mesh.spd <- SpatialPolygonsDataFrame(
        mesh.polys, data = data.frame(z = z),
        match.ID = F)
      velox.rast$rasterize(
        spdf = mesh.spd, field = 'z', background = NA)
      rast.mesh <- velox.rast$as.RasterLayer()
    } else {
      rast.mesh <- rasterize(mesh.polys, rast.origin)
      values(rast.mesh) <- z[values(rast.mesh)] 
    }
    return(rast.mesh)
  }
}

map.to.colors <- function(cols, vals) {
  # map each vals to a color based on the cols from cols
  out.cols <- vector(length = length(vals))
  
  cuts <- cut(vals, length(cols))
  out.cols <- cols[as.numeric(cuts)]
  
  return(out.cols)
}

interpolate.mesh <- function(
  coords, rast.from, method,
  nmax = stop('Set the number of neighbors nmax'),
  use.velox = T, return.type = 'raster') {
  # interpolate a new rasters according to the coords and values
  # from rast.from. Method can be idw, voronoi, or delaunay.
  # coords is should have cols named 'x', 'y', 'row', and 'col'.
  #
  # nmax is the number of neighbors for 'idw'.
  # use.velox specifies whether to use faster raster in 'voronoi'.
  #
  if (method == 'idw') {
    require(gstat)
    require(raster)
    
    if (!all(c('x', 'y', 'row', 'col') %in% colnames(coords))) {
      stop('coords should have x, y, row, and col as column names')
    }
    
    selected.vals <- values(rast.from)[
      ((coords[, 'row'] - 1) * ncol(rast.from) + coords[, 'col'])]
    
    df <- data.frame(
      x = coords[, 'x'], y = coords[, 'y'], z = selected.vals)
    
    mod <- gstat(id = 'z', formula = z~1, locations = ~x+y,
                 data = df, nmax = nmax)
    rast.int <- raster::interpolate(rast.from, mod, debug.level = 0)
    
  } else if (method == 'voronoi') {
    require(rgeos)
    require(deldir)
    
    if (!all(c('x', 'y') %in% colnames(coords))) {
      stop('coords should have x and y as column names')
    }
    
    offset <- rep(res(rast.from), each = 2) * .1 * c(1, -1, 1, -1)
    extent <- as(extent(rast.from) + offset, 'SpatialPolygons')
    
    polys <- voronoipolygons(coords)
    polys <- gIntersection(polys, extent, byid = T)
    
    rast.int <- meshing(rast.origin = rast.from,
                        mesh.polys = polys,
                        return.type = return.type,
                        use.velox = use.velox)
    
  } else if (method == 'delaunay') {
    require(rgeos)
    require(deldir)
    require(spatstat)
    require(maptools)
    
    if (!all(c('x', 'y') %in% colnames(coords))) {
      stop('coords should have x and y as column names')
    }
    
    W <- owin(xrange = range(coords[, 'x']),
              yrange = range(coords[, 'y']))
    
    polys <- as(delaunay(as.ppp(coords[, c('x', 'y')], W = W)),
                'SpatialPolygons')
    
    rast.int <- meshing(rast.origin = rast.from,
                        mesh.polys = polys,
                        return.type = return.type,
                        use.velox = use.velox)
  } else {
    stop(paste('Invalid method:', method))
  }
  
  return(rast.int)
}