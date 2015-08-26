setwd("~/Desktop/PhD/Thesis/2ndChapter/Simulation_Data/")
# examples for latitude
# Uses the latitude of the midpoint of a cell, e.g. 55.5 for a cell with Denmark in it
cellsurfarea <- function(latitude, cellsize.lat = 1, cellsize.long = cellsize.lat) {
  degtorad <- function(deg) 
    return(deg*pi/180)
  surfarea <- function(R,lambda1, lambda2,phi1, phi2) 
    return(R^2 * (lambda2-lambda1) * (sin(phi2) - sin(phi1)))
  R <- 6378.1
  return(surfarea(R, degtorad(1 * cellsize.long), degtorad(0), degtorad(latitude + 0.5 * cellsize.lat), degtorad(latitude - 0.5 * cellsize.lat))) 
}


# An example giving the area of all 0.01x1 degree grid cell bands, and then summing up:
cellsize.lat <- 0.01
midpoints <- seq(0, 90, by = cellsize.lat)
midpoints <- midpoints[- 9001] + 0.5 * cellsize.lat
areas <- sapply(midpoints, cellsurfarea, cellsize.lat = cellsize.lat, cellsize.long = 1)
cellareas <- tapply(areas, rep(0:89 + 0.5, each = 100), sum)
cols<- 360
rows<- 180
area_matrix <- matrix(ncol=4, nrow=(cols*rows))
colnames(area_matrix) <- c("Cell", "Lat", "Lon", "Area")
area_matrix[,c(1,2,3)] <- cbind(seq(1,cols*rows, by=1), rep(seq((cols/2*-1)+.5,cols/2-.5, by=1), times=rows), rep(seq((rows/2*-1)+.5,(rows/2)-.5, by=1), each=cols))
for (i in seq_along(rownames(as.matrix(cellareas)))){
  area_matrix[which(abs(area_matrix[,3]) == as.numeric(rownames(as.matrix(cellareas))[i])),4] <- as.matrix(cellareas)[i]
}
write.table(x=area_matrix, file="surface_area.txt", sep="\t")
#create a matrix to map the area as an image
image<- matrix(ncol=cols, nrow=rows)
image[,1] <- unique(area_matrix[,3])
for (j in seq_along(unique(area_matrix[,3]))){
image[j,] <- area_matrix[which(area_matrix[,3] == unique(area_matrix[,3])[j]),4]
}                       
image(t(image))
write.table(x=image, file="image_area.txt", sep="\t")
library(raster)
raster_area <- raster(file.choose())
plot(raster_area)
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-180, 180), ylim = c(0,90), asp=1, add=T)

