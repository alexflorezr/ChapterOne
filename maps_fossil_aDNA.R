install.packages("rworldmap")
library(rworldmap)
points_colors <- c("#8B4513","#8B5A00","#D2691E","#EE7600","#FFA500","#CDCD00","#9ACD32","#43CD80","#008B45", "#006400")
points_time  <- seq(50000, 5000, by=-5000)
points_table <- cbind(points_colors, points_time)
for (i in seq_along(species_vector)){
  temp_points <- read.delim(paste(species_vector[i], "_db_seq_clean.txt", sep=""), header=T, stringsAsFactors=F)
  newmap <- getMap(resolution = "low")
  par(mar=c(1,1,5,1))
  plot(newmap, xlim = c(-180, 180), ylim = c(0,90), asp=1, main="Cualquier")
  for (k in seq_along(temp_points[,1])){
    for (j in seq_along(points_table[,2])){
      if (temp_points$median_OxCal[k] <= as.numeric(points_table[j,2]) && temp_points$median_OxCal[k] >= as.numeric(points_table[j,2])-5000){
        temp_points$Map_color[k] <- points_table[j,1]
      }
    }
  }
  points(temp_points$Longitude, temp_points$Latitude, col=temp_points$Map_color, cex=1.5, pch=21, bg=paste(temp_points$Map_color, 90, sep=""))
}