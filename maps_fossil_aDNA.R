##############################
### Read the full database ###
##############################

Full_DB <- read.delim(file.choose(), header=T, stringsAsFactors=F)
Full_DB$Longitude <- as.numeric(Full_DB$Longitude)
Full_DB$Latitude <- as.numeric(Full_DB$Latitude)
full_vector <- as.vector(NULL)
for (i in seq_along(Full_DB$Latitude)){
  if (is.na(Full_DB$Longitude[i]) || is.na(Full_DB$Latitude[i])){
    full_vector <- c(full_vector,i)
  }
}
# remove the records without longitude and or latitude
Full_DB_LL <- Full_DB[-full_vector,]
# add two empty colums to assign color and type of point in the map
For_map <- matrix(NA,nrow=length(Full_DB_LL$Latitude), ncol=2)
colnames(For_map) <- c("Map_color", "Map_type")
Full_DB_map <- cbind(Full_DB_LL, For_map)

##############################
### maps: species together ###
##############################

library(rworldmap)
points_colors <- c("#8B4513","#8B5A00","#D2691E","#EE7600","#FFA500","#CDCD00","#9ACD32","#43CD80","#008B45", "#006400")
points_time  <- seq(50000, 5000, by=-5000)
points_table <- cbind(points_colors, points_time)
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-180, 180), ylim = c(0,90), asp=1, main=paste("all the species, n = ", length(Full_DB_LL$Latitude), sep=""))
for (k in seq_along(Full_DB_map[,1])){
  for (j in seq_along(points_table[,2])){
    if (Full_DB_map$Median_Age[k] <= as.numeric(points_table[j,2]) && Full_DB_map$Median_Age[k] >= as.numeric(points_table[j,2])-5000){
        Full_DB_map$Map_color[k] <- points_table[j,1]
    }
    if (length(strsplit(Full_DB_map$Sequence[k], split="")[[1]]) > 1){
        Full_DB_map$Map_type[k] <- 16
    }
    if (length(strsplit(Full_DB_map$Sequence[k], split="")[[1]]) <= 1){
        Full_DB_map$Map_type[k] <- 21
    }
  }
}
points(Full_DB_map$Longitude, Full_DB_map$Latitude, col=Full_DB_map$Map_color, cex=1.5, pch=Full_DB_map$Map_type, bg=paste(Full_DB_map$Map_color, 90, sep=""))


############Delete above this line #####
