##Analysis for Sheep Golgi. 
##Remember to change the working directory as appropriate

##read backbone table for formatting
table_backbone <- read.table(file = "C:/Users/Sebastian/Desktop/Golgi_test/Cortical Plate/Sheep/GA 90/19222/table_backbone.csv", header = F, sep = ",", colClasses = "character", as.is = "!stringsAsFactors")
colnames(table_backbone)
table_backbone[1,]
colnames(table_backbone) <- table_backbone[1,]

##Start from here
##Set working directory
setwd("C:/Users/Sebastian/Desktop/Golgi_test/Cortical Plate/Sheep/GA 90/19226/Sholl output/")

#Create file lists for both gyri and sulci tables
file.names.gyri <- dir(path = "../Sholl output/", pattern = "*-G.traces*")
file.names.sulci <- dir(path = "../Sholl output/", pattern = "*-S.traces*")

#Convert table backbone into final table
final_table <- table_backbone[2:22,]
final_table[,2] <- final_table[,3]
final_table[,1] <- seq(0,300, by = 15)
final_table <- final_table[,c(1,2)]
final_table

library("dplyr")
##join gyri files into final table
for (x in 1:length(file.names.gyri)) {
  files <- read.table(file = file.names.gyri[x], header = T, sep = ",")[,c(1,2)]
  final_table <- left_join(final_table, files, by = "Radius") 
}

##join sulci files into final table
for (x in 1:length(file.names.sulci)) {
  files <- read.table(file = file.names.sulci[x], header = T, sep = ",")[,c(1,2)]
  final_table <- left_join(final_table, files, by = "Radius") 
}


#Explore final table
final_table
dim(final_table)
str(final_table)

#Format final table for output
final_table <- final_table[,c(1,3:42)]
colnames(final_table) <- c("Radius",colnames(table_backbone)[2:41])
final_table[is.na(final_table)] <- 0

#Write out the final table
write.csv(final_table, file = "final_table.csv")
