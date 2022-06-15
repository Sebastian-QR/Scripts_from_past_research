##read final table backup for formatting
final_table_backup <- tail(read.csv(file = "C:/Users/Sebastian/Desktop/MAP2 co-localisations - (ROI)/CDK5 MAP2/110d/Analysis output/final_table_backup.csv", header = F, sep = ","),28)
rownames(final_table_backup) <- final_table_backup[,1]

##Set working directory
setwd("C:/Users/Sebastian/Desktop/MAP2 co-localisations - (ROI)/CDK5 MAP2/110d/Analysis output/2144.2/")
path <- dir(path = "C:/Users/Sebastian/Desktop/MAP2 co-localisations - (ROI)/CDK5 MAP2/110d/Analysis output/2144.2/", pattern = "*.txt")
##create blank output table
final_table <- final_table_backup
##Create file list
file.names <- dir(path = "C:/Users/Sebastian/Desktop/MAP2 co-localisations - (ROI)/CDK5 MAP2/110d/Analysis output/2144.2/", pattern = "*.txt")

##cbind all files into final table
for (x in 1:length(file.names)) {
  files <- tail(read.csv(file.names[x], header = F, sep = ","),28)
  final_table <- cbind(final_table, files[,2])
}
tail(final_table)
dim(final_table)
rownames(final_table) <- rownames(final_table_backup)
colnames(final_table) <- c("NILL","NILL",as.character(file.names))
write.csv(final_table, file="2144.2output.csv")
