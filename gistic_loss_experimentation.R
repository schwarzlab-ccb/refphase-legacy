#Fiffling with gistic losses t see if we can make them make sense


a <- read.table("LUSC_cl_loss_LOH_GISTIC_input.txt", sep = "\t", header = F , stringsAsFactors =F)
 orig <- a
saveRDS(a, file="orig_segmentation")

#so we don't want any gains:

a[as.numeric(a$V6) > 0,]$V6 <- as.numeric(a[as.numeric(a$V6) > 0,]$V6) - floor(as.numeric(a[as.numeric(a$V6) > 0,]$V6))

a[as.numeric(a$V6) > 0.4,]$V6 <- as.numeric(a[as.numeric(a$V6) > 0.4,]$V6) - 0.4


a[as.numeric(a$V6) > 0,]$V6

write.table(a, file = "LUSC_cl_loss_LOH_GISTIC_input.txt", quote = FALSE, row.names=FALSE, col.names =FALSE, sep ="\t")


a[as.numeric(a$V6) < -1,]$V6 <- a[as.numeric(a$V6) < -1.5,]$V6  + ((a[as.numeric(a$V6) < -1.5,]$V6 * -1) - 1.5)


a[as.numeric(a$V6) < -1.2,]$V6 <- -1.4


