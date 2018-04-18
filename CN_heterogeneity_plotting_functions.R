
#These small functions should bve quick to rewrite:

fun.chrom.length<- function(seg){
  # This function returns the length of chromosome 1,...,22
  # seg: first three columns of the minimum conistent region matrix
  print("Called fun.chrom.length")
  chrom.length <- c()
  for(i in 1:22){
    sub.chrom <- subset(seg, seg[, 1] == i)
    chrom.length <- c(chrom.length, sub.chrom[nrow(sub.chrom), 3])
    
  }
  print("Returning from fun.chrom.length")
  return(chrom.length)
  
}


fun.add.chrom<-function(seg, chrom.length){
  print("Called fun.add.chrom")
  # This function adds the length of chromosome 1,...,(i - 1) to chromsomes 2,...,22
  # seg: first three columns of the minimum conistent region matrix
  # chrom.length: output of fun.chrom.length(...)
  if(1 %in% seg[, 1])
    seg.tmp <-  subset(seg, seg[, 1] == 1)
  else
    seg.tmp <- matrix(NA, nr = 0, nc = ncol(seg))  
  for(i in 2:22){
    if(!(i %in% seg[, 1]))
      next()
    sub.chrom <- subset(seg, seg[, 1] == i)
    sub.chrom[, 2] <- sub.chrom[, 2] + sum(chrom.length[1:(i - 1)])
    sub.chrom[, 3] <- sub.chrom[, 3] + sum(chrom.length[1:(i - 1)])
    seg.tmp <- rbind(seg.tmp, sub.chrom)
  }
  print("Returning from fun.add.chrom")
  return(seg.tmp)
  
}

fun.chrom.mean<-function(seg){
  print("Called fun.chrom.mean")
  # This function finds the middle of the positions for each chromosome.
  # seg: Matrix with three columns, containing the positions.
  chrom.mean <- c()
  for(i in 1:22){
    sub.chrom <- subset(seg, seg[, 1] == i)
    chrom.mean <- c(chrom.mean, (min(sub.chrom[, 2]) + max(sub.chrom[, 3]))/2)
    
  }
  print("Returning from fun.chrom.mean")
  return(chrom.mean)
}

