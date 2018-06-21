

##
##Params:
##
##Returns:
##
readInMEDICCInputFormatParallel <- function(snp.file.df, processed.orig.cpn.segs.df, num.cores = 8) {
  # TODO add check that every single all.samp.SNP.list has the same SNPs..
  # Note error checking of removing NAs from the min.cons.seg.df should already have been performed
  registerDoParallel(cores = num.cores)
  print(foreach::getDoParWorkers())
  print(getDoParName())
  print(getDoParVersion())

  seg.name.vec <- paste(processed.orig.cpn.segs.df[, 1], processed.orig.cpn.segs.df[, 2], processed.orig.cpn.segs.df[, 3], sep=":")
  all.samp.vec <- processed.orig.cpn.segs.df$sample
  sample.vec <- unique(all.samp.vec)

  ASCAT.summary.df.list <- list()
  for (samp in sample.vec) {

  	samp.ASCAT.df <- data.frame(chr = as.numeric(snp.file.df$Chr),
                  				pos = as.numeric(snp.file.df$Pos),
     							BAF = as.numeric(snp.file.df[, paste0(gsub("-" , "\\.", samp), "_BAF")]),
     							BAFsegmented = rep(NA, length(snp.file.df$Chr)),
     							logR = as.numeric(snp.file.df[, paste0(gsub("-", "\\.", samp), "_LogR")]),
     							logRsegmented = rep(NA, length(snp.file.df$Chr)),
     							stringsAsFactors = FALSE)
  	samp.ASCAT.df$chr_pos <- paste0(samp.ASCAT.df$chr, "_", samp.ASCAT.df$pos)
  	ASCAT.summary.df.list[[samp]] <- samp.ASCAT.df

  }#end ASCAT.summary.df.list

  iterable.orig.cpn.segs.df <- cbind(sample = all.samp.vec, 
  								seg_name = processed.orig.cpn.segs.df$seg.name.vec,
  								chr = processed.orig.cpn.segs.df$chr,
  								startpos = processed.orig.cpn.segs.df$startpos,
  								endpos =  processed.orig.cpn.segs.df$endpos,
  								nMinor = processed.orig.cpn.segs.df$nMinor,
  								nMajor = processed.orig.cpn.segs.df$nMajor)

  itx <- iter(iterable.min.cons.seg.df, by="row")

	#results.list <- foreach(iterable.orig.cpn.segs.df = itx, .combine = c) %dopar% {
		for(i in 1:nrow(iterable.orig.cpn.segs.df))

  		seg <- processed.orig.cpn.segs.df[processed.orig.cpn.segs.df$sample == iterable.orig.cpn.segs.df$sample,]

  		ASCAT.summary <- ASCAT.summary.df.list[[iterable.orig.cpn.segs.df$sample]]
  		#So now we need to work out how many rows there will be in df corresponding to the sample segment that we return should be:
  		seg.ASCAT.summary <- ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
  					  						& ASCAT.summary$pos >= seg[j ,]$startpos 
  					  						& ASCAT.summary$pos <= seg[j ,]$endpos , , drop = FALSE]

  		seg.ASCAT.df <- data.frame(chr = as.numeric(seg.ASCAT.summary$Chr),
  		            				pos = as.numeric(seg.ASCAT.summary$Pos),
  		     						BAF = as.numeric(seg.ASCAT.summary[, paste0(gsub("-" , "\\.", samp), "_BAF")]),
  		     						BAFsegmented = rep(NA, length(seg.ASCAT.summary$Chr)),
  		     						logR = as.numeric(seg.ASCAT.summary[, paste0(gsub("-", "\\.", samp), "_LogR")]),
  		     						logRsegmented = rep(NA, length(seg.ASCAT.summary$Chr)),
  		     						sample = rep(processed.orig.cpn.segs.df$sample, length(seg.ASCAT.summary$Chr)), 
  		     						stringsAsFactors = FALSE)

  		if (nrow(seg.ASCAT.df) > 0) {

  			if (nrow(seg.ASCAT.df[seg.ASCAT.df$chr ==  iterable.orig.cpn.segs.df$chr
  									& seg.ASCAT.df$pos >=  iterable.orig.cpn.segs.df$startpos 
  									& seg.ASCAT.df$pos <=  iterable.orig.cpn.segs.df$endpos 
  									& seg.ASCAT.df$BAF != 2 ,]) > 0) {
  				#Is this segment balanced, if so set the BAFsegmented to 0.5
  			  if (iterable.orig.cpn.segs.df$nMajor == iterable.orig.cpn.segs.df$nMinor) {

  			    seg.ASCAT.df[seg.ASCAT.df$chr == iterable.orig.cpn.segs.df$chr
  			    			  & seg.ASCAT.df$pos >= iterable.orig.cpn.segs.df$startpos
  			    			  & seg.ASCAT.df$pos <= iterable.orig.cpn.segs.df$endpos,]$BAFsegmented <- 0.5
  			    
  			  } else {#Segment is balanced so we need to calculate the value of the BAF segmented and apply it to all these values 
  			

  			    avg.seg.baf <- 0.5 - mean(abs(as.numeric(seg.ASCAT.df[seg.ASCAT.df$chr == iterable.orig.cpn.segs.df$chr
  			    	 													& seg.ASCAT.df$pos >= iterable.orig.cpn.segs.df$startpos
  			    	 													& seg.ASCAT.df$pos <= iterable.orig.cpn.segs.df$endpos
  			    	 													& seg.ASCAT.df$BAF != 2 ,]$BAF)          - 0.5))

  			    seg.ASCAT.df[seg.ASCAT.df$chr == iterable.orig.cpn.segs.df$chr 
  			    			  & seg.ASCAT.df$pos >= iterable.orig.cpn.segs.df$startpos 
  			    			  & seg.ASCAT.df$pos <= iterable.orig.cpn.segs.df$endpos ,]$BAFsegmented <- avg.seg.baf

  			  }
  			
  			}#end if(nrow(ASCAT.summary[ASCAT.summary$chrom==seg[j,]$chr & ASCAT.summary$pos...

  		}#if (nrow(seg.ASCAT.df) > 0) {

  		seg.ASCAT.df
  	}#end results <- foreach(iterable.orig.cpn.segs.df = itx, .combine = c) %dopar% {

	all.samp.seg.ASCAT.summary.df <- do.call(rbind, results)

	#Now let's remake the ASCAT.summary.df.list 
	final.ASCAT.summary.df.list <- list()
	for (samp in sample.vec) {

		orig.samp.ASCAT.summary.df <- ASCAT.summary.df.list[[samp]]
		orig.samp.ASCAT.summary.df$chr_pos <- paste0(orig.samp.ASCAT.summary.df$chr, "_", orig.samp.ASCAT.summary.df$pos)
		unsorted.new.samp.summary.df <- ASCAT.all.samp.seg.ASCAT.summary.df[all.samp.seg.ASCAT.summary.df$sample == samp,]
		missing.samp.ASCAT.summary.df <- orig.samp.ASCAT.summary.df[ ! orig.samp.ASCAT.summary.df$chr_pos %in% unsorted.new.samp.summary.df$chr_pos]
		all.samp.positions.ASCAT.summary.df <- rbind(unsorted.new.samp.summary.df, missing.samp.ASCAT.summary.df)
		sorted.all.samp.positions.ASCAT.summary.df <- all.samp.positions.ASCAT.summary.df[with(all.samp.positions.ASCAT.summary.df, order(chr, pos)), ]
		final.ASCAT.summary.df.list[[samp]] <- sorted.all.samp.positions.ASCAT.summary.df

	}

  return(final.ASCAT.summary.df.list)

}#end cleanAndErrorCheckSNPsParallel <- function(all.samp.SNP.list, min.cons.seg.df,



