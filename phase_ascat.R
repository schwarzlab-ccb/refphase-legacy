# Apply phasing to an existing ascat.bc object and recalculate all necessary
# values within the object
phase_ascat <- function(ascat_bc, phased_bafs, segs) {

		#So let's get the phased BAF:
		samp.pha.baf.df <- patient.phased.BAF.df.list[[samp]]
		#So we're going to need the indices of the BAF snps...
		samp.pha.baf.df$snp_pos <- phased.baf.info.df$snp_pos
		samp.pha.baf.df$snp_idx <- phased.baf.info.df$snp_idx
		#So now we need the segments that refer to this sample (as these are joint segments
		#these should all be same for each sample regardless.)
		samp.seg.df <- segs

		#Let's do the LogR
		tum.logr.vec <- ascat.bc[["Tumor_LogR"]][,1]
		if (any(is.na(tum.logr.vec))) {
			message("There are NA values in the ascat.bc Tumor_LogR field which need to be removed. Revisit generation of this input. For changing the NAs to 0.")
			na.idx.vec <- which(is.na(tum.logr.vec))
			tum.logr.vec[na.idx.vec] <- 0
		}

		tum.logr.segmented.vec <- ascat.bc[["Tumor_LogR_segmented"]][,1]
		if(any(is.na(tum.logr.segmented.vec))) {
			message("There are NA values in the ascat.bc Tumor_LogR_segmented field which need to be removed. Revisit generation of this input. For changing the NAs to 0.")
			na.segmented.idx.vec <- which(is.na(tum.logr.segmented.vec))
			tum.logr.segmented.vec[na.segmented.idx.vec] <- 0
		}

		for (i in 1:nrow(samp.seg.df)) {
      tum.logr.segmented.vec[samp.seg.df[i,]$start_snp_idx:samp.seg.df[i,]$end_snp_idx] <- mean(tum.logr.vec[samp.seg.df[i,]$start_snp_idx:samp.seg.df[i,]$end_snp_idx])
		}

		# Let's do the BAF recalculation
		tum.baf.segmented.vec <- ascat.bc[["Tumor_BAF_segmented"]][[1]][,1]
    names(tum.baf.segmented.vec) <- rownames(ascat.bc[["Tumor_BAF_segmented"]][[1]])

    warning("******************** FIXME: Disabling recalculation of BAFsegmented for now")

# Temporarily disable and so something more simple
# 
# 
 		#So now this is a little more complicated as the BAF segmented is only in the SNPs
 		#Would be nice to make the compaison with the germline here as well. We have the data so could...
 		for (i in seq_len(nrow(samp.seg.df))) { 
 
 		  #So first get the indices of the BAF indices within the segment
 		  seg.hz.snp.df <- gl.hz.snp.df[gl.hz.snp.df$snp_idx <= samp.seg.df[i,]$end_snp_idx & gl.hz.snp.df$snp_idx >= samp.seg.df[i,]$start_snp_idx,]
 		  if (nrow(seg.hz.snp.df) == 0) {
 		  	next
 		  }
 
 		  #Get indices of hz SNPs that are present ... so we need some form of mapping here to the smaller
 		  #These data contain all samples from this patient:
 		  seg.min.hz.snp.idx<- seg.hz.snp.df[1,]$hz_snp_idx
 		  seg.max.hz.snp.idx<- seg.hz.snp.df[nrow(seg.hz.snp.df),]$hz_snp_idx
 		  #Let's get the germline baf data specfically so that we can perform a crude correction
 		  gl.col.name <- grep("BAF", colnames(seg.hz.snp.df))
 		  gl.seg.hz.snp.df <- seg.hz.snp.df[, c(1, 2, gl.col.name)]
 		  colnames(gl.seg.hz.snp.df)[3] <- "BAF"
 		  #Sample specific!
 		  #So let's get the SNPs that are within the phased baf df
 		  seg.ph.hz.snp.df <- samp.pha.baf.df[samp.pha.baf.df$snp_idx <= samp.seg.df[i,]$end_snp_idx & samp.pha.baf.df$snp_idx >= samp.seg.df[i,]$start_snp_idx,]
 
 		  if(nrow(seg.ph.hz.snp.df) != nrow(gl.seg.hz.snp.df)) {
 		  	message("Number of phased SNPs is different to the actual number hz snps present in seg, perhaps to few to phase by MPhase?")
 		  	message("nrow(seg.ph.hz.snp.df) = ", nrow(seg.ph.hz.snp.df), " nrow(gl.seg.hz.snp.df) = ", nrow(gl.seg.hz.snp.df))
 		  	next
 		  }
 
 		  #GL specific
 		  gl.seg.hz.snp.df <- cbind(gl.seg.hz.snp.df, seg.ph.hz.snp.df[, c("maj_is_cpnA_logic_vec", "maj_is_cpnB_logic_vec", "snp_pos", "snp_idx")])
 
 		  if (nrow(seg.hz.snp.df) == 0 | nrow(seg.ph.hz.snp.df) == 0) {
 		  	message("No heterozygous SNPs present within segment.")
 		  }
 
 		  if (nrow(seg.hz.snp.df) != nrow(seg.ph.hz.snp.df)) {
 		  	message("Differing numbers of hz snps in the original hz snp df and the phased. Not sure why.")
 		  }
 
 		  #Ok so...
 		  #Let's get the phasing and calculate difference for the sample and the germline
 		  #We're assuming that these are in the same order at this point:
 		  #germline segmentation value:
 		  #Get groups of indices...
 		  haplotype.A.snp.pos.vec <- seg.ph.hz.snp.df[seg.ph.hz.snp.df$maj_is_cpnA_logic_vec == 1,]$snp_pos
 		  m.hap.A.BAF <- suppressWarnings(mean(seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.A.snp.pos.vec,]$BAF, na.rm=TRUE))
 		  haplotype.B.snp.pos.vec <- seg.ph.hz.snp.df[seg.ph.hz.snp.df$maj_is_cpnA_logic_vec == 0,]$snp_pos
 		  m.hap.B.BAF <- suppressWarnings(mean(seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.B.snp.pos.vec,]$BAF, na.rm=TRUE))
 
 		  #Now let's do the same for the germline data.. this actually isn't necessary..
 		  gl.m.hap.A.BAF <- suppressWarnings(mean(gl.seg.hz.snp.df[gl.seg.hz.snp.df$snp_pos %in% haplotype.A.snp.pos.vec,]$BAF, na.rm=TRUE)) - 0.5
 		  gl.m.hap.B.BAF <- suppressWarnings(mean(gl.seg.hz.snp.df[gl.seg.hz.snp.df$snp_pos %in% haplotype.B.snp.pos.vec,]$BAF, na.rm=TRUE)) - 0.5
 
 		  #So we actually
 		  #So we'll actually work out the gl correction factor here?
 
 		  #Let's get the equivalent germline snps 
 		  #Default should be to make the BAF segmented 0.5 if we don't have any info w/regards to hz snps
       seg.md.ph.BAF.vec <- c()
       # We could attempt to add in a germline correction here...
       if(! is.na(m.hap.A.BAF) && m.hap.A.BAF > 0.5) {
         seg.md.ph.BAF.vec <- c(seg.md.ph.BAF.vec, abs(seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.A.snp.pos.vec,]$BAF - 1) - gl.m.hap.A.BAF)
       } else if (! is.na(m.hap.A.BAF) && m.hap.A.BAF < 0.5) {
         seg.md.ph.BAF.vec <- c(seg.md.ph.BAF.vec, seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.A.snp.pos.vec,]$BAF - gl.m.hap.A.BAF)
       }
 
 		  #We could attempt to add in a germline correction here...
 		  if(! is.na(m.hap.B.BAF) && m.hap.B.BAF > 0.5) {
 		  	seg.md.ph.BAF.vec <- c(seg.md.ph.BAF.vec, abs(seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.B.snp.pos.vec,]$BAF - 1) - gl.m.hap.B.BAF)
 		  } else if (! is.na(m.hap.B.BAF) && m.hap.B.BAF < 0.5) {
 		    seg.md.ph.BAF.vec <- c(seg.md.ph.BAF.vec, seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.B.snp.pos.vec,]$BAF - gl.m.hap.B.BAF)
 		  }
 
 		  #Format
 		  if (length(seg.md.ph.BAF.vec) > 0) {
 		  	seg.BAF.segmented.val <- rep(mean(seg.md.ph.BAF.vec, na.rm=TRUE), length(seg.min.hz.snp.idx:seg.max.hz.snp.idx))
 		  } else {
 		  	seg.BAF.segmented.val <- rep(0.5, length(seg.min.hz.snp.idx:seg.max.hz.snp.idx)) 
 		  }
 
 		  #Now let's update the segmented baf values with the correct indices in the vector
 		  #So to do this we need to check the actual
 		  tum.baf.segmented.vec[seg.min.hz.snp.idx:seg.max.hz.snp.idx] <- seg.BAF.segmented.val
       #rownames(tum.baf.segmented.vec)
 
 		  #Now we could edit this to compare to the germline and work out whether we should actually call this as separate....
 
 		}#end 	for (i in 1:nrow(samp.seg.df)){ 

		if(any(is.na(tum.logr.segmented.vec))) {
			message("There are NA values in the tum.logr.segmented.vec after applying phasing which need to be removed. Revisit generation of this input. For changing the NAs to 0.")
			na.segmented.idx.vec <- which(is.na(tum.logr.segmented.vec))
			tum.logr.segmented.vec[na.segmented.idx.vec] <- 0
		}

		#First let's do the free runs:
    message("Setting LogR segmented")
		ascat.bc[["Tumor_LogR_segmented"]][,1] <- tum.logr.segmented.vec
    message("tum.baf.segmented.vec length =", length(tum.baf.segmented.vec))
    message("ascat.bc[[Tumor_BAF_segmented]][[1]] length =", nrow(ascat.bc[["Tumor_BAF_segmented"]][[1]]))
		#ascat.bc[["Tumor_BAF_segmented"]][[1]][,1] <- tum.baf.segmented.vec
		#ascat.bc[["Tumor_BAF_segmented"]][[1]] <- matrix(tum.baf.segmented.vec, nrow=1, dimnames(list(as.character(seq_along(snps$BAFsegmented)))))
    #print(tum.baf.segmented.vec)
    message("Finished setting BAF segmented")
    message("Finished setting segmented")

    ascat_phased <- ascat.bc

		#End of section generating the spoofed ascat.bc object that utilises the phasing information. This should be made into a function.

		# #TODO there may be issues when no soution that passes ASCAT's goodness of fit thresholds.
		# #This will prevent any output in the runs that don't have a specified purity and ploidy. 
		# print("Doing the free run ASCAT steps....")

		# if (purity.est.vec[samp]  < 0.15) {

		# 	MINRHO <- 0.04#This is normally set within the ASCAT script but here I'm modifying based on the input...

		# } else {

		# 	MINRHO <- 0.10#These should probably be constants set at the to of the script.
		# }

		# setwd(free.run.dir)
		# ascat.plotSegmentedData(ascat.bc)
		# ascat.output = ascat.runAscat(ascat.bc,gamma=1)#, rho_manual = rho_solution, psi_manual = psi_solution)
		# save(ascat.bc, file=paste(samp, 'ascat.bc.RData',sep='.'))
		# save(ascat.output, file=paste(samp, 'ascat.output.RData',sep='.'))
		# #Calculate number of heterozygous positions per segment
		# tmp <- ascat.output$segments_raw
		# n.het <- 0#FIXME should really ad
		# acf <- ascat.output$aberrantcellfraction
		# ploidy <- ascat.output$ploidy
		
		# tmp <- cbind(tmp[,c(1:4)], n.het, tmp[,5] + tmp[,6], tmp[,5:6], rep(ploidy, nrow(tmp)), rep(acf, nrow(tmp)),tmp[,7],tmp[,8])
		# colnames(tmp)[c(6,9:12)] <- c('cnTotal', 'Ploidy', 'ACF', 'nAraw','nBraw')
		# chr <- sub('chr','', tmp[,2])
		# chr[chr %in% c('x','X')] <- 23
		# chr[chr %in% c('y', 'Y')] <- 24
		# tmp[,2] <- as.numeric(chr)
		# #add a column defining gain and losslength
		# callState <- rep('neutral', nrow(tmp))
		# callState[log2((tmp[,'nAraw']+tmp[,'nBraw'])/tmp[,'Ploidy']) >= log2(2.5/2)] <- 'Gain'
		# callState[log2((tmp[,'nAraw']+tmp[,'nBraw'])/tmp[,'Ploidy']) <= log2(1.5/2)] <- 'Loss'
		# callState[(tmp[,'nAraw']+tmp[,'nBraw']) >= ((2*tmp[,'Ploidy'])+1)] <- 'Amp'
		# callState[(tmp[,'nAraw']+tmp[,'nBraw']) <= 0.5] <- 'Del'

		# tmp <- cbind(tmp,callState)
		# seg <- tmp
		# save(seg, file=paste(samp,'ascat.seg.RData',sep='.'))

		# rm(ascat.output)
		# rm(seg)

		#So what we want to do here is to generate all of the solutions that ASCAT proposed using the phased data... 
		#TODO - pull out the code that saves the solutions
		#and get a loop that reruns ASCAT using all of them, each in its own directory.

    return(ascat_phased)

}
