
joinSegAndReturnNewArmGranges <- function(arm.idx.vec, curr.gr.idx, adj.idx, adj.segs.mat, sorted.arm.grs) {

	print("called joinSegAndReturnNewArmGranges")
	print("arm.idx.vec")
	print(arm.idx.vec)
	print("curr.gr.idx")
	print(curr.gr.idx)
	print("adj.segs.mat")
	print(adj.segs.mat)
	print("sorted.arm.grs")
	print(sorted.arm.grs)

	kept.arm.idxs <- arm.idx.vec[ ! arm.idx.vec %in% c(curr.gr.idx, adj.segs.mat[adj.idx, "gr_idx"] ) ]

	print("sorted.arm.grs[curr.gr.idx]")
	print(sorted.arm.grs[curr.gr.idx])

	print("adj.segs.mat[adj.idx,]")
	print(adj.segs.mat[adj.idx,])

	if (length(kept.arm.idxs) > 0) {

		kept.arm.chr.vec <- as.numeric(as.character(seqnames(sorted.arm.grs[kept.arm.idxs])))
		kept.arm.start.vec <- start(sorted.arm.grs[kept.arm.idxs])
		kept.arm.end.vec <- end(sorted.arm.grs[kept.arm.idxs])
		
	} else {

		kept.arm.chr.vec <- c()
		kept.arm.start.vec <- c()
		kept.arm.end.vec <- c()

	}

	new.seg.chr <- adj.segs.mat[adj.idx, "chrom"]
	new.seg.start <- min(adj.segs.mat[adj.idx, "startpos"], start(sorted.arm.grs[curr.gr.idx]))
	new.seg.end <- max(adj.segs.mat[adj.idx, "endpos"], end(sorted.arm.grs[curr.gr.idx]))

	print("new.seg.chr")
	print(new.seg.chr)
	print("new.seg.start")
	print(new.seg.start)
	print("new.seg.end")
	print(new.seg.end)


	new.arm.chr.vec <- c(new.seg.chr, kept.arm.chr.vec)
	new.arm.start.vec <- c(new.seg.start, kept.arm.start.vec)
	new.arm.end.vec <- c(new.seg.end, kept.arm.end.vec)

	new.arm.grs <- GRanges(seqnames=new.arm.chr.vec,
							ranges=IRanges(new.arm.start.vec, new.arm.end.vec),
							strand=rep("*", length(new.arm.chr.vec)))
	print("new.arm.grs")
	print(new.arm.grs)

	return(new.arm.grs)

}#end joinSegAndReturnNewArmGranges

#Try to join using phasing 
attemptToJoinSegbyPhasedCPN <- function(curr.seg.mat, adj.segs.mat, samp.arm.sites.df, arm.idx.vec, curr.gr.idx, sorted.arm.grs, MIN_PH_CPN) {
print("called attemptToJoinSegbyPhasedCPN")
	new.unsorted.arm.grs <- sorted.arm.grs

	# print("curr.seg.mat")
	# print(curr.seg.mat)
	# print("adj.segs.mat")
	# print(adj.segs.mat)

	for (k in 1:nrow(adj.segs.mat)) {

		if ( curr.seg.mat[, "num_ph_cpn"] >= MIN_PH_CPN & adj.segs.mat[k,  "num_ph_cpn"] >= MIN_PH_CPN) {

			cur.ph.cpn.vec <- samp.arm.sites.df[samp.arm.sites.df$pos >= curr.seg.mat[,"startpos"]
											  & samp.arm.sites.df$pos <= curr.seg.mat[, "endpos"], cpnX.string]

			adj.ph.cpn.vec <- samp.arm.sites.df[samp.arm.sites.df$pos >= adj.segs.mat[k, "startpos"]
											  & samp.arm.sites.df$pos <= adj.segs.mat[k, "endpos"], cpnX.string]

			print("!!!!!!!!!PHASED COPY NUMBER JOINING ATTEMPTED!!!!!!!!!!!!!!!")
			print(curr.seg.mat[, "num_ph_cpn"])
			print(adj.segs.mat[k,  "num_ph_cpn"])
			ph.p.val <- wilcoxTestWithErrorHandling(cur.ph.cpn.vec, adj.ph.cpn.vec)$p.val

			if (ph.p.val > PVAL_THRESH_TO_JOIN_PH_SEGS) {
				print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
				print("!!!!!!!!!!!!!!!!!!! SEGS ARE JOINED !!!!!!!!!!!!!!!!!!!!!!!!")
				print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
					#Join segment and make a new set of granges for the arm 
				new.unsorted.arm.grs <- joinSegAndReturnNewArmGranges(arm.idx.vec = arm.idx.vec,
												 					  curr.gr.idx = curr.gr.idx,
												 					  adj.idx = k,
																	  adj.segs.mat = adj.segs.mat,
																	  sorted.arm.grs = sorted.arm.grs)
				return(new.unsorted.arm.grs)


			} else {#if (ph.p.val > PVAL_THRESH_TO_JOIN_PH_SEGS) {=
				print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
				print("!!!!!!!!!!!!!!!!!! SEGS STAY SEPARATE !!!!!!!!!!!!!!!!!!!!!!")
				print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
			}
		}# if ( curr.seg.mat[,"num_ph_cpn"] >= MIN_PH_CPN & adj.segs.mat[ k,"num_ph_cpn"] >= MIN_PH_CPN)){
	}# for (k in 1:nrow(adj.segs.mat)) {

	return(new.unsorted.arm.grs)

}#end attemptToJoinSegbyPhasedCPN

#If join return
attemptToJoinSegbyTotalCPN <- function(curr.seg.mat, adj.segs.mat, samp.arm.sites.df, arm.idx.vec, curr.gr.idx, sorted.arm.grs, NUM_TOT_CPN) {
print("called attemptToJoinSegbyTotalCPN")
	new.unsorted.arm.grs <- sorted.arm.grs

	# print("curr.seg.mat")
	# print(curr.seg.mat)
	# print("adj.segs.mat")
	# print(adj.segs.mat)


	for (k in 1:nrow(adj.segs.mat)) {

		#If joining by phased copy number is possible don't join

		if (curr.seg.mat[,"num_tot_cpn"] >= NUM_TOT_CPN & adj.segs.mat[k, "num_tot_cpn"] >= NUM_TOT_CPN) {

			join.by.phased.cpn.not.possible <- curr.seg.mat[,"num_ph_cpn"] < MIN_PH_CPN | adj.segs.mat[k, "num_ph_cpn"] < MIN_PH_CPN 
			if (join.by.phased.cpn.not.possible == FALSE) {

				print("NOT ATTEMPTING TO JOIN BY TOT CPN AS ALREADY ATTEMPTED PHASED.")

			} else {

				cur.tot.cpn.vec <- rle(samp.arm.sites.df[samp.arm.sites.df$pos >= curr.seg.mat[,"startpos"]
													  & samp.arm.sites.df$pos <= curr.seg.mat[, "endpos"], "tot_raw_un_ph_cpn"])$values

				adj.tot.cpn.vec <- rle(samp.arm.sites.df[samp.arm.sites.df$pos >= adj.segs.mat[k, "startpos"]
													  & samp.arm.sites.df$pos <= adj.segs.mat[k, "endpos"], "tot_raw_un_ph_cpn"])$values

				print("!!!!!!!!!PHASED TOTAL NUMBER JOINING ATTEMPTED!!!!!!!!!!!!!!")
				print(curr.seg.mat[, "num_tot_cpn"])
				print(adj.segs.mat[k,  "num_tot_cpn"])
				tot.p.val <- wilcoxTestWithErrorHandling(cur.tot.cpn.vec, adj.tot.cpn.vec)$p.val
					
				if (tot.p.val > PVAL_THRESH_TO_JOIN_TOT_SEGS) {

					print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
					print("!!!!!!!!!!!!!!!!!!! SEGS ARE JOINED !!!!!!!!!!!!!!!!!!!!!!!!")
					print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
					new.unsorted.arm.grs <- joinSegAndReturnNewArmGranges(arm.idx.vec = arm.idx.vec,
													 					  curr.gr.idx = curr.gr.idx,
													 					  adj.idx = k,
																		  adj.segs.mat = adj.segs.mat,
																		  sorted.arm.grs = sorted.arm.grs)

					return(new.unsorted.arm.grs)

				} else {

					print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
					print("!!!!!!!!!!!!!!!!!! SEGS STAY SEPARATE !!!!!!!!!!!!!!!!!!!!!!")
					print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

				}#if (ph.p.val > PVAL_THRESH_TO_JOIN_PH_SEGS) {
			}#end 	if (join.by.phased.cpn.not.possible == TRUE) {
		}#end for (k in 1:nrow(adj.segs.mat)) {
	}#end for (k in 1:nrow(adj.segs.mat)) {

	return(new.unsorted.arm.grs)
}#end attemptToJoinSegbyTotalCPN


attemptToJoinSegbyMeanTotalCPNAlone <- function(curr.seg.mat, adj.segs.mat, samp.arm.sites.df, arm.idx.vec, curr.gr.idx, sorted.arm.grs, NUM_TOT_CPN, CPN_AVG_THRESH) {
print("called attemptToJoinSegbyMeanTotalCPNAlone")
	new.unsorted.arm.grs <- sorted.arm.grs

	# print("curr.seg.mat")
	# print(curr.seg.mat)
	# print("adj.segs.mat")
	# print(adj.segs.mat)

	#Allow joining by means...
	for (k in 1:nrow(adj.segs.mat)) {

		if (curr.seg.mat[, "num_tot_cpn"] > 0 & curr.seg.mat[, "num_tot_cpn"] < NUM_TOT_CPN & adj.segs.mat[k, "num_tot_cpn"] > 0) {

			curr.seg.mean <- curr.seg.mat[, "mean_tot_cpn"]
			curr.seg.median <- curr.seg.mat[, "median_tot_cpn"]

			adj.seg.mean <- adj.segs.mat[k, "mean_tot_cpn"]
			adj.seg.median <- adj.segs.mat[k, "median_tot_cpn"]

			if ( abs(curr.seg.mean - adj.seg.mean) <= CPN_AVG_THRESH | abs(curr.seg.median - adj.seg.median) <= CPN_AVG_THRESH ) {

				new.unsorted.arm.grs <- joinSegAndReturnNewArmGranges(arm.idx.vec = arm.idx.vec,
											 						  curr.gr.idx = curr.gr.idx,
											 						  adj.idx = k,
																	  adj.segs.mat = adj.segs.mat,
																	  sorted.arm.grs = sorted.arm.grs)

				return(new.unsorted.arm.grs)
			}#end if ( abs(curr.seg.mean - adj.seg.mean) <= 0.1 | abs(curr.seg.median - adj.seg.median) <= 0.1 ) {
		} else {

			print("Pair of segs do not qualify for mean joining.")
		}#if (curr.seg.mat[1 ,"num_tot_cpn"] > 0 & curr.seg.mat[1 ,"num_tot_cpn"] < NUM_TOT_CPN & any(adj.segs.mat[, "num_tot_cpn"] > 0)) {
	}#end for (k in 1:nrow(adj.segs.mat)) {

	return(new.unsorted.arm.grs)
	
}#end attemptToJoinSegbyMeanTotalCPNAlone


joinBlankSegsToMostConfidentCall <- function(curr.seg.mat, adj.segs.mat, samp.arm.sites.df, arm.idx.vec, curr.gr.idx, sorted.arm.grs) {

	print("called joinBlankSegsToMostConfidentCall")
	print("arm.idx.vec")
	print(arm.idx.vec)
	print("curr.gr.idx")
	print(curr.gr.idx)
	print("adj.segs.mat")
	print(adj.segs.mat)
	print("sorted.arm.grs")
	print(sorted.arm.grs)

	new.unsorted.arm.grs <- sorted.arm.grs

	adj.segs.mat <- adj.segs.mat[order(adj.segs.mat[,"num_ph_cpn"],adj.segs.mat[,"num_tot_cpn"]), ,drop = FALSE]
	print("adj.segs.mat")
	print(adj.segs.mat)

	adj.segs.mat <- adj.segs.mat[1, , drop = FALSE]

	print("adj.segs.mat")
	print(adj.segs.mat)
	#Allow joining by means...
	for (k in 1:nrow(adj.segs.mat)) {

		if (curr.seg.mat[, "num_tot_cpn"] == 0) {


				new.unsorted.arm.grs <- joinSegAndReturnNewArmGranges(arm.idx.vec = arm.idx.vec,
											 						  curr.gr.idx = curr.gr.idx,
											 						  adj.idx = k,
																	  adj.segs.mat = adj.segs.mat,
																	  sorted.arm.grs = sorted.arm.grs)

				return(new.unsorted.arm.grs)

		} else {

			print("No blank segs to join.")

		}#if (curr.seg.mat[1 ,"num_tot_cpn"] > 0 & curr.seg.mat[1 ,"num_tot_cpn"] < NUM_TOT_CPN & any(adj.segs.mat[, "num_tot_cpn"] > 0)) {
	}#end for (k in 1:nrow(adj.segs.mat)) {

	return(new.unsorted.arm.grs)
  

}



creatingArmJoiningInfoMatrixAndGranges <- function(arm.grs, samp.arm.sites.df, cpnX.string) {

	proc.arm.grs <- arm.grs
	proc.arm.grs$width <- width(proc.arm.grs)
	joining.vals.mat.col.names.vec <- c("chrom", "startpos", "endpos", "gr_idx",
									"num_ph_cpn", "mean_ph_cpn", "median_ph_cpn",
									"num_tot_cpn", "mean_tot_cpn", "median_tot_cpn")
	joining.vals.mat <- matrix(ncol = 10, nrow=length(proc.arm.grs), 0)
	colnames(joining.vals.mat) <- joining.vals.mat.col.names.vec

	for ( idx in 1:length(proc.arm.grs) ) {
		#Get the idx of the grange in the sorted by size arm grs list
		joining.vals.mat[idx, "chrom"] <- as.numeric(as.character(seqnames(proc.arm.grs[idx])))
		joining.vals.mat[idx, "startpos"] <- start(proc.arm.grs[idx])
		joining.vals.mat[idx, "endpos"] <- end(proc.arm.grs[idx])

		seg.sites.info.df <- samp.arm.sites.df[samp.arm.sites.df$chrom == as.numeric(as.character(seqnames(proc.arm.grs[idx])))
											&	samp.arm.sites.df$pos >= start(proc.arm.grs[idx])
											&  samp.arm.sites.df$pos <= end(proc.arm.grs[idx]) ,]
		seg.ph.cpn.sites.df <- seg.sites.info.df[! is.na(seg.sites.info.df[, cpnX.string]), ]						
		num.ph.cpn <- length(seg.ph.cpn.sites.df[,cpnX.string])
		if (num.ph.cpn > 0) {

			joining.vals.mat[idx, "num_ph_cpn"]    <- num.ph.cpn
			joining.vals.mat[idx, "mean_ph_cpn"]   <- mean(seg.ph.cpn.sites.df[,cpnX.string])
			joining.vals.mat[idx, "median_ph_cpn"] <- median(seg.ph.cpn.sites.df[,cpnX.string])

		}#end if (num.ph.cpn > 0) {
		#using run length encodings here as we have the issue that many logr values, due to being 1000G positiosn
		#are often very close together and thus are, in reality, a single datapoint...
		#Could be a problem with very low coverage samples...?todo, check this
		num.tot.cpn <- length(rle(seg.sites.info.df$tot_raw_un_ph_cpn)$values)
		if (num.tot.cpn > 0) {

			joining.vals.mat[idx, "num_tot_cpn"]    <- num.tot.cpn
			joining.vals.mat[idx, "mean_tot_cpn"]   <- mean(rle(seg.sites.info.df$tot_raw_un_ph_cpn)$values)
			joining.vals.mat[idx, "median_tot_cpn"] <- median(rle(seg.sites.info.df$tot_raw_un_ph_cpn)$values)

		}#end if (num.tot.cpn > 0) {

	}#end for ( ord.idx in 1:length(all.segs.considered.idx.vec) ) {

	proc.arm.grs$num_ph_cpn <- joining.vals.mat[,"num_ph_cpn"]
	proc.arm.grs$num_tot_cpn <- joining.vals.mat[,"num_tot_cpn"]
	sorted.arm.grs <- proc.arm.grs[ order(proc.arm.grs$num_ph_cpn, proc.arm.grs$num_tot_cpn, proc.arm.grs$width) ]
	joining.vals.mat <- joining.vals.mat[match(start(sorted.arm.grs), joining.vals.mat[, "startpos"]) , , drop = FALSE]
	joining.vals.mat[, "gr_idx"] <- 1:nrow(joining.vals.mat)

	return(list(joining.vals.mat, sorted.arm.grs))

}#end creatingArmJoiningInfoMatrix <- function(arm.grs, samp.arm.sites.df) {


attemptToJoinOneSegmentWithinChromArm <- function(current.arm.grs, joining.vals.mat, samp.arm.sites.df, sorted.arm.grs, cpnX.string,MIN_PH_CPN, NUM_TOT_CPN, CPN_AVG_THRESH) {

	#Edge case where this only one segment:
	if(length(current.arm.grs)  == 1) {

		return(current.arm.grs)

	}

	arm.idx.vec <- 1:length(current.arm.grs)
	old.length.current.arm.grs <- length(current.arm.grs)
	curr.ord.idx <- 1
	while ( curr.ord.idx <= length(current.arm.grs) & length(current.arm.grs) == old.length.current.arm.grs ) {

		curr.gr <- current.arm.grs[curr.ord.idx]
		precede.and.follow.gr.idx.vec <- c(follow(curr.gr, current.arm.grs), precede(curr.gr, current.arm.grs))
		sort.non.na.gr.pre.and.fol.idx.vec <- sort(precede.and.follow.gr.idx.vec[ ! is.na(precede.and.follow.gr.idx.vec)])
		curr.seg.mat <- joining.vals.mat[joining.vals.mat[,"gr_idx"] == curr.ord.idx, ,drop = FALSE]

		print("curr.seg.mat")
		print(curr.seg.mat)
		adj.segs.mat <- joining.vals.mat[sort.non.na.gr.pre.and.fol.idx.vec, ,drop = FALSE]
		print("adj.segs.mat")
		print(adj.segs.mat)

		new.unsorted.arm.grs <- attemptToJoinSegbyPhasedCPN(curr.seg.mat = curr.seg.mat,
																	adj.segs.mat = adj.segs.mat,
																	samp.arm.sites.df = samp.arm.sites.df,
																	arm.idx.vec = arm.idx.vec,
																	curr.gr.idx = curr.ord.idx,
																	sorted.arm.grs = current.arm.grs,
																	MIN_PH_CPN = MIN_PH_CPN)

		if( length(new.unsorted.arm.grs) != length(current.arm.grs)) {

			print("Joined segs using phased cpn.")
			current.arm.grs <- new.unsorted.arm.grs
			next

		}
		print("length(current.arm.grs)")
		print(length(current.arm.grs))
		new.unsorted.arm.grs <- attemptToJoinSegbyTotalCPN(curr.seg.mat = curr.seg.mat,
																	adj.segs.mat = adj.segs.mat,
																	samp.arm.sites.df = samp.arm.sites.df,
																	arm.idx.vec = arm.idx.vec,
																	curr.gr.idx = curr.ord.idx,
																	sorted.arm.grs = current.arm.grs,
																	NUM_TOT_CPN = NUM_TOT_CPN)
		print("length(new.unsorted.arm.grs)")
		print(length(new.unsorted.arm.grs))
		print("Returning from attemptToJoinSegbyTotalCPN")
		if( length(new.unsorted.arm.grs) != length(current.arm.grs)) {

			print("Joined segs using tot cpn.")
			current.arm.grs <- new.unsorted.arm.grs
			next

		}

		new.unsorted.arm.grs <- attemptToJoinSegbyMeanTotalCPNAlone(curr.seg.mat = curr.seg.mat,
																	adj.segs.mat = adj.segs.mat,
																	samp.arm.sites.df = samp.arm.sites.df,
																	arm.idx.vec = arm.idx.vec,
																	curr.gr.idx = curr.ord.idx,
																	sorted.arm.grs = current.arm.grs,
																	NUM_TOT_CPN = NUM_TOT_CPN,
																	CPN_AVG_THRESH = CPN_AVG_THRESH)  
		print("Returning from current.arm.grs")
		if( length(new.unsorted.arm.grs) != length(current.arm.grs)) {

			print("Joined segs using mean tot cpn.")
			current.arm.grs <- new.unsorted.arm.grs
			next

		}
		#joinBlankSegsToMostConfidentCall <- function(curr.seg.mat, adj.segs.mat, samp.arm.sites.df, arm.idx.vec, curr.gr.idx, sorted.arm.grs) {
		new.unsorted.arm.grs <- joinBlankSegsToMostConfidentCall(curr.seg.mat = curr.seg.mat,
																adj.segs.mat  = adj.segs.mat,
																samp.arm.sites.df = samp.arm.sites.df,
																arm.idx.vec = arm.idx.vec,
																curr.gr.idx = curr.ord.idx,
																sorted.arm.grs = current.arm.grs)
		if( length(new.unsorted.arm.grs) != length(current.arm.grs)) {

			print("Joined segs using mean tot cpn.")
			current.arm.grs <- new.unsorted.arm.grs
			next

		}




		if( length(new.unsorted.arm.grs) != length(current.arm.grs)) {

			current.arm.grs <- new.unsorted.arm.grs
			next

		}
		
		curr.ord.idx <- curr.ord.idx + 1

	}#end while ( curr.gr.idx <= length(sorted.arm.grs) ) {

	return(current.arm.grs)

}#end attemptToJoinOneSegmentsWithinChromArm











# #So now we want to know how many 
# #full_sample <- sapply(rownames(all.bp.all.samps.segs.df), function(x) strsplit(x, split = ":")[[1]][1])
# 	joined.samp.arms <- list()

# 	for( samp in sample.vec) {

# 		for ( chr.arm in chr.arms.present.vec) {

# 			samp.arm.sites.df <- samp.arms[[samp]][[chr.arm]]
# 			print("head(samp.arm.sites.df)")
# 			print(head(samp.arm.sites.df))
# 			samp.chr.arm.segs.df <- full.joint.segs.df[full.joint.segs.df$chr_arm %in% chr.arm
# 							                         & full.joint.segs.df$sample %in% samp, ]
# 			samp.chr.arm.seg.name.vec <- samp.chr.arm.segs.df$seg_name

# 			#samp.arm.sites.df <- joinable.df

# 			chr.vec <- as.numeric(sapply(samp.chr.arm.seg.name.vec, function(x) strsplit(x, split = ":")[[1]][1]))
# 			start.vec <- as.numeric(sapply(samp.chr.arm.seg.name.vec, function(x) strsplit(x, split = ":")[[1]][2]))
# 			end.vec <- as.numeric(sapply(samp.chr.arm.seg.name.vec, function(x) strsplit(x, split = ":")[[1]][3]))
# 			arm.grs <- GRanges(seqnames=chr.vec,
# 							ranges=IRanges(start.vec, end.vec),
# 							strand=rep("*", length(chr.vec)))


# 			print("arm.grs")
# 			print(arm.grs)
# 			if(length(getOverlapGRanges(arm.grs, arm.grs)) != length(arm.grs)) {

# 				print(getOverlapGRanges(arm.grs, arm.grs))
# 				stop("Error, overlapping granges within an arm!")
# 			}

# 			current.arm.grs <- arm.grs
# 			old.length.current.arm.grs <- 0
# 			while (length(current.arm.grs) != old.length.current.arm.grs & length(current.arm.grs) != 1) {

# 				joining.vals.mat.and.arm.grs.list <- creatingArmJoiningInfoMatrixAndGranges(current.arm.grs, samp.arm.sites.df, cpnX.string)
# 				joining.vals.mat <- joining.vals.mat.and.arm.grs.list[[1]]
# 				current.arm.grs <- joining.vals.mat.and.arm.grs.list[[2]]
# 				print("############")
# 				print(current.arm.grs)
# 				print(joining.vals.mat)
# 				print("%%%%%%%%%%%%%")
# 				# readline(prompt="Press [enter] to continue")
# 				old.length.current.arm.grs <- length(current.arm.grs)
# 				current.arm.grs <- attemptToJoinOneSegmentWithinChromArm(current.arm.grs, joining.vals.mat, samp.arm.sites.df, sorted.arm.grs, cpnX.string, MIN_PH_CPN, NUM_TOT_CPN, CPN_AVG_THRESH) 

# 				if(length(getOverlapGRanges(current.arm.grs, current.arm.grs)) != length(current.arm.grs)) {
# 					print(getOverlapGRanges(current.arm.grs, current.arm.grs))
# 					stop("Error, overlapping granges within an arm!")
# 				}
# 				print(length(current.arm.grs))
# 				print(old.length.current.arm.grs)
# 			} #end while (length(current.arm.grs) == old.length.current.arm.grs ) {


# 			print("Original arm grs")
# 			print(arm.grs)
# 			print("Final current.arm.grs")
# 			print(current.arm.grs)

# 		# readline(prompt="Press [enter] to continue")

# 		sorted.joining.vals.mat <- joining.vals.mat[order(joining.vals.mat[,"startpos"]),]
# 		joined.samp.arms[[samp]][[chr.arm]] <- as.data.frame(sorted.joining.vals.mat)



# 		}#end for ( chr.arm in chr.arms.present.vec) {

# 	} # end for( samp in sample.vec) {

# 	#Cool so now we can do both of these for cpn A and cpn B






