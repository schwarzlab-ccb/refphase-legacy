#Error checking look

head(a[, (ncol(a)-5):ncol(a)])

focal_gain_id@22:16000001:51304566
pj.mcs.df[pj.mcs.df$chrom %in% 22, (ncol(pj.mcs.df)-8):ncol(pj.mcs.df)]

#Original function looks like this:

getSubclonalFocalEvents <- function(possible.samp.segs.vec,
									mcs.sc.cpn.event.to.map.df,
									event.type) {


	possible.focal.segs.alone.vec <- sapply(possible.samp.segs.vec, function(x) strsplit(x, split="@")[[1]][2])
	possible.focal.segs.alone.df <- data.frame(chrom = sapply(possible.focal.segs.alone.vec, function(x) strsplit(x, split=":")[[1]][1]),
										 start = as.numeric(sapply(possible.focal.segs.alone.vec, function(x) strsplit(x, split=":")[[1]][2])),
										 end = as.numeric(sapply(possible.focal.segs.alone.vec, function(x) strsplit(x, split=":")[[1]][3])),
										stringsAsFactors=FALSE)
	possible.focal.segs.alone.grs <- makeGRangesFromDataFrame(possible.focal.segs.alone.df,
		                         keep.extra.columns=FALSE,
		                         ignore.strand=TRUE,
		                         seqinfo=NULL,
		                         seqnames.field=c("chrom"),
		                         start.field=c("start"),
		                         end.field=c("end"),
		                         starts.in.df.are.0based=FALSE)
	reduced.poss.f.sc.gr <- reduce(possible.focal.segs.alone.grs)
	reduced.poss.f.sc.gr$focal_event_id <- paste0("focal_", event.type, "_id@", as.character(seqnames(reduced.poss.f.sc.gr)), ":", start(reduced.poss.f.sc.gr), ":", end(reduced.poss.f.sc.gr))

	mcs.sc.cpn.event.to.map.gr <- makeGRangesFromDataFrame(mcs.sc.cpn.event.to.map.df,
		                         keep.extra.columns=TRUE,
		                         ignore.strand=TRUE,
		                         seqinfo=NULL,
		                         seqnames.field=c("chrom", "chr"),
		                         start.field=c("startpos"),
		                         end.field=c("endpos"),
		                         starts.in.df.are.0based=FALSE)

	print(mcs.sc.cpn.event.to.map.gr)
	print(reduced.poss.f.sc.gr)
	overlapping.f.sc.grs <- suppressWarnings(getOverlapGRanges(mcs.sc.cpn.event.to.map.gr, reduced.poss.f.sc.gr))
	print(overlapping.f.sc.grs)
	mcs.sc.cpn.event.to.map.df$focal_event_id <- "none"
	mcs.sc.cpn.event.to.map.df[overlapping.f.sc.grs$query_hit_idx, ]$focal_event_id <- overlapping.f.sc.grs$focal_event_id
	print(mcs.sc.cpn.event.to.map.df)
	return(mcs.sc.cpn.event.to.map.df)

	}#end getSubclonalFocalEvents




#So input is samp_seg_name vector that looks something like this:
#"MCS/L_LTX064_SU_T1-R7@18:53177938:54270110"

#Each focal type needs to be done separately


#N - 1 samples

#Check if first attempt or no change since last joining:
	#go to next n 

#Identify longest sequence (in terms of size) for all combinations of n-1 on each chromosome?

#Get ID from intersection/reduced total size

#Mark these segments/the largest area of them off limits now

#repeat 

cols.needed.vec <- c("sample", "uniq_joined_seg_name", "seg_name", "samp_seg_name","cpn_event_vs_ploidy","pj_focal_by_size_and_prop" ,"cl_focal_mcs_seg_cpn_vs_ploidy", "cl_focal_mcs_seg_LOH", "cl_focal_mcs_seg_LOH_or_loss","sc_focal_mcs_seg_cpn_vs_ploidy","sc_focal_mcs_seg_LOH")
pj.mcs.df[pj.mcs.df$chrom %in% 7,(ncol(pj.mcs.df)-7):ncol(pj.mcs.df)]
pj.mcs.df[pj.mcs.df$chrom %in% 22, "cpn_event_vs_ploidy", ]



event.type <- "gain"
qualifying.mcs.seg.vec <- pj.mcs.df[pj.mcs.df$arm %in% "11q" & pj.mcs.df$pj_focal_by_size_and_prop %in% c("focal")
							& pj.mcs.df$cpn_event_vs_ploidy %in% c("gain", "amp")
							& pj.mcs.df$across_samp_cl_cpn_vs_ploidy_present %in% c("subclonal_events_present")
							, c("chr", "startpos", "endpos","uniq_joined_seg_name", "samp_seg_name")]$samp_seg_name

#So we're going to do this by looking within arms...
#So let's cut down the datafram

arm.pj.mcs.df <- pj.mcs.df[pj.mcs.df$arm %in% "11q",]

#So we need to know which columns we're actually dealing with so that we can check them for pre-existing focal definitions:
if(event.type == "gain"){

	focal.col <- "sc_focal_mcs_seg_cpn_vs_ploidy"
	cpn.vs.ploidy.allowed <- c("gain", "amp")

} else if (event.type == "loss") {

	focal.col <- "sc_focal_mcs_seg_cpn_vs_ploidy"
	cpn.vs.ploidy.allowed <- c("loss", "deep")

} else if (event.type == "LOH") {

	focal.col <- "sc_focal_mcs_seg_cpn_vs_ploidy"
	cpn.vs.ploidy.allowed <- c("loss", "deep")


} else if (event.type == "loss_or_LOH") {


}





#Segs to consider

pj.mcs.focal.loss.df <- pj.mcs.df[pj.mcs.df$pj_focal_by_size_and_prop %in% c("focal") 
			 					& pj.mcs.df$cpn_event_vs_ploidy %in% c("loss", "deep_loss")
			 					& pj.mcs.df$across_samp_cl_cpn_vs_ploidy_present %in% c("subclonal_events_present"),
			  							c("chr", "startpos", "endpos", "uniq_joined_seg_name", "samp_seg_name")]






#potential.cl.foc.mcs.segs.vec all the samp_seg_names that match the condition of the subclonal focal we're looking for
#all chr/arm post joining minimum consistent segmentation segmentation dataframe
#focal column that will be updated by this function
#What is the type of the event we're looking to capture and the string will be placed in the focal id name
markClonalFocalEvents <- function(potential.cl.foc.mcs.segs.vec, pj.mcs.df, col.to.update, event.type) {

	arm.vec <- unique(pj.mcs.df$arm)
	num.samples  <- length(unique(pj.mcs.df$sample))

	for (arm in arm.vec ) {

		arm.pj.mcs.df <- pj.mcs.df[pj.mcs.df$arm %in% arm,]
		arm.pj.mcs.df$old_focal_gain_definition <- arm.pj.mcs.df[, col.to.update]
		arm.pj.mcs.df[, col.to.update] <- "none"
		#So now we need to create the column that tell sus whether we actually can look for runs:
		arm.pj.mcs.df$potential_new_focal <- arm.pj.mcs.df$samp_seg_name %in% potential.cl.foc.mcs.segs.vec
														
		arm.unclassified.focal.mat <- matrix(FALSE, nrow=length(unique(arm.pj.mcs.df$sample)), ncol=length(unique(arm.pj.mcs.df$seg_name)))
		colnames(arm.unclassified.focal.mat) <- unique(arm.pj.mcs.df$seg_name)
		rownames(arm.unclassified.focal.mat) <- unique(arm.pj.mcs.df$sample)
		#arm.unclassified.focal.mat.allele <- matrix(FALSE, nrow=length(unique(arm.pj.mcs.df$sample)), ncol=length(unique(arm.pj.mcs.df$seg_name)))
		#colnames(arm.unclassified.focal.mat.allele)  <- colnames(arm.unclassified.focal.mat) <- unique(arm.pj.mcs.df$seg_name)
		#rownames(arm.unclassified.focal.mat.allele)  <- rownames(arm.unclassified.focal.mat) <- unique(arm.pj.mcs.df$sample)
		for(samp in unique(arm.pj.mcs.df$sample)) {
			for(seg in unique(arm.pj.mcs.df$seg_name)){
				arm.unclassified.focal.mat[samp, seg] <- arm.pj.mcs.df[arm.pj.mcs.df$sample == samp & arm.pj.mcs.df$seg_name == seg,]$potential_new_focal
				#arm.unclassified.focal.mat.allele[samp, seg] <- arm.pj.mcs.df[arm.pj.mcs.df$sample == samp & arm.pj.mcs.df$seg_name == seg,]$dominant_allele
			}#end for(seg in unique(arm.pj.mcs.df$seg_name)){
		}#end for(samp in unique(arm.pj.mcs.df$sample)) {

		#The number of unique 
		f.acr.samp.vec <- colSums(arm.unclassified.focal.mat) == num.samples
		unclassified.cl.focals.remain <- any(f.acr.samp.vec)

		while (unclassified.cl.focals.remain) {
			print(unclassified.cl.focals.remain)

			cl.foc.rle.df <- with(rle(f.acr.samp.vec), 
								data.frame(number = values,
											start = cumsum(lengths) - lengths + 1,
											end = cumsum(lengths))[order(values),])

			curr.cl.fcl.rle.idx <- min(which(cl.foc.rle.df$number == TRUE))
			curr.start.cl.fcl.idx <- cl.foc.rle.df[curr.cl.fcl.rle.idx,]$start
			curr.end.cl.fcl.idx <- cl.foc.rle.df[curr.cl.fcl.rle.idx,]$end
			#So now we want the middle
			mid.cl.fcl.idx <- curr.start.cl.fcl.idx + ceiling((curr.end.cl.fcl.idx - curr.start.cl.fcl.idx) / 2)

			#So now we need to know whether there are contiguous segments
			centre.seg.name <- colnames(arm.unclassified.focal.mat)[mid.cl.fcl.idx]
			samps.w.centre.seg <- rownames(arm.unclassified.focal.mat[arm.unclassified.focal.mat[,max.seg.name]==TRUE,max.seg.name, drop = FALSE])
			joined.seg.size.vec <- arm.pj.mcs.df[arm.pj.mcs.df$sample %in% samps.w.max.seg & arm.pj.mcs.df$seg_name %in% max.seg.name,]$MPhase_joined_seg_width
			print("joined.seg.size.vec")
			print(joined.seg.size.vec)
			names(joined.seg.size.vec) <- arm.pj.mcs.df[arm.pj.mcs.df$sample %in% samps.w.max.seg & arm.pj.mcs.df$seg_name %in% max.seg.name,]$uniq_joined_seg_name
			joined.focal.seg.matching.back.to <- names(which.max(joined.seg.size.vec))
			print(joined.seg.size.vec)
			#So now we have the segment we're matching back to, we need its start and end so we can give these segments this id
			chr.start.end.vec <- as.numeric(strsplit(strsplit(joined.focal.seg.matching.back.to, split="@")[[1]][2], split = ":")[[1]])
			#So now we can make the new ID
			new.focal.id <- paste0("focal_", event.type,"_id","@", chr.start.end.vec[1],":",chr.start.end.vec[2],":",chr.start.end.vec[3])
			#We can also rejig the matrix to eliminate this from a potential new focal.
			start.col.idx <- grep(chr.start.end.vec[2], colnames(arm.unclassified.focal.mat))
			end.col.idx <- grep(chr.start.end.vec[3], colnames(arm.unclassified.focal.mat))
			seg.names.to.update <- colnames(arm.unclassified.focal.mat)[start.col.idx:end.col.idx]
			arm.unclassified.focal.mat[,start.col.idx:end.col.idx] <- FALSE
			#Update the original matrix with the new focal ids
			arm.pj.mcs.df[arm.pj.mcs.df$sample%in%samps.w.centre.seg & arm.pj.mcs.df$seg_name%in%seg.names.to.update, col.to.update] <- new.focal.id
			f.acr.samp.vec <- colSums(arm.unclassified.focal.mat) == length(unique(arm.pj.mcs.df$sample))
			unclassified.cl.focals.remain <- any(f.acr.samp.vec)
			print(new.focal.id)
			print("##########3")

		}#end while (unclassified.focals.remain) {

			#Update dataframe with this arms results
			pj.mcs.df[pj.mcs.df$arm %in% arm, col.to.update] <- arm.pj.mcs.df[, col.to.update]

	}#end for (arm in arm.vec ) {

	#Return fully updated dataframe
	return(pj.mcs.df)

}#end markClonalFocalEvents <- function(potential.cl.foc.mcs.segs.vec, pj.mcs.df, col.to.update, event.type) {



#Note this function assumes that all clonal focal samples of this type have already been removed from the analysis.

#potential.sc.foc.mcs.segs.vec all the samp_seg_names that match the condition of the subclonal focal we're looking for
#all chr/arm post joining minimum consistent segmentation segmentation dataframe
#focal column that will be updated by this function
#What is the type of the event we're looking to capture and the string will be placed in the focal id name
markSubclonalFocalEvents <- function(potential.sc.foc.mcs.segs.vec, pj.mcs.df, col.to.update, event.type) {

	arm.vec <- unique(pj.mcs.df$arm)

	for (arm in arm.vec ) {

		arm.pj.mcs.df <- pj.mcs.df[pj.mcs.df$arm %in% arm,]
		arm.pj.mcs.df$old_focal_gain_definition <- arm.pj.mcs.df[, col.to.update]
		arm.pj.mcs.df[, col.to.update] <- "none"
		#So now we need to create the column that tell sus whether we actually can look for runs:
		arm.pj.mcs.df$potential_new_focal <- arm.pj.mcs.df$samp_seg_name %in% potential.sc.foc.mcs.segs.vec
														
		arm.unclassified.focal.mat <- matrix(FALSE, nrow=num.samples, ncol=length(unique(arm.pj.mcs.df$seg_name)))
		colnames(arm.unclassified.focal.mat) <- unique(arm.pj.mcs.df$seg_name)
		rownames(arm.unclassified.focal.mat) <- unique(arm.pj.mcs.df$sample)
		#arm.unclassified.focal.mat.allele <- matrix(FALSE, nrow=length(unique(arm.pj.mcs.df$sample)), ncol=length(unique(arm.pj.mcs.df$seg_name)))
		#colnames(arm.unclassified.focal.mat.allele)  <- colnames(arm.unclassified.focal.mat) <- unique(arm.pj.mcs.df$seg_name)
		#rownames(arm.unclassified.focal.mat.allele)  <- rownames(arm.unclassified.focal.mat) <- unique(arm.pj.mcs.df$sample)
		for(samp in unique(arm.pj.mcs.df$sample)) {
			for(seg in unique(arm.pj.mcs.df$seg_name)){
				arm.unclassified.focal.mat[samp, seg] <- arm.pj.mcs.df[arm.pj.mcs.df$sample == samp & arm.pj.mcs.df$seg_name == seg,]$potential_new_focal
				#arm.unclassified.focal.mat.allele[samp, seg] <- arm.pj.mcs.df[arm.pj.mcs.df$sample == samp & arm.pj.mcs.df$seg_name == seg,]$dominant_allele
			}#end for(seg in unique(arm.pj.mcs.df$seg_name)){	
		}#end for(samp in unique(arm.pj.mcs.df$sample)) {

		unclassified.focals.remain <- any(arm.unclassified.focal.mat==TRUE)
		while (unclassified.focals.remain) {
			print(arm.unclassified.focal.mat)
			max.seg.name <- names(which.max(colSums(arm.unclassified.focal.mat)))
			samps.w.max.seg <-rownames(arm.unclassified.focal.mat[arm.unclassified.focal.mat[,max.seg.name]==TRUE,max.seg.name, drop = FALSE])
			joined.seg.size.vec <- arm.pj.mcs.df[arm.pj.mcs.df$sample %in% samps.w.max.seg & arm.pj.mcs.df$seg_name %in% max.seg.name,]$MPhase_joined_seg_width
			print("joined.seg.size.vec")
			print(joined.seg.size.vec)
			names(joined.seg.size.vec) <- arm.pj.mcs.df[arm.pj.mcs.df$sample %in% samps.w.max.seg & arm.pj.mcs.df$seg_name %in% max.seg.name,]$uniq_joined_seg_name
			joined.focal.seg.matching.back.to <- names(which.max(joined.seg.size.vec))
			print(joined.seg.size.vec)
			#So now we have the segment we're matching back to, we need its start and end so we can give these segments this id
			chr.start.end.vec <- as.numeric(strsplit(strsplit(joined.focal.seg.matching.back.to, split="@")[[1]][2], split = ":")[[1]])
			#So now we can make the new ID
			new.focal.id <- paste0("focal_", event.type,"_id","@", chr.start.end.vec[1],":",chr.start.end.vec[2],":",chr.start.end.vec[3])
			#We can also rejig the matrix to eliminate this from a potential new focal.
			start.col.idx <- grep(chr.start.end.vec[2], colnames(arm.unclassified.focal.mat))
			end.col.idx <- grep(chr.start.end.vec[3], colnames(arm.unclassified.focal.mat))
			seg.names.to.update <- colnames(arm.unclassified.focal.mat)[start.col.idx:end.col.idx]
			arm.unclassified.focal.mat[,start.col.idx:end.col.idx] <- FALSE
			#Update the original matrix with the new focal ids
			arpj.mcs.df[arm.pj.mcs.df$sample%in%samps.w.max.seg & arm.pj.mcs.df$seg_name%in%seg.names.to.update, col.to.update] <- new.focal.id
			unclassified.focals.remain <- any(arm.unclassified.focal.mat==TRUE)
			print(new.focal.id)
			print("##########3")

		}#end while (unclassified.focals.remain) {

			#Update the dataframe with this arm's results
			pj.mcs.df[pj.mcs.df$arm %in% arm, col.to.update] <- arm.pj.mcs.df[, col.to.update]

	}#end for (arm in arm.vec ) {

	#return updated dataframe
	return(pj.mcs.df)

}#end markClonalFocalEvents <- function() {


























arm.unclassified.focal.mat <- as.logical(arm.unclassified.focal.mat)
new.arm.unclassified.focal.mat <- arm.unclassified.focal.mat

new.arm.unclassified.focal.mat[new.arm.unclassified.focal.mat==TRUE] <- "gain"
new.arm.unclassified.focal.mat[new.arm.unclassified.focal.mat=="FALSE"] <- "none"
	
test2 <- matrix(paste(new.arm.unclassified.focal.mat,arm.unclassified.focal.mat.allele,sep="_"),nrow(arm.unclassified.focal.mat),ncol(arm.unclassified.focal.mat))
test3 <- test2 == "gain_cpnA"
colSums(test3)
arm.pj.mcs.df[arm.pj.mcs.df$seg_name%in%"22:38035915:39908248", c(focal.col, "cpn_event_vs_ploidy", "pj_focal_by_size_and_prop", "potential_new_focal")]



#So now we want to essentially repeat what we did for the clonal stuff here




#Olf clonal focal definitions and sizes:
# mcs.seg.name.vec <- unique(pj.mcs.df$seg_name)

# for (mcs.seg in mcs.seg.name.vec) {

# 	if (all(pj.mcs.df[pj.mcs.df$seg_name %in% mcs.seg,]$pj_focal_by_size_and_prop %in% "focal")) {

# 		if(all(pj.mcs.df[pj.mcs.df$seg_name %in% mcs.seg, ]$cpn_event_vs_ploidy %in% c("gain", "amp"))) {

# 			pj.mcs.df[pj.mcs.df$seg_name %in% mcs.seg ,]$cl_focal_mcs_seg_cpn_vs_ploidy <- "gain"

# 		}

# 		if(all(pj.mcs.df[pj.mcs.df$seg_name %in% mcs.seg, ]$cpn_event_vs_ploidy %in% c("loss", "deep_loss"))) {

# 			pj.mcs.df[pj.mcs.df$seg_name %in% mcs.seg ,]$cl_focal_mcs_seg_cpn_vs_ploidy <- "loss"

# 		}

# 		if(all(pj.mcs.df[pj.mcs.df$seg_name %in% mcs.seg, ]$within_samp_cl_LOH_present %in% TRUE)) {

# 			pj.mcs.df[pj.mcs.df$seg_name %in% mcs.seg ,]$cl_focal_mcs_seg_LOH <- "LOH"

# 		}

# 		if(all(pj.mcs.df[pj.mcs.df$seg_name %in% mcs.seg, ]$within_samp_cl_LOH_present %in% TRUE |
# 			pj.mcs.df[pj.mcs.df$seg_name %in% mcs.seg, ]$cpn_event_vs_ploidy %in% c("loss", "deep_loss"))) {

# 			pj.mcs.df[pj.mcs.df$seg_name %in% mcs.seg ,]$cl_focal_mcs_seg_LOH_or_loss <- "loss_and_any_LOH"

# 		}

# 	}# if (all(pj.mcs.df[pj.mcs.df$seg_name %in% mcs.seg,]$pj_focal_by_size_and_prop %in% "focal")) {

# }#end for (mcs.seg in mcs.seg.name.vec) {

# #Now we have identified clonal focal events in their various categories, we need to give them a unique ID that also has the size
# addUniqueIDsToClonalFocalSegVec <- function (cl.focal.cpn.vs.ploidy.vec, arm.vec) {

# 	require("stringr")
# 	cl.focal.cpn.vs.ploidy.rle <- rle(paste0(cl.focal.cpn.vs.ploidy.vec, "_", arm.vec))
# 	cl.focal.cpn.vs.ploidy.end.vec <- cumsum(cl.focal.cpn.vs.ploidy.rle$lengths)
# 	cl.focal.cpn.vs.ploidy.start.vec <- c(1, (cl.focal.cpn.vs.ploidy.end.vec[-length(cl.focal.cpn.vs.ploidy.end.vec)]+1))
# 	cl.focal.cpn.vs.ploidy.id.vec <- c(1:length(cl.focal.cpn.vs.ploidy.end.vec))
# 	cl.focal.cpn.vs.ploidy.id.vec <- str_pad(cl.focal.cpn.vs.ploidy.id.vec, 6, pad = "0")


# 	for (j in 1:length(cl.focal.cpn.vs.ploidy.start.vec)) {

# 		cl.focal.cpn.vs.ploidy.vec[cl.focal.cpn.vs.ploidy.start.vec[j]:cl.focal.cpn.vs.ploidy.end.vec[j]] <- paste0(cl.focal.cpn.vs.ploidy.vec[cl.focal.cpn.vs.ploidy.start.vec[j]:cl.focal.cpn.vs.ploidy.end.vec[j]], "_", cl.focal.cpn.vs.ploidy.id.vec[j])

# 	}
# 	#Make all the none_whatever_number back into nones:
# 	cl.focal.cpn.vs.ploidy.vec[grepl("none",cl.focal.cpn.vs.ploidy.vec)] <- "none"
# 	#So now return s
# 	return(cl.focal.cpn.vs.ploidy.vec)
# }

# sample.vec <- unique(pj.mcs.df$sample)
# samp <- sample.vec[1]
# arm.vec <- pj.mcs.df[pj.mcs.df$sample == samp,]$arm

# cl.focal.cpn.vs.ploidy.vec <- pj.mcs.df[pj.mcs.df$sample == samp,]$cl_focal_mcs_seg_cpn_vs_ploidy
# names(cl.focal.cpn.vs.ploidy.vec) <-  pj.mcs.df[pj.mcs.df$sample == samp,]$seg_name
# cl.focal.cpn.vs.ploidy.vec <- addUniqueIDsToClonalFocalSegVec(cl.focal.cpn.vs.ploidy.vec, arm.vec)

# cl.focal.mcs.seg.LOH.vec <- pj.mcs.df[pj.mcs.df$sample == samp,]$cl_focal_mcs_seg_LOH
# names(cl.focal.mcs.seg.LOH.vec) <-  pj.mcs.df[pj.mcs.df$sample == samp,]$seg_name
# cl.focal.mcs.seg.LOH.vec <- addUniqueIDsToClonalFocalSegVec(cl.focal.mcs.seg.LOH.vec, arm.vec)

# cl.focal.mcs.seg.LOH.or.loss.vec <- pj.mcs.df[pj.mcs.df$sample == samp,]$cl_focal_mcs_seg_LOH_or_loss
# names(cl.focal.mcs.seg.LOH.or.loss.vec) <-  pj.mcs.df[pj.mcs.df$sample == samp,]$seg_name
# cl.focal.mcs.seg.LOH.or.loss.vec <- addUniqueIDsToClonalFocalSegVec(cl.focal.mcs.seg.LOH.or.loss.vec, arm.vec)

# for (seg.name in names(cl.focal.cpn.vs.ploidy.vec)) {

# 	pj.mcs.df[pj.mcs.df$seg_name == seg.name,]$cl_focal_mcs_seg_cpn_vs_ploidy <- cl.focal.cpn.vs.ploidy.vec[seg.name]
# 	pj.mcs.df[pj.mcs.df$seg_name == seg.name,]$cl_focal_mcs_seg_LOH <- cl.focal.mcs.seg.LOH.vec[seg.name]
# 	pj.mcs.df[pj.mcs.df$seg_name == seg.name,]$cl_focal_mcs_seg_LOH_or_loss <- cl.focal.mcs.seg.LOH.or.loss.vec[seg.name]

# }

# #We need the column:
# #TODO:
# getSizeOfClonalFocals <- function(pj.mcs.df, focal.col.name.to.modify) {

# 	#If it's not none then we need to calculate the size of the event, for now we'll use the size of the run to be size of the event since it's clonal.....
# 	#Maybe revisit at a later date.
# 	new.pj.mcs.df <- pj.mcs.df
# 	#Get unique ids:
# 	uniq.cl.focal.id.vec <- unique(new.pj.mcs.df[, focal.col.name.to.modify])
# 	uniq.cl.focal.id.vec <- uniq.cl.focal.id.vec[uniq.cl.focal.id.vec != "none" ] 

# 	if( length(uniq.cl.focal.id.vec) == 0 ){

# 		return(new.pj.mcs.df[, focal.col.name.to.modify])

# 	} else {

# 		for (orig.uniq.id in uniq.cl.focal.id.vec) {

# 			focal.run.chr <- new.pj.mcs.df[new.pj.mcs.df[, focal.col.name.to.modify] %in% orig.uniq.id,]$chrom[1]
# 			focal.run.start.pos <- min(new.pj.mcs.df[new.pj.mcs.df[, focal.col.name.to.modify] %in% orig.uniq.id,]$startpos)
# 			focal.run.end.pos <- max(new.pj.mcs.df[new.pj.mcs.df[, focal.col.name.to.modify] %in% orig.uniq.id,]$endpos)
# 			new.focal.id <- paste0("focal_", strsplit(orig.uniq.id,"_")[[1]][1], "_id@", focal.run.chr, ":", focal.run.start.pos, ":", focal.run.end.pos)
# 			idxs.to.change <- new.pj.mcs.df[, focal.col.name.to.modify] %in% orig.uniq.id

# 			new.pj.mcs.df[idxs.to.change, focal.col.name.to.modify] <- new.focal.id

# 		}#for (orig.id in uniq.cl.focal.ids) {

# 		return(new.pj.mcs.df[, focal.col.name.to.modify])

# 	}#end if( length(uniq.cl.focal.ids) == 0 ){

# }#getMeanSizeAndStandardDevFocal




	cl.focal.cpn.vs.ploidy.new.id.vec <- getSizeOfClonalFocals(pj.mcs.df = pj.mcs.df, focal.col.name.to.modify = "cl_focal_mcs_seg_cpn_vs_ploidy")
	pj.mcs.df$cl_focal_mcs_seg_cpn_vs_ploidy <- cl.focal.cpn.vs.ploidy.new.id.vec
	cl.focal.mcs.seg.LOH.new.id.vec <- getSizeOfClonalFocals(pj.mcs.df = pj.mcs.df, focal.col.name.to.modify = "cl_focal_mcs_seg_LOH")
	pj.mcs.df$cl_focal_mcs_seg_LOH <- cl.focal.mcs.seg.LOH.new.id.vec
	cl.focal.mcs.seg.LOH.or.loss.new.id.vec <- getSizeOfClonalFocals(pj.mcs.df = pj.mcs.df, focal.col.name.to.modify = "cl_focal_mcs_seg_LOH_or_loss")
	pj.mcs.df$cl_focal_mcs_seg_LOH_or_loss <- cl.focal.mcs.seg.LOH.or.loss.new.id.vec




#Old subclonal stuff

# #So the function works like this:
# #Exclude any clonal focal segments of the same type from the vector
# #Make a dataframe and then reduced genomic ranges of the remaining potentially qualifying focal events
# #Overlap the
# #Make sure to check both dataframes have any data in them before calling function.
# getSubclonalFocalEvents <- function(possible.samp.segs.vec,
# 									mcs.sc.cpn.event.to.map.df,
# 									event.type) {


# 	possible.focal.segs.alone.vec <- sapply(possible.samp.segs.vec, function(x) strsplit(x, split="@")[[1]][2])
# 	possible.focal.segs.alone.df <- data.frame(chrom = sapply(possible.focal.segs.alone.vec, function(x) strsplit(x, split=":")[[1]][1]),
# 										 start = as.numeric(sapply(possible.focal.segs.alone.vec, function(x) strsplit(x, split=":")[[1]][2])),
# 										 end = as.numeric(sapply(possible.focal.segs.alone.vec, function(x) strsplit(x, split=":")[[1]][3])),
# 										stringsAsFactors=FALSE)
# 	possible.focal.segs.alone.grs <- makeGRangesFromDataFrame(possible.focal.segs.alone.df,
# 		                         keep.extra.columns=FALSE,
# 		                         ignore.strand=TRUE,
# 		                         seqinfo=NULL,
# 		                         seqnames.field=c("chrom"),
# 		                         start.field=c("start"),
# 		                         end.field=c("end"),
# 		                         starts.in.df.are.0based=FALSE)
# 	reduced.poss.f.sc.gr <- reduce(possible.focal.segs.alone.grs)
# 	reduced.poss.f.sc.gr$focal_event_id <- paste0("focal_", event.type, "_id@", as.character(seqnames(reduced.poss.f.sc.gr)), ":", start(reduced.poss.f.sc.gr), ":", end(reduced.poss.f.sc.gr))

# 	mcs.sc.cpn.event.to.map.gr <- makeGRangesFromDataFrame(mcs.sc.cpn.event.to.map.df,
# 		                         keep.extra.columns=TRUE,
# 		                         ignore.strand=TRUE,
# 		                         seqinfo=NULL,
# 		                         seqnames.field=c("chrom", "chr"),
# 		                         start.field=c("startpos"),
# 		                         end.field=c("endpos"),
# 		                         starts.in.df.are.0based=FALSE)

# 	print(mcs.sc.cpn.event.to.map.gr)
# 	print(reduced.poss.f.sc.gr)
# 	overlapping.f.sc.grs <- suppressWarnings(getOverlapGRanges(mcs.sc.cpn.event.to.map.gr, reduced.poss.f.sc.gr))
# 	print(overlapping.f.sc.grs)
# 	mcs.sc.cpn.event.to.map.df$focal_event_id <- "none"
# 	mcs.sc.cpn.event.to.map.df[overlapping.f.sc.grs$query_hit_idx, ]$focal_event_id <- overlapping.f.sc.grs$focal_event_id
# 	print(mcs.sc.cpn.event.to.map.df)
# 	return(mcs.sc.cpn.event.to.map.df)

# }#end getSubclonalFocalEvents