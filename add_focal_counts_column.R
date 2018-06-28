


## Overlap types
## upstream_overlap: startup & endin
## Q --------------
## S       -------------
##
## downstream_overlap: startin & enddown
## Q       -------------
## S --------------
##
## query_contained_in_subject: startin & endin 
## Q     -----
## S --------------
##
## subject_inside_query: startup & enddown
## Q --------------
## S     -----
#Params:
# query: one or more granges/iranges that you wish to compare to the subject
# subject: one or more grangaes/ranges that you wish to be compared with the query
# output: specifies the output type either "gr" for genomic range or "df" for dataframe 
#Returns a dataframe/grange with following columns/fields:
# seqnames: this is just Grange name for chrom of the input GRange
# start: start of the input GRange
# end: end of the input GRange
# width: width of the input GRange
# strand: strand (mostly useless)
# chrom: seqname of input renamed as chrom
# query_hit_idx: idx of the GRange from input query GRange for this row
# subject_hit_idx: matching idx of the overlapping output Subject segment found for this row
# query_start: genomic start of query
# query_end: genomic end of query
# subject_start: genomic start of matching segment
# subject_end: genomic end of matching segment
# overlap_start: genomic start of overlap
# overlap_end: genomic end of overlap
# overlap_length: length of overlap
# overlap_prop_query: proportion of query Grange covered by overlap
# overlap_prop_subject: proportion of subject/result grange covered by overlap
# overlap_type: what kind of overlap it is (edited)
# 
#Also returns all the metadata in the subject, ie original columns holding non-genomic position information
#Todo add widths for everything
getOverlapGRanges <- function(query, subject, output="gr", ...) {
        require(GenomicRanges); require(IRanges)
        
        ## Input check
        if(!((class(query) == "GRanges" & class(subject) == "GRanges") | (class(query) == "IRanges" & class(subject) == "IRanges"))) {
                stop("Query and subject need to be of same class, either GRanges or IRanges!")
        }
        
        ## Find overlapping ranges
        if(class(query) == "GRanges") {
          seqlengths(query) <- rep(NA, length(seqlengths(query)))
          seqlengths(subject) <- rep(NA, length(seqlengths(subject)))
        }
  overlap.index.mat <- as.matrix(findOverlaps(query, subject, ...))

        query <- query[overlap.index.mat[, 1]]
        subject <- subject[overlap.index.mat[, 2]]
        overlap.mat <- cbind(query_start = start(query), query_end = end(query), subject_start = start(subject), subject_end = end(subject))
  
        ## Pre-queries for overlaps
  startup <- overlap.mat[,"subject_start"] < overlap.mat[,"query_start"]
  enddown <- overlap.mat[,"subject_end"] > overlap.mat[,"query_end"]
  startin <- overlap.mat[,"subject_start"] >= overlap.mat[,"query_start"] & overlap.mat[,"subject_start"] <= overlap.mat[, "query_end"]
  endin <- overlap.mat[,"subject_end"] >= overlap.mat[,"query_start"] & overlap.mat[,"subject_end"] <=  overlap.mat[, "query_end"]

  ## Overlap types
  olup <- startup & endin
  oldown <- startin & enddown
  inside <- startin & endin 
  contained <- startup & enddown
  
        ## Overlap types in one vector
  overlap_type <- rep("", length(overlap.mat[,"query_start"]))
  overlap_type[olup] <- "upstream_overlap"
  overlap_type[oldown] <- "downstream_overlap"
  overlap_type[inside] <- "subject_inside_query" 
  overlap_type[contained] <- "query_contained_in_subject"
  
        ## Overlap positions
  overlap_start <- rep(0, length(overlap.mat[, "query_start"]))
  overlap_end <- rep(0, length(overlap.mat[, "query_start"]))
  overlap_start[olup] <- overlap.mat[, "query_start"][olup]
  overlap_end[olup] <- overlap.mat[, "subject_end"][olup]
  overlap_start[oldown] <- overlap.mat[, "subject_start"][oldown]
  overlap_end[oldown] <- overlap.mat[, "query_end"][oldown]
  overlap_start[inside] <- overlap.mat[, "subject_start"][inside]
  overlap_end[inside] <- overlap.mat[, "subject_end"][inside]
  overlap_start[contained] <- overlap.mat[, "query_start"][contained]
  overlap_end[contained] <- overlap.mat[, "query_end"][contained]

  ## Absolute and relative length of overlaps
  overlap_length <- (overlap_end - overlap_start) + 1
        overlap_prop_query <- overlap_length / width(query) * 100
        overlap_prop_subject <- overlap_length / width(subject) * 100
  
  ## Output type
        oldf <- data.frame(query_hit_idx=overlap.index.mat[,1], subject_hit_idx=overlap.index.mat[,2], overlap.mat, overlap_start, overlap_end, overlap_length, overlap_prop_query, overlap_prop_subject, overlap_type)
        if(class(query) == "GRanges") {
    oldf <- cbind(chrom=as.character(seqnames(query)), oldf)
  }
        if(output=="df") {
                return(cbind(oldf, as.data.frame(mcols(subject))))
        }
        if(output=="gr") {
                if(class(query)=="GRanges") {
                        elementMetadata(query) <- cbind(as.data.frame(elementMetadata(query)), oldf, mcols(subject))
                }
                if(class(query)=="IRanges") {#Todo add the metadata in th eiranges 
                        query <- GRanges(seqnames = Rle(rep("dummy", length(query))), ranges = IRanges(start=oldf[,"query_start"], end=oldf[,"query_end"]), strand = Rle(strand(rep("+", length(query)))), oldf)  
                }
                return(query)
        }
}


 rle_w_na_handling<-function (x)
 {
      if (!is.vector(x)&&  !is.list(x))
          stop("'x' must be an atomic vector")
      n<- length(x)
      if (n == 0L)
          return(structure(list(lengths = integer(), values = x),
              class = "rle"))

      #### BEGIN NEW SECTION PART 1 ####
      naRepFlag<-F
      if(any(is.na(x))){
          naRepFlag<-T
          IS_LOGIC<-ifelse(typeof(x)=="logical",T,F)

          if(typeof(x)=="logical"){
              x<-as.integer(x)
              naMaskVal<-2
          }else if(typeof(x)=="character"){
              naMaskVal<-paste(sample(c(letters,LETTERS,0:9),32,replace=T),collapse="")
          }else{
              naMaskVal<-max(0,abs(x[!is.infinite(x)]),na.rm=T)+1
          }	

          x[which(is.na(x))]<-naMaskVal
      }
      #### END NEW SECTION PART 1 ####

      y<- x[-1L] != x[-n]
      i<- c(which(y), n)

      #### BEGIN NEW SECTION PART 2 ####
      if(naRepFlag)
          x[which(x==naMaskVal)]<-NA

      if(IS_LOGIC)
          x<-as.logical(x)
      #### END NEW SECTION PART 2 ####

      structure(list(lengths = diff(c(0L, i)), values = x[i]),
          class = "rle")
 }



#This script is ignoring that there is (as of 26/06/2018) an issue with the orig_seg_is_arm_or_focal meaning that some events are left without arm or focal annotation and have none left.
fjs.df <- full.joint.segs.df
fjs.df$samp_seg_name <- paste0(fjs.df$sample,"@",fjs.df$seg_name)

# fjs.df$focal_cpn_event_id <- NA
# fjs.df$focal_LOH_id <- NA
# fjs.df$focal_LOH_or_loss_id <- NA
#Ok so now say we want to find out how many clonal and subclonal focal events there are, we need to:

#Make focal original id column:
fjs.df$orig_focal_samp_seg_name <- NA
fjs.df[ ! is.na(fjs.df$orig_seg_is_arm_or_focal),]$orig_focal_samp_seg_name <- paste0(fjs.df[ ! is.na(fjs.df$orig_seg_is_arm_or_focal),]$orig_seg_is_arm_or_focal, "_", fjs.df[ ! is.na(fjs.df$orig_seg_is_arm_or_focal),]$orig_samp_seg_name)
#We could remove the arm defined ones here? we're not interested them anyway?
fjs.df[grepl("arm", fjs.df$orig_focal_samp_seg_name),]$orig_focal_samp_seg_name <- NA
#fjs.df[ ! is.na(fjs.df$orig_seg_is_arm_or_focal),]$focal_samp_seg_name <- paste0(fjs.df[ ! is.na(fjs.df$orig_seg_is_arm_or_focal),]$orig_seg_is_arm_or_focal, "_", fjs.df


#Event is focal in all regions:
#Get all the unique seg names 
uniq.seg.name.vec <- unique(fjs.df$seg_name)
fjs.df$cl_focal_mcs_seg_cpn_vs_ploidy <- "none"
fjs.df$cl_focal_mcs_seg_LOH <- "none"
fjs.df$cl_focal_mcs_seg_LOH_or_loss <- "none"
for (i in 1:length(uniq.seg.name.vec)) {

	if (all(fjs.df[fjs.df$seg_name %in% uniq.seg.name.vec[i] ,]$orig_seg_is_arm_or_focal %in% "focal")) {

		#So now we also want to know if they all go in the same direction in terms of cpn vs ploidy
		if(all(fjs.df[fjs.df$seg_name %in% uniq.seg.name.vec[i] ,]$cpn_event_vs_ploidy %in% c("gain", "amp"))) {

			fjs.df[fjs.df$seg_name %in% uniq.seg.name.vec[i] ,]$cl_focal_mcs_seg_cpn_vs_ploidy <- "gain"

		}
		if(all(fjs.df[fjs.df$seg_name %in% uniq.seg.name.vec[i] ,]$cpn_event_vs_ploidy %in% c("loss", "deep_loss"))) {

			fjs.df[fjs.df$seg_name %in% uniq.seg.name.vec[i] ,]$cl_focal_mcs_seg_cpn_vs_ploidy <- "loss"

		}

		if(all(fjs.df[fjs.df$seg_name %in% uniq.seg.name.vec[i] ,]$within_samp_cl_LOH_present %in% TRUE  )) {

			fjs.df[fjs.df$seg_name %in% uniq.seg.name.vec[i] ,]$cl_focal_mcs_seg_LOH <- "LOH"

		}

		if(all(fjs.df[fjs.df$seg_name %in% uniq.seg.name.vec[i] ,]$within_samp_cl_LOH_present %in% TRUE | fjs.df[fjs.df$seg_name %in% uniq.seg.name.vec[i] ,]$cpn_event_vs_ploidy %in% c("loss", "deep_loss") )) {

			fjs.df[fjs.df$seg_name %in% uniq.seg.name.vec[i] ,]$cl_focal_mcs_seg_LOH_or_loss <- "loss_and_any_LOH"

		}

	}

}


#So now we want to get all these coulmns and give them a unique id to all the clonal ones, simple task would be to label all the runs in the columns?
test.vec <- c("none", "none", "loss", "loss", "loss", "none", "loss", "none")
#test.vec <-  rev(rep(6:10, 1:5))

#Ok so let's get the run lengths and starts/ends 


#So let's annotate these clonal ones... we sort of want to do it by sample




addUniqueIDsToClonalFocalSegVec <- function (cl.focal.cpn.vs.ploidy.vec) {

	require("stringr")
	cl.focal.cpn.vs.ploidy.rle <- rle(cl.focal.cpn.vs.ploidy.vec)
	cl.focal.cpn.vs.ploidy.end.vec <- cumsum(cl.focal.cpn.vs.ploidy.rle$lengths)
	cl.focal.cpn.vs.ploidy.start.vec <- c(1, (cl.focal.cpn.vs.ploidy.end.vec[-length(cl.focal.cpn.vs.ploidy.end.vec)]+1))
	cl.focal.cpn.vs.ploidy.id.vec <- c(1:length(cl.focal.cpn.vs.ploidy.end.vec))
	cl.focal.cpn.vs.ploidy.id.vec <- str_pad(cl.focal.cpn.vs.ploidy.id.vec, 6, pad = "0")


	for (j in 1:length(cl.focal.cpn.vs.ploidy.start.vec)) {

		cl.focal.cpn.vs.ploidy.vec[cl.focal.cpn.vs.ploidy.start.vec[j]:cl.focal.cpn.vs.ploidy.end.vec[j]] <- paste0(cl.focal.cpn.vs.ploidy.vec[cl.focal.cpn.vs.ploidy.start.vec[j]:cl.focal.cpn.vs.ploidy.end.vec[j]], "_", cl.focal.cpn.vs.ploidy.id.vec[j])

	}
	#Make all the none_whatever_number back into nones:

	cl.focal.cpn.vs.ploidy.vec[grepl("none",cl.focal.cpn.vs.ploidy.vec)] <- "none"
	#So now return s
	return(cl.focal.cpn.vs.ploidy.vec)
}


sample.vec <- unique(fjs.df$sample)
samp <- sample.vec[1]

cl.focal.cpn.vs.ploidy.vec <- fjs.df[fjs.df$sample == samp,]$cl_focal_mcs_seg_cpn_vs_ploidy
names(cl.focal.cpn.vs.ploidy.vec) <-  fjs.df[fjs.df$sample == samp,]$seg_name
cl.focal.cpn.vs.ploidy.vec <- addUniqueIDsToClonalFocalSegVec(cl.focal.cpn.vs.ploidy.vec)

cl.focal.mcs.seg.LOH.vec <- fjs.df[fjs.df$sample == samp,]$cl_focal_mcs_seg_LOH
names(cl.focal.mcs.seg.LOH.vec) <-  fjs.df[fjs.df$sample == samp,]$seg_name
cl.focal.mcs.seg.LOH.vec <- addUniqueIDsToClonalFocalSegVec(cl.focal.mcs.seg.LOH.vec)

cl.focal.mcs.seg.LOH.or.loss.vec <- fjs.df[fjs.df$sample == samp,]$cl_focal_mcs_seg_LOH_or_loss
names(cl.focal.mcs.seg.LOH.or.loss.vec) <-  fjs.df[fjs.df$sample == samp,]$seg_name
cl.focal.mcs.seg.LOH.or.loss.vec <- addUniqueIDsToClonalFocalSegVec(cl.focal.mcs.seg.LOH.or.loss.vec)

for (seg.name in names(cl.focal.cpn.vs.ploidy.vec)) {

	fjs.df[fjs.df$seg_name == seg.name,]$cl_focal_mcs_seg_cpn_vs_ploidy <- cl.focal.cpn.vs.ploidy.vec[seg.name]
	fjs.df[fjs.df$seg_name == seg.name,]$cl_focal_mcs_seg_LOH <- cl.focal.mcs.seg.LOH.vec[seg.name]
	fjs.df[fjs.df$seg_name == seg.name,]$cl_focal_mcs_seg_LOH_or_loss <- cl.focal.mcs.seg.LOH.or.loss.vec[seg.name]

}


#Make the subclonal ones:
fjs.df$sc_focal_mcs_seg_cpn_vs_ploidy <- "none"
fjs.df$sc_focal_mcs_seg_LOH <- "none"
fjs.df$sc_focal_mcs_seg_LOH_or_loss <- "none"


#So now we have all the names of the focal ones that are clonal we want to get the counts of all the events that that 
#are subclonal focal events
#So now get all the focal samp seg names
all.focal.orig.seg.names.vec <-  fjs.df[fjs.df$orig_seg_is_arm_or_focal %in% c("focal")  ,]$orig_samp_seg_name

#So get those focals that are  in the clonal loss or the clonal gain
cl.gain.and.loss.focal.orig.seg.names.vec <- fjs.df[grepl("gain|loss",fjs.df$cl_focal_mcs_seg_cpn_vs_ploidy) ,]$orig_samp_seg_name
#also get those that are in clonal LOH
cl.loh.focal.orig.seg.names.vec <- fjs.df[grepl("LOH",fjs.df$cl_focal_mcs_seg_LOH) ,]$orig_samp_seg_name
#We want to mark these segments?
#get those athat are loss or  LOH
cl.loh.and.or.loss.focal.orig.seg.names.vec <- fjs.df[grepl("loss_and_any_LOH", fjs.df$cl_focal_mcs_seg_LOH_or_loss) ,]$orig_samp_seg_name

#So now if we want to know the actual number of unique focals we need to minus the vector of interested from all the focal ids
#Then make granges from the remaining IDs
#Collapse and take union from these, count the length of the final GRanges object.

#Ok so how many subclonal Granges are there
possible.focal.sc.gains.losses.vs.ploidy.vec <- all.focal.orig.seg.names.vec[ ! all.focal.orig.seg.names.vec %in% cl.gain.and.loss.focal.orig.seg.names.vec]
#Make into dataframe to make a granges object


focal.sc.gains.losses.segs.alone.vec <- sapply(possible.focal.sc.gains.losses.vs.ploidy.vec, function(x) strsplit(x, split="@")[[1]][2])
poss.f.sc.gain.loss.df <- data.frame(chrom = sapply(focal.sc.gains.losses.segs.alone.vec, function(x) strsplit(x, split=":")[[1]][1]),
									 start = as.numeric(sapply(focal.sc.gains.losses.segs.alone.vec, function(x) strsplit(x, split=":")[[1]][2])),
									 end = as.numeric(sapply(focal.sc.gains.losses.segs.alone.vec, function(x) strsplit(x, split=":")[[1]][3])),
									stringsAsFactors=FALSE)

poss.f.sc.gain.loss.gr <- makeGRangesFromDataFrame(poss.f.sc.gain.loss.df,
	                         keep.extra.columns=FALSE,
	                         ignore.strand=FALSE,
	                         seqinfo=NULL,
	                         seqnames.field=c("seqnames", "seqname",
	                                          "chromosome", "chrom",
	                                          "chr", "chromosome_name",
	                                          "seqid"),
	                         start.field=c("Start_Position","start", "startpos"),
	                         end.field=c("End_Position", "end", "stop", "endpos"),
	                         strand.field="strand",
	                         starts.in.df.are.0based=FALSE)
reduced.poss.f.sc.gain.loss.gr <- reduce(poss.f.sc.gain.loss.gr)
reduced.poss.f.sc.gain.loss.gr$focal_gain_id <- paste0("focal_gain_id@", as.character(seqnames(reduced.poss.f.sc.gain.loss.gr)), ":", start(reduced.poss.f.sc.gain.loss.gr), ":", end(reduced.poss.f.sc.gain.loss.gr))
reduced.poss.f.sc.gain.loss.gr$focal_loss_id <- paste0("focal_loss_id@", as.character(seqnames(reduced.poss.f.sc.gain.loss.gr)), ":", start(reduced.poss.f.sc.gain.loss.gr), ":", end(reduced.poss.f.sc.gain.loss.gr))



#FOCAL GAINS
subclonal.focal.gains.df <- fjs.df[fjs.df$orig_seg_is_arm_or_focal %in% c("focal") & fjs.df$cpn_event_vs_ploidy %in% c("gain", "amp") & fjs.df$across_samp_cl_cpn_vs_ploidy_present %in% c("subclonal_events_present") , c("chr", "startpos", "endpos","orig_samp_seg_name", "samp_seg_name")]
subclonal.focal.gains.df$focal_gain_df_id <- "none"
subclonal.focal.gains.gr <- makeGRangesFromDataFrame(subclonal.focal.gains.df,
	                         keep.extra.columns=TRUE,
	                         ignore.strand=FALSE,
	                         seqinfo=NULL,
	                         seqnames.field=c("seqnames", "seqname",
	                                          "chromosome", "chrom",
	                                          "chr", "chromosome_name",
	                                          "seqid"),
	                         start.field=c("Start_Position","start", "startpos"),
	                         end.field=c("End_Position", "end", "stop", "endpos"),
	                         strand.field="strand",
	                         starts.in.df.are.0based=FALSE)


overlapping.fg.grs <- getOverlapGRanges(subclonal.focal.gains.gr, reduced.poss.f.sc.gain.loss.gr)
print(overlapping.fg.grs)
print(subclonal.focal.gains.df)
subclonal.focal.gains.df[overlapping.fg.grs$query_hit_idx, ]$focal_gain_df_id <- overlapping.fg.grs$focal_gain_id
fjs.df[match(subclonal.focal.gains.df$samp_seg_name, fjs.df$samp_seg_name),]$sc_focal_mcs_seg_cpn_vs_ploidy <- subclonal.focal.gains.df$focal_gain_df_id

#FOCAL LOSSES 
subclonal.focal.losses.df <- fjs.df[fjs.df$orig_seg_is_arm_or_focal %in% c("focal") & fjs.df$cpn_event_vs_ploidy %in% c("loss", "deep") & fjs.df$across_samp_cl_cpn_vs_ploidy_present %in% c("subclonal_events_present") , c("chr", "startpos", "endpos", "orig_samp_seg_name", "samp_seg_name")]
subclonal.focal.losses.df$focal_loss_df_id <- "none"

subclonal.focal.losses.gr <- makeGRangesFromDataFrame(subclonal.focal.losses.df,
	                         keep.extra.columns=TRUE,
	                         ignore.strand=FALSE,
	                         seqinfo=NULL,
	                         seqnames.field=c("seqnames", "seqname",
	                                          "chromosome", "chrom",
	                                          "chr", "chromosome_name",
	                                          "seqid"),
	                         start.field=c("Start_Position","start", "startpos"),
	                         end.field=c("End_Position", "end", "stop", "endpos"),
	                         strand.field="strand",
	                         starts.in.df.are.0based=FALSE)


overlapping.fl.grs <- getOverlapGRanges(subclonal.focal.losses.gr, reduced.poss.f.sc.gain.loss.gr)
subclonal.focal.losses.df[overlapping.fl.grs$query_hit_idx, ]$focal_loss_df_id <- overlapping.fl.grs$focal_loss_id
fjs.df[match(subclonal.focal.losses.df$samp_seg_name, fjs.df$samp_seg_name),]$sc_focal_mcs_seg_cpn_vs_ploidy <- subclonal.focal.losses.df$focal_loss_df_id
 


#FOCAL LOH

#Ok so how many subclonal Granges are there
possible.focal.loh.alone.vec <- all.focal.orig.seg.names.vec[ ! all.focal.orig.seg.names.vec %in% cl.loh.focal.orig.seg.names.vec ]
#Make into dataframe to make a granges object


possible.focal.loh.seg.alone.vec <- sapply(possible.focal.loh.alone.vec, function(x) strsplit(x, split="@")[[1]][2])
poss.f.sc.gain.loh.df <- data.frame(chrom = sapply(possible.focal.loh.seg.alone.vec, function(x) strsplit(x, split=":")[[1]][1]),
									 start = as.numeric(sapply(possible.focal.loh.seg.alone.vec, function(x) strsplit(x, split=":")[[1]][2])),
									 end = as.numeric(sapply(possible.focal.loh.seg.alone.vec, function(x) strsplit(x, split=":")[[1]][3])),
									stringsAsFactors=FALSE)
poss.f.sc.gain.loh.df$focal_loh_id <- "none"
poss.f.sc.gain.loh.gr <- makeGRangesFromDataFrame(poss.f.sc.gain.loh.df,
	                         keep.extra.columns=FALSE,
	                         ignore.strand=FALSE,
	                         seqinfo=NULL,
	                         seqnames.field=c("seqnames", "seqname",
	                                          "chromosome", "chrom",
	                                          "chr", "chromosome_name",
	                                          "seqid"),
	                         start.field=c("Start_Position","start", "startpos"),
	                         end.field=c("End_Position", "end", "stop", "endpos"),
	                         strand.field="strand",
	                         starts.in.df.are.0based=FALSE)
reduced.poss.f.sc.loh.gr <- reduce(poss.f.sc.gain.loh.gr)
reduced.poss.f.sc.loh.gr$focal_loh_id <- paste0("focal_LOH_id@", as.character(seqnames(reduced.poss.f.sc.loh.gr)), ":", start(reduced.poss.f.sc.loh.gr), ":", end(reduced.poss.f.sc.loh.gr))


subclonal.focal.loh.df <- fjs.df[fjs.df$orig_seg_is_arm_or_focal %in% c("focal") & fjs.df$within_samp_cl_LOH_present %in% TRUE & fjs.df$across_samp_cl_LOH %in% FALSE, c("chr", "startpos", "endpos","orig_samp_seg_name", "samp_seg_name")]
subclonal.focal.loh.df$focal_loh_df_id <- "none"
subclonal.focal.loh.gr <- makeGRangesFromDataFrame(subclonal.focal.loh.df,
	                         keep.extra.columns=TRUE,
	                         ignore.strand=FALSE,
	                         seqinfo=NULL,
	                         seqnames.field=c("seqnames", "seqname",
	                                          "chromosome", "chrom",
	                                          "chr", "chromosome_name",
	                                          "seqid"),
	                         start.field=c("Start_Position","start", "startpos"),
	                         end.field=c("End_Position", "end", "stop", "endpos"),
	                         strand.field="strand",
	                         starts.in.df.are.0based=FALSE)


overlapping.loh.grs <- getOverlapGRanges(subclonal.focal.loh.gr, reduced.poss.f.sc.loh.gr)
subclonal.focal.loh.df[overlapping.loh.grs$query_hit_idx, ]$focal_loh_df_id <- overlapping.loh.grs$focal_loh_id
fjs.df[match(subclonal.focal.loh.df$samp_seg_name, fjs.df$samp_seg_name),]$sc_focal_mcs_seg_LOH <- subclonal.focal.loh.df$focal_loh_df_id


#LOH and/or loss alone


#Ok so how many subclonal Granges are there
possible.focal.loss.or.loh.vec <-  all.focal.orig.seg.names.vec[ ! all.focal.orig.seg.names.vec %in% cl.loh.and.or.loss.focal.orig.seg.names.vec ]
#Make into dataframe to make a granges object


possible.focal.loh.or.loss.seg.alone.vec <- sapply(possible.focal.loss.or.loh.vec, function(x) strsplit(x, split="@")[[1]][2])
poss.f.sc.loss.or.loh.df <- data.frame(chrom = sapply(possible.focal.loh.or.loss.seg.alone.vec, function(x) strsplit(x, split=":")[[1]][1]),
									 start = as.numeric(sapply(possible.focal.loh.or.loss.seg.alone.vec, function(x) strsplit(x, split=":")[[1]][2])),
									 end = as.numeric(sapply(possible.focal.loh.or.loss.seg.alone.vec, function(x) strsplit(x, split=":")[[1]][3])),
									stringsAsFactors=FALSE)
poss.f.sc.loss.or.loh.df$focal_loss_or_loh_id <- "none"
poss.f.sc.loss.or.loh.gr <- makeGRangesFromDataFrame(poss.f.sc.loss.or.loh.df,
	                         keep.extra.columns=FALSE,
	                         ignore.strand=FALSE,
	                         seqinfo=NULL,
	                         seqnames.field=c("seqnames", "seqname",
	                                          "chromosome", "chrom",
	                                          "chr", "chromosome_name",
	                                          "seqid"),
	                         start.field=c("Start_Position","start", "startpos"),
	                         end.field=c("End_Position", "end", "stop", "endpos"),
	                         strand.field="strand",
	                         starts.in.df.are.0based=FALSE)
reduced.poss.f.sc.loss.or.loh.gr <- reduce(poss.f.sc.loss.or.loh.gr)
reduced.poss.f.sc.loss.or.loh.gr$focal_loss_or_loh_id <- paste0("focal_LOH_id@", as.character(seqnames(reduced.poss.f.sc.loss.or.loh.gr)), ":", start(reduced.poss.f.sc.loss.or.loh.gr), ":", end(reduced.poss.f.sc.loss.or.loh.gr))


subclonal.focal.loss.or.loh.df <- fjs.df[fjs.df$orig_seg_is_arm_or_focal %in% c("focal") 
										 & ((fjs.df$within_samp_cl_LOH_present %in% TRUE & fjs.df$across_samp_cl_LOH %in% FALSE) 
										 | (fjs.df$cpn_event_vs_ploidy %in% c("loss", "deep") & fjs.df$across_samp_cl_cpn_vs_ploidy_present %in% c("subclonal_events_present"))), c("chr", "startpos", "endpos","orig_samp_seg_name", "samp_seg_name")]

subclonal.focal.loss.or.loh.df$focal_loss_or_loh_df_id <- "none"
subclonal.focal.loss.or.loh.gr <- makeGRangesFromDataFrame(subclonal.focal.loss.or.loh.df,
	                         keep.extra.columns=TRUE,
	                         ignore.strand=FALSE,
	                         seqinfo=NULL,
	                         seqnames.field=c("seqnames", "seqname",
	                                          "chromosome", "chrom",
	                                          "chr", "chromosome_name",
	                                          "seqid"),
	                         start.field=c("Start_Position","start", "startpos"),
	                         end.field=c("End_Position", "end", "stop", "endpos"),
	                         strand.field="strand",
	                         starts.in.df.are.0based=FALSE)


overlapping.loss.or.loh.grs <- getOverlapGRanges(subclonal.focal.loss.or.loh.gr, reduced.poss.f.sc.loss.or.loh.gr)
subclonal.focal.loss.or.loh.df[overlapping.loss.or.loh.grs$query_hit_idx, ]$focal_loss_or_loh_df_id <- overlapping.loss.or.loh.grs$focal_loss_or_loh_id
fjs.df[match(subclonal.focal.loss.or.loh.df$samp_seg_name, fjs.df$samp_seg_name),]$sc_focal_mcs_seg_LOH_or_loss <- subclonal.focal.loss.or.loh.df$focal_loss_or_loh_df_id



#So in summary the 


