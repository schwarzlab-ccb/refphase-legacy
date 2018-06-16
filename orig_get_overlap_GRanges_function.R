## Utility: identify overlaps in range data sets, such as annotation or alignment positions defined
## by two IRanges/GRanges objects.  

## Overlap types
	## olup: startup & endin
	## Q --------------
	## S       -------------

	## oldown: startin & enddown
	## Q       -------------
	## S --------------

	## inside: startin & endin 
	## Q     -----
	## S --------------

	## contained: startup & enddown
	## Q --------------
	## S     -----
olRanges <- function(query, subject, output="gr", ...) {
        require(GenomicRanges); require(IRanges)
        
        ## Input check
        if(!((class(query)=="GRanges" & class(subject)=="GRanges") | (class(query)=="IRanges" & class(subject)=="IRanges"))) {
                stop("Query and subject need to be of same class, either GRanges or IRanges!")
        }
        
        ## Find overlapping ranges
        if(class(query)=="GRanges") {
        	seqlengths(query) <- rep(NA, length(seqlengths(query)))
        	seqlengths(subject) <- rep(NA, length(seqlengths(subject)))
        }
	olindex <- as.matrix(findOverlaps(query, subject, ...))

        query <- query[olindex[,1]]
        subject <- subject[olindex[,2]]
        olma <- cbind(Qstart=start(query), Qend=end(query), Sstart=start(subject), Send=end(subject))
	
        ## Pre-queries for overlaps
	startup <- olma[,"Sstart"] < olma[,"Qstart"]
	enddown <- olma[,"Send"] > olma[,"Qend"]
	startin <- olma[,"Sstart"] >= olma[,"Qstart"] & olma[,"Sstart"] <= olma[,"Qend"]
	endin <- olma[,"Send"] >= olma[,"Qstart"] & olma[,"Send"] <=  olma[,"Qend"]

	## Overlap types
	olup <- startup & endin
	oldown <- startin & enddown
	inside <- startin & endin 
	contained <- startup & enddown
	
        ## Overlap types in one vector
	OLtype <- rep("", length(olma[,"Qstart"]))
	OLtype[olup] <- "olup"
	OLtype[oldown] <- "oldown"
	OLtype[inside] <- "inside" 
	OLtype[contained] <- "contained"
	
        ## Overlap positions
	OLstart <- rep(0, length(olma[,"Qstart"]))
	OLend <- rep(0, length(olma[,"Qstart"]))
	OLstart[olup] <- olma[,"Qstart"][olup]
	OLend[olup] <- olma[,"Send"][olup]
	OLstart[oldown] <- olma[,"Sstart"][oldown]
	OLend[oldown] <- olma[,"Qend"][oldown]
	OLstart[inside] <- olma[,"Sstart"][inside]
	OLend[inside] <- olma[,"Send"][inside]
	OLstart[contained] <- olma[,"Qstart"][contained]
	OLend[contained] <- olma[,"Qend"][contained]

	## Absolute and relative length of overlaps
	OLlength <- (OLend - OLstart) + 1
        OLpercQ <- OLlength/width(query)*100
        OLpercS <- OLlength/width(subject)*100
	
	## Output type
        oldf <- data.frame(Qindex=olindex[,1], Sindex=olindex[,2], olma, OLstart, OLend, OLlength, OLpercQ, OLpercS, OLtype)
        if(class(query) == "GRanges") {
		oldf <- cbind(space=as.character(seqnames(query)), oldf)
	}
        if(output=="df") {
                return(cbind(oldf, as.data.frame(mcols(subject))))
        }
        if(output=="gr") {
                if(class(query)=="GRanges") {
                        elementMetadata(query) <- cbind(as.data.frame(elementMetadata(query)), oldf, mcols(subject))
                }
                if(class(query)=="IRanges") {#Todo add the metadata in th eiranges 
                        query <- GRanges(seqnames = Rle(rep("dummy", length(query))), ranges = IRanges(start=oldf[,"Qstart"], end=oldf[,"Qend"]), strand = Rle(strand(rep("+", length(query)))), oldf)  
                }
                return(query)
        }
}

