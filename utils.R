

##OVERLOADING FUNCTION FOR NOW... 
#Refactor to pass sample.ph.hz.ascat.df.list as a variable
##Params:
##
##Returns:
##
getPhasedASCATRawCopyNumberAB <- function(ref.ss.lr.baf.df, sample, sample.ph.hz.ascat.df.list) {
  ph.raw.cpn.df <- sample.ph.hz.ascat.df.list[[sample]]
  return(ph.raw.cpn.df[ph.raw.cpn.df$snp_pos %in% paste0(ref.ss.lr.baf.df$chrom, "_", ref.ss.lr.baf.df$pos), c("ph_ascat_nAraw","ph_ascat_nBraw")])

}#end calculateSNPCopyNumberAB <- function(

##
##Params:
##
##Returns:
##
calculateSNPsRawCopyNumberAB <- function(logr.vec, baf.vec, tumor.ploidy, tumor.purity, gamma.val) {

  cpnA <- (tumor.purity - 1 + baf.vec * 2 ^ (logr.vec / gamma.val) * ((1 - tumor.purity) * 2 + tumor.purity * tumor.ploidy)) / tumor.purity
  cpnB <- (tumor.purity - 1 - (baf.vec - 1) * 2 ^ (logr.vec / gamma.val) * ((1 - tumor.purity) * 2 + tumor.purity * tumor.ploidy)) / tumor.purity

  return(cbind(cpnA, cpnB))

}#end calculateSNPCopyNumberAB <- function(



wilcoxTestWithErrorHandling <- function(x, y = NULL,
            alternative = c("two.sided", "less", "greater"),
            mu = 0, paired = FALSE, exact = FALSE, correct = TRUE,
            conf.int = FALSE, conf.level = 0.95, errors.not.permitted="")
             {
    test.results <- tryCatch({

      test.results <- wilcox.test(x = x,
                                  y = y,
                                  alternative = alternative,
                                  mu = mu,
                                  paired = paired,
                                  exact = exact,
                                  correct = correct,
                                  conf.int = conf.int,
                                  conf.level = conf.level)

      if ( is.na(test.results[["p.value"]]) | is.nan(test.results[["p.value"]])) {

        test.results[["stat"]] <- 1
        test.results[["conf.int"]] <- c(NA, NA)
        test.results[["p.value"]] <- 1
        test.results[["p.val"]] <- 1

      }
      
      return(test.results)

      },
       error = function(err) {
          message(err)
          write(paste(err$message, "\n"), stderr())
          write("\n", stderr())
          if(grepl("NULL", err$message) & "NULL" %in% errors.not.permitted){
            stop("t.test ERROR with a NULL")
          } else if (grepl("enough", err$message) & "enough" %in% errors.not.permitted){
            stop("t.test ERROR with not enough")
          } else if (grepl("constant", err$message) & "constant" %in% errors.not.permitted){
            stop("t.test ERROR with constant data")
          } else if(grepl("enough|constant|NULL", err$message) & "other" %in% errors.not.permitted){
            stop(err$message)
          }
          test.results <- list(stat = 1,
                               conf.int = c(NA, NA),
                               p.val = 1)
          return(test.results)
        })

        if (is.na(test.results$p.val)) {

            test.results$p.val <- 1
        }

        return(test.results)
}#end wilcoxTestWithErrorHandling
  



#
tTestWithErrorHandling <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
                                   mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95, errors.not.permitted=""){
  
  test.results <- tryCatch({
    test.results <- t.test(x = x, y = y, alternative = alternative, mu = mu, paired = paired, var.equal = var.equal, conf.level = conf.level)

    if ( is.na(test.results[["p.value"]]) | is.nan(test.results[["p.value"]])) {

      test.results[["stat"]] <- 1
      test.results[["conf.int"]] <- c(NA, NA)
      test.results[["p.value"]] <- 1
      test.results[["p.val"]] <- 1

    }
    return(test.results)

    },
    error = function(err) {
      message(err)
      write(paste(err$message, "\n"), stderr())
      write("\n", stderr())
      if(grepl("NULL", err$message) & "NULL"%in% errors.not.permitted){
        stop("t.test ERROR with a NULL")
      } else if (grepl("enough", err$message) & "enough"%in% errors.not.permitted){
        stop("t.test ERROR with not enough")
      } else if (grepl("constant", err$message) & "constant" %in% errors.not.permitted){
        stop("t.test ERROR with constant data")
      } else if(grepl("enough|constant|NULL", err$message) & "other" %in% errors.not.permitted){
        stop(err$message)
      }
      test.results <- list()
      test.results[["stat"]] <- 1
      test.results[["conf.int"]] <- c(NA, NA)
      test.results[["p.value"]] <- 1
      test.results[["p.val"]] <- 1
      return(test.results)
    })
  return(test.results)
}



checkSegGainOrLossSignificance <- function(logr.vec, samp.to.check, ploidy.est.vec, purity.est.vec, threshold, cpn.change.type = "gain", min.num.for.logr.test = 3) {

  min.thresh.raw.exp.samp       <-  ploidy.est.vec[samp.to.check] * (2 ^ (threshold))
  min.thresh.logr.exp.for.samp   <-  calcExpectedLogR(min.thresh.raw.exp.samp, 0, ploidy.est.vec[samp.to.check], purity.est.vec[samp.to.check])
  median.seg.logr                  <- median(logr.vec)
  
  if (length(logr.vec) >= min.num.for.logr.test) {

     gen.t.test.results        <- tTestWithErrorHandling(logr.vec, mu = min.thresh.logr.exp.for.samp)

    if (cpn.change.type == "gain") {

          specific.t.test.alt <- "greater"

    } else if (cpn.change.type == "loss") {

          specific.t.test.alt <- "less"

    }

    specific.t.test.result <- tTestWithErrorHandling(logr.vec, mu = min.thresh.logr.exp.for.samp, alternative = specific.t.test.alt)

    low.obs.logr.conf.int  <- specific.t.test.result$conf.int[1]
    high.obs.logr.conf.int <- specific.t.test.result$conf.int[2]
    gen.t.test.pval     <- gen.t.test.results$p.value
    alt.t.test.pval     <- specific.t.test.result$p.value
    
  }#end if (nrow(ss.baf.lr.baf.df) >= 3) {

  #Todo, refactor this magic number...
  if (length(logr.vec) < 3) {

    low.obs.logr.conf.int  <- NA
    high.obs.logr.conf.int <- NA
    gen.t.test.pval         <- NA
    alt.t.test.pval     <- NA

  }#end f (nrow(ss.baf.lr.baf.df) < 3) {
  
  out <- cbind(min.thresh.raw.exp.samp,
               min.thresh.logr.exp.for.samp,
               median.seg.logr,
               low.obs.logr.conf.int,
               high.obs.logr.conf.int,
               gen.t.test.pval,
               alt.t.test.pval)  
  return(out)
  
}#end checkSegGainOrLossSignificance <- function(ss.baf.lr.baf.df, samp.to.check, ploidy.est.vec, purity.est.vec, threshold) {



calcExpectedLogR <- function(exp.intA, exp.intB, tumor.ploidy, tumor.purity) {

  exp.logr <- log2((2 * (1 - tumor.purity) + tumor.purity * (exp.intA + exp.intB)) / (2 * (1 - tumor.purity) + (tumor.purity * tumor.ploidy)))

  return(exp.logr)

}#end calcExpectedLogR <- function(exp.intA, exp.intB, tumor.ploidy, tumor.purity) {


calcExpectedBAF <- function(exp.intA, exp.intB, tumor.purity) {

  exp.BAF <- (1 - tumor.purity + (tumor.purity * exp.intB)) / (2 - 2 * tumor.purity + tumor.purity * (exp.intA + exp.intB))

  return(exp.BAF)

}#end calcExpectedBAF <- function(exp.intA, exp.intB, tumor.purity)

# Convert "chr1"-style chromosome names to integers
# chr1 -> 1, ..., chrX -> 23, chrY -> 24
chrom2integer <- function(chrom) {
  return(as.integer(gsub("Y", "24",
                    gsub("X", "23",
                    gsub("chr", "", chrom, fixed = TRUE), fixed = TRUE), fixed = TRUE)))
}

