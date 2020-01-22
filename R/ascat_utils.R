# Convert the internal refphase snps object into an ascat.bc object for use in
# internal rerun of ASCAT
#
# > str(ascat.bc)
# List of 14
# $ Tumor_LogR          :'data.frame':   938946 obs. of  1 variable:
# ..$ A_LTX121_SU_T1-R1: num [1:938946] -0.0542 -0.2103 -0.2103 -0.4974 -0.3978 ...
# $ Tumor_BAF           :'data.frame':   938946 obs. of  1 variable:
# ..$ BAF: num [1:938946] 1 NA NA NA NA ...
# $ Tumor_LogR_segmented: num [1:938946, 1] -0.335 -0.335 -0.335 -0.335 -0.335 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:938946] "1" "2" "3" "4" ...
# .. ..$ : chr "A_LTX121_SU_T1-R1"
# $ Tumor_BAF_segmented :List of 1
# ..$ : num [1:33666, 1] 0.23 0.23 0.23 0.23 0.23 ...
# .. ..- attr(*, "dimnames")=List of 2
# .. .. ..$ : chr [1:33666] "41" "105" "123" "135" ...
# .. .. ..$ : NULL
# $ Germline_LogR       :'data.frame':   938946 obs. of  1 variable:
# ..$ A_LTX121_SU_T1-R1: num [1:938946] 0 NA NA NA NA NA NA NA NA NA ...
# $ Germline_BAF        :'data.frame':   938946 obs. of  1 variable:
# ..$ BAF: num [1:938946] 1.01 NA NA NA NA ...
# $ SNPpos              :'data.frame':   938946 obs. of  2 variables:
# ..$ chrom: chr [1:938946] "1" "1" "1" "1" ...
# ..$ pos  : int [1:938946] 762273 762496 762497 777237 777318 777335 777361 777456 783071 783147 ...
# $ ch                  :List of 23
# ..$ : int [1:92902] 1 2 3 4 5 6 7 8 9 10 ...
# ..$ : int [1:68365] 92903 92904 92905 92906 92907 92908 92909 92910 92911 92912 ...
# ..$ : int [1:53503] 161268 161269 161270 161271 161272 161273 161274 161275 161276 161277 ...
# ..$ : int [1:22578] 899048 899049 899050 899051 899052 899053 899054 899055 899056 899057 ...
# ..$ : int [1:17321] 921626 921627 921628 921629 921630 921631 921632 921633 921634 921635 ...
# $ chr                 :List of 28
# ..$ : int [1:50594] 1 2 3 4 5 6 7 8 9 10 ...
# ..$ : int [1:42308] 50595 50596 50597 50598 50599 50600 50601 50602 50603 50604 ...
# ..$ : int [1:68365] 92903 92904 92905 92906 92907 92908 92909 92910 92911 92912 ...
# ..$ : int [1:53503] 161268 161269 161270 161271 161272 161273 161274 161275 161276 161277 ...
# ..$ : int [1:17321] 921626 921627 921628 921629 921630 921631 921632 921633 921634 921635 ...
# $ chrs                : chr [1:23] "1" "2" "3" "4" ...
# $ samples             : chr "A_LTX121_SU_T1-R1"
# $ gender              : chr "XX"
# $ sexchromosomes      : chr [1:2] "X" "Y"
# $ failedarrays        : NULL
snps_to_ascat_bc <- function(snps, sample_name) {
  bafseg <- matrix(snps$BAFsegmented, ncol = 1, dimnames = list(as.character(seq_along(snps$BAFsegmented))))
  # Remove positions with NA in BAF
  bafseg <- bafseg[!is.na(snps$BAF), , drop = FALSE]

  ascat_bc <- list(
    Tumor_LogR = data.frame(sample = snps$logR),
    Tumor_BAF = data.frame(sample = snps$BAF),
    Tumor_LogR_segmented = matrix(snps$logRsegmented,
      ncol = 1,
      dimnames = list(as.character(seq_along(snps$logRsegmented)))
    ),
    Tumor_BAF_segmented = list(bafseg),
    Germline_LogR = data.frame(),
    Germline_BAF = data.frame(),
    SNPpos = data.frame(
      chrs = as.character(snps$chrom),
      pos = snps$pos
    ),
    ch = NA,
    chr = NA,
    chrs = as.character(unique(sort(snps$chrom))),
    samples = c(sample_name),
    gender = c("XX"),
    sexchromosomes = c("23", "24"),
    failedarrays = NULL
  )

  ascat_bc$ch <- lapply(ascat_bc$chrs, function(x) {
    which(ascat_bc$SNPpos$chrs == x)
  })
  ascat_bc$chr <- split_genome(ascat_bc$SNPpos) # From ascat.R

  # colnames(ascat_bc$Germline_BAF) <- "BAF"
  colnames(ascat_bc$Tumor_BAF) <- "BAF"
  colnames(ascat_bc$Tumor_LogR) <- sample_name

  ascat_bc
}
