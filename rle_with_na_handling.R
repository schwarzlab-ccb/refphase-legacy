 # I'm not sure if this should be posted here - so apologies if this is the case.
 # I've found a problem while using rle and am proposing a solution to the issue.

 # Description:
 # I ran into a niggle with rle today when working with vectors with NA values
 # (using R 2.31.0 on Windows 7 x64). It transpires that a run of NA values
 # is not encoded in the same way as a run of other values. See the following
 # example as an illustration:

 # Example:
 # The example
 #          rv<-c(1,1,NA,NA,3,3,3);rle(rv)
 # Returns
 #          Run Length Encoding
 #            lengths: int [1:4] 2 1 1 3
 #            values : num [1:4] 1 NA NA 3
 # not
 #          Run Length Encoding
 #            lengths: int [1:3] 2 2 3
 #            values : num [1:3] 1 NA 3
 # as I expected. This caused my code to fail later (unsurprising).

 # Analysis:
 # The problem stems from the test
 #           y<- x[-1L] != x[-n]
 # in line 7 of the rle function body. In this test, NA values return logical NA
 # values, not TRUE/FALSE (again, unsurprising).

 # Resolution:
 # I modified the rle function code as included below. As far as I tested, this
 # modification appears safe. The convoluted construction of naMaskVal
 # should guarantee that the NA masking value is always different from
 # any value in the vector and should be safe regardless of the input vector
 # form (a raw vector is not handled since the NA values do not apply here).

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

#  Conclusion:
#  I think that the proposed code modification is an improvement on the existing
#  implementation of rle. Is it impertinent to suggest this R-modification to the
#  gurus at R?

# > Best wishes (in flame-war trepidation),