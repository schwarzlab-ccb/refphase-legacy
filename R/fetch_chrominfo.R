#!/usr/bin/env Rscript
#
# Quick hack to put together a chromosome information table that can be used
# for mphase.R
#
# Needs:
# - chrom
# - chr_chrom
# - p_start
# - p_end
# - q_start
# - q_end
# - q_name
# - p_name
# - acrocentric

#library(optparse)
#library(digest)
#
#option_list <- list(
#  make_option(c("-n", "--genome"),
#    default = "hg19",
#    help = "Goldenpath genome name (default=hg19)"
#  ),
#  make_option(c("-d", "--seed"),
#    default = 531,
#    help = "Set the random seed at the beginning of the script"
#  )
#)
#args <- parse_args(OptionParser(option_list = option_list),
#  convert_hyphens_to_underscores = TRUE
#)
#conf <- args
#conf$id <- substr(digest(args), 1, 8)
#set.seed(args$seed)

fetch_chrominfo <- function(genome) {

  genome <- conf$genome
  chromosomes <- paste0("chr", c(1:22, "X", "Y"))
  is_acrocentric <- c("chr13", "chr14", "chr15", "chr21", "chr22", "chrY")
  arms <- c("p", "q")

  chrominfo_file <- paste0("http://hgdownload.cse.ucsc.edu/goldenPath/", genome, "/database/chromInfo.txt.gz")
  cytobands_file <- paste0("http://hgdownload.cse.ucsc.edu/goldenPath/", genome, "/database/cytoBand.txt.gz")

  read_tsv_url <- function(path, header = FALSE) {
    con <- gzcon(url(path))
    txt <- readLines(con)
    csv <- read.delim(textConnection(txt), header = header)
    return(csv)
  }

  chrom_lengths <- read_tsv_url(chrominfo_file)
  colnames(chrom_lengths) <- c("chrom", "length", "file")

  cytobands <- read_tsv_url(cytobands_file)
  colnames(cytobands) <- c("chrom", "start", "end", "cytoband", "stain")
  cytobands$arm <- substr(cytobands$cytoband, 0, 1)

  chrom_data <- chrom_lengths[chrom_lengths$chrom %in% chromosomes, ]
  row.names(chrom_data) <- chrom_data$chrom
  # Sort
  chrom_data <- chrom_data[chromosomes, ]
  chrom_data <- cbind(
    chrom_data,
    data.frame(
      p_start = NA,
      p_end = NA,
      p_size = NA,
      q_start = NA,
      q_end = NA,
      q_size = NA,
      p_name = NA,
      q_name = NA,
      acrocentric = NA
    )
  )

  for (chrom in chromosomes) {
    chrom_data[chrom, "acrocentric"] <- chrom %in% is_acrocentric
    for (arm in arms) {
      # cat("chrom", chrom, "arm", arm, "\n")
      arm_info <- cytobands[cytobands$chrom == chrom & cytobands$arm == arm, ]
      arm_start <- min(arm_info$start)
      arm_end <- max(arm_info$end)

      chrom_data[chrom, paste0(arm, "_start")] <- arm_start + 1
      chrom_data[chrom, paste0(arm, "_end")] <- arm_end
      chrom_data[chrom, paste0(arm, "_size")] <- arm_end - (arm_start + 1)
      chrom_data[chrom, paste0(arm, "_name")] <- paste0(gsub("chr", "", chrom), arm)
    }
  }

  row.names(chrom_data) <- gsub("chr", "", row.names(chrom_data))
  chrom_data$chr_chrom <- chrom_data$chrom
  chrom_data$chrom <- row.names(chrom_data)
  chrom_data <- chrom_data[, c(
    "chrom", "chr_chrom", "p_start", "p_end", "p_size",
    "q_start", "q_end", "q_size", "p_name", "q_name",
    "acrocentric"
  )]
  #write.table(chrom_data, file = paste0("ucsc_", genome, "_chrominfo.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

  return(chrom_data)
}
