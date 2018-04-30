#We need to write a script that will allow running of mphase from the command line
#For now it only needs to take the following input: file containing segmentation for all samples to be investigated with 
#!/usr/bin/env Rscript
# Copyright 2018-2013 Trevor L Davis <trevor.l.davis@gmail.com>
# Copyright 2008 Allen Day
#  
#  This file is free software: you may copy, redistribute and/or modify it  
#  under the terms of the GNU General Public License as published by the  
#  Free Software Foundation, either version 2 of the License, or (at your  
#  option) any later version.  
#  
#  This file is distributed in the hope that it will be useful, but  
#  WITHOUT ANY WARRANTY; without even the implied warranty of  
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
#  General Public License for more details.  
#  
#  You should have received a copy of the GNU General Public License  
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.  
# options(scipen=999);
# suppressPackageStartupMessages(library("optparse"))
# suppressPackageStartupMessages(library("stats"))
# suppressPackageStartupMessages(library(compiler)
# suppressPackageStartupMessages(library(parallel)
# suppressPackageStartupMessages(library(data.table)
# suppressPackageStartupMessages(library(mclust))
# suppressPackageStartupMessages(library(GenomicRanges))
# suppressPackageStartupMessages(library(Rsamtools))
# suppressPackageStartupMessages(library(gtools))
# suppressPackageStartupMessages(library(diptest))
# suppressPackageStartupMessages(library(gplots))


#What we need to run it
#Full info in SNP file
#Full info in seg file
#Directory to write out to
#we might want to add more options later.


# # specify our desired options in a list
# # by default OptionParser will add an help option equivalent to 
# # make_option(c("-h", "--help"), action="store_true", default=FALSE, 
# #               help="Show this help message and exit")
# option_list <- list( 
#     make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
#         help="Print extra output [default]"),
#     make_option(c("-q", "--quietly"), action="store_false", 
#         dest="verbose", help="Print little output"),
#     make_option(c("-s", "--segments"), type="integer", default=5, 
#         help="Number of random normals to generate [default %default]",
#         metavar="number"),
#     make_option("--generator", default="rnorm", 
#         help = "Function to generate random deviates [default \"%default\"]"),
#     make_option("--mean", default=0, 
#         help="Mean if generator == \"rnorm\" [default %default]"),
#     make_option("--sd", default=1, metavar="standard deviation",
#         help="Standard deviation if generator == \"rnorm\" [default %default]")
#     )
                                        
# # get command line options, if help option encountered print help and exit,
# # otherwise if options not found on command line then set defaults, 
# parser <- parse_args(OptionParser(option_list=option_list))


# # check date is provided
# if (!is.na(parser$date)) {
#   mydate <- as.Date(parser$date)
# }
# else {
#   stop("date parameter must be provided. See script usage (--help)")
# }


# # print some progress messages to stderr if "quietly" wasn't requested
# if ( parser$verbose ) { 
#     write("writing some verbose output to standard error...\n", stderr()) 
# }


# # do some operations based on user input
# if( parser$generator == "rnorm") {
#     cat(paste(rnorm(parser$count, mean=parser$mean, sd=parser$sd), collapse="\n"))
# } else {
#     cat(paste(do.call(parser$generator, list(parser$count)), collapse="\n"))
# }
# cat("\n")