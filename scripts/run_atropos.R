#!/usr/bin/env Rscript

rm(list = ls())

cat("\n############################################################\n")
cat("Starting \n")
cat("############################################################\n\n")

## ------------------------------------------------------------------------

require("optparse")

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(c("-i", "--input_directory"), type="character", default=NULL, 
              help="Path of the input directory containing reads in their respective run sub-directories \n
              e.g., -i raw [contains raw/run1 and raw/run2]\n
              N.B.: sample name is extracted from .fastq.gz samples  before the first '_' e.g., XXX-XX \n
              sample names must be unambigious and unique e.g., sample-1 = sample-11, sample-01 != sample-10", metavar="character"),
  make_option(c("--atropos_bin"), type="character", default = "atropos", 
              help="atropos binary", metavar="character"),
   make_option(c("-r", "--primer_removed_directory"), type="character", default = "00_atropos_primer_removed", 
              help="Directory to store primer trimmed/demultiplexed samples", metavar="character"),
  make_option(c("-o", "--output_directory"), type="character", default= "dada2", 
              help="Name of the output directory", metavar="character"),
  make_option(c("--fwd_primer"), type="character", default = NULL, 
              help="Sequence of the Fwd primer to remove", metavar="character"),
  make_option(c("--rev_primer"), type="character", default = NULL, 
              help="Sequence of the Rev primer to remove", metavar="character"),
  make_option(c("-c", "--raw_file_pattern"), type="character", default= c("_R1_.fastq.gz","_R2_.fastq.gz"), 
              help="Pattern of raw files", metavar="character"),
  make_option(c("-T", "--threads"), type="character", default= NULL, 
              help="Number of threads", metavar="character"),
  make_option(c("-f", "--fun_dir"), type="character", default= NULL, 
              help="Directory containing the R functions", metavar="character")
); 
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser)

## ------------------------------------------------------------------------
# parse_args(opt_parser, args = c("--help"))
# 
# if(is.null(opt$input_directory)) {
#   print_help(opt_parser)
#   stop("You must give an input directory (-i)")
# }
# 
# if(is.null(opt$atropos_binary)) {
#   print_help(opt_parser)
#   stop("You must provide path to atropos program (-a)")
# }

source(opt$fun_dir)

## ------------------------------------------------------------------------
unlist(lapply(strsplit(opt$raw_file_pattern, ","), as.character)) -> opt$raw_file_pattern


cat(paste0('\n# Input directory: ',opt$input_directory,'.\n'))
cat(paste0('\n# Directory for primer removed files: ',opt$primer_removed_directory,'.\n'))
cat(paste0('\n# Output directory: ',opt$output_directory,'.\n'))
cat(paste0('\n# Fwd Primer sequence: ',opt$fwd_primer,'.\n'))
cat(paste0('\n# Rev Primer sequence: ',opt$rev_primer,'.\n'))
cat(paste0('\n# Number of threads: ',opt$threads,'.\n'))
cat(paste0('\n# Number of threads: ',opt$threads,'.\n'))
cat(paste0('\n# Raw file pattern: ',opt$raw_file_pattern,'.\n'))

## ------------------------------------------------------------------------

run_atropos(raw_files_path = opt$input_directory,
            atropos = opt$atropos_bin,
            cut_dir = opt$primer_removed_directory,
            raw_file_pattern = opt$raw_file_pattern,
            output  = opt$output_directory,
            PRIMER_F = opt$fwd_primer,
            PRIMER_R = opt$rev_primer,
            NSLOTS = opt$threads)

## ------------------------------------------------------------------------


cat("\n\nHere some Info for the record keeping:\n\n")
print(sessionInfo())

cat("Parameters Used:\n\n")

cat(paste0('\n# Input directory: ',opt$input_directory,'.\n'))
cat(paste0('\n# Directory for primer removed files: ',opt$primer_removed_directory,'.\n'))
cat(paste0('\n# Output directory: ',opt$output_directory,'.\n'))
cat(paste0('\n# Fwd Primer sequence: ',opt$fwd_primer,'.\n'))
cat(paste0('\n# Rev Primer sequence: ',opt$rev_primer,'.\n'))
cat(paste0('\n# Number of threads: ',opt$threads,'.\n'))
cat(paste0('\n# Number of threads: ',opt$threads,'.\n'))
cat(paste0('\n# Raw file pattern: ',opt$raw_file_pattern,'.\n'))



cat("\n############################################################\n")
cat("All done \n")
cat("############################################################\n\n")

print(Sys.time())

## ------------------------------------------------------------------------
