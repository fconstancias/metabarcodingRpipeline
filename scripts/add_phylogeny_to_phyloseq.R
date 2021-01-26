#!/usr/bin/env Rscript

rm(list = ls())

cat("\n############################################################\n")
cat("Starting \n")
cat("############################################################\n\n")

## ------------------------------------------------------------------------

require("optparse")
source("https://raw.githubusercontent.com/fconstancias/metabarcodingRpipeline/master/scripts/functions.R")

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(c("-p", "--phyloseq_path"), type="character", default=NULL, 
              help="Path of the input phyloseq object", metavar="character"),
  make_option(c("-m", "--method"), type="character", default = "R", 
              help="Method for phylogenetic reconstruction [default DECIPHER phangorn R packages see: <https://f1000research.com/articles/5-1492>]", metavar="character"),
  make_option(c("-o", "--output_phyloseq"), type="character", default= "phyloseq_tree", 
              help="Output path of the phyloseq object", metavar="character"),
  make_option(c("-f", "--fun_dir"), type="character", default= NULL, 
              help="Directory containing the R functions", metavar="character")
  
); 
opt_parser = OptionParser(option_list=option_list);
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


## ------------------------------------------------------------------------
source(opt$fun_dir)

add_phylogeny_to_phyloseq(phyloseq_path = opt$phyloseq_path,
                          method = opt$method,
                          output_phyloseq = opt$output_phyloseq)

## ------------------------------------------------------------------------


cat("\n\nHere some Info for the record keeping:\n\n")
print(sessionInfo())


cat("\n############################################################\n")
cat("All done \n")
cat("############################################################\n\n")

print(Sys.time())

## ------------------------------------------------------------------------
