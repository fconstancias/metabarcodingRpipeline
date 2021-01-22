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
  make_option(c("-o", "--output_directory"), type="character", default= "dada2", 
              help="Name of the output directory", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default= "dada", 
              help="Path of the Excel metadata file - ", metavar="character"),
  make_option(c("--phylo"), type="character", default= FALSE, 
              help="was phylognetic analysis performed", metavar="character"),
  make_option(c("--rooted_tree"), type="character", default= NULL, 
              help="select rooted tree", metavar="character"),
  make_option(c("--collapse_no_mis"), type="character", default= TRUE, 
              help="Did you perform 100% similarity clustering for ASV (collapse_no_mis) ?", metavar="character"),
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
cat(paste0('\n# Input directory: ',opt$input_directory,'.\n'))
cat(paste0('\n# Output directory: ',opt$output_directory,'.\n'))

cat(paste0('\n# Metadata excel file path',opt$metadata,'.\n'))
cat(paste0('\n# Was phylognetic analysis performed ? ',opt$phylo,'.\n'))
cat(paste0('\n# If so, would you include a rooted phylogenetic tree ? ',opt$rooted_tree,'.\n'))
cat(paste0('\n# Did you perform 100% similarity clustering for ASV (collapse_no_mis) ? ',opt$collapse_no_mis,'.\n'))


## ------------------------------------------------------------------------

run_merge_phyloseq(raw_files_path = opt$input_directory,
                   metadata = opt$metadata ,#"/Users/fconstan/Documents/GitHub/metabarcodingRpipeline/test-data/metadata.xlsx",
                   phylo = as.logical(opt$phylo),
                   taxa_dir = "04_dada2_taxonomy",
                   phylo_dir = "05_phylo",
                   merged_run_dir = "03_dada2_merged_runs_chimera_removed",
                   rooted_tree = opt$rooted_tree,
                   output = opt$output_directory,
                   collapseNoMis = opt$collapse_no_mis) -> physeq

physeq
## ------------------------------------------------------------------------


cat("\n\nHere some Info for the record keeping:\n\n")
print(sessionInfo())

cat("Parameters Used:\n\n")

## ------------------------------------------------------------------------
cat(paste0('\n# Input directory: ',opt$input_directory,'.\n'))
cat(paste0('\n# Output directory: ',opt$output_directory,'.\n'))

cat(paste0('\n# Metadata excel file path',opt$metadata,'.\n'))
cat(paste0('\n# Was phylognetic analysis performed ? ',opt$phylo,'.\n'))
cat(paste0('\n# If so, would you include a rooted phylogenetic tree ? ',opt$rooted_tree,'.\n'))
cat(paste0('\n# Did you perform 100% similarity clustering for ASV (collapse_no_mis) ? ',opt$collapse_no_mis,'.\n'))

cat("\n############################################################\n")
cat("All done \n")
cat("############################################################\n\n")

print(Sys.time())

## ------------------------------------------------------------------------



