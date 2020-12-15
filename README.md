# metabarcodingRpipeline

## Clone the repository:

Change the directory where you would like to clone the repository.

	$ cd my_directory

Use ``git clone`` to clone on your computer the repository including the functions and test data.

	$ git clone https://github.com/fconstancias/metabarcodingRpipeline.git


## Install the proper conda envirionment:
### Create a dedicated conda environment:
	$ conda create -n metabarcodingRpipeline -y
### Activate the conda environment:
	$ conda activate metabarcodingRpipeline
### Install R and atropos:
	(metabarcodingRpipeline)$ conda install -c bioconda atropos -y; conda install -c conda-forge r-base -y
### run R:
	(metabarcodingRpipeline)$ R
### install necessary R packages - from the R session:

Running the following commands will take a bit of time.

	install.packages("devtools");install.packages("optparse");devtools::install_github("tidyverse/tidyverse");devtools::install_github("KlausVigo/phangorn");devtools::install_github("benjjneb/dada2")
	
	if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
	BiocManager::install("ShortRead");BiocManager::install("DECIPHER");BiocManager::install("phyloseq")
	
	
### quit R:
	quit(save = "no")
	

## Running the pipeline:

Activate the dedicated conda environment:

	$ conda activate metabarcodingRpipeline



Use ``Rscript`` to run the pipeline and specify some necessary parameters e.g.: *databases*

	(metabarcodingRpipeline)$ Rscript scripts/dada2_metabar_pipe.R \
		-i test-data \
		-o dada2 \
		-V V3V4 \
		--metadata test-data/metadata.xlsx \
		--database ~/db/DADA2/silva_nr99_v138_train_set.fa.gz \
		--database_for_species_assignments ~/db/DADA2/silva_species_assignment_v138.fa.gz > run_logs.txt 2>&1
		
The ``> mylogs.txt 2>&1`` trick will redirect what is printed on the screen to a file including potential errors and also parameters that you used.

## Constructing and adding a phylogenetic tree from a phyloseq object:
default using <https://f1000research.com/articles/5-1492>

Depending on your machine and richness of your phyloseq object it might take quite some time. Hence, phylo is set to FALSE by default on the main pipeline.

OTU/ASV sequences must be stored in the @refseq dataset of the corresponding phyloseq object.

	(metabarcodingRpipeline)$  Rscript scripts/add_phylogeny_phyloseq.R \
	-p dada2/physeq.RDS \
	-o dada2/physeq_phylo


## TO DO:

- <s>add phylogenetic tree to a phyloseq object</s>
- replace taxonomic assignments of a phyloseq object using alternative method/ database

## Help:


activate the dedicated conda environment:

	$ conda activate metabarcodingRpipeline

print help:
	
	(metabarcodingRpipeline)$ Rscript scripts/dada2_metabarcoding_pipeline.R --help


	Usage: scripts/dada2_metabarcoding_pipeline.R [options]


	Options:
	
	-i CHARACTER, --input_directory=CHARACTER
	Path of the input directory containing raw _R1_ and _R2_ raw reads in their respective run sub-directories 
 
              e.g., -i raw [contains raw/run1 and raw/run2 

              sample names are taken from .fastq.gz files names after the first '_' ]

	-a CHARACTER, --atropos_binary=CHARACTER
		Path of atropos program [used for primer removal]

	-o CHARACTER, --output_directory=CHARACTER
		Name of the output directory

	-V CHARACTER, --pipeline=CHARACTER
		V4 or V3V4 will use default primers and parameters as used in the FBT lab [primers, trunc, maxee, overlap, expected length, ...]

	-t CHARACTER, --tax_method=CHARACTER
		User can specify using dada2 or DECIPHER for taxonomic assignments [default dada]

	--tax_threshold=NUMERIC
		Thershold for taxonomic assignments [if --tax_metod = dada: he minimum bootstrap confidence for assigning a taxonomic level. if --tax_method DECIPHER: Numeric specifying the confidence at which to truncate the output taxonomic classifications. ]

	--metadata=CHARACTER
		Path to excel document containing metadata [Sample identifier column should be sample_name]

	--database=CHARACTER
		Path to the taxonomic database

	--database_for_species_assignments=CHARACTER
		Path to the speies-level taxonomic database [only for -t dada]

	--phylo=CHARACTER
		Compute phylogenetic tree from the ASV sequence ?

	--PRIMER_F=CHARACTER
		Sequence of the gene specific Fwd primer to be removed with atropos [if using -V V4 or V3V4, this parameter is already set]

	--PRIMER_R=CHARACTER
		Sequence of the gene specific Rev primer to be removed with atropos [if using -V V4 or V3V4, this parameter is already set]

	--minover=NUMERIC
		Minimum overlap for merginf R1 and R2 reads [if using -V V3 or V3V4, this parameter is already set]

	--trunclen=NUMERIC
		Nucleotide position to truncate the Fwd and Rev reads at [if using -V V4 or V3V4, this parameter is already set]

	--trim_length=NUMERIC
		ASV of length outside the range will be discarded [i.e., insilco size exclusion of ASV - if using -V V4 or V3V4, this parameter is already set]

	--maxee=NUMERIC
		Maximum expected error for Fwd and Rev reads [if using -V V4 or V3V4, this parameter is already set]

	--minLen=NUMERIC
		Minimul read length [if using -V V4 or V3V4, this parameter is already set]

	-T NUMERIC, --slots=NUMERIC
		Number of threads to perform the analyses

	-h, --help
		Show this help message and exit

