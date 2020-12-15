#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'
#'

run_atropos <- function(raw_files_path,
                        atropos,
                        cut_dir = "00_atropos_primer_removed",
                        output = "dada2",
                        PRIMER_F,
                        PRIMER_R,
                        NSLOTS = 4,
                        raw_file_pattern = c("*_R1_*.gz","*_R2_*.gz"),
                        cut_file_pattern = c("_primersout_R1_.fastq.gz","_primersout_R2_.fastq.gz"),
                        MIN_L = 100){
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2); require(ShortRead); require(Biostrings)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using ShortRead version ", packageVersion('ShortRead'),'\n\n'))
  cat(paste0('\n##',"You are using Biostrings version ", packageVersion('Biostrings'),'\n\n'))
  cat(paste0('\n##',"You are using atropo version "));  system2(atropos, args = "trim --version")
  
  cat('################################\n\n')
  
  ## ------------------------------------------------------------------------
  ## get run directories (under the raw_files_path)
  list.dirs(path = raw_files_path, 
            full.names = FALSE, 
            recursive = FALSE) %>% as.vector() -> run_list  
  
  ## ------------------------------------------------------------------------
  ## get run directories (under the raw_files_path)
  setwd(raw_files_path)
  setwd("./..")
  
  for(i in seq_along(run_list)) {
    
    
    cut_path <- file.path(output, 
                          cut_dir, 
                          run_list[i])
    
    dir.create(cut_path, showWarnings = TRUE, recursive = TRUE)
    
    cat(paste0('\n# output dir :  ',cut_path,'\n'))
    
    fnFs <- sort(list.files(file.path(raw_files_path,run_list[i]), pattern = glob2rx(as.character(raw_file_pattern[1])), full.names = TRUE))
    fnRs <- sort(list.files(file.path(raw_files_path,run_list[i]), pattern = glob2rx(as.character(raw_file_pattern[2])), full.names = TRUE))
    
    exists <- file.exists(fnFs) & file.exists(fnRs)
    fnFs <- fnFs[exists]
    fnRs <- fnRs[exists]
    
    if(length(fnRs) != length(fnFs)) stop("Forward and reverse files do not match.")
    
    sample.names <- basename(fnFs) %>%
      str_extract("[^_]+")
    
    cat(paste0('\n# sample names list starts with : \n'))
    head(sample.names)
    
    fnFs_cut <- file.path(cut_path, paste0(sample.names, cut_file_pattern[1]))
    fnRs_cut <- file.path(cut_path, paste0(sample.names, cut_file_pattern[2]))
    
    sum(nchar(PRIMER_F),nchar(PRIMER_R))/2 * 2/3 -> MIN_F
    
    ## ------------------------------------------------------------------------
    ## run atropos
    for(i in seq_along(fnFs)) {
      system2(atropos, args = c("trim", "--pair-filter any",
                                "--no-indels", "--discard-untrimmed", "--max-n 0", " -T ", NSLOTS,
                                paste("-g", PRIMER_F) , paste("-G", PRIMER_R), "-n", 2, #paste("-g", PRIMER_F, "-a", dada2:::rc(PRIMER_R)) , paste("-G", PRIMER_R, "-A", dada2:::rc(PRIMER_F)), "-n", 4,
                                "-O", MIN_F %>% round(0), #primer match is >= 2/3 of primer length, from Fred Mahé's swarm pipeline
                                "-o", fnFs_cut[i], "-p", fnRs_cut[i], # output files
                                "-pe1", fnFs[i], "-pe2", fnRs[i],
                                "--minimum-length ", MIN_L)) # input files
    }
  }
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'
#'


run_dada2_qplot <- function(raw_files_path,
                            prop.sample = 20,
                            aggregate = TRUE,
                            cut_dir = "00_atropos_primer_removed",
                            qplot_dir = "01_dada2_quality_profiles",
                            output = "dada2",
                            cut_file_pattern = c("_primersout_R1_.fastq.gz","_primersout_R2_.fastq.gz"),
                            seed_value = 123){
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  
  cat('################################\n\n')
  ## ------------------------------------------------------------------------
  ## get run directories (under the raw_files_path)
  list.dirs(path = raw_files_path, 
            full.names = FALSE, 
            recursive = FALSE) %>% as.vector() -> run_list  
  
  ## ------------------------------------------------------------------------
  ## get run directories (under the raw_files_path)
  setwd(raw_files_path)
  setwd("./..")
  
  # out <- vector("list", length(run_list)*2) # Fwd and Rev
  
  for(i in seq_along(run_list)) {
    
    cut_path <- file.path(output, 
                          cut_dir, 
                          run_list[i])
    
    out_qplot<- file.path(output, 
                          qplot_dir, 
                          run_list[i])
    
    
    dir.create(out_qplot, showWarnings = TRUE, recursive = TRUE)
    
    cat(paste0('\n# output dir :  ',out_qplot,'\n'))
    
    fnFs_cut <- sort(list.files(cut_path, pattern = glob2rx(paste0("*", cut_file_pattern[1])), full.names = TRUE))
    fnRs_cut <- sort(list.files(cut_path, pattern = glob2rx(paste0("*", cut_file_pattern[2])), full.names = TRUE))
    
    
    exists <- file.exists(fnFs_cut) & file.exists(fnRs_cut)
    
    fnFs <- fnFs_cut[exists]
    fnRs <- fnRs_cut[exists]
    
    # readFastq(fnRs) %>% #idea to filter based on length
    #   ShortRead::sread()
    
    if(length(fnRs) != length(fnFs)) stop ("Forward and reverse files do not match.")
    
    sample.names <- basename(fnFs) %>%
      str_extract("[^_]+")
    
    cat(paste0('\n# sample names list starts with : \n'))
    head(sample.names)
    
    ## ------------------------------------------------------------------------
    ### Helper functions to improve  plots
    lines_seq <- function(x) {
      y <-  x%%10
      if (y > 5) {
        upper <- x+(10-y)
      } else {
        upper <- x-y
      }
      lower <- upper - 50
      return(seq(upper, lower, -10))
    } # from https://github.com/adriaaula/dada2_guidelines
    
    improve_headers <- theme(strip.background = element_blank(),
                             strip.text.x =element_text(margin = margin(.05, 0, .1, 0, "cm")),
                             strip.text = element_text(size = 10),
                             axis.text.y = element_text(size = 10))
    
    ## ------------------------------------------------------------------------
    ### Forward
    set.seed(seed_value) #random  generator necessary for reproducibility
    
    ii <- sample(length(sample.names),round(length(sample.names) * (prop.sample/100),0)+ 1 ) 
    
    qplot_f <- plotQualityProfile(fnFs[ii], aggregate = as.logical(aggregate))
    
    qplot_f_lines <- qplot_f +
      geom_vline(xintercept = lines_seq(max(qplot_f$data$Cycle)), alpha = 0.3, size = 0.3, linetype = 'dashed') +
      improve_headers +
      ggtitle(paste0(run_list[i]," - ","Forward reads"))
    
    ggsave(plot = qplot_f_lines,
           path= out_qplot,
           device="pdf",
           filename = paste0(run_list[i],"_forward.pdf"),
           width = 410,
           height = 300,
           units = 'mm')
    
    ## ------------------------------------------------------------------------
    ### Reverse
    
    qplot_r <- plotQualityProfile(fnRs[ii], aggregate = aggregate)
    
    qplot_r_lines <- qplot_r +
      geom_vline(xintercept = lines_seq(max(qplot_r$data$Cycle)), alpha = 0.3, size = 0.3, linetype = 'dashed') +
      improve_headers +
      ggtitle(paste0(run_list[i]," - ","Reverse reads"))
    
    ggsave(plot = qplot_r_lines,
           path= out_qplot,
           device="pdf",
           filename = paste0(run_list[i],"_reverse.pdf"),
           width = 410,
           height = 300,
           units = 'mm')
    
  }
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'

run_dada2_filter_denoise_merge_reads <- function(raw_files_path,
                                                 trunclen,
                                                 maxee,
                                                 truncQ = 6,
                                                 minLen,
                                                 nthreads = 8,
                                                 nbases = 20000000,
                                                 pool = "pseudo",
                                                 minover,
                                                 cut_dir = "00_atropos_primer_removed",
                                                 filt_dir = "02_dada2_filtered_denoised_merged",
                                                 output = "dada2",
                                                 cut_file_pattern = c("_primersout_R1_.fastq.gz","_primersout_R2_.fastq.gz"),
                                                 filt_pattern = c("_R1_filtered.fastq.gz","_R2_filtered.fastq.gz"),
                                                 seed_value = 123)
{
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  
  cat('################################\n\n')
  ## ------------------------------------------------------------------------
  ## get run directories (under the raw_files_path)
  list.dirs(path = raw_files_path, 
            full.names = FALSE, 
            recursive = FALSE) %>% as.vector() -> run_list  
  
  ## ------------------------------------------------------------------------
  ## get run directories (under the raw_files_path)
  setwd(raw_files_path)
  setwd("./..")
  
  for(i in seq_along(run_list)) {
    
    
    cut_path <- file.path(output, 
                          cut_dir, 
                          run_list[i])
    
    filt_path <- file.path(output, 
                           filt_dir, 
                           run_list[i])
    
    
    dir.create(filt_path, showWarnings = TRUE, recursive = TRUE)
    
    cat(paste0('\n# output dir :  ',filt_path,'\n'))
    
    fnFs_cut <- sort(list.files(cut_path, pattern = glob2rx(paste0("*", cut_file_pattern[1])), full.names = TRUE))
    fnRs_cut <- sort(list.files(cut_path, pattern = glob2rx(paste0("*", cut_file_pattern[2])), full.names = TRUE))
    
    exists <- file.exists(fnFs_cut) & file.exists(fnRs_cut)
    
    fnFs <- fnFs_cut[exists]
    fnRs <- fnRs_cut[exists]
    
    # readFastq(fnRs) %>% #idea to filter based on length
    #   ShortRead::sread()
    
    if(length(fnRs) != length(fnFs)) stop ("Forward and reverse files do not match.")
    
    sample.names <- basename(fnFs) %>%
      str_extract("[^_]+")
    
    cat(paste0('\n# sample names list starts with : \n'))
    head(sample.names)
    
    if(sample.names %>% length()  <= 1) stop("There is something wrong with the provided samples or you just provided one sample.")
    
    ## ------------------------------------------------------------------------
    filtFs <- file.path(filt_path, paste0(sample.names, filt_pattern[1]))
    filtRs <- file.path(filt_path, paste0(sample.names, filt_pattern[2]))
    
    ## ------------------------------------------------------------------------
    
    cat(paste0('\n# filterAndTrim \n'))
    
    # error when multi = T and loads of samples ?
    # https://github.com/benjjneb/dada2/issues/273
    # ?filterAndTrim
    # If memory is an issue, execute in a clean environment and reduce the chunk size n and/or the number of threads.
    # now OMP = TRUE and multi = FALSE
    
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=trunclen, trimLeft = 0, trimRight = 0,
                         maxEE=maxee,  rm.phix=TRUE, maxN=0, minLen=minLen, verbose = T,
                         compress=TRUE, multithread= nthreads, truncQ = truncQ, n = 1e5, OMP = FALSE)
    
    out2 <- data.frame(row.names = sample.names,
                       out)
    
    cat('\n# Filtering and trimming done with the following parameters:')
    cat(str_c('\n# Forward pair: trimming at ',trunclen[1],' nts and max expected error ',maxee[1]))
    cat(str_c('\n# Reverse pair: trimming at ',trunclen[2],' nts and max expected error ',maxee[2],'\n\n'))
    
    ## ------------------------------------------------------------------------
    
    cat(str_c('\n# Filtered fastq files were generated in : "', filt_path,'" \n'))
    
    ## ------------------------------------------------------------------------
    
    # get filtered reads if they still exist after trimming
    filtFs <- file.path(filt_path, paste0(sample.names, filt_pattern[1]))
    filtRs <- file.path(filt_path, paste0(sample.names, filt_pattern[2]))
    
    exists <- file.exists(filtFs) & file.exists(filtRs)
    filtFs <- filtFs[exists]
    filtRs <- filtRs[exists]
    
    if(length(filtFs) != length(filtRs)) stop("Forward and reverse filtered files do not match.")
    
    ## ------------------------------------------------------------------------
    cat(str_c('\n# derepFastq \n'))
    
    derepFs <- derepFastq(filtFs, verbose=FALSE)
    derepRs <- derepFastq(filtRs, verbose=FALSE)
    
    # Name the derep-class objects by the sample names
    names(derepFs) <- sample.names[exists]
    names(derepRs) <- sample.names[exists]
    
    cat('\n# Dereplication done\n')
    
    ## ------------------------------------------------------------------------
    set.seed(seed_value) #random  generator necessary for reproducibility
    cat(str_c('\n# learnErrors  \n'))
    
    errF <- learnErrors(derepFs, 
                        multithread=TRUE, 
                        MAX_CONSIST = 12, 
                        randomize = TRUE, 
                        nbases = nbases)
    
    errR <- learnErrors(derepRs,
                        multithread=TRUE, 
                        MAX_CONSIST = 12, 
                        randomize = TRUE, 
                        nbases = nbases)
    
    ## ------------------------------------------------------------------------
    cat(str_c('\n# plotErrors  \n'))
    
    err.plotf <- plotErrors(errF, nominalQ=TRUE) + 
      ggtitle(str_c(run_list[i], " Run - error rates - Forward reads"))
    
    ggsave(str_c(filt_path, "/" ,"errors_",run_list[i],"_fwd.pdf"),
           plot=err.plotf, width = 9, height = 8)
    
    err.plotr <- plotErrors(errR, nominalQ=TRUE) + 
      ggtitle(str_c(run_list[i], " Run - error rates - Reverse reads"))
    
    ggsave(str_c(filt_path, "/" ,"errors_",run_list[i],"_rev.pdf"),
           plot=err.plotr, width = 9, height = 8)
    
    cat('\n# Errors learnt and ploted \n')
    
    ## ------------------------------------------------------------------------
    cat(paste0('\n# dada \n'))
    set.seed(seed_value) #random  generator necessary for reproducibility
    
    if(identical(names(derepFs), names(derepRs)) != TRUE ) stop("Samples names are not consistent between Forward and Reverse dereplicated samples")
    
    dadaFs <- dada(derepFs, 
                   err=errF, 
                   multithread=TRUE, 
                   pool=ifelse(pool == "FALSE", as.logical(pool), pool))
    
    dadaRs <- dada(derepRs, 
                   err=errR, 
                   multithread=TRUE, 
                   pool=ifelse(pool == "FALSE", as.logical(pool), pool))
    
    cat('\n# DADA2 algorithm performed \n')
    
    # https://github.com/benjjneb/dada2/issues/77
    # dada2:::checkConvergence(dadaFs[[1]])
    # dada2:::checkConvergence(dadaRs[[1]])
    
    ## ------------------------------------------------------------------------
    cat(str_c('\n# mergePairs with defined ',minover ,' nt overlap \n'))
    set.seed(seed_value) #random  generator necessary for reproducibility
    
    if(sample.names %>% length()  > 1){
      if(identical(names(dadaFs), names(dadaRs)) != TRUE ) stop("Samples names are not consistent between Forward and Reverse dadas")
      if(identical(names(derepFs), names(derepRs)) != TRUE ) stop("Samples names are not consistent between Forward and Reverse dereplicated sequences")
    }
    
    
    mergers <- mergePairs(dadaFs, derepFs,
                          dadaRs, derepRs,
                          minOverlap = minover,
                          justConcatenate = FALSE,
                          maxMismatch = 0)
    
    cat('# Pairs were merged\n')
    
    ## ------------------------------------------------------------------------
    seqtab <- makeSequenceTable(mergers)
    
    cat(paste0('# Number of samples: ',dim(seqtab)[1], '\n'))
    cat(paste0('# Number of detected variants (ASVs): ',dim(seqtab)[2]))
    cat("# The variants (ASVs) have the following length distribution:")
    table(nchar(getSequences(seqtab))) # plot in the future ?
    
    cat(str_c('\n# saving seqtab as ',str_c(filt_path,"/",run_list[i],"_seqtab.rds") ,'\n'))
    saveRDS(seqtab, str_c(filt_path,"/",run_list[i],"_seqtab.rds"))
    
    ## ------------------------------------------------------------------------
    cat(str_c('\n# Plotting Sequences/ASV distribution to ',str_c(filt_path, "/", "seq_distrib_",run_list[i],".pdf") ,'\n\n'))
    
    plotLengthDistro <- function(st) {
      tot.svs <- table(nchar(colnames(st)))
      tot.reads <- tapply(colSums(st), nchar(colnames(st)), sum)
      df <- data.frame(Length=as.integer(c(names(tot.svs), names(tot.reads))),
                       Count=c(tot.svs, tot.reads),
                       Type=rep(c("ASVs", "Reads"), times=c(length(tot.svs), length(tot.reads))))
      p <- ggplot(data=df, aes(x=Length, y=Count, color=Type)) + geom_point() + 
        facet_wrap(~Type, scales="free_y") + theme_bw() + xlab("Amplicon Length")
    }
    
    plotLengthDistro(seqtab) + scale_y_log10() + 
      ggtitle(str_c("Sequence / ASV length distribution : ",run_list[i], " Run")) -> p
    
    ggsave(str_c(filt_path, "/","seq_distrib_",run_list[i],".pdf"), plot=p, width = 9, height = 8)
    
    ## ------------------------------------------------------------------------
    cat(str_c('\n# Generating summary \n'))
    
    getN <- function(x) sum(getUniques(x))
    
    data.frame(sample=sample.names,
               input=as.data.frame(out2)$reads.in,
               filtered=as.data.frame(out2)$reads.out,
               denoisedF=sapply(dadaFs, getN),
               denoisedR=sapply(dadaRs, getN),
               merged=sapply(mergers, getN),
               tabled=rowSums(seqtab)) %>%
      mutate(filtered_pc = round(filtered/input, 3)) %>%
      mutate(denoisedF_pc = round(denoisedF/filtered, 3)) %>% 
      mutate(denoisedR_pc = round(denoisedR/filtered, 3)) %>%     
      mutate(merged_pc = round(merged/denoisedF, 3)) %>% 
      mutate(filtered_merged_pc = round(merged/filtered, 3)) %>%
      mutate(input_merged_pc = round(merged/input, 3)) -> track
    
    write_tsv(data.frame(track),str_c(filt_path,"/",run_list[i],"_track_analysis.tsv"))
    
    cat("\n# The distribution of merged kept is the following:\n")
    summary(track$merged_pc)
    
    # assign(paste0(name.run,"_track"), track)
    # assign(paste0(name.run,"_seqtab"), seqtab)
    
    save(track, 
         seqtab,
         file=paste0(filt_path,"/",run_list[i],".RData"))
    #load(paste0(output,"/",name.run,".RData"))
    
    file.remove(filtFs, filtRs, fnFs_cut, fnRs_cut, showWarnings = TRUE)
    
  }
}



#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'

run_dada2_mergeRuns_removeBimeraDenovo <- function(raw_files_path,
                                                   merged_run_dir = "03_dada2_merged_runs_chimera_removed",
                                                   chimera_method = "consensus",
                                                   trim_length,
                                                   collapseNoMis = TRUE,
                                                   minOverlap = 20, # https://github.com/benjjneb/dada2/issues/518
                                                   filt_dir = "02_dada2_filtered_denoised_merged",
                                                   output = "dada2",
                                                   seed_value = 123)
{
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  
  cat('################################\n\n')
  
  
  ## ------------------------------------------------------------------------
  setwd(raw_files_path)
  setwd("./..")
  
  merged_run_path <- file.path(output, 
                               merged_run_dir)
  
  dir.create(merged_run_path, showWarnings = TRUE, recursive = TRUE)
  
  cat(paste0('\n# output dir :  ',merged_run_path,'\n'))
  
  ## ------------------------------------------------------------------------
  ## get .rds seqtables and tsv summary files:
  seqtables <- sort(list.files(paste0(output,"/",filt_dir),
                               pattern = glob2rx("*_seqtab.rds*"),
                               full.names = TRUE,
                               recursive = TRUE))
  
  track <- sort(list.files(paste0(output,"/",filt_dir),
                           pattern = glob2rx("*_track_analysis.tsv*"),
                           full.names = TRUE,
                           recursive = TRUE))
  
  summary <- map(track, read_tsv) %>% bind_rows()
  ## ------------------------------------------------------------------------
  
  cat(str_c('\n# removeBimeraDenovo running on ',seqtables, ' file '))
  
  ## ------------------------------------------------------------------------
  # https://github.com/benjjneb/dada2/issues/345
  
  list.df <- map(seqtables, readRDS)
  
  st.all <- mergeSequenceTables(tables = list.df)
  
  cat('\n# mergeSequenceTables done\n')
  
  ## ------------------------------------------------------------------------
  cat('\n# removeBimeraDenovo start\n')
  seqtab.raw <- removeBimeraDenovo(st.all, method = chimera_method,
                                   multithread = TRUE, 
                                   verbose = FALSE)
  
  cat('\n# removeBimeraDenovo done\n')
  
  num.chimera.removed <- ncol(st.all) - ncol(seqtab.raw)
  perc.num.chimera.removed <- round(100*num.chimera.removed/ncol(st.all) %>% round(2),2)
  reads.chimera.removed <- sum(colSums(st.all)) - sum(colSums(seqtab.raw))
  perc.reads.chimera.removed <- round(100*reads.chimera.removed/sum(colSums(st.all)) %>% round(2),2)
  
  cat(paste0('# ',num.chimera.removed," chimera were found and removed\n"))
  cat(paste0('# These represent ',perc.num.chimera.removed,'% of total ASVs and ',perc.reads.chimera.removed,'% of total reads\n'))
  
  ## ------------------------------------------------------------------------
  # Distribution of variants
  cat("\n# The variants (ASVs) have the following length distribution:\n")
  table(nchar(getSequences(seqtab.raw)))
  
  plotLengthDistro <- function(st) {
    tot.svs <- table(nchar(colnames(st)))
    tot.reads <- tapply(colSums(st), nchar(colnames(st)), sum)
    df <- data.frame(Length=as.integer(c(names(tot.svs), names(tot.reads))),
                     Count=c(tot.svs, tot.reads),
                     Type=rep(c("ASVs", "Reads"), times=c(length(tot.svs), length(tot.reads))))
    p <- ggplot(data=df, aes(x=Length, y=Count, color=Type)) + geom_point() + facet_wrap(~Type, scales="free_y") + theme_bw() + xlab("Amplicon Length")
  }
  plot <- plotLengthDistro(seqtab.raw) + scale_y_log10() + 
    ggtitle(str_c("Overall Sequence / ASV length distribution ")) +
    geom_vline(xintercept = trim_length[1], size = 0.1, colour = "red", alpha = 0.8, linetype = 2, show.legend = "min") + geom_vline(xintercept = trim_length[2], size = 0.1, colour = "red", alpha = 0.8, linetype = 2 ) -> p
  ggsave(str_c(merged_run_path,"/","seqtab_distrib",".pdf"),plot=plot, width = 9, height = 8)
  
  ## ------------------------------------------------------------------------
  # Trim the unespecific amplifications from our dataset
  seqtab <- seqtab.raw[,nchar(colnames(seqtab.raw)) %in% seq(trim_length[1],
                                                             trim_length[2])]
  
  cat(paste0('\n# Reads shorter than ',trim_length[1],'bp and longer than ',trim_length[2], 'bp were removed.\n'))
  cat("\n# The variants (ASVs) after length filtering have the following length distribution:\n")
  table(nchar(getSequences(seqtab)))
  cat(paste0("\n# A total of ", round((sum(colSums(seqtab)) * 100) / sum(colSums(seqtab.raw)), digits = 2), "% reads were kept after length filtering.\n\n"))
  cat(paste0("\n# A total of ", round((dim(seqtab)[2] * 100) / dim(seqtab.raw)[2], digits = 2), "% ASVs were kept after length filtering.\n\n"))
  
  ## ------------------------------------------------------------------------
  full_join(summary, 
            data.frame(sample = rownames(st.all),
                       tabled_joined = rowSums(st.all),
                       chimera_out = rowSums(seqtab.raw),
                       length_filtered = rowSums(seqtab)), by='sample') %>% 
    mutate(tabled_pc = round(tabled_joined /merged, 2)) %>%
    mutate(chimera_out_pc = round(chimera_out/tabled, 2)) %>% 
    mutate(length_filtered_pc = round(length_filtered/chimera_out, 2)) -> track
  
  #write_tsv(data.frame(track),str_c(output,"/",name,"_track_analysis.tsv"))
  
  cat("\n# The distribution of chimera reads kept is the following:\n")
  summary(track$chimera_out_pc)
  
  ## ------------------------------------------------------------------------
  if(collapseNoMis==TRUE){
    cat(str_c('\n# Saving uncollapsed .rds and fasta files as well as summary .tsv \n'))
    
    saveRDS(seqtab, str_c(merged_run_path,"/uncollapsed_no-chim-seqtab.rds"))
    
    uniquesToFasta(seqtab,
                   str_c(merged_run_path,"/uncollapsed_no-chim-seqtab.fasta"),
                   ids= str_c("asv",c(1:ncol(seqtab)), ";size=", colSums(seqtab)))
    
    write_tsv(track, str_c(merged_run_path,"/uncollapsed_track_analysis.tsv"))
    
    cat('\n# You have decided to run collapseNoMismatch on your dataset. Please note that it is only helpful IF you are working with several sequencing runs and it might take a long time to run. You might want to go further (taxonomy, ...) on your uncollapsed seqtable while it is running \n')
    
    collapsed_100 <- collapseNoMismatch(seqtab,  
                                        minOverlap = minOverlap,
                                        identicalOnly = FALSE)
    
    cat(str_c('\n# Saving collapsed .rds and fasta files as well as summary .tsv  \n'))
    
    saveRDS(collapsed_100, str_c(merged_run_path,"/minOverlap_",minOverlap,"_collapse_no_mismatch_no-chim-seqtab.rds"))
    
    uniquesToFasta(collapsed_100,
                   str_c(merged_run_path,"/minOverlap_",minOverlap,"_collapse_no_mismatch_no-chim-seqtab.fasta"),
                   ids= str_c("asv",c(1:ncol(collapsed_100)), ";size=", colSums(collapsed_100)))
    
    track %>% 
      mutate(collapsed_100 = rowSums(collapsed_100)) %>%
      mutate(collapsed_100_pc = round(collapsed_100 / length_filtered, digits = 3)) -> track.final
    
    write_tsv(track.final, str_c(merged_run_path,"/track_analysis.tsv"))
    
    cat(str_c('# Your final 100% clustered ASV table can be found in "', paste0(merged_run_path,"/seqtab.rds"),'"\n'))
    cat(str_c('# A FASTA file with your final ASVs was written in "',paste0(merged_run_path,"/seqtab.fasta"), '"\n'))
    
    cat(str_c('# In "',paste0(merged_run_path,"/track_analysis.tsv"),"\" you will find a table where you can check the loss of reads in each step. Check it out to see if everything's correct!",'\n'))
    # cat(str_c('# You have to copy them to your local computer using "scp [your.user.id]@euler.ethz.ch:',str_c(output,"/",name,"_track_analysis.tsv "),'." and go further.\n'))
    
  }
  ## ------------------------------------------------------------------------
  if(collapseNoMis==FALSE){
    cat(str_c('\n# Saving .rds and fasta files as well as summary .tsv \n'))
    
    saveRDS(seqtab, str_c(merged_run_path,"/no-chim-seqtab.rds"))
    
    uniquesToFasta(seqtab,
                   str_c(merged_run_path,"/no-chim-seqtab.fasta"),
                   ids= str_c("asv",c(1:ncol(seqtab)), ";size=", colSums(seqtab)))
    
    write_tsv(track, str_c(merged_run_path,"/track_analysis.tsv"))
    
    cat(str_c('# Your final ASV table can be found in "', paste0(merged_run_path,"/seqtab.rds"),'"\n'))
    cat(str_c('# A FASTA file with your final ASVs was written in "',paste0(merged_run_path,"/seqtab.fasta"), '"\n'))
    
    cat(str_c('# In "',paste0(merged_run_path,"/track_analysis.tsv"),"\" you will find a table where you can check the loss of reads in each step. Check it out to see if everything's correct!",'\n'))
    # cat(str_c('# You have to copy them to your local computer using "scp [your.user.id]@euler.ethz.ch:',str_c(output,"/",name,"_track_analysis.tsv "),'." and go further.\n'))
    
    cat('\n# chimera removal step is done. You can go further into the analysis (taxonomy, phylogeny) or explore collapseNoMismatch IF you are dealing with multiple runs ... it might take very long \n\n')
  }
  ## ------------------------------------------------------------------------
  # return(list("unrooted_tree" = fitGTR$tree,
  #           "rooted_tree" = ape::multi2di(fitGTR$tree))) #https://github.com/joey711/phyloseq/issues/936
  
}



#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'



run_dada_DECIPHER_taxonomy <- function(raw_files_path,
                                       taxa_dir = "04_dada2_taxonomy",
                                       method = "dada", # method = "DECIPHER" or "dada" 
                                       threshold = 60,  # used for DECIPHER and dada2 if outputBootstraps = FALSE
                                       tryRC = FALSE,
                                       collapseNoMis = TRUE,
                                       db,
                                       db_species,
                                       outputBootstraps = TRUE, #3 only for dada2 method# outputBootstraps <- TRUE 
                                       allowMultiple = TRUE,
                                       merged_run_dir = "03_dada2_merged_runs_chimera_removed",
                                       output = "dada2",
                                       seed_value = 123)
{
  
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2); require(DECIPHER)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using DECIPHER version ", packageVersion('DECIPHER'),'\n\n'))
  cat('################################\n\n')
  
  ## ------------------------------------------------------------------------
  setwd(raw_files_path)
  setwd("./..")
  
  taxa_path <- file.path(output, 
                         taxa_dir)
  
  dir.create(taxa_path, showWarnings = TRUE, recursive = TRUE)
  
  cat(paste0('\n# output dir :  ',taxa_path,'\n'))
  
  ## ------------------------------------------------------------------------
  
  merged_run_path <- file.path(output, 
                               merged_run_dir)
  
  if(collapseNoMis==FALSE){
    list.files(merged_run_path,
               pattern = glob2rx("*uncollapsed_no-chim-seqtab.rds"),
               full.names = TRUE,
               recursive = TRUE) -> seqtab.nochim
  }else{
    list.files(merged_run_path,
               pattern = glob2rx("*collapse_no_mismatch_no-chim-seqtab.rds"),
               full.names = TRUE,
               recursive = TRUE) -> seqtab.nochim
  }
  
  ## ------------------------------------------------------------------------
  
  dbname <- str_extract(basename(db), "[^.]+")
  
  set.seed(seed_value) #random  generator necessary for reproducibility
  
  seqtab.nochim <- readRDS(seqtab.nochim)
  # seqtab.nochim <- seqtab.nochim[,1:10] # for testing purpose only
  ## ------------------------------------------------------------------------
  cat(paste0('\n# You have decided to use : ',method ,' method against : ', dbname), ' database \n')
  
  ## ------------------------------------------------------------------------
  if(method=="DECIPHER"){
    
    dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
    
    load(db)
    
    ids <- IdTaxa(dna,
                  trainingSet,
                  strand = if_else(tryRC == TRUE, "both", "top"),
                  processors = NULL,
                  verbose = TRUE,
                  threshold = threshold)
    
    ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") # ranks of interest
    # need to check if this work for GTDB also...
    
    # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
    taxid <- t(sapply(ids, function(x) {
      m <- match(ranks, x$rank)
      taxa <- x$taxon[m]
      taxa[startsWith(taxa, "unclassified_")] <- NA
      taxa
    }))
    
    colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
    
    # Saving rds object
    saveRDS(taxid, paste0(taxa_path,"/", dbname,"_assignation.rds"))
    
    # Create a merged table with counts and tax
    taxid <- as_tibble(taxid, rownames = 'ASV')
    merged_table <- as_tibble(t(seqtab.nochim), rownames = 'ASV') %>%
      left_join(taxid, by = 'ASV') %>%
      mutate(ASV_id = paste0("asv",c(1:nrow(.)))) %>%
      select(ASV_id, everything())
    
    write_tsv(x = merged_table,
              path = paste0(taxa_path,"/",dbname,"_table.tsv"))
    
    cat(paste0('# The obtained taxonomy file can be found in "', paste0(output,"/", name,"_", dbname,"_assignation.rds"), '"\n\n'))
    cat(paste0('# Although we always recommend you to work directly in R with .rds files, we created a .tsv in "',paste0(output,"/", name,"_", dbname,"_table.tsv"),'" with tax and counts tables merged\n\n'))
    # cat(paste0('# You have to copy them to your local computer using "scp [your.user.id]@euler.ethz.ch:',paste0(output,"/", name,"_", dbname,"*"),'." and go further .\n\n'))
    
  }
  
  if(method=="dada"){
    
    taxa <- assignTaxonomy(seqtab.nochim, db,
                           outputBootstraps = as.logical(outputBootstraps),
                           multithread = TRUE,
                           verbose = TRUE,
                           minBoot = threshold,
                           tryRC = as.logical(tryRC))
    
    if(outputBootstraps==FALSE)
    {
      
      if(!file.exists(db_species))
      {
        if(str_extract(basename(db), "[^.]+")=="hitdb_v1")
        {  
          taxa %>% mutate(Kingdom = "Bacteria")
          colnames(taxa) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
          
        }
        # Create a merged table with counts and tax
        taxaid <- as_tibble(taxa, rownames = 'ASV')
        
        merged_table <- as_tibble(t(seqtab.nochim), rownames = 'ASV') %>%
          left_join(taxaid, by = 'ASV') %>%
          mutate(ASV_id = paste0("asv",c(1:nrow(.)))) %>%
          select(ASV_id, everything())
        
        write_tsv(x = merged_table,
                  path = paste0(taxa_path,"/", dbname,"_table.tsv"))
        
        saveRDS(taxaid, paste0(taxa_path,"/", dbname,"_assignation.rds"))
        
        cat(paste0('\n# The obtained taxonomy file can be found in "', paste0(taxa_path,"/", dbname,"_assignation.rds"), '"\n'))
        cat(paste0('# You have to copy them to your local computer using "scp [your.user.id]@euler.ethz.ch:',paste0(taxa_path,"/", dbname,"*"),'." and go further ... \n\n'))
      }
      
      if(file.exists(db_species))
      {
        taxa_Species <- addSpecies(taxa, db_species,
                                   verbose = TRUE,
                                   allowMultiple = as.logical(allowMultiple),
                                   tryRC = as.logical(tryRC))
        
        
        # Create a merged table with counts and tax
        taxa_full <- as_tibble(taxa_Species, rownames = 'ASV')  %>% 
          unite("Species",Species:tail(names(.), 1), na.rm = TRUE, remove = TRUE, sep = "|")
        
        merged_table <- as_tibble(t(seqtab.nochim), rownames = 'ASV') %>%
          left_join(taxa_full, by = 'ASV') %>%
          mutate(ASV_id = paste0("asv",c(1:nrow(.)))) %>%
          select(ASV_id, everything())
        
        write_tsv(x = merged_table,
                  path = paste0(taxa_path,"/", dbname,"_table.tsv"))
        
        saveRDS(as_tibble(taxa_Species, rownames = 'ASV')
                , paste0(taxa_path,"/", dbname,"_assignation.rds"))
        
        cat(paste0('# The obtained taxonomy file can be found in "', paste0(taxa_path,"/", dbname,"_assignation.rds"), '"\n'))
        cat(paste0('# Mutliple species assignments can be returned if you set allowMultiple = TRUE. \n'))
        # cat(paste0('# You have to copy them to your local computer using "scp [your.user.id]@euler.ethz.ch:',paste0(output,"/", name,"_", dbname,"*"),'." and go further ... \n\n'))
        
      }
      
    }
    if(outputBootstraps==TRUE)
    {
      
      if(file.exists(db_species))
      {
        boot_taxa <- taxa$boot
        
        taxa_Species <- addSpecies(taxa$tax, db_species,
                                   verbose = TRUE,
                                   allowMultiple = as.logical(allowMultiple),
                                   tryRC = as.logical(tryRC))
        
        # Create a merged table with counts and tax
        taxa_full <- left_join(as_tibble(taxa_Species, rownames = 'ASV') %>% 
                                 unite("Species",Species:tail(names(.), 1), na.rm = TRUE, sep = "|"), # because / sometimes already Pseudomonas_koreensis(AF468452)/koreensis
                               (as_tibble(taxa$boot, rownames = 'ASV')),
                               by = 'ASV', suffix = c("", "_Boot"))
        
        merged_table <- as_tibble(t(seqtab.nochim), rownames = 'ASV') %>%
          left_join(taxa_full, by = 'ASV') %>%
          mutate(ASV_id = paste0("asv",c(1:nrow(.)))) %>%
          select(ASV_id, everything())
        
        write_tsv(x = merged_table,
                  path = paste0(taxa_path,"/", dbname,"_table.tsv"))
        
        # saveRDS(as_tibble(boot_taxa, rownames = 'ASV'), 
        #         paste0(output,"/", name,"_", dbname,"_boot.rds"))
        
        # saveRDS(as_tibble(taxa_Species, rownames = 'ASV'), 
        #         paste0(output,"/", name,"_", dbname,"_assignation.rds"))
        
        saveRDS(list(as_tibble(taxa_Species, rownames = 'ASV'),
                     as_tibble(boot_taxa, rownames = 'ASV')),
                paste0(taxa_path,"/", dbname,"_assignation.rds"))
        
        cat(paste0('# The obtained taxonomy file can be found in "', paste0(taxa_path,"/", dbname,"_assignation.rds"), '"\n'))
        # cat(paste0('# You have to copy them to your local computer using "scp [your.user.id]@euler.ethz.ch:',paste0(output,"/", name,"_", dbname,"*"),'." and go further ... \n\n'))
      }
      
      if(!file.exists(db_species))
      {
        if(str_extract(basename(db), "[^.]+") == "hitdb_v1")
        {  
          taxa$tax %>% mutate(Kingdom = "Bacteria")
          colnames(taxa$tax) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
          
          boot_taxa %>% mutate(Kingdom_Boot = 100)
        }
        # Create a merged table with counts and tax
        taxaid <- left_join(as_tibble(taxa$tax, rownames = 'ASV'),(as_tibble(taxa$boot, rownames = 'ASV')),by = 'ASV', suffix = c("", "_Boot"))
        
        merged_table <- as_tibble(t(seqtab.nochim), rownames = 'ASV') %>%
          left_join(taxaid, by = 'ASV') %>%
          mutate(ASV_id = paste0("asv",c(1:nrow(.)))) %>%
          select(ASV_id, everything())
        
        write_tsv(x = merged_table,
                  path = paste0(taxa_path,"/", dbname,"_table.tsv"))
        
        saveRDS(as_tibble(taxa, rownames = 'ASV'), 
                paste0(taxa_path,"/", dbname,"_assignation.rds"))
        
        cat(paste0('# The obtained taxonomy file can be found in "', paste0(taxa_path,"/", dbname,"_assignation.rds"), '"\n'))
        # cat(paste0('# You have to copy them to your local computer using "scp [your.user.id]@euler.ethz.ch:',paste0(output,"/", name,"_", dbname,"*"),'." and go further ... \n\n'))
      }
    }
  }
  
  #return(phyloseq_object)
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'

run_DECIPHER_phangorn_phylogeny <- function(raw_files_path,
                                            method = "R",
                                            output = "dada2",
                                            phylo_dir = "05_phylo",
                                            merged_run_dir = "03_dada2_merged_runs_chimera_removed",
                                            collapseNoMis = TRUE)
{
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2); require(DECIPHER); require(phangorn)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using DECIPHER version ", packageVersion('DECIPHER'),'\n\n'))
  cat(paste0('\n##',"You are using phangorn version ", packageVersion('phangorn'),'\n\n'))
  
  cat('################################\n\n')
  
  ## ------------------------------------------------------------------------
  setwd(raw_files_path)
  setwd("./..")
  
  phylo_path <- file.path(output, 
                          phylo_dir)
  
  dir.create(phylo_path, showWarnings = TRUE, recursive = TRUE)
  
  cat(paste0('\n# output dir :  ',phylo_path,'\n'))
  
  ## ------------------------------------------------------------------------
  
  merged_run_path <- file.path(output, 
                               merged_run_dir)
  
  if(collapseNoMis==FALSE){
    list.files(merged_run_path,
               pattern = glob2rx("*uncollapsed_no-chim-seqtab.rds"),
               full.names = TRUE,
               recursive = TRUE) %>% readRDS() -> seqtab.nochim
  }else{
    list.files(merged_run_path,
               pattern = glob2rx("*collapse_no_mismatch_no-chim-seqtab.rds"),
               full.names = TRUE,
               recursive = TRUE) %>% readRDS() -> seqtab.nochim
  }
  
  ## ------------------------------------------------------------------------
  if(method=="R"){
    
    sequences <- DNAStringSet(getSequences(seqtab.nochim))
    names(sequences) <- sequences  # this propagates to the tip labels of the tree
    
    alignment <- AlignSeqs(DNAStringSet(sequences),
                           anchor=NA)
    
    phang_align <- phyDat(as(alignment, 'matrix'), type='DNA')
    
    dm <- dist.ml(phang_align)
    
    treeNJ <- NJ(dm)  # note, tip order != sequence order
    
    fit = pml(treeNJ, data=phang_align)
    
    ## negative edges length changed to 0!
    fitGTR <- update(fit, k = 4, inv = 0.2)
    
    fitGTR <- optim.pml(fitGTR, model = 'GTR', optInv = TRUE, optGamma = TRUE,
                        rearrangement = 'stochastic',
                        control = pml.control(trace = 0))
    
    detach('package:phangorn', unload = TRUE)
    detach('package:DECIPHER', unload = TRUE)
    
  }
  ## ------------------------------------------------------------------------
  if(method=="qiime2_MAFFT_Fastree"){
    # TODO
  }
  # physeq@refseq = Biostrings::DNAStringSet(taxa_names(physeq))# https://github.com/benjjneb/dada2/issues/613
  # 
  # taxa_names(physeq)  <- paste0("ASV", str_pad(seq(ntaxa(physeq)), 
  #                       nchar(ntaxa(physeq)), 
  #                       pad = "0"))
  
  ## ------------------------------------------------------------------------
  saveRDS(fitGTR$tree, paste0(phylo_path,"/unrooted_tree.rds"))
  saveRDS(phangorn::midpoint(fitGTR$tree), paste0(phylo_path,"/rooted_tree.rds"))
  
  return(list("unrooted_tree" = fitGTR$tree,
              "rooted_tree" = phangorn::midpoint(fitGTR$tree))) #https://github.com/joey711/phyloseq/issues/936
  
  
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'
#'


run_merge_phyloseq <- function(raw_files_path,
                               metadata,
                               taxa_dir = "04_dada2_taxonomy",
                               phylo = FALSE,
                               phylo_dir = "05_phylo",
                               merged_run_dir = "03_dada2_merged_runs_chimera_removed",
                               rooted_tree = TRUE,
                               output = "dada2",
                               collapseNoMis = TRUE,
                               clean = FALSE)
{
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2); require(phyloseq)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using phyloseq version ", packageVersion('phyloseq'),'\n\n'))
  
  cat('################################\n\n')
  
  ## ------------------------------------------------------------------------
  setwd(raw_files_path)
  setwd("./..")
  
  list.files(file.path(output, taxa_dir),
             pattern = glob2rx("*.tsv"),
             full.names = TRUE,
             recursive = TRUE) %>% read_tsv() -> tmp
  
  tmp %>%
    column_to_rownames("ASV") %>%
    select_if(is.numeric) %>%
    select(-contains("_boot")) %>%
    as.matrix() -> table2
  
  tmp %>%
    select(Kingdom:Species, ASV) %>%
    mutate(Species = str_replace_all(Species, "[/]", replacement = NA_character_)) %>%
    mutate_all(funs(str_replace_all(., c("unclassified|unidentified|Ambiguous_taxa"), replacement = NA_character_))) %>% # not shure this is working but that's the idea
    mutate_all(funs(str_replace_all(., c("metagenome|unidentified|uncultured"), replacement = "uncultured"))) %>%
    replace(is.na(.), "unknown") %>%
    column_to_rownames("ASV") %>%
    as.matrix() -> tax_table
  # fo increase genericity need to ignore case for taxonomy and ignore if rank is not present.
  
  phyloseq(tax_table(tax_table),
           otu_table(table2,
                     taxa_are_rows = TRUE))  %>%
    filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> physeq
  
  prune_samples(sample_sums(physeq) > 0, physeq) -> physeq
  
  ## ------------------------------------------------------------------------
  ## add ASV as refseq part of the phyloseq object
  
  physeq@refseq = Biostrings::DNAStringSet(taxa_names(physeq)) # https://github.com/benjjneb/dada2/issues/613
  
  ## ------------------------------------------------------------------------
  if(phylo == TRUE){
    if(rooted_tree == TRUE){
      list.files(file.path(output, phylo_dir),
                 pattern = glob2rx("rooted_tree.rds"),
                 full.names = TRUE,
                 recursive = TRUE) %>% readRDS() -> tree
      physeq@phy_tree = tree
    }
    if(rooted_tree == FALSE){
      list.files(file.path(output, phylo_dir),
                 pattern = glob2rx("*unrooted_tree.rds"),
                 full.names = TRUE,
                 recursive = TRUE) %>% readRDS() -> tree
      physeq@phy_tree = tree
    }
    # physeq <- merge_phyloseq(physeq,
    #                          tree)
  }
  
  taxa_names(physeq) <- paste0("ASV", str_pad(seq(ntaxa(physeq)),
                                              nchar(ntaxa(physeq)),
                                              pad = "0"))
  
  # ifelse(collapseNoMis == TRUE, track_file = "*track_analysis.tsv", track_file = "uncollapsed_track_analysis.tsv")
  
  list.files(file.path(output, merged_run_dir),
             pattern = glob2rx("track_analysis.tsv"),
             full.names = TRUE,
             recursive = TRUE) %>% read_tsv() -> track
  
  if(file.exists(metadata))
  {
    full_join(track,
              readxl::read_xlsx(metadata),
              by = c("sample" = "sample_name")) %>%
      column_to_rownames("sample") -> meta
    
    physeq <- merge_phyloseq(physeq,
                             meta %>% sample_data())
  }
  
  if(!file.exists(metadata))
  {
    
    physeq <- merge_phyloseq(physeq,
                             track %>% column_to_rownames("sample_name") %>% sample_data())
  }
  
  
  saveRDS(physeq, 
          paste0(output,"/","physeq.RDS"))
  
  return(physeq)
}

#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'
#'

run_pipe <- function(raw_files_path,
                     atropos_bin = "atropos",
                     out_dir = "dada2",
                     V = "V4",
                     PRIMER_F,
                     PRIMER_R,
                     tax_threshold = 60,
                     trim_length = c(240,400),
                     trunclen = c(260,250),
                     maxee = c(3,4),
                     minLen = 100,
                     minover = 15,
                     phylo = FALSE,
                     tryRC = FALSE,
                     tax_method = "dada",
                     metadata = NULL,
                     db = "~/db/DADA2/silva_nr_v138_train_set.fa.gz",
                     db_species = "~/db/DADA2/silva_species_assignment_v138.fa.gz",
                     SLOTS = 6){
  if(V == "V4") {
    
    PRIMER_F = "GTGCCAGCMGCCGCGGTAA"
    PRIMER_R = "GGACTACHVGGGTWTCTAAT"
    trim_length = c(220,280)
    trunclen =  c(170,160)
    maxee = c(3,4)
    minLen = 120
    minover = 40
  } 
  if(V == "V3V4"){
    
    PRIMER_F = "CCTAYGGGRBGCASCAG"
    PRIMER_R = "GGACTACNNGGGTATCTAAT"
    trim_length = c(240,500)
    trunclen =  c(260,250)
    maxee = c(4,5)
    minLen = 160
    minover = 10
    
  }
  cat(paste0('\n##',"running run_atropos() '\n\n'"))
  
  run_atropos(raw_files_path = raw_files_path,
              atropos = atropos_bin,
              PRIMER_F = PRIMER_F,
              PRIMER_R = PRIMER_R,
              output = out_dir,
              NSLOTS = SLOTS)   
  
  cat(paste0('\n##',"running run_dada2_qplot() '\n\n'"))
  
  run_dada2_qplot(raw_files_path = raw_files_path,
                  output = out_dir)
  
  cat(paste0('\n##',"running run_dada2_filter_denoise_merge_reads() '\n\n'"))
  
  run_dada2_filter_denoise_merge_reads(raw_files_path = raw_files_path,
                                       trunclen = trunclen,
                                       maxee = maxee,
                                       minLen = minLen,
                                       minover = minover,
                                       output = out_dir)
  
  cat(paste0('\n##',"running run_dada2_mergeRuns_removeBimeraDenovo() '\n\n'"))
  
  run_dada2_mergeRuns_removeBimeraDenovo(raw_files_path = raw_files_path,
                                         trim_length = trim_length,
                                         output = out_dir)
  
  cat(paste0('\n##',"running run_dada_DECIPHER_taxonomy() '\n\n'"))
  
  run_dada_DECIPHER_taxonomy(raw_files_path = raw_files_path,
                             method = tax_method, # "DECIPHER" or "dada" 
                             threshold = tax_threshold,  # used for DECIPHER and dada2 if outputBootstraps = FALSE
                             tryRC = tryRC,
                             collapseNoMis = TRUE,
                             db = db, # db = "~/db/DADA2/GTDB_r89-mod_June2019.RData"
                             db_species = db_species, # only for dada2 method
                             output = out_dir
  )
  
  if(phylo == TRUE) {
    cat(paste0('\n##',"running run_DECIPHER_phangorn_phylogeny() '\n\n'"))
    
    run_DECIPHER_phangorn_phylogeny(raw_files_path = raw_files_path,
                                    method = "R",
                                    output = out_dir) -> phylo
    
    run_merge_phyloseq(raw_files_path = raw_files_path,
                       metadata = metadata,
                       phylo = TRUE,
                       output = out_dir)
  }
  
  if(phylo == FALSE) {
    run_merge_phyloseq(raw_files_path = raw_files_path,
                       metadata = metadata,
                       phylo = FALSE,
                       output = out_dir)
  }
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'
#'


add_phylogeny_to_phyloseq <- function(phyloseq_path,
                                      method = "R",
                                      output_phyloseq = "dada2_phylo"){
  
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2); require(DECIPHER); require(phangorn)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using DECIPHER version ", packageVersion('DECIPHER'),'\n\n'))
  cat(paste0('\n##',"You are using phangorn version ", packageVersion('phangorn'),'\n\n'))
  
  cat('################################\n\n')
  
  ## ------------------------------------------------------------------------
  
  phyloseq_path %>%
    readRDS() -> physeq
  
  ## ------------------------------------------------------------------------
  if(method=="R"){
    
    sequences <- DNAStringSet(physeq@refseq)
    names(sequences) <- sequences  # this propagates to the tip labels of the tree
    
    alignment <- AlignSeqs(DNAStringSet(sequences),
                           anchor=NA)
    
    phang_align <- phyDat(as(alignment, 'matrix'), type='DNA')
    
    dm <- dist.ml(phang_align)
    
    treeNJ <- NJ(dm)  # note, tip order != sequence order
    
    fit = pml(treeNJ, data=phang_align)
    
    ## negative edges length changed to 0!
    fitGTR <- update(fit, k = 4, inv = 0.2)
    
    fitGTR <- optim.pml(fitGTR, model = 'GTR', optInv = TRUE, optGamma = TRUE,
                        rearrangement = 'stochastic',
                        control = pml.control(trace = 0))
    
    detach('package:phangorn', unload = TRUE)
    detach('package:DECIPHER', unload = TRUE)
    
  }
  ## ------------------------------------------------------------------------
  if(method=="qiime2_MAFFT_Fastree"){
    # TODO
  }
  # physeq@refseq = Biostrings::DNAStringSet(taxa_names(physeq))# https://github.com/benjjneb/dada2/issues/613
  # 
  # taxa_names(physeq)  <- paste0("ASV", str_pad(seq(ntaxa(physeq)), 
  #                       nchar(ntaxa(physeq)), 
  #                       pad = "0"))
  
  ## ------------------------------------------------------------------------
  physeq@phy_tree <- phangorn::midpoint(fitGTR$tree)
  
  physeq %>%
    saveRDS(file = paste0(output_phyloseq, ".RDS"))
  
  return(physeq)
  
}

