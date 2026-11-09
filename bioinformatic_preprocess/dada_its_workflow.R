library("dada2")
library("tidyverse")
library("ShortRead")
library("Biostrings")
library("microViz")
library("phyloseq")
packageVersion("dada2")

# prepare data information
path <- "/mnt/data/chomic/seqme2025/tatry_gradient/its/dada_tatry_gradient_its/"

# fastq filenames 
# input files should be unzipped because of cutadapt
fnR1s <- sort(list.files("/mnt/data/chomic/seqme2025/tatry_gradient/its/raw", pattern="_R1.fq", full.names = TRUE))
fnR2s <- sort(list.files("/mnt/data/chomic/seqme2025/tatry_gradient/its/raw/", pattern="_R2.fq", full.names = TRUE))

count_reads <- function(file) {
  con <- gzfile(file, "r")
  lines <- 0
  while (length(chunk <- readLines(con, n = 100000)) > 0) {
    lines <- lines + length(chunk)
  }
  close(con)
  return(lines / 4)
}

read_counts <- sapply(fnR1s, count_reads)
data.frame(file = basename(fnR1s), reads = read_counts)



# Extract sample names
extract_sample_ids <- function(paths) {
  filenames <- basename(paths)
  nums <- sub("^([0-9]+).*_R[12]\\.fq$", "\\1", filenames)
  formatted <- sprintf("s%03d", as.numeric(nums))
  return(formatted)
}

sample_names <- extract_sample_ids(fnR1s)


# PRIMERS
FWD <- "GTGARTCATCGARTCTTTG"  ## gITS7
REV <- "TCCTCCGCTTATTGATATGC"    ## ITS4


# check presence and orientation of the primers in the reads

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)


# discard all reads with Ns prior to primer matching
fnR1sfiltN <- file.path(path, "filtN", basename(fnR1s)) # Put N-filterd files in filtN/ subdirectory
fnR2sfiltN <- file.path(path, "filtN", basename(fnR2s))
filterAndTrim(fnR1s, fnR1sfiltN, fnR2s, fnR2sfiltN, maxN = 0, multithread = TRUE, compress = F)

# check frequency of primers in the first pair of read files
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnR1sfiltN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnR2sfiltN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnR1sfiltN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnR2sfiltN[[1]]))

# Remove primers
cutadapt <- "/usr/bin/cutadapt" # path to cutadapt directory
system2(cutadapt, args = "--version") # check functionality

path.cut <- file.path(path, "cutadapt") # path to cutdapt-ed folder
if(!dir.exists(path.cut)) dir.create(path.cut)
fnR1s.cut <- file.path(path.cut, basename(fnR1s)) # path to cutadapt-ed R1 files
fnR2s.cut <- file.path(path.cut, basename(fnR2s)) # path to cutadapt-ed R2 files

# make primer reverse complement versions
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run cutadapt
for(i in seq_along(fnR1s)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads,
                             "--cores=190", # multithreading (190 cores)
                             "--minimum-length 1", # 
                             "-o", fnR1s.cut[i], "-p", fnR2s.cut[i], # output files
                             fnR1sfiltN[i], fnR2sfiltN[i])) # input files
}

# check presence of primers in cutadapt-ed files (first pair of files)
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnR1s.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnR2s.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnR1s.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnR2s.cut[[1]]))



# Forward and reverse fastq filenames have the format:
cutR1s <- sort(list.files(path.cut, pattern = "_R1.fq", full.names = TRUE))
cutR2s <- sort(list.files(path.cut, pattern = "_R2.fq", full.names = TRUE))

cutR1_fq <- readFastq("cutadapt", pattern ="R1.fq")
cutR2_fq <- readFastq("cutadapt", pattern ="R2.fq")

# inspect quality profiles
cutR1_qual_plot <- plotQualityProfile(cutR1s, n = 1e5, aggregate = TRUE) + ggtitle("cutR1")
cutR2_qual_plot <- plotQualityProfile(cutR2s, n = 1e5, aggregate = TRUE) + ggtitle("cutR2")


# length distribution cut R1
width(sread(cutR1_fq)) %>%
  hist(breaks = 250, xlim = c(200,250), main = "Read Length Distribution", xlab = "Read Length (bp)", 
       col = "khaki4", border = "white")

# length distribution cut R2
width(sread(cutR2_fq)) %>%
  hist(breaks = 250, xlim = c(200,250), main = "Read Length Distribution", xlab = "Read Length (bp)", 
       col = "khaki4", border = "white")

cutR1_cplx_plot <- plotComplexity(cutR1s, kmerSize = 2, # Size of kmers (default is 2)
                                  window = 50,
                                  by = 10,   # Step size for sliding window
                                  n = 1e5,  # Number of reads to sample
                                  bins = 100,     # Histogram bins
                                  aggregate = TRUE  # Aggregate across multiple files
) + ggtitle("cutR1")

cutR2_cplx_plot <- plotComplexity(cutR2s, kmerSize = 2, window = 50, by = 10, n = 1e5, bins = 100,
                                  aggregate = TRUE) + ggtitle("cutR2") 

# remove reads with poly-G tail 
cutR1_polyGcheck <- checkPolyG(cutR1s)
cutR2_polyGcheck <- checkPolyG(cutR2s)

removePolyG_paired(cutR1s, cutR2s)

nopolyR1s <- sort(list.files("nopolyG", pattern="_R1_nopolyG.fq.gz", full.names = TRUE))
nopolyR2s <- sort(list.files("nopolyG", pattern="_R2_nopolyG.fq.gz", full.names = TRUE))

nopolyR1_polyGcheck <- checkPolyG(nopolyR1s)
nopolyR2_polyGcheck <- checkPolyG(nopolyR2s)

# discard reads with  >1 expected errors, length >230, low complexity (<8 acc. to complexity plots)
#  truncate at length 230 bp
filtR1s <- file.path(path, "filtered", paste0(sample_names, "_R1_filt.fq.gz"))
filtR2s <- file.path(path, "filtered", paste0(sample_names, "_R2_filt.fq.gz"))
names(filtR1s) <- sample_names
names(filtR2s) <- sample_names

out <- filterAndTrim(nopolyR1s, filtR1s, nopolyR2s, filtR2s, maxLen = 230,
                     maxEE=1, truncQ=2, rm.phix=TRUE, rm.lowcomplex = 8,
                     compress=TRUE, multithread=TRUE, verbose = TRUE) # On Windows set multithread=FALSE

head(out)

filtR1_fq <- readFastq("filtered", pattern ="R1_filt.fq.gz")
filtR2_fq <- readFastq("filtered", pattern ="R2_filt.fq.gz")

# inspect quality profiles
filtR1_qual_plot <- plotQualityProfile(filtR1s, n = 1e5, aggregate = TRUE) + ggtitle("filtR1")
filtR2_qual_plot <- plotQualityProfile(filtR2s, n = 1e5, aggregate = TRUE) + ggtitle("filtR2")

# length distribution filt R1
width(sread(filtR1_fq)) %>%
  hist(breaks = 250, xlim = c(200,250), main = "Read Length Distribution", xlab = "Read Length (bp)", 
       col = "khaki4", border = "white")

width(sread(filtR2_fq)) %>%
  hist(breaks = 250, xlim = c(200,250), main = "Read Length Distribution", xlab = "Read Length (bp)", 
       col = "khaki4", border = "white")


filtR1_cplx_plot <- plotComplexity(filtR1s, kmerSize = 2, window = 50, by = 10, n = 1e5, bins = 100,
                                   aggregate = TRUE) + ggtitle("filtR1")
filtR2_cplx_plot <- plotComplexity(filtR2s, kmerSize = 2, window = 50, by = 10, n = 1e5, bins = 100,
                                   aggregate = TRUE) + ggtitle("filtR2")


#select representative subset of samples
library(dplyr)
metadata <- read.table("../../map_tatry_gradient_250806.txt")

# 1 sample from each tzpe of sample and sampling time
selected_samples <- metadata %>%
  group_by(site, sampling_time) %>%
  slice_sample(n = 1) %>%  # nebo n = 2 pro větší pokrytí
  ungroup()

selected_ids <- selected_samples$sample_id
filtR1s_sub <- filtR1s[names(filtR1s) %in% selected_ids]
filtR2s_sub <- filtR2s[names(filtR2s) %in% selected_ids]


# learn error rates
# errR1_old <- learnErrors(filtR1s_sub, nbases = Inf, multithread=TRUE, verbose = TRUE, MAX_CONSIST = 20)
# 
# plotErrors(errR1_old, nominalQ=TRUE)

binnedQs <- c(2, 11, 25, 37)
binnedQualErrfun <- makeBinnedQualErrfun(binnedQs)

errR1 <- learnErrors(filtR1s_sub, errorEstimationFunction=binnedQualErrfun, MAX_CONSIST = 20,
                     nbases = Inf, multithread=TRUE, verbose = TRUE)
errR2 <- learnErrors(filtR2s_sub, errorEstimationFunction=binnedQualErrfun, MAX_CONSIST = 20,
                     nbases = Inf, multithread=TRUE, verbose = TRUE)

plotErrors(errR1, nominalQ=FALSE)
plotErrors(errR2, nominalQ=FALSE)


# construct ASV table, remove chimeras
# pooling FALSE, TRUE nebo pseudo
dadaR1s <- dada(filtR1s, err=errR1, multithread=TRUE, errorEstimationFunction=binnedQualErrfun,
                pool = "pseudo", OMEGA_A = 1e-40, verbose = TRUE)


dadaR2s <- dada(filtR2s, err=errR1, multithread=TRUE, errorEstimationFunction=binnedQualErrfun,
                pool = "pseudo", OMEGA_A = 1e-40, verbose = TRUE)

# merge paired-end reads
merged <- mergePairs(dadaR1s, filtR1s, dadaR2s, filtR2s, minOverlap = 150,
                     maxMismatch = 0, verbose = T, trimOverhang = T)

# ASV table
seqtab <- makeSequenceTable(merged, orderBy = "abundance")

seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)



asv_sums <- colSums(seqtab_nochim)
asv_prevalence <- colSums(seqtab_nochim > 0)

hist(asv_prevalence, breaks = 105, xlim = c(0,10))

keep_asvs <- (asv_total >= 10) & (asv_prevalence >= 2)
seqtab_final <- seqtab_nochim[, keep_asvs]






# ASV counts
ncol(seqtab_nochim) # 5061
ncol(seqtab_final) # 1807


depths <- rowSums(seqtab_nochim)
summary(depths)
hist(depths, breaks = 100)


# write ASVs FASTA file for further taxonomic assignment
write.fasta.dada<-function(dada2, file){
  seqs<-dada2::getSequences(dada2)
  hash<-paste0(">ASV_",seq(1,length(colnames(dada2))))
  write(c(rbind(hash, seqs)),file)
}

write.fasta.dada(seqtab_final, paste0(path,"asvs.fa"))

seqtab_phy <- t(seqtab_final)
rownames(seqtab_phy) <- paste0("ASV_",seq(1,length(rownames(seqtab_phy))))

library("openxlsx")

phy <- phyloseq(otu_table(seqtab_phy, taxa_are_rows = T),
                sample_data(metadata),
                tax_table(as.matrix(read.xlsx("blast/tatry_gradient_BLAST_10_best_hits.xlsx", 
                                              sheet = "tax_table", rowNames = T))))

ntaxa(phy) # 1807



  
  
phy_sh <-
  phy %>% # ntaxa 1807
  subset_samples(catchment != "ctrl") %>%
  subset_taxa(Kingdom == "Fungi") %>%   #ntaxa()# 1736
  tax_glom("specSH")  #%>% ntaxa() # 1506

tax_fix_interactive(phy_sh)

ord_explore(phy_sh)

saveRDS(phy_sh, "tatry_gradient_its_phy_sh.rds")





######## FUNCTIONS ########
checkPolyG <- function(files, pattern = "G{8,}$") {
  results <- data.frame(
    sample = character(),
    polyG_count = integer(),
    total_reads = integer(),
    polyG_percent = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (file in files) {
    fq <- ShortRead::readFastq(file)
    seqs <- as.character(ShortRead::sread(fq))
    
    polyG_hits <- grepl(pattern, seqs)
    
    sample_name <- tools::file_path_sans_ext(basename(file))
    polyG_count <- sum(polyG_hits)
    total_reads <- length(seqs)
    polyG_percent <- round(100 * polyG_count / total_reads, 2)
    
    results <- rbind(
      results,
      data.frame(
        sample = sample_name,
        polyG_count = polyG_count,
        total_reads = total_reads,
        polyG_percent = polyG_percent,
        stringsAsFactors = FALSE
      )
    )
    
    # Úklid
    rm(fq, seqs, polyG_hits)
    gc(verbose = FALSE)
  }
  
  rownames(results) <- results$sample
  results$sample <- NULL
  return(results)
}



removePolyG <- function(infiles, outdir = "nopolyG", pattern = "G{8,}$") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  for (infile in infiles) {
    fq <- readFastq(infile)
    seqs <- sread(fq)
    keep <- !grepl(pattern, as.character(seqs))
    fq_filtered <- fq[keep]
    
    # Získání názvu souboru
    fname <- basename(infile)
    
    # Detekce R1 nebo R2
    if (grepl("_R1", fname)) {
      suffix <- "_R1_nopolyG.fq.gz"
      sample <- sub("_R1.*", "", fname)
    } else if (grepl("_R2", fname)) {
      suffix <- "_R2_nopolyG.fq.gz"
      sample <- sub("_R2.*", "", fname)
    } else {
      warning(paste("Soubor neobsahuje R1/R2:", fname))
      next
    }
    
    # Vytvoření výstupní cesty
    outfile <- file.path(outdir, paste0(sample, suffix))
    
    # Přepis souboru, pokud existuje
    if (file.exists(outfile)) file.remove(outfile)
    writeFastq(fq_filtered, outfile, compress = TRUE)
    
    cat(fname, ": removed", sum(!keep), "poly-G sequences out of", length(seqs), "\n")
  }
}

removePolyG_paired <- function(R1_files, R2_files, outdir = "nopolyG", pattern = "G{10,}$") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  for (i in seq_along(R1_files)) {
    fq1 <- ShortRead::readFastq(R1_files[i])
    fq2 <- ShortRead::readFastq(R2_files[i])
    
    seqs1 <- as.character(ShortRead::sread(fq1))
    seqs2 <- as.character(ShortRead::sread(fq2))
    
    polyG1 <- grepl(pattern, seqs1)
    polyG2 <- grepl(pattern, seqs2)
    
    # Odstranit čtení, pokud poly-G je v R1 nebo R2
    keep <- !(polyG1 | polyG2)
    
    fq1_filtered <- fq1[keep]
    fq2_filtered <- fq2[keep]
    
    sample <- sub("_R1.*", "", basename(R1_files[i]))
    out1 <- file.path(outdir, paste0(sample, "_R1_nopolyG.fq.gz"))
    out2 <- file.path(outdir, paste0(sample, "_R2_nopolyG.fq.gz"))
    
    if (file.exists(out1)) file.remove(out1)
    if (file.exists(out2)) file.remove(out2)
    
    ShortRead::writeFastq(fq1_filtered, out1, compress = TRUE)
    ShortRead::writeFastq(fq2_filtered, out2, compress = TRUE)
    
    cat(sample, ": removed", sum(!keep), "paired reads with poly-G\n")
  }
}


#  BLAST commands
# ======================
# blastn -query asvs.fa -db /mnt/data/chomic/blast/sh_general_release_04.04.2024/sh_general_release_dynamic_04.04.2024.fasta -num_threads 190 -out tatry_10_best_hits.xml -outfmt=5 2>&1 -max_target_seqs 10
# 
# 
# conda activate bioconda
# 
# python /mnt/data/chomic/blast/Blast_xml_parse.py tatry_10_best_hits.xml
# 
# gawk 'BEGIN{FS=OFS="\t"}{print $1,$4,$6,$7,$9,$10,$12,$13,$15,$16,$18,$19,$21,$22,$24,$25,$27,$28,$30,$31}' < 10hits.txt | sed '1 i\SeqID+first_hit+score+e-value+query len+query start+query end+target len+target start+target end+align len+identities+gaps+coverage%+id%+next_hits' > tatry_BLAST_10_best_hits.txt
