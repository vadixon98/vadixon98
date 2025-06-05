############################################################
# AOS PHYLOGENY & FUNCTIONAL ANALYSIS PIPELINE
# Script: Final_Complete_Annotated.R
# Description: This pipeline retrieves Allene Oxide Synthase (AOS; CYP74A) sequences,
#              extracts promoter regions, performs multiple sequence alignment,
#              scans for motifs, and builds phylogenetic trees with bootstrap support.
# Author: Victoria Dixon
# Date: 2025-05-02 
############################################################

# === Setup Directories ===
# Create 'data' and 'results' directories if they do not exist to organize outputs
if (!dir.exists('data')) dir.create('data', recursive = TRUE)
if (!dir.exists('results')) dir.create('results', recursive = TRUE)

# === 0. Install & Load Required Packages ===
# Define CRAN and Bioconductor package lists
cran_pkgs <- c("biomartr", "ape", "phangorn", "tidyverse", "seqinr")
bioc_pkgs <- c("rentrez", "msa", "Biostrings", "ggtree", "biomaRt", "rtracklayer", "IRanges")

# Install missing CRAN packages from specified repository
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Ensure BiocManager is available to install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Install missing Bioconductor packages using BiocManager
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load all required libraries into the session
invisible(lapply(c(cran_pkgs, bioc_pkgs), library, character.only = TRUE))

# === 1. Sequence Retrieval ===
# Download coding sequences (CDS) and protein sequences for tomato
message('Downloading CDS and proteome for Solanum lycopersicum...')
cds_file  <- getCDS(organism = 'Solanum lycopersicum', path = 'data')
all_cds   <- Biostrings::readDNAStringSet(cds_file)
prot_file <- getProteome(organism = 'Solanum lycopersicum', path = 'data')
all_prot  <- Biostrings::readAAStringSet(prot_file)

# Define pattern to identify AOS sequences in headers (case-insensitive)
pattern   <- 'CYP74A|allene oxide synthase|AOS'
# Subset sequences matching the AOS pattern
aos_cds   <- all_cds[grep(pattern, names(all_cds), ignore.case = TRUE)]
aos_prot  <- all_prot[grep(pattern, names(all_prot), ignore.case = TRUE)]

# Check that sequences were found, else stop with informative error
if (length(aos_cds) == 0) stop('No AOS CDS sequences found; check header patterns')
if (length(aos_prot) == 0) stop('No AOS protein sequences found; check header patterns')

# Write selected AOS sequences to FASTA files for downstream analysis
Biostrings::writeXStringSet(aos_cds,  'data/Sl_AOS_CDS.fasta')
Biostrings::writeXStringSet(aos_prot, 'data/Sl_AOS_protein.fasta')

# === 2. Promoter Retrieval ===
message('Retrieving ~2kb promoter sequences using GFF and genome FASTA...')
# Download GFF3 annotation and genome FASTA for tomato
gff3_file <- biomartr::getGFF(organism = 'Solanum lycopersicum', path = 'data')
genome_fa <- biomartr::getGenome(organism = 'Solanum lycopersicum', path = 'data')
# Import GFF3 and genome for processing
gff       <- rtracklayer::import.gff3(gff3_file)
genome    <- Biostrings::readDNAStringSet(genome_fa)

# Extract gene features matching the AOS pattern in annotations
attrs     <- mcols(gff)$attributes
aos_genes <- gff[gff$type == 'gene' & grepl('allene oxide synthase|CYP74A|AOS', attrs, ignore.case = TRUE)]
if (length(aos_genes) == 0) stop('No AOS gene found; check GFF attributes')

# Function to extract ~2kb upstream promoter sequence for each gene
promoter_seqs <- lapply(aos_genes, function(g) {
  chr    <- as.character(seqnames(g))
  strand <- as.character(strand(g))
  # Determine start and end positions based on gene orientation
  if (strand == '+') {
    st <- start(g) - 2000; ed <- start(g) - 1
  } else {
    st <- end(g) + 1; ed <- end(g) + 2000
  }
  st <- max(st, 1)  # Ensure coordinates do not fall below 1
  seq <- Biostrings::getSeq(genome, names = chr, start = st, end = ed)
  seq_set <- Biostrings::DNAStringSet(as.character(seq))
  # Reverse-complement if gene is on the minus strand
  if (strand == '-') seq_set <- Biostrings::reverseComplement(seq_set)
  # Name sequence by gene ID or Name attribute
  nm <- if (!is.null(g$Name) && nzchar(g$Name)) g$Name else g$ID
  names(seq_set) <- nm
  seq_set
})
# Combine promoter sequences into a single DNAStringSet
promoter_chars <- sapply(promoter_seqs, function(x) as.character(x)[1])
promoters_set  <- Biostrings::DNAStringSet(promoter_chars)
names(promoters_set) <- names(promoter_seqs)
# Save promoters to FASTA for motif analysis
Biostrings::writeXStringSet(promoters_set, 'data/AOS_promoters_2kb.fasta')

# === 3. Multiple Sequence Alignment & Motif Scanning ===
# Align AOS protein sequences using ClustalW via the msa package
prot_seqs <- Biostrings::readAAStringSet('data/Sl_AOS_protein.fasta')
aln       <- msa::msa(prot_seqs, method = 'ClustalW')
msa_aa    <- as(aln, "AAStringSet")
# Save alignment for record
Biostrings::writeXStringSet(msa_aa, 'results/AOS_protein_aln.fasta')

# Example motif scan: heme-binding motif F-GGPRC in protein sequences
pattern_heme <- 'F.GGPRC'
heme_hits    <- Biostrings::vmatchPattern(pattern_heme, prot_seqs)
print(heme_hits)  # Display motif hit positions for user review

# Run FIMO from the MEME suite to scan promoter sequences for jasmonate-responsive elements
system('fimo --oc results/fimo_JA_elements motifs/JAElems.meme data/AOS_promoters_2kb.fasta')

# === 4. Phylogenetic Reconstruction ===
message('Building phylogeny from protein alignment...')
# Convert MSA to AAbin format and then to phylogenetic data format
alignment_aabin <- msa::msaConvert(aln, type = "ape::AAbin")
phydat          <- phangorn::phyDat(alignment_aabin, type = 'AA')

# Calculate maximum likelihood tree with JTT model and gamma distribution
dist_mat <- phangorn::dist.ml(phydat)
NJ_tree  <- ape::nj(dist_mat)                        # Neighbor-Joining starting tree
fit_nj   <- phangorn::pml(NJ_tree, phydat)           # Initial ML fit
fit_ml   <- phangorn::optim.pml(fit_nj, model = 'JTT', optGamma = TRUE, optInv = TRUE)

# Bootstrap analysis to assess node support (1000 replicates)
trees_bs <- phangorn::bootstrap.pml(fit_ml, bs = 1000, optNni = TRUE)
bs_tree   <- phangorn::plotBS(fit_ml$tree, trees_bs, quiet = TRUE)

# Plot and save Neighbor-Joining tree with tip labels
gg_nj <- ggtree::ggtree(NJ_tree) + ggtree::geom_tiplab(size = 3)
ggsave('results/AOS_NJ_tree.png', gg_nj, width = 6, height = 8)

# Plot and save ML tree with bootstrap support labels on nodes
p_ml <- ggtree::ggtree(bs_tree) + 
  ggtree::geom_tiplab(size = 3) +
  ggtree::geom_nodelab(aes(label = label), hjust = -0.3, size = 2)
ggsave('results/AOS_ML_tree.png', p_ml, width = 6, height = 8)

# === 5. Pipeline Completion ===
message('Pipeline complete! Results are saved in the results directory.')
