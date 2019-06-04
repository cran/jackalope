## ----setup, echo = FALSE, message = FALSE--------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)
library(jackalope)
library(scrm)
library(ape)
set.seed(65456156)

## ----examples-create-assembly--------------------------------------------
ref <- create_genome(10, 1e6)

## ----examples-print-assembly, echo = FALSE-------------------------------
print(ref)

## ----examples-mevo-objects-----------------------------------------------
sub <- sub_TN93(pi_tcag = c(0.1, 0.2, 0.3, 0.4),
                alpha_1 = 0.0001, alpha_2 = 0.0002,
                beta = 0.00015)
ins <- indels(rate = 2e-5, max_length = 10)
del <- indels(rate = 1e-5, max_length = 40)

## ----examples-reads-for-assembly-pacbio, eval = FALSE--------------------
#  pacbio(ref, out_prefix = "pacbio", n_reads = 2 * 500e3)

## ----examples-reads-for-assembly-hybrid, eval = FALSE--------------------
#  pacbio(ref, out_prefix = "pacbio", n_reads = 500e3)
#  illumina(ref, out_prefix = "illumina", n_reads = 500e6, paired = TRUE,
#           read_length = 100)

## ----examples-reads-for-assembly-illumina, eval = FALSE------------------
#  illumina(ref, out_prefix = "ill_pe", n_reads = 500e6, paired = TRUE,
#           read_length = 100)
#  illumina(ref, out_prefix = "ill_mp", seq_sys = "MSv3",
#           read_length = 250, n_reads = 50e6, matepair = TRUE,
#           frag_mean = 3000, frag_sd = 500)

## ----examples-divergence-scrm--------------------------------------------
ssites <- scrm(paste("40", ref$n_seqs(), "-t 100 -I 2 20 20 100"))

## ----examples-divergence-create------------------------------------------
vars <- create_variants(ref, vars_ssites(ssites), sub, ins, del)

## ----examples-divergence-create-print, echo = FALSE----------------------
print(vars)

## ----examples-divergence-write-vcf, eval = FALSE-------------------------
#  write_vcf(vars, "variants")

## ----examples-divergence-illumina-pool, eval = FALSE---------------------
#  illumina(vars, out_prefix = "vars_illumina", n_reads = 500e6, paired = TRUE,
#           read_length = 100, barcodes = c(rep("AACCGCGG", 20),
#                                           rep("GGTTATAA", 20)))

## ----examples-divergence-illumina-individual, eval = FALSE---------------
#  illumina(vars, out_prefix = "vars_illumina", n_reads = 500e6, paired = TRUE,
#           read_length = 100, sep_files = TRUE)

## ----examples-phylogeny-tree---------------------------------------------
tree <- rcoal(10)

## ----examples-phylogeny-tree-create-for-show-----------------------------
vars <- create_variants(ref, vars_phylo(tree), sub, ins, del)

## ----examples-phylogeny-tree-variants-print, echo = FALSE----------------
print(vars)

## ----examples-phylogeny-tree-illumina, eval = FALSE----------------------
#  illumina(vars, out_prefix = "phylo_tree", seq_sys = "MSv3",
#           read_length = 250, n_reads = 50e6)

## ----examples-phylogeny-gtrees-scrm--------------------------------------
# Run scrm for one sequence size:
one_seq <- function(size) {
    sims <- scrm(
        paste("24 1",
              # Output gene trees:
              "-T",
              # Recombination:
              "-r 1", size,
              # 3 species with no ongoing migration:
              "-I 3", paste(rep("8", 3), collapse = " "), "0",
              # Species 2 derived from 1 at time 1.0:
              "-ej 0.5 2 1",  
              # Species 3 derived from 2 at time 0.5:
              "-ej 0.25 3 2"
        ))
    return(sims$trees[[1]])
}
# For all sequences:
gtrees <- list(trees = lapply(ref$sizes(), one_seq))

## ----examples-phylogeny-gtrees-create-variants---------------------------
vars <- create_variants(ref, vars_gtrees(gtrees),
                        sub, ins, del)

## ----examples-phylogeny-gtrees-variants-print, echo = FALSE--------------
print(vars)

## ----examples-phylogeny-write-vcf, eval = FALSE--------------------------
#  write_vcf(vars, out_prefix = "var_gtrees",
#            sample_matrix = matrix(1:vars$n_vars(), ncol = 2, byrow = TRUE))

## ----examples-phylogeny-gtrees-illumina, eval = FALSE--------------------
#  illumina(ref, out_prefix = "phylo_gtrees",
#           seq_sys = "MSv3",
#           read_length = 250,
#           n_reads = 50e6)

