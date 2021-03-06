#' jackalope: An efficient, flexible molecular evolution and sequencing simulator.
#'
#' `jackalope` simply and efficiently
#' simulates (i) haplotypes from reference genomes and (ii) reads from both Illumina
#' and Pacific Biosciences (PacBio) platforms.
#' It can either read reference genomes from FASTA files or simulate new ones.
#' Variant haplotypes can be simulated using summary statistics, phylogenies,
#' Variant Call Format (VCF) files, and coalescent simulations—the latter of which
#' can include selection, recombination, and demographic fluctuations.
#' `jackalope` can simulate single, paired-end, or mate-pair Illumina reads,
#' as well as reads from Pacific Biosciences
#' These simulations include sequencing errors, mapping qualities, multiplexing,
#' and optical/PCR duplicates. All outputs can be written to standard file formats.
#'
#'
#' @importFrom Rcpp evalCpp
#' @import zlibbioc
#' @useDynLib jackalope, .registration = TRUE
#'
#' @docType package
#' @name jackalope
NULL
