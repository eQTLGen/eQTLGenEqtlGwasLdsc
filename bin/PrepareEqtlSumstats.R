#!/usr/bin/env Rscript

library(argparse)
library(arrow)
library(data.table)
library(tidyverse)
library(IGUtilityPackage)

parser <- ArgumentParser(description = 'Prepare eQTL sumstats to the format usable for LDSC.')

parser$add_argument('--eqtl_folder', metavar = 'file', type = 'character', 
help = 'GWAS parquet folder format with per phenotype GWAS sumstats files.')
parser$add_argument('--snplist', metavar = 'file', type = 'character', 
help = 'File with matching HapMap3 variants.')
parser$add_argument('--snpref', metavar = 'file', type = 'character', 
help = 'File with eQTLGen SNP reference.')
parser$add_argument('--i2_thresh', type = 'numeric', default = 40, help = 'Heterogeneity threshold. Defaults to <40%.')
parser$add_argument('--gene', type = 'character', 
help = 'Name of gene.')
parser$add_argument('--remove_eqtl', type = 'character', default = "no",
help = '"no": eQTL regions are kept in; "yes": remove eQTL regions from LDSC.')
parser$add_argument('--remove_hla', type = 'character', default = "no",
help = '"no": HLA region is kept in; "yes": remove HLA region from LDSC.')
parser$add_argument('--window', type = 'numeric', default = 1000000,
help = 'What is the genomic window to remove around each eQTL locus.')
parser$add_argument('--p_thresh', type = 'numeric', default = 5e-8,
help = 'What is the P-value threshold to declare significant eQTLs.')

args <- parser$parse_args()

# Functions
ZtoP <- function(Z, largeZ = FALSE, log10P = TRUE){
    if (!is.numeric(Z)) {
        message("Some of the Z-scores are not numbers! Please check why!")
        message("Converting the non-numeric vector to numeric vector.")
        Z <- as.numeric(Z)
    }
    if (largeZ == TRUE) {
        P <- log(2) + pnorm(abs(Z), lower.tail = FALSE, log.p = TRUE)
        if (largeZ == TRUE & log10P == TRUE) {
            P <- -(P * log10(exp(1)))
        }
    }
    else {
        P <- 2 * pnorm(abs(Z), lower.tail = FALSE)
        if (min(P) == 0) {
            P[P == 0] <- .Machine$double.xmin
            message("Some Z-score indicates very significant effect and P-value is truncated on 2.22e-308. If relevant, consider using largeZ = TRUE argument and logarithmed P-values instead.")
        }
    }
    return(P)
}

# HapMap3 SNP list

snplist <- fread(args$snplist, header = FALSE)

message("eQTL file loading...")
eqtl <- open_dataset(args$eqtl_folder) %>%
  filter(phenotype %in% args$gene & i_squared <= args$i2_thresh) %>% collect() %>% as.data.table()
message("eQTL file loading...done!")

custom_schema <- schema(
  variant_index = int64(),
  variant =  string(),
  bp = int32(),
  chromosome = int32(),
  non_eff_allele = string(),
  eff_allele = string()
 )

message("SNP reference loading...")
snpref <- open_dataset(args$snpref, schema = custom_schema)
print(snpref)
snpref <- snpref %>% filter(variant_index %in% eqtl$variant_index) %>% collect() %>% as.data.table()
message("SNP reference loading...done!")

message("eQTL and ref merging...")
eqtl <- merge(eqtl, snpref, by = "variant_index")
message("eQTL and ref merging...done!")

eqtl$p <- ZtoP(eqtl$beta / eqtl$standard_error)

colnames(eqtl)[9] <- "chr"

message(paste(nrow(eqtl), "variants in data."))

if (args$remove_eqtl == "yes"){

message("Removing eQTLs from the calculation")

lead_variants <- IdentifyLeadSNPs(
eqtl, 
window = args$window, 
Pthresh = args$p_thresh,
snp_id_col = "variant",
snp_chr_col = "chr",
snp_pos_col = "bp",
eff_all_col = "eff_allele",
other_all_col = "non_eff_allele",
beta_col = "beta",
se_col = "standard_error"
)

lead_variants[, start := pos - args$window]
lead_variants[, end := pos + args$window]

eqtl[, start := bp]
eqtl[, end := bp]
setkey(eqtl, chr, start, end)
setkey(lead_variants, chr, start, end)
overlap_results <- foverlaps(eqtl, lead_variants, nomatch = 0)
eqtl <- eqtl[!variant_index %in% overlap_results$variant_index]

message(paste(nrow(eqtl), "variants in data after removal of eQTLs."))

}

eqtl <- eqtl[variant %in% snplist$V1]
message(paste(nrow(eqtl), "variants in data after filtering on HapMap3 variants."))

if (args$remove_hla == "yes"){
    message("Removing HLA region!")
    eqtl$chr <- as.numeric(eqtl$chr)

    eqtl <- eqtl[!(chr == 6 & bp > 25000000 & bp < 34000000)]
    message(paste(nrow(eqtl), "variants kept after HLA removal."))
}

eqtl <- data.table(snpid = eqtl$variant,
                   A1 = eqtl$eff_allele,
                   A2 = eqtl$non_eff_allele,
                   beta = eqtl$beta,
                   se = eqtl$standard_error,
                   P = eqtl$p,
                   N = eqtl$sample_size
)

message("Writing output...")
fwrite(eqtl, file = paste0(args$gene, "_InputToMunge.txt.gz"), sep = "\t")
message("Writing output...done!")
