#!/usr/bin/env Rscript

library(argparse)
library(arrow)
library(data.table)
library(tidyverse)

parser <- ArgumentParser(description = 'Prepare GWAS sumstats to the format usable for LDSC.')

parser$add_argument('--gwas_folder', metavar = 'file', type = 'character', 
help = 'GWAS parquet folder format with per phenotype GWAS sumstats files.')
parser$add_argument('--pheno_manifest', metavar = 'file', type = 'character', 
help = 'GWAS manifest file, specifying if phenoytype is continuous or case-control and sample size.')
parser$add_argument('--snplist', metavar = 'file', type = 'character', 
help = 'File with matching HapMap3 variants.')
parser$add_argument('--snpref', metavar = 'file', type = 'character', 
help = 'eQTLGen SNP reference file.')
parser$add_argument('--phenotype', type = 'character', 
help = 'Name of GWAS phenotype.')
parser$add_argument('--output', metavar = 'file', type = 'character', 
help = 'Output file.')

args <- parser$parse_args()

# HapMap3 SNP list
snplist <- fread(args$snplist, header = FALSE)
snpref <- open_dataset(args$snpref)
message("SNP ref read!")
snpref <- snpref %>% filter(variant %in% snplist$V1) %>% as.data.table()
message("SNP ref filtered!")

gwas <- open_dataset(args$gwas_folder) %>% 
filter(variant_index %in% snpref$variant_index & gwas_id %in% args$phenotype) %>%
collect() %>% as.data.table()

message("GWAS filtered!")

gwas <- data.table(snpid = gwas$variant,
                   A1 = gwas$str_allele2,
                   A2 = gwas$str_allele1,
                   beta = gwas$beta,
                   se = gwas$se,
                   P = as.numeric(gwas$p)
)

# Manifest file
manifest <- fread(args$pheno_manifest)
manifest <- manifest[manifest$phenotype %in% args$phenotype]

gwas$N <- manifest$N

# Merge with manifest file
message("Writing output...")
fwrite(gwas, file = paste0(args$phenotype, "_InputToMunge.txt"), sep = "\t")
message("Writing output...done!")
