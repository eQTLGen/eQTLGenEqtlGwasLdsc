#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name="eQtlGwasLdsc"

# Here load needed system tools (Java 1.8 is required, one of singularity or anaconda - python 2.7 are needed,
# depending on the method for dependancy management)

module load jdk/16.0.1
module load openjdk/11.0.2
module load squashfs
module load singularity

set -f

nextflow_path=[path to Nextflow executable]

eqtl=[eQTLGen eQTL sumstats folder in parquet hive format]
sig=[File with significant eQTL effects]

gwas=[Folder with GWAS sumstats in parquet hive format]
gwas_manifest=[File with GWAS file annotations]

ref=[eQTLGen p2 SNP reference file in parquet format]
gtf=[ENSEMBL .gtf file]
output_folder=[Optional output folder name]

NXF_VER=23.04.1 ${nextflow_path}/nextflow run main.nf \
--sig_eqtls ${sig} \
--eqtl_files ${eqtl} \
--gwas_files ${gwas} \
--gwas_manifest ${gwas_manifest} \
--allele_info ${ref} \
--gtf ${gtf} \
--p_thresh 5e-8 \
--OutputDir ${output_folder} \
--minN_thresh 0 \
--maxN_thresh 0.5 \
--i2_thresh 40 \
--leadvar_window  1000000 \
--cis_window 1000000 \
--trans_window 1000000 \
-profile slurm,singularity \
-resume
