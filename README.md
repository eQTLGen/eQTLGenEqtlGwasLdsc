# Pipeline for running colocalisation analysis between eQTLGen eQTL and UKBB GWAS summary statistics

**TBA**

This pipeline runs eQTL summary statistics from eQTLGen phase 2 project and runs genetic correlation (LDSC) analysis for every gene against set of prepared GWAS summary statistics.

## Usage information

### Requirements for the system

- Have access to HPC with multiple cores.
- Have Bash >=3.2 installed.
- Have Java >=11 installed.
- Have Slurm scheduler managing the jobs in the HPC.
- HPC has Singularity installed and running.

### Setup of the pipeline

You can either clone it by using git (if available in HPC):

`git clone https://github.com/eQTLGen/eQTLGenEqtlGwasLdsc.git`

Or just download this from the github download link and unzip.

### Required inputs

#### eQTLs

`--eqtl_files` Folder of genome-wide eQTL meta-analysis parquet files. Output of [MetaAnalysis](https://github.com/eQTLGen/MetaAnalysis) pipeline.

`--allele_info` eQTLGen SNP reference file in .parquet format, containing SNP ID, chr, pos, and alleles.

#### GWASs

`--gwas_files` Folder with formatted GWAS summary statistics. These need to be in parquet format (hive structure) and harmonised to match eQTLGen eQTL summary statistics (IDs matched, genomic positions in hg38 and effect directions harmonised to match eQTLGen alternative allele).

`--gwas_manifest` .

#### LDSC

`--ldsc_folder` Folder with LDSC, as pulled from official LDSC repo (https://github.com/bulik/ldsc).

#### Annotation

`--gtf` ENSEMBL .gtf file for annotating gene symbols and cis/trans effects. In eQTLGen p2 project we use [v106 version](https://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz).

`--ldsc_ref` Folder with LDSC reference files, as shared by LDSC authors.

### Optional inputs

`--OutputDir` Optional output folder path. Defaults to `results` inside the pipeline directory.

### Additional settings

`--minN_thresh` Minimal sample size threshold for including genes and variants to the analysis. Defaults to 0 (no filter).

`--maxN_thresh` Variant filter for including eQTL per-locus variants into the analysis. Defaults to 0.5, meaning that variants having sample size <50% of the maximal sample size in this locus are excluded from the analysis. This is because eQTLGen data entials meta-analysis: some variants are tested in only part of the cohorts. This filter is meant to partly address the potential issues with emerging from different power for different variants. 

`--i2_thresh` Meta-analysis heterogeneity I2 threshold. Defaults 40 (<40%).

### Running the pipeline

TBA

### Outputs

TBA

## Acknowledgements

This pipeline was written by Urmo Võsa and Robert Warmerdam.

LDSC is work of Brendan Bulik-Sullivan and Hilary Finucane

Repo: https://github.com/bulik/ldsc

Citations: 
[Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies. Nature Genetics, 2015](https://www.nature.com/articles/ng.3211)

[Bulik-Sullivan, B., et al. An Atlas of Genetic Correlations across Human Diseases and Traits. Nature Genetics, 2015](https://www.nature.com/articles/ng.3406)