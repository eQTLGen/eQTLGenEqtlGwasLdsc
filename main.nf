#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2


def helpmessage() {

log.info"""

EqtlGwasLdsc v${workflow.manifest.version}"
===========================================================
Pipeline for running LDSC genetic correlation analyses analyses (https://www.nature.com/articles/s41467-020-20885-8) between eQTLGen eQTL summary statistics and GWAS summary statistics.

Usage:

nextflow run main.nf 
--eqtl_files \
--allele_info \
--ldsc_ref \
--gtf \
--OutputDir

Mandatory arguments:
--eqtl_files                eQTLGen parquet dataset.
--gwas_files                GWAS parquet dataset. Needs to be harmonised to follow same format as eQTL dataset.
--allele_info               Parquet file with alleles and SNP positions for eQTL dataset.
--snplist                   File with HapMap3 variants matching to eQTL dataset.
--ldsc_folder               Folder with LDSC scripts (https://github.com/bulik/ldsc).
--ldsc_ref                  Folder with LDSC reference files, as shared by LDSC authors.
--gwas_manifest             GWAS manifest file.
--remove_eqtl               Whether to remove significant eQTLs from the input eQTL file (yes/no). Default is no.
--window                    Window for removing variants for eQTL peak. Defaults to 1000000.
--pthresh                   P-value threshold for declaring significant eQTL. Defaults to 5e-8.

Optional arguments:
--OutputDir                 Output directory. Defaults to "results".
--i2_thresh                 Heterogeneity threshold. Defaults to 100 (<=100%).
--gene_filter               File to filter the genes included to the analysis. File with ENSG IDs, no header. Defaults that no filtering done.

""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

// Default parameters
params.OutputDir = 'results'
params.i2_thresh = 100
params.gene_filter = 'data/help_input.txt'
params.remove_eqtl = 'no'
params.window = 1000000
params.pthresh = 5e-8

//Show parameter values
log.info """=======================================================
eQTLGen eQTL-GWAS colocalisation pipeline v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Version']                         = workflow.manifest.version
summary['Current user']                             = "$USER"
summary['Current home']                             = "$HOME"
summary['Current path']                             = "$PWD"
summary['Working dir']                              = workflow.workDir
summary['Script dir']                               = workflow.projectDir
summary['Config Profile']                           = workflow.profile
summary['Container Engine']                         = workflow.containerEngine
if(workflow.containerEngine) summary['Container']   = workflow.container
summary['Output directory']                         = params.OutputDir
summary['eQTL folder']                              = params.eqtl_files
summary['GWAS folder']                              = params.gwas_files
summary['GWAS manifest file']                       = params.gwas_manifest
summary['LDSC folder']                              = params.ldsc_folder
summary['LDSC reference']                           = params.ldsc_ref
summary['Allele info file']                         = params.allele_info
summary['HapMap3 SNP list file']                    = params.snplist
summary['I2 threshold']                             = params.i2_thresh
summary['Gene filter']                              = params.gene_filter
summary['Remove eQTL']                              = params.remove_eqtl
summary['eQTL window']                              = params.window
summary['P threshold']                              = params.pthresh
summary['Gene filter']                              = params.gene_filter

// import modules
include { PREPAREGWAS; MUNGEGWAS; LDSC; COLLECTGWAS; PREPAREEQTL; PrepareGwas; MungeGwas; PrepareEqtl; Ldsc; CollectGwas } from './modules/EqtlGwasLdsc.nf'

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

// Define argument channels
// Get eQTL channel
eqtl_ch = Channel.fromPath(params.eqtl_files, type: 'dir').ifEmpty { exit 1, "eQTL files not found!" }
gwas_ch = Channel.fromPath(params.gwas_files, type: 'dir').ifEmpty { exit 1, "GWAS files not found!" }
gwas_manifest_ch = Channel.fromPath(params.gwas_manifest, type: 'file').ifEmpty { exit 1, "GWAS manifest file not found!" }
allele_ch = Channel.fromPath(params.allele_info).ifEmpty { exit 1, "eQTLGen reference not found!" }
snplist_ch = Channel.fromPath(params.snplist).ifEmpty { exit 1, "SNP list not found!" }
ldsc_folder_ch = Channel.fromPath(params.ldsc_folder).ifEmpty { exit 1, "LDSC scripts folder not found!" }
ldsc_ref_ch = Channel.fromPath(params.ldsc_ref).ifEmpty { exit 1, "LDSC reference folder not found!" }
filter_ch = Channel.fromPath(params.gene_filter)

// Get phenotype names
gwas_id_ch = Channel
    .fromPath("${params.gwas_files}/gwas_id*/*")
    .map { path -> path.toString().replaceAll(/.*gwas_id=/, '').replaceAll(/\/.*/, '') }.unique()

// Get gene names
gene_id_ch = Channel
    .fromPath("${params.eqtl_files}/phenotype*/*")
    .map { path -> path.toString().replaceAll(/.*phenotype=/, '').replaceAll(/\/.*/, '') }.unique()

gwas_input_ch = gwas_ch.combine(gwas_manifest_ch).combine(allele_ch).combine(snplist_ch).combine(ldsc_folder_ch).combine(ldsc_ref_ch).combine(gwas_id_ch)

i2_thresh = Channel.value(params.i2_thresh)
rm_eqtl = Channel.value(params.remove_eqtl)
window = Channel.value(params.window)
pthresh = Channel.value(params.pthresh)

workflow {
        PREPAREGWAS(gwas_input_ch)
        MUNGEGWAS(PREPAREGWAS.out.combine(ldsc_folder_ch).combine(ldsc_ref_ch))
        
        collectgwas_input_ch = MUNGEGWAS.out.map { it[0] }.collect()
        COLLECTGWAS(collectgwas_input_ch)

        prepareeqtl_input_ch = eqtl_ch.combine(gene_id_ch).combine(allele_ch).combine(snplist_ch).combine(i2_thresh)

        PREPAREEQTL(prepareeqtl_input_ch.combine(rm_eqtl).combine(window).combine(pthresh))
        LDSC(PREPAREEQTL.out.combine(COLLECTGWAS.out).combine(ldsc_folder_ch).combine(ldsc_ref_ch))
        LDSC.out.collectFile(name: 'EqtlGwasLdscResults.txt', keepHeader: true, sort: true, storeDir: "${params.OutputDir}")
        }


workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
