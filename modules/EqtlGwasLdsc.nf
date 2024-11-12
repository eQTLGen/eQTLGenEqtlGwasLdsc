#!/bin/bash nextflow


process PrepareGwas {
    scratch true

    publishDir "${params.OutputDir}", mode: 'copy', overwrite: true, pattern: "*.log"
    container 'quay.io/urmovosa/pqtlvseqtl:v0.1'

    input:
        tuple path(gwas), path(gwas_manifest), path(ref), path(snplist), path(ldsc_folder), path(ldsc_ref), val(phenotype)

    output:
        tuple file("*_InputToMunge.txt"), val(phenotype)

    script:
        """
        ##############################################################
        # Convert to LDSC .txt format and filter to HapMap3 variants #
        ##############################################################

        PrepareGwasSumstats.R \
        --gwas_folder ${gwas} \
        --pheno_manifest ${gwas_manifest} \
        --snplist ${snplist} \
        --snpref ${ref} \
        --phenotype ${phenotype} \
        --output ${phenotype}_InputToMunge.txt
        """
}

process MungeGwas {
    scratch true

    publishDir "${params.OutputDir}", mode: 'copy', overwrite: true, pattern: "*.log"
    container 'manninglab/ldsc'
    input:
        tuple path(gwas), val(phenotype), path(ldsc_folder), path(ldsc_ref)

    output:
        tuple path("*.sumstats.gz"), path("*.log")

    script:
        """
        #########
        # Munge #
        #########

        python2 ${ldsc_folder}/munge_sumstats.py \
        --sumstats ${phenotype}_InputToMunge.txt \
        --out ${phenotype} \
        --merge-alleles ${ldsc_ref}/w_hm3.snplist

        """
}

process CollectGwas {
    scratch true

    input:
        path(gwas)

    output:
        path("gwas_folder")

    script:
        """
        mkdir gwas_folder
        mv *sumstats.gz gwas_folder
        """
}


process PrepareEqtl {
    scratch true
    container 'quay.io/urmovosa/pqtlvseqtl:v0.1'

    input:
        tuple path(eqtl), val(gene), path(ref), path(snplist), val(i2), val(analysis_type), val(window), val(pthresh)

    output:
        tuple path("*_InputToMunge.txt.gz"), val(gene)

    script:
        """

        ###########################################################################
        # Convert eQTL sumstat to LDSC .txt format and filter to HapMap3 variants #
        ###########################################################################

        PrepareEqtlSumstats.R \
            --eqtl_folder ${eqtl} \
            --snplist ${snplist} \
            --snpref ${ref} \
            --i2_thresh ${i2} \
            --gene ${gene} \
            --remove_eqtl ${analysis_type} \
            --window ${window}\
            --p_thresh ${pthresh}

        """
}



process Ldsc {

    container 'manninglab/ldsc'
    input:
        tuple path(eqtl), val(gene), path(gwas_sumstats), path(ldsc_folder), path(ldsc_ref)

    output:
        path("*vsAll.txt")

    script:
        """
        ###################################
        # Get the string of GWAS sumstats #
        ###################################
        
        gwas_sumstats_string=\$(ls -d -1 ${gwas_sumstats}/*sumstats.gz | tr '\\n' ',' | sed 's/,\$//')

        #########
        # Munge #
        #########

        zcat ${gene}_InputToMunge.txt.gz > ${gene}_InputToMunge.txt

        python2 ${ldsc_folder}/munge_sumstats.py \
            --sumstats ${gene}_InputToMunge.txt \
            --out ${gene} \
            --merge-alleles ${ldsc_ref}/w_hm3.snplist

        rm ${gene}_InputToMunge.txt

        ################################
        # Run LDSC genetic correlation #
        ################################

        ${ldsc_folder}/ldsc.py \
            --rg ${gene}.sumstats.gz,\${gwas_sumstats_string} \
            --ref-ld-chr ${ldsc_ref}/eur_w_ld_chr/ \
            --w-ld-chr ${ldsc_ref}/eur_w_ld_chr/ \
            --out ${gene}vsAll
        
        ##########################################
        # Extract results from the LDSC log file #
        ##########################################

        # Define start of table
        start_line=\$(grep -n "Summary of Genetic Correlation Results" ${gene}vsAll.log | cut -d: -f1)
        
        # Extract table and save to new file
        awk -v start="\$start_line" '
        NR > start {
            # Stop if we encounter a blank line or "Analysis finished at"
            if (/^[[:space:]]*\$/ || /Analysis finished at/) exit
            # Only print non-blank lines
            if (NF) print
        }
    ' ${gene}vsAll.log > ${gene}vsAll.txt

        """
}

workflow PREPAREGWAS {
    take:
        data

    main:
        loci_ch = PrepareGwas(data)
        
    emit:
        loci_ch

}

workflow MUNGEGWAS {
    take:
        data

    main:
        loci_ch = MungeGwas(data)
        
    emit:
        loci_ch

}

workflow COLLECTGWAS {
    take:
        data

    main:
        loci_ch = CollectGwas(data)
        
    emit:
        loci_ch

}

workflow PREPAREEQTL {
    take:
        data

    main:
        loci_ch = PrepareEqtl(data)
        
    emit:
        loci_ch

}

workflow LDSC {
    take:
        data

    main:
        ldsc_output_ch = Ldsc(data)

    emit:
        ldsc_output_ch
}
