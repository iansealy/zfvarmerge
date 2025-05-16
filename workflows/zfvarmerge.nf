/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_GVCF            } from '../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_VCF             } from '../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_CONCAT        } from '../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_STATS         } from '../modules/nf-core/bcftools/stats/main'
include { GATK4_GENOMICSDBIMPORT } from '../modules/local/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS    } from '../modules/nf-core/gatk4/genotypegvcfs/main'
include { BCFTOOLS_MAKEBEDS as BCFTOOLS_MAKEBEDS_FREEBAYES  } from '../modules/local/bcftools/makebeds/main'
include { BCFTOOLS_MAKEBEDS as BCFTOOLS_MAKEBEDS_BCFTOOLS   } from '../modules/local/bcftools/makebeds/main'
include { GNU_SORT as GNU_SORT_FREEBAYES_SNPS               } from '../modules/nf-core/gnu/sort/main'
include { GNU_SORT as GNU_SORT_FREEBAYES_INDELS             } from '../modules/nf-core/gnu/sort/main'
include { GNU_SORT as GNU_SORT_BCFTOOLS_SNPS                } from '../modules/nf-core/gnu/sort/main'
include { GNU_SORT as GNU_SORT_BCFTOOLS_INDELS              } from '../modules/nf-core/gnu/sort/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_FREEBAYES_SNPS   } from '../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_FREEBAYES_INDELS } from '../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_BCFTOOLS_SNPS    } from '../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_BCFTOOLS_INDELS  } from '../modules/nf-core/bedtools/merge/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_zfvarmerge_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ZFVARMERGE {

    take:
    ch_gatk_vcf           // channel: samplesheet read in from --input
    ch_freebayes_vcf      // channel: samplesheet read in from --input
    ch_bcftools_vcf       // channel: samplesheet read in from --input
    ch_fasta              // channel: [ val(meta), [ fasta ] ]
    ch_fasta_fai          // channel: [ val(meta), [ fai ] ]
    ch_fasta_dict         // channel: [ val(meta), [ dict ] ]
    ch_genome_bed         // channel: [ val(meta), [ bed ] ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: BCFtools index gVCFs
    //
    BCFTOOLS_INDEX_GVCF (
        ch_gatk_vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_GVCF.out.versions)

    //
    // MODULE: GATK GenomicsDBImport
    //
    ch_gatk_vcfs = ch_gatk_vcf.map { meta, vcf -> vcf }.collect()
    ch_gatk_tbis = BCFTOOLS_INDEX_GVCF.out.tbi.map { meta, tbi -> tbi }.collect()
    ch_genomicsdbimport_input = ch_genome_bed
        .map{ interval -> [ [ id:'genomicsdb', order:interval[1].baseName ], interval[1], [], "${params.genomicsdb}/genomicsdb.${interval[1].baseName}" ] }
    GATK4_GENOMICSDBIMPORT (
        ch_genomicsdbimport_input,
        ch_gatk_vcfs,
        ch_gatk_tbis,
        false // not providing sample name map file
    )
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)

    //
    // MODULE: GATK GenotypeGVCFs
    //
    ch_genomicsdb = GATK4_GENOMICSDBIMPORT.out.genomicsdb
        .map{ meta, genomicsdb -> [ meta, genomicsdb, [], [], [] ] }
    GATK4_GENOTYPEGVCFS (
        ch_genomicsdb,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        [[], []], // not providing dbSNP VCF file
        [[], []]  // not providing dbSNP VCF index file
    )
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)

    //
    // MODULE: BCFtools concat
    //
    ch_vcfs_tbis = GATK4_GENOTYPEGVCFS.out.vcf.map{ meta, vcf -> [meta, vcf, meta.order] }
        .map { it[0].remove("order"); it }
        .groupTuple()
        .map{ meta, vcfs, order -> [[id:"gatk"], order.withIndex().sort().collect { vcfs[it[1]] }, []] }
    BCFTOOLS_CONCAT (
        ch_vcfs_tbis
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    //
    // MODULE: BCFtools index VCF
    //
    BCFTOOLS_INDEX_VCF (
        BCFTOOLS_CONCAT.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_VCF.out.versions)

    //
    // MODULE: BCFtools stats
    //
    ch_vcf_tbi = BCFTOOLS_CONCAT.out.vcf.join(BCFTOOLS_INDEX_VCF.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .map{ meta, vcf, tbi -> [meta, vcf, tbi] }
    BCFTOOLS_STATS (
        ch_vcf_tbi,
        [[], []], // no need for regions file
        [[], []], // no need for targets file
        [[], []], // no need for samples file
        [[], []], // no need for exons file
        [[], []]  // no need for fasta file
    )
    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS.out.stats.collect{it[1]})
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)

    //
    // MODULE: Make freebayes BED files
    //
    BCFTOOLS_MAKEBEDS_FREEBAYES (
        ch_freebayes_vcf,
        ch_fasta,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MAKEBEDS_FREEBAYES.out.versions)

    //
    // MODULE: Sort freebayes SNP BED files
    //
    ch_freebayes_snp_beds = BCFTOOLS_MAKEBEDS_FREEBAYES.out.snpbed
        .map{ bed -> [ [id:"freebayes.snps.unmerged"], bed ] }
        .groupTuple()
    GNU_SORT_FREEBAYES_SNPS (
        ch_freebayes_snp_beds
    )
    ch_versions = ch_versions.mix(GNU_SORT_FREEBAYES_SNPS.out.versions)

    //
    // MODULE: Merge freebayes SNP BED files
    //
    ch_freebayes_snp_bed = GNU_SORT_FREEBAYES_SNPS.out.sorted
        .map{ meta, bed -> [ [id:"freebayes.snps"], bed ] }
    BEDTOOLS_MERGE_FREEBAYES_SNPS (
        ch_freebayes_snp_bed
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE_FREEBAYES_SNPS.out.versions)

    //
    // MODULE: Sort freebayes indel BED files
    //
    ch_freebayes_indel_beds = BCFTOOLS_MAKEBEDS_FREEBAYES.out.indelbed
        .map{ bed -> [ [id:"freebayes.indels.unmerged"], bed ] }
        .groupTuple()
    GNU_SORT_FREEBAYES_INDELS (
        ch_freebayes_indel_beds
    )
    ch_versions = ch_versions.mix(GNU_SORT_FREEBAYES_INDELS.out.versions)

    //
    // MODULE: Merge freebayes indel BED files
    //
    ch_freebayes_indel_bed = GNU_SORT_FREEBAYES_INDELS.out.sorted
        .map{ meta, bed -> [ [id:"freebayes.indels"], bed ] }
    BEDTOOLS_MERGE_FREEBAYES_INDELS (
        ch_freebayes_indel_bed
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE_FREEBAYES_INDELS.out.versions)

    //
    // MODULE: Make BCFtools BED files
    //
    BCFTOOLS_MAKEBEDS_BCFTOOLS (
        ch_bcftools_vcf,
        ch_fasta,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MAKEBEDS_BCFTOOLS.out.versions)

    //
    // MODULE: Sort BCFtools SNP BED files
    //
    ch_bcftools_snp_beds = BCFTOOLS_MAKEBEDS_BCFTOOLS.out.snpbed
        .map{ bed -> [ [id:"bcftools.snps.unmerged"], bed ] }
        .groupTuple()
    GNU_SORT_BCFTOOLS_SNPS (
        ch_bcftools_snp_beds
    )
    ch_versions = ch_versions.mix(GNU_SORT_BCFTOOLS_SNPS.out.versions)

    //
    // MODULE: Merge BCFtools SNP BED files
    //
    ch_bcftools_snp_bed = GNU_SORT_BCFTOOLS_SNPS.out.sorted
        .map{ meta, bed -> [ [id:"bcftools.snps"], bed ] }
    BEDTOOLS_MERGE_BCFTOOLS_SNPS (
        ch_bcftools_snp_bed
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE_BCFTOOLS_SNPS.out.versions)

    //
    // MODULE: Sort BCFtools indel BED files
    //
    ch_bcftools_indel_beds = BCFTOOLS_MAKEBEDS_BCFTOOLS.out.indelbed
        .map{ bed -> [ [id:"bcftools.indels.unmerged"], bed ] }
        .groupTuple()
    GNU_SORT_BCFTOOLS_INDELS (
        ch_bcftools_indel_beds
    )
    ch_versions = ch_versions.mix(GNU_SORT_BCFTOOLS_INDELS.out.versions)

    //
    // MODULE: Merge BCFtools indel BED files
    //
    ch_bcftools_indel_bed = GNU_SORT_BCFTOOLS_INDELS.out.sorted
        .map{ meta, bed -> [ [id:"bcftools.indels"], bed ] }
    BEDTOOLS_MERGE_BCFTOOLS_INDELS (
        ch_bcftools_indel_bed
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE_BCFTOOLS_INDELS.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'zfvarmerge_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_replace_names = params.multiqc_replace_names ?
        Channel.fromPath(params.multiqc_replace_names, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_sample_names  = params.multiqc_sample_names ?
        Channel.fromPath(params.multiqc_sample_names, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        ch_multiqc_replace_names.toList(),
        ch_multiqc_sample_names.toList()
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
