/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GATK4_GENOMICSDBIMPORT } from '../modules/nf-core/gatk4/genomicsdbimport/main'
include { BCFTOOLS_INDEX         } from '../modules/nf-core/bcftools/index/main'
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
    // MODULE: BCFtools index
    //
    BCFTOOLS_INDEX (
        ch_gatk_vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions)

    ch_gatk_vcf_tbi = ch_gatk_vcf.join(BCFTOOLS_INDEX.out.tbi, failOnDuplicate: true, failOnMismatch: true)
    ch_gatk_vcfs = []
    ch_gatk_vcf_tbi.collect {it[1]}.subscribe(onNext: { it -> ch_gatk_vcfs = it })
    ch_gatk_tbis = []
    ch_gatk_vcf_tbi.collect {it[2]}.subscribe(onNext: { it -> ch_gatk_tbis = it })
    ch_genomicsdbimport_input = ch_genome_bed
        .map{ interval -> [ [ id:'genomicsdb', order:interval[1].baseName ], ch_gatk_vcfs, ch_gatk_tbis, interval[1], [], params.genomicsdb ? "${params.genomicsdb}/genomicsdb.${interval[1].baseName}" : [] ] }

    //
    // MODULE: GATK GenomicsDBImport
    //
    GATK4_GENOMICSDBIMPORT (
        ch_genomicsdbimport_input,
        false,                            // not getting interval list
        params.genomicsdb ? true : false, // whether updating existing workspace
        false                             // not providing sample name map file
    )
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)

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
