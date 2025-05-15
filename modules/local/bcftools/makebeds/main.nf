process BCFTOOLS_MAKEBEDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5acacb55c52bec97c61fd34ffa8721fce82ce823005793592e2a80bf71632cd0/data':
        'community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11' }"

    input:
    tuple val(meta), path(vcf)    // channel: [ val(meta), [ vcf ] ]
    tuple val(meta2), path(fasta) // channel: [ val(meta), [ fasta ] ]
    tuple val(meta3), path(fai)   // channel: [ val(meta), [ fai ] ]

    output:
    path "*.snps.bed"   , emit: snpbed
    path "*.indels.bed" , emit: indelbed
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools +fill-tags $vcf -- -t AC | \\
    bcftools view -V indels --min-ac=1 - | \\
    bcftools norm -c x -f $fasta - | \\
    bcftools query -f '%CHROM\t%POS0\t%END\\n' - > ${prefix}.snps.bed

    bcftools +fill-tags $vcf -- -t AC | \\
    bcftools view -v indels - | \\
    bcftools norm -c x -f $fasta - | \\
    bcftools query -f '%CHROM\t%POS0\t%END\\n' - > ${prefix}.indels.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.snps.bed ${prefix}.indels.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
