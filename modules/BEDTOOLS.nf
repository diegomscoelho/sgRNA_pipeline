process SAMTOBED {
    label 'BEDTOOLS'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path sam

    output:
    path "${sam.baseName}.bed"

    script:
    """
    set -euxo pipefail

    bedtools bamtobed -i ${sam} > ${sam.baseName}.bed

    """
}