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

    

    """
}