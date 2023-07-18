process SAMTOBED {
    label 'BEDTOOLS'

    input:
    path sam

    output:
    path "${sam.baseName}.bed", emit: bed

    script:
    """
    set -euxo pipefail

    bedtools bamtobed -i ${sam} > ${sam.baseName}.bed

    """
}