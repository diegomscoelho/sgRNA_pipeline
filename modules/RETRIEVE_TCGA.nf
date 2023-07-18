process RETRIEVE_TCGA {
    label 'R'
    publishDir "${params.outdir}/${name}", mode: 'copy'

    input:
    path gene_names
    val name

    output:
    path "TCGA_retrieved.csv"

    script:
    template 'retrieve_TCGA.R'

}