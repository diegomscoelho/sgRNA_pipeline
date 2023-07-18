process RETRIEVE_TCGA {
    debug true
    label 'R'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path gene_names

    output:
    path "TCGA_retrieved.csv"

    script:
    template 'retrieve_TCGA.R'

}