process ASSIGN_GENES_MATCH {
    debug true
    label 'R'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path gtf_rds
    path bed

    output:
    path "gene_name_match.rds", emit: rds

    script:
    template 'assign_genes_match.R'

}