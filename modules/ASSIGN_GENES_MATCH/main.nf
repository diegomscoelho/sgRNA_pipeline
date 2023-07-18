process ASSIGN_GENES_MATCH {
    label 'R'
    publishDir "${params.outdir}/${bed.baseName}", mode: 'copy'

    input:
    path gtf_rds
    path bed

    output:
    path "gene_name_match.rds", emit: rds
    path "annotated.csv"
    val "${bed.baseName}", emit: name

    script:
    template 'assign_genes_match.R'

}