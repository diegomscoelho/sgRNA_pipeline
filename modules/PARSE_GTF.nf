process PARSE_GTF {
    label 'R'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path gtf

    output:
    path "${gtf.simpleName}.rds"

    script:
    template 'parse_gtf.R'

}