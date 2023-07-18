process PARSE_GTF {
    label 'R'

    input:
    path gtf

    output:
    path "${gtf.simpleName}.rds", emit: gtf

    script:
    template 'parse_gtf.R'

}