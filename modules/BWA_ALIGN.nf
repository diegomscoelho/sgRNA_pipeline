process BWA_ALIGN {
    label 'BWA'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path ref_fasta
    path fasta

    output:
    path "${fasta.baseName}.sam", emit: sam

    script:
    """
    set -euxo pipefail

    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`;

    bwa aln \$INDEX ${fasta} > reads.sai;
    
    bwa samse \$INDEX reads.sai ${fasta} > ${fasta.baseName}.sam;

    """
}