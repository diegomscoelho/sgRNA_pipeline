process BWA_INDEX {
    label 'BWA'

    input:
    path fasta

    output:
    path bwa, emit: index

    script:
    """
    set -euxo pipefail
    
    mkdir bwa
    bwa index -p bwa/${fasta.baseName} $fasta

    """

}