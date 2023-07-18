/*
 * pipeline input parameters
 */

 // print usage
if (params.help) {
    log.info """\
        ===================================
        S  G  R  N  A   P I P E L I N E
        ===================================
        Usage:
            nextflow run diegomscoelho/sgRNA
        Options:
            --fasta FASTA,FASTQ                 'FASTA,FASTQ files inside `inputs` folder. Currently supports bgzipped file or not'
            --bwa_index DIRNAME                 'BWA index folder.'
            --ref_genome_fasta FASTA            'FASTA used to create BWA index.'
            --ref_gtf GTF,GFF3                  'GTF/GFF3 used to get gene_name'
        """
        .stripIndent()
    exit 1
}

params.fasta = "$projectDir/inputs/*{fa,fq,fastq,fasta}"
bwa_index_folder = params.bwa_index? file(params.bwa_index) : file("$projectDir/ref/bwa/")
params.outdir = "results"

log.info """\
    ===================================
      S  G  R  N  A   P I P E L I N E
    ===================================
    Fasta files            : ${params.fasta}
    """
    .stripIndent()


// Include all modules

include { BWA_INDEX } from "$baseDir/modules/BWA_INDEX.nf"
include { BWA_ALIGN } from "$baseDir/modules/BWA_ALIGN.nf"
include { SAMTOBED } from "$baseDir/modules/BEDTOOLS.nf"
include { PARSE_GTF } from "$baseDir/modules/PARSE_GTF.nf"

// Get all necessary files

if (!params.ref_gtf) {
    params.ref_gtf = file("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz")
}

workflow {

    fasta = Channel.fromPath(params.fasta)
    bwa_index = Channel.fromPath(bwa_index_folder)
    
    // If BWA_INDEX does not exist, automatically downloads Reference fasta and create index
    if (!bwa_index_folder.exists()) {
        ref_genome_fasta = params.ref_genome_fasta? file(params.ref_genome_fasta) : file("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz")
        bwa_index_ch = BWA_INDEX(ref_genome_fasta)
        bwa_index = bwa_index_ch.index
    }


    align_ch = BWA_ALIGN(bwa_index, fasta)
    samtobed_ch = SAMTOBED(align_ch.sam)

    parsed_gtf_ch = PARSE_GTF(params.ref_gtf)

}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir\n" : "Oops .. something went wrong" )
}