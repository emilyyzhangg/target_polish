#!/usr/bin/env nextflow
params.fasta = "$projectDir/data/draft_assembly.fa"
params.reads = "$projectDir/data/long_reads.fq"
params.bed = "$projectDir/data/polishing_coords.bed"
params.prefix = "polished_assembly"
params.output_file = "polished_assembly.fa"
params.length = 15
params.s = 100
params.x = 150
params.max_threads = 12
params.benchmark = false
params.sensitive = true
params.script_path = "$projectDir/scripts"
params.outputDir = "output"

log.info """\
    P I P E L I N E
    ===================================
    fasta	: ${params.fasta}
    reads	: ${params.reads}
    bed	: ${params.bed}
    """
    .stripIndent()


process EXTRACT_SEQ {
    errorStrategy 'retry', maxRetries: 2
    tag "Extract regions to polish from FASTA file"

    input:
    path fasta
    path bed 

    output:
    path "${params.prefix}.gaps.fa"

    script:
    """
    python ${params.script_path}/extract-seq.py -f $fasta -b $bed -l ${params.length} -o ${params.prefix}.gaps.fa 
    """
}

process RUN_MINIMAP2 {
    errorStrategy 'retry', maxRetries: 2
    tag "Run minimap2 on FASTA file"

    input:
    path fasta
    path reads

    output:
    path "${params.prefix}.unpolished.minimap2.paf"

    script:
    """
    minimap2 -t ${params.max_threads} $fasta $reads > ${params.prefix}.unpolished.minimap2.paf
    """
}

process RUN_GOLDPOLISH {
    errorStrategy 'retry', maxRetries: 2
    tag "Running GoldPolish on extracted seqs"

    input:
    path reads
    path mapping 
    path gaps 

    output:
    path "${params.prefix}.gaps.goldpolished.fa" 

    script:
    """
    goldpolish --mappings $mapping -s ${params.s} -x ${params.x} -t ${params.max_threads} $gaps $reads ${params.prefix}.gaps.goldpolished.fa
    """
}

process RUN_POST_PROCESSING {
    errorStrategy 'retry', maxRetries: 2
    tag "Re-inserting polished subseqs into FASTA file"
    publishDir params.outputDir, mode:'copy'
 
    input:
    path fasta
    path polished_gaps 

    output:
    path "${params.output_file}"

    script:
    """
    python ${params.script_path}/post-processing.py -f $fasta -g $polished_gaps -o ${params.output_file}
    """
}

workflow {
    fasta = file(params.fasta)
    reads = file(params.reads)
    bed = file(params.bed)

    mapping = RUN_MINIMAP2(fasta, reads)

    gaps = EXTRACT_SEQ(fasta, bed)

    polished_gaps = RUN_GOLDPOLISH(reads, mapping, gaps)

    RUN_POST_PROCESSING(fasta, polished_gaps)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! The ouput is in --> $projectDir/$params.outputDir\n" : "Oops .. something went wrong" )
}
