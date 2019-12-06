#!/usr/bin/env nextflow
 
/* 
 * Proof of concept Nextflow based fqscreen pipeline
 * 
 */ 

 
/*
 * Defines some parameters in order to specify 
 * read pairs by using the command line options
 */
params.reads = "./testNf/*_R{1,2}*.fastq.gz"
params.conf = "./fastq_screen.ht.conf"
params.outdir = 'results'
log.info """
         FASTQSCREEN   P I P E L I N E    
         =============================
         reads : ${params.reads}
         outdir: ${params.outdir}
         """
         .stripIndent()

 
/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 */
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs } 
 
 
/*
 * Step 1. remove optical duplicates from NovaSeq 
 */ 
process dedup {
    tag "$pair_id"
     
    input:
    set pair_id, file(reads) from read_pairs
 
    output:
    set pair_id, "${pair_id}.clump*.fq.gz" into totrim
 
    """
    export APP=/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app
    \$APP/BBmap/v38.34/clumpify.sh in1=${reads[0]} in2=${reads[1]} out1=${pair_id}.clump1.fq.gz out2=${pair_id}.clump2.fq.gz dupedist=12000 dedupe=t optical=t
    """
}
  
/*
 * Step 2. Trimming with cutadapt for nextera
 */
process trim {
    cpus 28
    tag "$pair_id"
    input:
    set pair_id, file("${pair_id}.clump*.fq.gz") from totrim
     
    output:
    set pair_id, file("${pair_id}.fq.gz") into tofqscreen
    set pair_id, file("${pair_id}.trimmed_1.fastq.gz") into tofqc
    file("*cutadapt.txt") into forMqc1
    """
    module load cutadapt
    cutadapt -A CTGTCTCTTATACACA -a CTGTCTCTTATACACA -G GATGTGTATAAGAGACAG -g GATGTGTATAAGAGACAG -o ${pair_id}.trimmed_1.fastq.gz -p ${pair_id}.trimmed_2.fastq.gz -j ${task.cpus} -O 6 -m 20 ${pair_id}.clump*.fq.gz > ${pair_id}.cutadapt.txt 
    cat ${pair_id}.trimmed_*.fastq.gz > ${pair_id}.fq.gz
    """
}
 
/*
 * Step 3. QC with fastqc
 */
process fastqc {
       
    input:
     set pair_id, file("${pair_id}.trimmed_1.fastq.gz") from tofqc
     
    output:
     file("fastqc_${pair_id}_logs") into forMqc2
 
    """
    mkdir fastqc_${pair_id}_logs
    module load fastqc
    fastqc -o fastqc_${pair_id}_logs -f fastq -q ${pair_id}.trimmed_1.fastq.gz
    """
}
/*
 * Step 4. QC and remove human reads with fastqscreen
 */
process fastqscreen {
    cpus 28
    publishDir params.outdir, mode: 'copy', pattern: '*tagged_filter*'  
    input:
     set pair_id, file("${pair_id}.fq.gz") from tofqscreen
    file config from file(params.conf) 
    output:
     file("${pair_id}_screen*") into forMqc3
     file("${pair_id}*tagged_filter*")  
 
    """
    module load bowtie2
    export APP=/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app
    \$APP/fastq_screen/v0.14.0/fastq_screen --conf ${config} --threads ${task.cpus} --tag --filter 0--------- ${pair_id}.fq.gz 
    """
}

/*
 * Step 5. QC summary with multiqc
 */


process multiqc {
    publishDir params.outdir, mode:'copy'

    input:
    file('./*') from forMqc3.mix(forMqc1,forMqc2).collect()

    output:
    file('multiqc_report.html') optional true

    script:
    """
    module load multiqc
    multiqc -v .
    """
}


workflow.onComplete { 
	println ( workflow.success ? "Done!Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
