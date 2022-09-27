#!/usr/bin/env nextflow
 
/*
 * The following pipeline parameters specify the refence genomes
 * and read pairs and can be provided as command line options
 */
params.cram = "/HG00512.alt_bwamem_GRCh38DH.20150724.CHS.sv_7kb_mate.cram"
params.ref = "/tools/GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.outdir = "results"
params.design = "/workspaces/ConsensuSV-nextflow-MEI/design.csv"
workflow {
    files = Channel.fromPath(params.design).splitCsv()
    UNCRAM(files)
    INDEX(UNCRAM.out.bam)
    DELLY(UNCRAM.out.bam, UNCRAM.out.sample, INDEX.out)
    BreakDancer(UNCRAM.out.bam, UNCRAM.out.sample, INDEX.out)
    CNVNator(UNCRAM.out.bam, UNCRAM.out.sample, INDEX.out)
}
all_chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
all_chromosomes_space = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
all_chromosomes_num = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
process UNCRAM {
    tag "Uncraming files"
 
    input:
    path cram
 
    output:
    path 'output.bam', emit: bam
    val cram.simpleName, emit:sample
 
    script:
    """
    samtools view -b  -T ${params.ref} -o output.bam $cram
    """
}
 

process INDEX {
    tag "Indexing files"
 
    input:
    path bam
 
    output:
    path "${bam}.bai"
 
    script:
    """
    samtools index $bam
    """
} 

process DELLY {
    tag "Calling DELLY"
 
    publishDir "${params.outdir}/${sample}/vcf"

    input:
    path bam
    val sample
    path bai
 
    output:
    path "delly.vcf"

    script:
    """
    delly_v0.8.7_linux_x86_64bit call -o inter_file -g ${params.ref} $bam && bcftools view inter_file > delly.vcf
    """
}

process BreakDancer {
    conda '/tools/anaconda/envs/breakseq'

    tag "Calling BreakDancer"
 
    publishDir "${params.outdir}/${sample}/vcf"

    input:
    path bam
    val sample
    path bai
 
    output:
    path "breakseq.vcf"

    script:
    """
    bam2cfg.pl $bam > config
    breakdancer-max config > inter_file
    python /tools/breakdancer-master/bin/breakdancer2vcf.py < inter_file > breakseq.vcf
    """
}

process TARDIS {
    tag "Calling Tardis"
 
    publishDir "${params.outdir}/${sample}/vcf"

    input:
    path bam
    val sample
    path bai
 
    output:
    path "delly.vcf"

    script:
    """
    tardis -i ${bam} --ref ${params.ref} --sonic /tools/GRCh38_1kg.sonic --out tardis.vcf --first-chr 0 --last-chr 24
    """
    
}

process CNVNator {
    tag "Calling CNVNator"
 
    publishDir "${params.outdir}/${sample}/vcf"

    input:
    path bam
    val sample
    path bai
 
    output:
    path "cnvnator.vcf"

    script:
    """
    cnvnator -root cnv.root -tree $bam -chrom $all_chromosomes_num X Y
    cnvnator -root cnv.root -his 1000 -fasta ${params.ref}
    cnvnator -root cnv.root -stat 1000
    cnvnator -root cnv.root -partition 1000
    cnvnator -root cnv.root -call 1000 > cnv.out
    cnvnator2VCF.pl -prefix $sample -reference GRCh38 cnv.out . > cnvnator.vcf
    """
}

process BreakSeq {
    tag "Calling BreakSeq"
 
    publishDir "${params.outdir}/${sample}/vcf"

    input:
    path bam
    val sample
    path bai
 
    output:
    path "breakseq.vcf"

    script:
    """
    run_breakseq2.py --reference ${params.ref} --bams $bam --work breakseq/ --bwa /tools/bwa-0.7.17/bwa --samtools /tools/samtools-0.1.19/samtools --bplib_gff /tools/breakseq2_bplib_20150129_chr.gff --nthreads 4 --sample $sample --chromosomes $all_chromosomes_space
    gunzip breakseq/breakseq.vcf.gz
    """
}