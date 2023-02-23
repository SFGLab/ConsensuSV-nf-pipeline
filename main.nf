#!/usr/bin/env nextflow
 
/*
 * The following pipeline parameters specify the refence genomes
 * and read pairs and can be provided as command line options
 */
 
params.ref = "/tools/GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.outdir = "results"
params.design = "/workspaces/ConsensuSV-nextflow-MEI/design.csv"
params.threads = 4
params.mem = 4
workflow {
    files = Channel.fromPath(params.design).splitCsv()
    UNCRAM(files)
    INDEX(UNCRAM.out.bam)
    DELLY(UNCRAM.out.bam, UNCRAM.out.sample, INDEX.out)
    BreakDancer(UNCRAM.out.bam, UNCRAM.out.sample, INDEX.out)
    TARDIS(UNCRAM.out.bam, UNCRAM.out.sample, INDEX.out)
    CNVNator(UNCRAM.out.bam, UNCRAM.out.sample, INDEX.out)
    BreakSeq(UNCRAM.out.bam, UNCRAM.out.sample, INDEX.out)
    Manta(UNCRAM.out.bam, UNCRAM.out.sample, INDEX.out)
    Lumpy(UNCRAM.out.bam, UNCRAM.out.sample, INDEX.out)
    Whamg(UNCRAM.out.bam, UNCRAM.out.sample, INDEX.out)
    //FAKE(UNCRAM.out.sample, DELLY.out, BreakDancer.out, TARDIS.out, CNVNator.out, Manta.out, Lumpy.out, Whamg.out)
    ConsensuSV(UNCRAM.out.sample, DELLY.out, BreakDancer.out, TARDIS.out, CNVNator.out, Manta.out, Lumpy.out, Whamg.out)
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
    if [[ $cram == *.cram ]]; then
        samtools view -b  -T ${params.ref} -o output.bam $cram
    elif [[ $cram == *.bam ]]; then
        cp $cram output.bam
    else
        exit 5
    fi
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
 
    publishDir "${params.outdir}/vcfs/${sample}"

    input:
    path bam
    val sample
    path bai
 
    output:
    path "delly.vcf"

    script:
    """
    delly_v1.1.6_linux_x86_64bit call -o delly.vcf -g ${params.ref} $bam
    """
}

process BreakDancer {
    conda '/tools/anaconda/envs/breakseq'

    tag "Calling BreakDancer"
 
    publishDir "${params.outdir}/vcfs/${sample}"

    input:
    path bam
    val sample
    path bai
 
    output:
    path "breakdancer.vcf"

    script:
    """
    bam2cfg.pl $bam > config
    breakdancer-max config > inter_file
    python /tools/breakdancer-master/bin/breakdancer2vcf.py < inter_file > breakdancer.vcf
    """
}

process TARDIS {
    tag "Calling Tardis"

    time { 1.minute }
    errorStrategy "retry"
    maxRetries 1 // total of 2 - if tardis is longer than 1h, then it has this weird thing that calculates very long
 
    publishDir "${params.outdir}/vcfs/${sample}"

    input:
    path bam
    val sample
    path bai
 
    output:
    path "tardis.vcf"

    script:
    """
    if [ ${task.attempt} -lt 2 ] 
    then
        tardis -i ${bam} --ref ${params.ref} --sonic /tools/GRCh38_1kg.sonic --out tardis.vcf --first-chr 0 --last-chr 24
    else
        cat /tools/ConsensuSV-core/header > tardis.vcf
    fi
    """
    
}

process CNVNator {
    tag "Calling CNVNator"
 
    publishDir "${params.outdir}/vcfs/${sample}"

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
    conda '/tools/anaconda/envs/breakseq'

    tag "Calling BreakSeq"
 
    publishDir "${params.outdir}/vcfs/${sample}"

    input:
    path bam
    val sample
    path bai
 
    output:
    path "breakseq.vcf"

    script:
    """
    run_breakseq2.py --reference ${params.ref} --bams $bam --work breakseq/ --bwa /tools/bwa-0.7.17/bwa --samtools /tools/samtools-0.1.19/samtools --bplib_gff /tools/breakseq2_bplib_20150129_chr.gff --nthreads ${params.threads} --sample $sample --chromosomes $all_chromosomes_space
    gunzip breakseq/breakseq.vcf.gz
    mv breakseq/breakseq.vcf .
    """
}

process Manta {
    conda '/tools/anaconda/envs/breakseq'
    
    errorStrategy "retry"
    maxRetries 2 // total of 3

    tag "Calling Manta"
 
    publishDir "${params.outdir}/vcfs/${sample}"

    input:
    path bam
    val sample
    path bai
 
    output:
    path "manta.vcf"

    script:
    """
    if [ ${task.attempt} -lt 3 ] 
    then
        configManta.py --bam $bam --referenceFasta ${params.ref} --runDir Manta --callRegions /tools/grch38.bed.gz
        Manta/runWorkflow.py -j ${params.threads} -g ${params.mem}
        gunzip Manta/results/variants/diploidSV.vcf.gz
        mv Manta/results/variants/diploidSV.vcf manta.vcf
    else
        cat /tools/ConsensuSV-core/header > manta.vcf
    fi
    """
}

process Lumpy {
    conda '/tools/anaconda/envs/breakseq'

    tag "Calling Lumpy"
 
    publishDir "${params.outdir}/vcfs/${sample}"

    input:
    path bam
    val sample
    path bai
 
    output:
    path "lumpy.vcf"

    script:
    """
    samtools view -b -F 1294 -@ ${params.threads} $bam > lumpy.discordants.unsorted.bam
    samtools sort -o lumpy.discordants.bam -@ ${params.threads} -m ${params.mem}G lumpy.discordants.unsorted.bam
    samtools view -h $bam | /tools/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > lumpy.splitters.unsorted.bam
    samtools sort -o lumpy.splitters.bam -@ ${params.threads} -m ${params.mem}G lumpy.splitters.unsorted.bam
    lumpyexpress -B $bam -S lumpy.splitters.bam -D lumpy.discordants.bam -R ${params.ref} -T lumpy -o lumpy.vcf
    """
}

process Whamg {

    tag "Calling Whamg"
 
    publishDir "${params.outdir}/vcfs/${sample}"

    input:
    path bam
    val sample
    path bai
 
    output:
    path "whamg.vcf"

    script:
    """
    whamg -c $all_chromosomes -a ${params.ref} -f $bam | perl /tools/wham/utils/filtWhamG.pl > whamg.vcf  2> wham.err
    """
}

process SVelter {
    conda '/tools/anaconda/envs/breakseq'

    tag "Calling SVelter"
 
    publishDir "${params.outdir}/vcfs/${sample}"

    input:
    path bam
    val sample
    path bai
 
    output:
    path "svelter.vcf"

    script:
    """
    svelter.py Setup --reference ${params.ref} --workdir SVelter --support /tools/svelter/Support/GRCh38/
    svelter.py --sample $bam --workdir SVelter --chromosome $all_chromosomes
    """
}

process ConsensuSV {

    tag "Calling ConsensuSV"
 
    publishDir "${params.outdir}/consensusv"

    input:
    val sample
    path vcf_DELLY
    path vcf_BreakDancer
    path vcf_TARDIS
    path vcf_CNVNator
    path vcf_Manta
    path vcf_Lumpy
    path vcf_Whamg

    output:
    path "${sample}.vcf"

    script:
    """
    mkdir -p consensus_working_results
    mkdir -p consensus_working_results/$sample
    cp $vcf_DELLY $vcf_BreakDancer $vcf_TARDIS $vcf_CNVNator $vcf_Manta $vcf_Lumpy $vcf_Whamg consensus_working_results/$sample
    python -u /tools/ConsensuSV-core/main.py -of \$PWD/ -f \$PWD/consensus_working_results/ -s $sample -c breakdancer,breakseq,cnvnator,delly,lumpy,manta,tardis,whamg -mod /tools/ConsensuSV-core/pretrained_1000g_illumina.model
    mv consensuSV__${sample}.vcf ${sample}.vcf
    """
}