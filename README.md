# ConsensuSV - Nextflow Pipeline [WORK IN PROGRESS]

Table of Contents
=================

* [What is ConsensuSV?](#what-is-consensusv)
* [Citation](#citation)
* [Testing scenarios](#testing-scenarios)
* [Preparation of your samples](#preparation-of-your-samples)
* [Running the pipeline](#running-the-pipeline)
* [Output location](#output-location)
* [Pipeline control webservice](#pipeline-control-webservice)
* [Pipeline details](#pipeline-details)
* [Setup on NVIDIA DGX A100 systems](#setup-on-nvidia-dgx-a100-systems)
* [Setup using singularity](#setup-using-singularity)
* [Benchmark](#benchmark)
## What is ConsensuSV?

Automatised pipeline of ConsensuSV workflow - special edition using NextFlow. Going from cram/bam data, to output vcf files (structural variants, indels and SNPs). Easy to run, scalable solution for variant discovery.

Docker image: https://hub.docker.com/repository/docker/mateuszchilinski/consensusv-nf-pipeline

## Citation

If you use ConsensuSV in your research, we kidnly as you to cite the following publication:

The citation will be updated upon publication.

```
@article{Chilinski_ConsensuSVfrom_the_whole-genome_2022,
author = {Chiliński, Mateusz and Plewczynski, Dariusz},
doi = {10.1093/bioinformatics/btac709},
journal = {Bioinformatics},
title = {{ConsensuSV—from the whole-genome sequencing data to the complete variant list}},
year = {2022}
}
```

## Testing scenario

For the ease of the verification of the algorithm, we created a testing scenarios. All of them process data of very shallow sequencing of HG00512, HG00513, HG00514, HG00731 done by 1000 Genomes consortium, that is why the obtained Structural Variants are really low in numbers. The main purpose of those tests, however, are to be sure all the systems work properly before setting up multi-scale analysis.

In all the testing scenarios, you need to get into the docker container having downloaded it using:

```shell
docker run -p 8082:8082 -it mateuszchilinski/consensusv-nf-pipeline:latest
```

The testing scenario uses all 4 samples (2 in form of CRAM, and 2 in form of BAM):

```shell
./test_nf.sh
```

## Preparation of your samples

Since the algorithm is easy-to-run, the preparation of the sample is minimal. The algorithm works with bam/cram files, and a csv file containing all the sample information is needed. There is only one column, that is the path to the bam/cram file.

An example csv file (provided with the software as well):

Path |
-------------- |
/HG00512.alt_bwamem_GRCh38DH.20150724.CHS.sv_7kb_mate.cram |
/HG00513.alt_bwamem_GRCh38DH.20150724.CHS.sv_7kb_mate.cram |
/HG00514.bam |
/HG00731.bam |

Bear in mind that the column headers are provided only for the ease of the example, and should not be present in the csv file.

## Running the pipeline

To get into the docker image for working with your data, it's best to mount local directories to the container:

```bash
docker run --mount type=bind,source=/mnt/,target=/mnt/ -p 8082:8082 -it mateuszchilinski/consensusv-nf-pipeline:latest
```

In that example, we bind folder /mnt/, where working directory and samples will be stored to in-container folder /mnt/. Once in the container, we can run the pipeline using the following command:

```bash
nextflow main.nf --design design.csv -with-conda
```
Remember to put your samples in the file samples.csv according to the [guidelines](#preparation-of-your-samples). 

All the parameters that can be used with the script are shown in the following table:

Parameter | Description
-------------- | ---------------
--design | File location of the csv file that described all the samples according to the [guidelines](#preparation-of-your-samples).
--threads | Max number of threads per task.
--mem | Max memory per thread.
--outdir | Output dir of ConsensuSV-core.
--ref | Reference genome

## Output location

The location of the output depends on your working directory, provided as the parameter. In that directory, two folder will be created:
* output - a folder where ConsensuSV calls of Structural Variants are stored
* vcfs - where you will find folders for each of the sample, containing all VCF files from the individual SV calling tools, along with file with SNPs and Indels (separately, SNPs.vcf and Indels.vcf)

## Pipeline details

The overall schema of the pipeline is shown on the following picture:

<p align="center">
<img src="https://github.com/SFGLab/ConsensuSV-nf-pipeline/blob/main/pipeline.png" />
</p>

If the input is cram file, it is unpacked to bam file. Then, the bam files are indexed, and multiple tools are run in parallel to obtain the Structural Variants, Indels and SNPs. We are using bcftools for the Indels and SNPs callings. The SV-callers used in this pipeline are: Delly, BreakDancer, Tardis, CNVNator, BreakSeq, Manta, Lumpy, and Whamg.

The final step is merging the Structural Variant calls into one unified file using our ConsensuSV-core algorithm. For the details on it, please refer to ConsensuSV-core github repository (https://github.com/SFGLab/ConsensuSV-core).

## Setup on NVIDIA DGX A100 systems

The consensusv-nf-pipeline can be run using HPC software. The parallelism nature of nextflow lets multiple samples to be processed at once. We have tested the pipeline on NVIDIA DGX A100 cluster, and the following requirements need to be met:
* slurm
* pyxis version at least 0.11.0
* enroot installed
* usage of --container-writable flag; most of the temp files are stored in the temporary direction provided by the user, however for testing scenarios & actual calculations some of the SV callers need to use space of the container - however, it is limited and there is no need for further adjustment

Preparation of the enroot image for DGX A100 systems can be done locally using the following commands:
```bash
enroot import docker://mateuszchilinski@mateuszchilinski/consensusv-nf-pipeline
```

Bear in mind that because of the size of the container, you might be required to change the temporary folder for enroot before executing the previous command (:
```bash
export TMPDIR=/home/dir_to_your_temp_folder/
```

After having created enroot image, you need to upload it (mateuszchilinski+consensusv-nf-pipeline.sqsh) to your DGX A100 system (using e.g. scp). After the upload is complete, you can run the container using the following command:

```bash
srun --pty --container-image ~/mateuszchilinski+consensusv-nf-pipeline.sqsh --container-writable /bin/bash
```

Sometimes, before running this commend you need to set the environmental variable:

```bash
XDG_RUNTIME_DIR=~
```

However, since the sample should be out-of-the-container, it is recommended to mount the folder containing samples to the container (in the following example, we mount /dir/to/data location on the cluster to /data folder in the container). You can also mount more folders to the container, depending on what you want to do:

```bash
srun --pty --container-mounts /dir/to/data:/data:rw --container-image ~/mateuszchilinski+consensusv-nf-pipeline.sqsh --container-writable /bin/bash
```

After that, you can start using the software, e.g. run:

```
./test_run_csv.sh
```

## Setup using singularity

To prepare singularity image use the following command (bear in mind, that some of the pieces of the software we are using use their own directory as temp file storage, so we need to create writable container - this we will run it in "sandbox" mode):

```bash
singularity build --sandbox DIR_WHERE_CONTAINER_WILL_BE/ docker://mateuszchilinski/consensusv-nf-pipeline
```

Then to run the container:

```bash
singularity run --pwd /workspace/ --writable DIR_WHERE_CONTAINER_WILL_BE/
```

Then you can simply start using the software as described in [Running the pipeline](#running-the-pipeline).

## Benchmark

We have used 9 samples for the benchmark provided by NYGC - HG00512, HG00513, HG00514, HG00731, HG00732, HG00733, NA19238, NA19239, NA19240. The comparisons were done using svbench (https://github.com/kcleal/svbench) and Venn diagrams of the common SVs. The results from svbench can be seen below:

<p align="center">
<img src="https://github.com/SFGLab/ConsensuSV-nf-pipeline/blob/main/benchmark.png" />
</p>

And the Venn diagrams can be seen there:

<p align="center">
<img src="https://github.com/SFGLab/ConsensuSV-nf-pipeline/blob/main/average.png" />
</p>

<p align="center">
<img src="https://github.com/SFGLab/ConsensuSV-nf-pipeline/blob/main/venns.png" />
</p>
