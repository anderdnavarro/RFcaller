# RFcaller
[![Support this project by running your production jobs at BatchX](https://images.batchx.io/gh-badge-logo.svg)](https://platform.batchx.io/uniovi/tools/labxa%2Frfcaller "Support this project by running your production jobs at BatchX")

A pipeline that uses read-level features and extra trees/random forest algorithms for accurate and fast detection of somatic mutations in next generation sequencing data.

## Index
1. [Installation](#Installation)
2. [Quick usage](#Quick-usage)
3. [Docker options](#Docker-options)
4. [RFcaller options](#RFcaller-options)
    - [Required inputs](#Required-inputs)
    - [Optional inputs](#Optional-inputs)
    - [Additional parameters](#Additional-parameters)
5. [Examples](#Examples)
	- [Writting input file](#Writting-input-file)
	- [Changing default paths](#Changing-default-paths)
	- [Changing TimeZone](#Changing-TimeZone)
6. [Outputs](#Outputs)

## Installation
We have created a docker image with all dependencies installed:

- If you don't have docker already installed in your system, please follow these [instructions](https://docs.docker.com/get-docker/).

```bash
docker pull labxa/rfcaller:1.0.0
```

The image has the following structure:

- `databases` directory contains dbSNP and ploidy files.
- `RFcaller` has the scripts used by the pipeline.
- `output` is the working directory.

```basic
/home
 |-- databases
      |-- dbSNP
           |-- UCSC_dbSNP153Common_hg19_combined.vcf.gz
           |-- UCSC_dbSNP153Common_hg19_combined.vcf.gz.tbi
           |-- UCSC_dbSNP153Common_hg38_combined.vcf.gz
           |-- UCSC_dbSNP153Common_hg38_combined.vcf.gz.tbi
      |-- ploidy_files
           |-- GRCm38.ploidy.file
           |-- GRCm39.ploidy.file
 |-- RFcaller
/output
```

## Quick usage
Here is basic configuration, for more information see [docker](#Docker-options) and [RFcaller](#RFcaller-options) options.

```bash
docker run --rm -v /BAMS_PATH/:/bams/ -v /GENOME_PATH/:/genome/ -v $(pwd):/output/ -u $(id -u):$(id -g) -it labxa/rfcaller:1.0.0 -@ INT -i INPUT.LIST --genome /genome/GENOME.FA --dbSNP hg19/hg38
```

## Docker options
Among all the options offered by docker (`docker run --help`), we recommend:
1. `--rm`: Automatically remove the container when it exits.
2. `-v, --volume`: Mount local volumes in the container.
    - With the option `-v $(pwd):/output/`, RFcaller results will be in your current directory.
3. `-u, --user`: Specify the user ID and its group ID. It's useful to not run the pipeline as root.
4. `-i, --interactive`: Keep STDIN open even if not attached.
5. `-t, --tty`: Allocate a pseudo-TTY. When combined with `-i` it allows you to connect your terminal with the container terminal.
6. `-e, --env`: Set environment variables.

## RFcaller options
### Required inputs
RFcaller has the following required inputs:
1. `-i, --input`: The input is a TSV (tab-separated values) file with five required columns  and an optional extra one. This last column is used to specify the file with the positions you want to analyze. 

|  |  |  |  |  |  |
|:---:|:---:|:---:|:---:|:---:|:---:|
| Normal name | Normal BAM<br>*(CRAM supported)* | Tumor name | Tumor BAM<br>*(CRAM supported)* | Output name | VCF/BED<br>(optional) |
|  |  |  |  |  |  |

2. `-g, --genome`: Reference genome in FASTA format *(bgzip supported)*.
3. `-p, --dbSNP`: VCF file provided by the image with common SNPs (MAF ≥ 1%) to eliminate these positions from the analysis. In case you want to use your own VCF *(bgzip supported)*, use also the `--ploidy_file` argument.
    - *Choices: `[hg19, hg38, your_dbSNP]`*

### Optional inputs
RFcaller provides additional configuration through the following optional inputs:
1. `-@`: Max number of threads to use.
    - *Default: 20*
2. `-b, --bamsDir`: Main directory where BAMs are located inside the container. By default, the program will look in the `/bams/` directory and that's why we set the docker option *`-v /BAMS_PATH/:/bams/`*. However, if you mount your BAMs in other directory inside the container, you have to use this option to specify it.
    - *Default: `/bams/`*
3. `-o, --outputDir`: Name for the output directory. By default, the results will be saved in your current directory, but with this option a new folder will be created to save them.
4. `--positions`: This argument refers to the file format of the `--input` optional column.
    - *Default: `VCF`*
    - *Choices: `[VCF, BED]`*
5. `--ploidy_file`: `bcftools call --ploidy_file` to use in variant calling. By default, if the user specifies *`--dbSNP hg19`* or *`--dbSNP hg38`*, `--ploidy_file` is GRCh37 or GRCh38 respectively. Otherwise, you must select between the following options: 
    - *Choices: `[GRCh37, GRCh38, GRCm38 (mouse), GRCm39 (mouse), None (consider all sites as diploid), your_ploidy_file]`*
6. `--includePatches`: Tag used to analyze both canonical chromosomes and patches.
7. `--keep`: Tag used to keep intermediate files.

### Additional parameters
RFcaller has some extra parameters to use as filters:
|      Parameter      |                                                                                                                                        Description                                                                                                                                       | Default |
|:-------------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:-------:|
|    contamination    | Percentage of tumor contamination in normal sample. For example,  if you expect that 10% of reads in normal comes from  tumor (10% contamination),you have to set `"contamination": 0.1`. This argument is used to increase the maximum number of mutated reads allowed in normal sample |   0.05  |
|      TD_cov_SNV     |                                                                                                                         Minimum coverage for tumor sample (SNVs)                                                                                                                         |    7    |
|     TD_cov_INDEL    |                                                                                                                        Minimum coverage for tumor sample (INDELs)                                                                                                                        |    7    |
|      ND_cov_SNV     |                                                                                                                         Minimum coverage for normal sample (SNVs)                                                                                                                        |    7    |
|     ND_cov_INDEL    |                                                                                                                        Minimum coverage for normal sample (INDELs)                                                                                                                       |    7    |
|      TD_mut_SNV     |                                                                                                 Minimum number of mutated reads for a position in tumor sample  to not discard it (SNVs)                                                                                                 |    3    |
|     TD_mut_INDEL    |                                                                                                Minimum number of mutated reads for a position in tumor sample  to not discard it (INDELs)                                                                                                |    4    |
|      ND_mut_SNV     |                                                                                                      Maximum number of mutated reads allowed for a position  in normal sample (SNVs)                                                                                                     |    3    |
|     ND_mut_INDEL    |                                                                                                     Maximum number of mutated reads allowed for a position  in normal sample (INDELs)                                                                                                    |    2    |
|      ND_window      |                                                                                                         Window size around a position to look for mutations  in normal (INDELs)                                                                                                         |    10   |
|    SNV_threshold    |                                                                                                                       Minimum QUAL value to consider a SNV as good                                                                                                                      | 10.1774 |
|   INDEL_threshold   |                                                                                                                     Minimum QUAL value to consider an INDEL as good                                                                                                                     | 32.1418 |
| polyINDEL_threshold |                                                                                                                Minimum QUAL value to consider an homopolymerINDEL (mononucleotide microsatellites) as good                                                                                                               |  0.7723 |

## Examples
### Writting input file
Imagine that you are working in the directory `/home/example` which has the following structure:

```basic
/home/example/
 |-- example.metadata
 |-- bams
      |-- normal.bam
	  |-- normal.bam.bai
	  |-- tumor.bam
	  |-- tumor.bam.bai
 |-- genome
 	  |-- hg19.fa
	  |-- hg19.fa.fai
```

If we use the following command to run RFcaller:

```bash
docker run --rm -v /home/example/bams/:/bams/ -v /home/example/genome/:/genome/ -v $(pwd):/output/ -u $(id -u):$(id -g) -it labxa/rfcaller:1.0.0 -i example.metadata --genome /genome/hg19.fa --dbSNP hg19
```

The example.metadata file should be:

```basic
normal	normal.bam	tumor	tumor.bam	example
```

This is because thanks to the `-v` command, the files inside the container are located in the following way:
|              Computer              |      Docker      |
| :--------------------------------: | :--------------: |
| /home/example/example.metadata | /output/example.metadata |
| /home/example/bams/normal.bam | /bams/normal.bam |
| /home/example/bams/tumor.bam  | /bams/tumor.bam  |
| /home/example/genome/hg19.fa  | /genome/hg19.fa  |

A more detailed explanation:
- We use `-i example.metadata` because inside the container, we are working in the `/output/` directory, which is where the `example.metadata` file is located. It would also be valid to write `-i /output/example.metadata`.
- We have said before that the default directory for the BAMs inside the container is `/bams/`. In this sense, as our BAMs are in this directory, we don't have to specify the full path in the `example.metadata` file, just write where these BAMs are located inside this folder.
- In the case of genome, we have to specify the full path inside the container, that is why we use `--genome /genome/hg19.fa`.

### Changing default paths
To make it clear how to use the Docker options together with RFcaller, we have prepared a more complex example. For this time, we start with the following directory structure:

```basic
/home/complex_example/
 |-- example.metadata
/data/
 |-- bams
      |-- normal_dir
      		|-- normal.bam
	  		|-- normal.bam.bai
      |-- tumor_dir
	  		|-- tumor.bam
	  		|-- tumor.bam.bai	
 |-- hg19_genome
 	  |-- hg19.fa
	  |-- hg19.fa.fai
```

In the event that we were working in the `/home/` directory:

```bash
docker run --rm -v /data/bams/:/example_bams/ -v /data/hg19_genome/:/genome/ -v $(pwd):/output/ -u $(id -u):$(id -g) -it labxa/rfcaller:1.0.0 -i complex_example/example.metadata --bamsDir /example_bams/ --genome /genome/hg19.fa --dbSNP hg19 
```

The distribution of the files will be:
|              Computer              |      Docker      |
| :--------------------------------: | :--------------: |
| /home/complex_example/example.metadata | /output/complex_example/example.metadata |
| /data/bams/normal_dir/normal.bam | /example_bams/normal_dir/normal.bam |
| /data/bams/tumor_dir/tumor.bam | /example_bams/tumor_dir/tumor.bam |
| /data/hg19_genome/genome/hg19.fa  | /genome/hg19.fa  |

Thus, the `metadata.example` file should be:

```basic
normal	normal_dir/normal.bam	tumor	tumor_dir/tumor.bam	example
```

A more detailed explanation:
- As we are using `-v $(pwd):/output/` and our working directory outside the container is `/home/`, the `example.metadata` file in the container is inside the folder `complex_example`. This is why we use `-i complex_example/example.metadata`.
- We have mounted the BAMs directory in `/example_bams/` instead of `/bams/`, so we also have to changed BAMs directory for RFcaller with the `-b, --bamsDir` option. Now, for the `example.metadata` we have to write the full BAMs path removing the main directory because it is defined with  `--bamsDir`.

### Changing TimeZone
The time zone of the docker image is `Europe/Madrid` and it's used to set the time for the log file. To change this, add the `-e` option to set environment variables:

```bash
docker run --rm -e "TZ=$(cat /etc/timezone)" -v /BAMS_PATH/:/bams/ -v /GENOME_PATH/:/genome/ -v $(pwd):/output/ -u $(id -u):$(id -g) -it labxa/rfcaller:1.0.0 -@ INT -i INPUT.LIST --genome /genome/GENOME.FA --dbSNP hg19/hg38
```

Or set it manually (to know your TZ visit: [TZ database](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones))

```bash
docker run --rm -e "TZ=America/Toronto" -v /BAMS_PATH/:/bams/ -v /GENOME_PATH/:/genome/ -v $(pwd):/output/ -u $(id -u):$(id -g) -it labxa/rfcaller:1.0.0 -@ INT -i INPUT.LIST --genome /genome/GENOME.FA --dbSNP hg19/hg38
```

## Outputs
For each case three directories and a final file will be created:

1. calling: bcftools basic call results
	- `${name}.vcf.gz` -> Original calling from bcftools
	- `${name}.norm.vcf.gz` -> After normalizing the indels
	- `${name}.filter.norm.vcf.gz` -> After removing SNPs and very low quality mutations
	- `${name}.SNVs.filter.norm.vcf` and `${name}.INDELs.filter.norm.vcf` -> Separate SNVs and INDELs into two files
	- `polyIndel.pos` -> BED file with those positions which contain an homopolymer INDEL
2. somaticSNV: RFcaller results for SNVs
	- `reduced.positions` -> All positions to analyze in `${name}.SNVs.filter.norm.vcf`
	- `reduced_*.bam` -> A small bam with just the reads overlapping positions in `reduced.positions` 
	- `*.mini.pileup` and `${name}.mini.pileup.vcf` -> Filtered mini.pileup and its VCF for the positions in `reduced.positions` 
	- `common.positions` -> File with the positions in `${name}.mini.pileup.vcf` with the format `CHR \t POS`
	- `${tumor}.read.names`, `${name}.mutations.interval`, `${name}.sequence`,  `${name}.stats`, `${name}.cigar`  and `${name}.mapq` -> Intermediate files with read-level features
	- `${name}.prepare_muts.csv` -> File generated from gathering the above features and filtering them
	- `regression_results_SNVs_${name}.txt`-> Regression algorithm result for SNVs
	- `mutations_SNVs_${name}.vcf.gz` -> Final set of SNVs
3. somaticINDEL: RFcaller results for INDELs
	- `reduced.positions` -> All positions to analyze in `${name}.INDELs.filter.norm.vcf`
	- `short_reduced.positions` and `long_reduced.positions` -> Positions separated by the size of the INDEL (short < 7 and long >= 7)
	- `reduced_*.bam` -> A small bam with just the reads overlapping positions in `reduced.positions` 
	- `*.mini.pileup` and `${name}.mini.pileup.vcf` -> Filtered mini.pileup and its VCF for the positions in `reduced.positions` 
	- `short_common.positions`, `long_common.positions` and `common.positions` -> File with the positions in `${name}.mini.pileup.vcf` with the format `CHR \t POS`
	- `${name}.features`, `${name}.mutations.interval`, `${name}.sequence`, `${name}.distribution`,  `${name}.stats`, `${name}.cigar`  and `${name}.mapq` -> Intermediate files with read-level features
	- `${name}.prepare_muts.csv` -> File generated from gathering the above features and filtering them
	- `regression_results_INDELs_${name}.txt`-> Regression algorithm result for INDELs
	- `mutations_INDELs_${name}.vcf.gz` -> Final set of INDELs
4. `mutations_${name}.vcf` -> Final VCF with the mutations from SNV and INDEL pipelines

If the tag `--keep` is set, all above directories are removed and it is only conserved the final VCF file `mutations_${name}.vcf`.

Two additional global files are created:

- RFcaller.log -> Log file with the commands that have been executed for each case and their results
- .stderr.log -> A file with errors and warnings that have occurred during the pipeline

------

**Authors:** Ander Díaz-Navarro, Pablo Bousquets-Muñoz & Xose S. Puente -- Universidad de Oviedo
