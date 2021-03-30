# RFcaller

A pipeline that uses read-level features and extra trees/random forest algorithms for accurate and fast detection of somatic mutations in next generation sequencing data.

## Index

1. [Installation](#Installation)
2. [Usage](#Usage)
	- [Docker options](#Docker-options)
	- [RFcaller options](#RFcaller-options)
	- [Examples](#Examples)
  3. [Outputs](#Outputs)

## Installation

We have created a docker image with all dependencies installed:

- If you don't have docker already installed in your system, please follow these [instructions](https://docs.docker.com/get-docker/)

```bash
docker pull labxa/rfcaller:latest
```

The image has the following structure:

- `databases` directory contains dbSNP and genome files
- `RFcaller` has the scripts used by the pipeline
- `output` is the default working directory

```basic
/home
 |-- databases
      |-- dbSNP
           |-- UCSC_dbSNP153Common_hg19_combined.vcf.gz
           |-- UCSC_dbSNP153Common_hg19_combined.vcf.gz.tbi
           |-- UCSC_dbSNP153Common_hg38_combined.vcf.gz
           |-- UCSC_dbSNP153Common_hg38_combined.vcf.gz.tbi
      |-- hg19
           |-- hs37d5.fa.gz
           |-- hs37d5.fa.gz.fai
           |-- hs37d5.fa.gz.gzi
      |-- hg38
           |-- hg38.fa.gz
           |-- hg38.fa.gz.fai
           |-- hg38.fa.gz.gzi
 |-- RFcaller
 |-- output
```



## Usage

Here is quick configuration, for more information see [docker](#Docker-options) and [RFcaller](#RFcaller-options) options.

```bash
docker run --rm --cpus INT -v /BAMS_PATH/:/bams/ -v $(pwd):/output/ -u $(id -u):$(id -g) -it labxa/rfcaller -@ INT -i INPUT.LIST
```

### Docker options

Among all the options offered by docker, we recommend:

```bash
docker run --rm --cpus INT -v /BAMS_PATH/:/bams/ -v $(pwd):/output/ -u $(id -u):$(id -g) -it labxa/rfcaller
```

###### --rm

Automatically remove the container when it exits.

###### --cpus

Docker has the option to regulate the number of CPUs that can be used by the container. Although this option is also inside pipeline by the argument `-@`, it is recommended to set the same number for both options.

###### -v, --volume

Mount local volumes in the container.

###### -u, --user

Specify the user ID and its group ID. It's useful to not run the pipeline as root.

### RFcaller options

###### -i, --input

The input is a TSV file with five required columns and an optional extra one:

1. Normal name
2. Normal BAM
3. Tumor name
4. Tumor BAM
5. Output name
6. VCF/BED (optional) -> A VCF/BED with specific mutations to test. If this field is provided, the pipeline will only analyze these positions rather than whole genome/exome

*Note: Using Docker, your local path it's not the same as the path inside the docker container. By default, the program will look for BAMs in the `/bams` folder, that's why we set the -v option `/BAMS_PATH/:/bams/`. However, if  BAMs are inside folders or in another path you have to specify it. See [Examples](#Set-BAMs'-path) for more information.*

###### -@

*Default: 24*

Max number of threads to use. 

###### -w, --workDir

*Default: /output/*

Directory where the pipeline will be run. The default option is ready to be used with the docker parameter *`-v $(pwd):/output/`*, so the results of the pipeline will be in your current directory.

###### -b, --bamsDir

*Default: /bams/*

Main directory where BAMs are located. By default, the program will look in the `/bams/` directory and that's why we set the -v option  `/BAMS_PATH/:/bams/`. However, if you mount your BAMs in other directory inside the container, you have to use this option to specify it.

###### -o, --outputDir

*Default: .*

Name for the output directory. By default, the results will be saved in your current directory, but with this option a new folder will be created to save them.

###### -g, --genome

*Default: /home/databases/hg19/hs37d5.fa.gz*

Reference genome in fasta format **(must be indexed)**. The default path points at the hg19 reference genome installed by default in the container. In case of you need hg38 it is also provided within the image *(/home/databases/hg38/hg38.fa.gz)*. If you want to use your own genome see the [Examples](#Provide-your-own-genome) section below.

*To index a fasta file use samtools faidx.*

###### -p, --dbSNP

*Default: /home/databases/dbSNP/UCSC_dbSNP153Common_hg19_combined.vcf.gz*

dbSNP database with common SNPs (MAF ≥ 1%) **(the VCF must be indexed)**. The default path is ready to be used inside the container. For hg38 dbSNP the path is: */home/databases/dbSNP/UCSC_dbSNP153Common_hg38_combined.vcf.gz*.

*To index a VCF file, first compress it with bgzip and then use tabix.*

###### --positions

*Default: VCF [VCF, BED]*
Format of the file with the positions to make the call.

###### --contamination

*Default: 0.05 [Range: 0-1]*

Percentage of tumor contamination in normal sample. For example, if you expect that 10% of reads in normal comes from tumor (10% contamination), you have to set `--contamination 0.1`. This argument is used to increase the maximum number of mutated reads allowed in normal sample.

###### --assembly

*Default: GRCh37 [GRCh37, GRCh38]* 

BCFtools tag for ploidy to use in variant calling. 

###### --TD_cov_SNV

*Default: 7*

Minimum coverage for tumor sample (SNVs).

###### --TD_cov_INDEL

*Default: 7*

Minimum coverage for tumor sample (INDELs).

###### --ND_cov_SNV

*Default: 7*

Minimum coverage for normal sample (SNVs).

###### --ND_cov_INDEL

*Default: 7*

Minimum coverage for normal sample (INDELs).

###### --TD_mut_SNV

*Default: 3*

Minimum number of mutated reads for a position in tumor sample to not discard it (SNVs).

###### --TD_mut_INDEL

*Default: 4*

Minimum number of mutated reads for a position in tumor sample to not discard it (INDELs).

###### --ND_mut_SNV

*Default: 3*

Maximum number of mutated reads allowed for a position in normal sample (SNVs).

###### --ND_mut_INDEL

*Default: 2*

Maximum number of mutated reads allowed for a position in normal sample (INDELs).

###### --ND_window

*Default: 10*

Window size around a position to look for mutations in normal (INDELs).

###### --SNV_threshold

*Default: 10.1774* 

Minimum QUAL value to consider a SNV as good.

###### --INDEL_threshold

*Default: 32.1418*

Minimum QUAL value to consider an INDEL as good.

###### --polyINDEL_threshold

*Default: 0.7723*

Minimum QUAL value to consider a polyINDEL as good.

###### --clean

*Default: False [True, False]*

Remove all intermediate files.

### Examples

#### Set BAMs' path

With these examples we'll try to explain how the connection between your computer and the docker works in order to write the BAMs' path in the `INPUT.list` file.

##### Basic example

If we use the following docker configuration:

```bash
docker run --rm --cpus INT -v /home/xalab/example/WGS/:/bams/ -v $(pwd):/output/ -u $(id -u):$(id -g) -it labxa/rfcaller -@ INT -i INPUT.LIST
```

|              Computer              |      Docker      |
| :--------------------------------: | :--------------: |
| /home/xalab/example/WGS/normal.bam | /bams/normal.bam |
| /home/xalab/example/WGS/tumor.bam  | /bams/tumor.bam  |

Because the default `bamsDir` is `/bams/` the `INPUT.list` file should be:

```basic
normal	normal.bam	tumor	tumor.bam	normalVStumor
```

##### Each BAM in a specific subfolder

It may happen that the BAMs are in different folders within the root directory:

```bash
docker run --rm --cpus INT -v /home/xalab/example/WGS/:/bams/ -v $(pwd):/output/ -u $(id -u):$(id -g) -it labxa/rfcaller -@ INT -i INPUT.LIST
```

|                     Computer                     |             Docker             |
| :----------------------------------------------: | :----------------------------: |
| /home/xalab/example/WGS/normal_sample/normal.bam | /bams/normal_sample/normal.bam |
|  /home/xalab/example/WGS/tumor_sample/tumor.bam  |  /bams/tumor_sample/tumor.bam  |

In this scenario, the `INPUT.list` file should be:

```basic
normal	normal_sample/normal.bam	tumor	tumor_sample/tumor.bam	normalVStumor
```

##### Changing the bamsDir option

In case you mount the BAMs in another directory instead of `/bams`, you have to change the `bamsDir` option:

```bash
docker run --rm --cpus INT -v /home/xalab/example/WGS/:/bams_dir2/ -v $(pwd):/output/ -u $(id -u):$(id -g) -it labxa/rfcaller -@ INT -i INPUT.LIST --bamsDir /bams_dir2
```

|              Computer              |        Docker         |
| :--------------------------------: | :-------------------: |
| /home/xalab/example/WGS/normal.bam | /bams_dir2/normal.bam |
| /home/xalab/example/WGS/tumor.bam  | /bams_dir2/tumor.bam  |

And the `INPUT.list` file:

```basic
normal	normal.bam	tumor	tumor.bam	normalVStumor
```

##### Several BAMs' folders

It is also possible that you have several folders containing the BAMs. In this case the `bamsDir`option must be `/`:

```bash
docker run --rm --cpus INT -v /home/xalab/example/WGS1/:/bams_dir1/ -v /home/xalab/example/WGS2/:/bams_dir2/ -v $(pwd):/output/ -u $(id -u):$(id -g) -it labxa/rfcaller -@ INT -i INPUT.LIST --bamsDir /
```

|              Computer               |        Docker         |
| :---------------------------------: | :-------------------: |
| /home/xalab/example/WGS1/normal.bam | /bams_dir1/normal.bam |
| /home/xalab/example/WGS2/tumor.bam  | /bams_dir2/tumor.bam  |

Thus the `INPUT.list` file must contain the absolute path:

```
normal	/bams_dir1/normal.bam	tumor	/bams_dir2/tumor.bam	normalVStumor
```

#### Provide your own genome

In case you want to use your own reference genome you have to mount also this volume and specify the new path. For example:

```bash
docker run --rm --cpus INT -v /BAMS_PATH/:/bams/ -v $(pwd):/output/ -v /PATH_TO_REF_DIR/:/genome/ -u $(id -u):$(id -g) -it labxa/rfcaller -@ INT -i INPUT.LIST -g /genome/REF_GENOME.FA
```

#### Change time zone

The time zone of the docker image is `Europe/Madrid` and it's used to set the time for the log file. To change this, add the `-e` option to set environment variables:

```bash
docker run --rm --cpus INT -e "TZ=$(cat /etc/timezone)" -v /BAMS_PATH/:/bams/ -v $(pwd):/output/ -u $(id -u):$(id -g) -it labxa/rfcaller -@ INT -i INPUT.LIST
```

Or set it manually (to know your TZ visit: [TZ database](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones))

```bash
docker run --rm --cpus INT -e "TZ=America/Toronto" -v /BAMS_PATH/:/bams/ -v $(pwd):/output/ -u $(id -u):$(id -g) -it labxa/rfcaller -@ INT -i INPUT.LIST
```

#### Run bash

By default, the entrypoint of the container is the *RFcaller* script, so it will be launched automatically at startup. However, there are times when you may want to go inside the container. To do that, replace the default behaviour with the `--entrypoint` option:

```bash
docker run --rm --cpus INT -v /BAMS_PATH/:/bams/ -v $(pwd):/output/ -u $(id -u):$(id -g) --entrypoint bash -it labxa/rfcaller
```

#### Show Help message

```bash
docker run --rm -it labxa/rfcaller
```

```bash
Usage: RFcaller [options]

Options:
  -i, --input           FILE    TSV FILE with the cases to analyze
                                It must contain the following columns:
                                NORMAL_NAME NORMAL_BAM TUMOR_NAME TUMOR_BAM OUTPUT_NAME (optional: VCF/BED)
  -@                    INT     Number of threads to use [24]
  -w, --workDir         DIR     Directory where run the pipeline [/output]
  -b, --bamsDir         DIR     Main directory where BAMs are located [/bams]
  -o, --outputDir       STR     Name for the output directory [.]
  -g, --genome          FILE    Reference genome in FASTA format [/home/databases/hg19/hs37d5.fa.gz]
  -p, --dbSNP           FILE    VCF with common SNPs [/home/databases/dbSNP/UCSC_dbSNP153Common_hg19_combined.vcf.gz]
  --positions           STR     Format of the file with the positions to make the call (choices: VCF/BED) [VCF]
  --contamination       FLOAT   Percentage of tumor contamination in normal sample [0-1] [0.05]
  --assembly            STR     BCFTOOLS tag for ploidy to use in variant calling (choices: GRCh37/GRCh38) [GRCh37]
  --TD_cov_SNV          INT     Minimum coverage for tumor (SNVs) [7]
  --TD_cov_INDEL        INT     Minimum coverage for tumor (INDELs) [7]
  --ND_cov_SNV          INT     Minimum coverage for normal (SNVs) [7]
  --ND_cov_INDEL        INT     Minimum coverage for normal (INDELs) [7]
  --TD_mut_SNV          INT     Minimum number of mutated reads in tumor (SNVs) [3]
  --TD_mut_INDEL        INT     Minimum number of mutated reads in tumor (INDELs) [4]
  --ND_mut_SNV          INT     Maximum number of mutated reads in normal (SNVs) [3]
  --ND_mut_INDEL        INT     Maximum number of mutated reads in normal (INDELs) [2]
  --ND_window           INT     Window size around a position to look for mutations in normal (INDELs) [10]
  --SNV_threshold       FLOAT   Minimum regression value to consider a SNV as good [10.726]
  --INDEL_threshold     FLOAT   Minimum regression value to consider an INDEL as good [32.1418]
  --polyINDEL_threshold FLOAT   Minimum regression value to consider a polyINDEL as good [0.7723]
  --clean               BOOLEAN Remove all intermediate files (choices: True/False) [False]
  -h, --help                    Show this help and exit
```



### Outputs

For each case, three directories and a final file will be created:

1. calling: bcftools basic call results
	- `${name}.vcf.gz` -> Original calling from bcftools
	- `${name}.norm.vcf.gz` -> After normalizing the indels
	- `${name}.filter.norm.vcf.gz` -> After removing SNPs
	- `${name}.SNVs.filter.norm.vcf` and `${name}.INDELs.filter.norm.vcf` -> Separate SNVs and INDELs into two files
	- `polyIndel.pos` -> BED file with those positions which contain an INDEL of the class TAAAA -> TAAA
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

If the option `--clean True` is set, all above directories are removed and it is only conserved the final VCF file `mutations_${name}.vcf`.

Two additional global files are created:

- RFcaller.log -> Log file with the commands that have been executed for each case and their results
- .stderr.log -> A file with errors and warnings that have occurred during the pipeline

------

**Authors:** Ander Díaz-Navarro & Xose S. Puente -- Universidad de Oviedo