[![license](https://img.shields.io/badge/license-MIT-yellow.svg)](https://github.com/xa-lab/RFcaller/tree/master/LICENSE) ![version](https://img.shields.io/badge/version-1.2.0-blue) [![zenodo](https://img.shields.io/badge/docs-zenodo-green)](https://zenodo.org/record/7113432#.YzG0Ay8RrSw) [![DOI](https://zenodo.org/badge/doi/10.1101/2022.05.11.491496.svg)](https://doi.org/10.1101/2022.05.11.491496)

[![License](https://images.batchx.io/gh-badge-logo.svg)](https://platform.batchx.io/uniovi/profile)

# Context

A pipeline that uses read-level features and [extra trees/random forest algorithms](https://en.wikipedia.org/wiki/Random_forest) for accurate and fast detection of somatic mutations in next generation sequencing data.

# Input

## Required inputs

This image has the following required inputs:

1. `normalBam/tumorBam`: Normal/tumor [BAM](https://genome.ucsc.edu/FAQ/FAQformat.html#format5.1) *(CRAM supported)*.
2. `normalName/tumorName`: Name of normal/tumor sample.
3. `outputName`: Name of the analysis.
4. `genomeFasta`: Reference genome in [FASTA format](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp) *(bgzip supported)*.
5. [`dbSNP`](https://en.wikipedia.org/wiki/DbSNP): [VCF](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/) files provided by the image with common SNPs (MAF ≥ 1%) to eliminate these positions from the analysis. In case you want to use your own VCFs, use the *`Other`* argument and complete the *`dbSNPfile`* option.
    - Default: `hg19`
    - Choices: `[hg19, hg38, Other]`

## Optional inputs

This image provides additional configuration through the following optional inputs:

1. `normalIndex/tumorIndex`: Index of normal/tumor BAM. Not required but recommended.
    - *To index a BAM file use `samtools index`*
2. `genomeIndexFai`: FAI index of FASTA genome. Not required but recommended.
    - *To index a FASTA file use `samtools faidx`*
3. `genomeIndexGzi`: GZI index of FASTA genome. Only when the FASTA is compressed with [bgzip](https://www.htslib.org/doc/bgzip.html). Not required but recommended.
    - *To index a FASTA file use `samtools faidx`*
4. `dbSNPfile`: VCF provided by the user with common SNPs to eliminate these positions from the analysis *(bgzip supported)*.
5. `dbSNPindex`: In case of upload a VCF.GZ (bgzip) file, provide also the index. Not required but recommended.
    - *To index a VCF (bgzip) file use `tabix`*
6. `PoN`: TSV provided by the image with positions that appear recurrently in a Panel of Normals to eliminate them from the analysis. In case you want to use your own PoN use the *`Other`* argument and complete the *`PoNuser`* option.
    - Default: `hg19`
    - Choices: `[hg19, hg38, None, Other]`
7. `PoNuser`: Same as *`PoN`* but provided by the user. Only when *`PoN`* argument is *`Other`*.
8. `regions`: File with the positions to be included in the analysis (VCF, [BED](https://www.ensembl.org/info/website/upload/bed.html) and TSV formats supported).
10. `assemblyTag`: [BCFTOOLS tag](http://samtools.github.io/bcftools/bcftools.html#call) for ploidy to use in variant calling. 
    - Default: `GRCh37`
    - Choices: `[GRCh37, GRCh38, GRCm38, GRCm39, Other, None]`
11. `assemblyFile`: File with ploidy information. For more information use the command `bcftools call --ploidy ?`.

## Additional parameters

The image has some extra parameters to use as filters:
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

# Output

At the end of RFcaller two files are generated:

1. `mutations_${outputName}.vcf`: Final VCF with the mutations from SNV and INDEL pipelines.
2. `RFcaller.log`: Log file with the commands that have been executed for each case and their results.


# Example

The files used in this example can be downloaded from [xa-lab/RFcaller](https://github.com/xa-lab/RFcaller/tree/master/example) repository and the reference genome (GRCh37) [here](https://dcc.icgc.org/releases/PCAWG/reference_data/pcawg-bwa-mem). Then, we can upload them to our BatchX filesystem:

```bash
bx cp RFcaller/example bx://example/
```

> Note: Before uploading the genome we have to compress it with `bgzip`.

## How to create BAMs and FASTA index files

RFcaller indexes files automatically if they are not provided by the user. However, it can be very useful to know how to create them from BatchX.

### *BAM indexation*

To index a BAM we will use [batchx@bioinformatics/samtools/index](https://platform.batchx.io/batchx/images/bioinformatics%2Fsamtools%2Findex):

```bash
bx submit -v=1 -m=4000 batchx@bioinformatics/samtools/index '{
  "bam": "bx://example/normal.bam"
}'

bx submit -v=1 -m=4000 batchx@bioinformatics/samtools/index '{
  "bam": "bx://example/tumor.bam"
}'
```

### *FASTA indexation*

As the downloaded genome is compressed with `gzip` instead of `bgzip`, first we need to change it. To do that, we will use these two commands in our computer. Finally, we can upload the genome:

```bash
gunzip genome.fa.gz
bgzip genome.fa
bx cp genome.fa.gz bx://reference/
```

To index the FASTA we will use the [batchx@bioinformatics/samtools/faidx](https://platform.batchx.io/batchx/images/bioinformatics%2Fsamtools%2Ffaidx):

```bash
bx submit -v=1 -m=3000 batchx@bioinformatics/samtools/faidx '{
  "fastaFile": "bx://reference/genome.fa.gz"
}'
```

##  Basic input
Once we have all the necessary files, we can continue with the example. This is a basic `input.json` to run RFcaller in BatchX.
```bash
bx submit -v=20 -m=20000 uniovi@labxa/rfcaller '{
  "normal": {
    "normalBam": "bx://example/normal.bam",
    "normalName": "normal",
    "normalIndex": "bx://example/normal.bam.bai",
},
  "tumor": {
    "tumorBam": "bx://example/tumor.bam",
    "tumorName": "tumor",
    "tumorIndex": "bx://example/tumor.bam.bai"
},
  "outputName": "test",
  "genome": {
    "genomeFasta": "bx://reference/genome.fa.gz",
    "genomeIndexFai": "bx://reference/genome.fa.gz.fai",
    "genomeIndexGzi": "bx://reference/genome.fa.gz.gzi"
},
  "dbSNP": "hg19"
}'
```

##  dbSNP provided by the user

In case you want to provide your own dbSNP, the `input.json` should be like this:

```bash
bx submit -v=20 -m=20000 uniovi@labxa/rfcaller '{
  "normal": {
    "normalBam": "bx://example/normal.bam",
    "normalName": "normal",
    "normalIndex": "bx://example/normal.bam.bai",
},
  "tumor": {
    "tumorBam": "bx://example/tumor.bam",
    "tumorName": "tumor",
    "tumorIndex": "bx://example/tumor.bam.bai"
},
  "outputName": "test",
  "genome": {
    "genomeFasta": "bx://reference/genome.fa.gz",
    "genomeIndexFai": "bx://reference/genome.fa.gz.fai",
    "genomeIndexGzi": "bx://reference/genome.fa.gz.gzi"
},
  "dbSNP": "Other",
  "dbSNPuser": {
    "dbSNPfile": "bx://custom_dbSNP.vcf.gz",
    "dbSNPindex": "bx://custom_dbSNP.vcf.gz.tbi"
  }
}'
```

# Tools Version

- [samtools](https://www.htslib.org/doc/1.10/samtools.html) (built with v1.10)
- [bcftools](https://www.htslib.org/doc/1.10/bcftools.html) (built with v1.10.2)
- [bedtools](https://bedtools.readthedocs.io/en/latest/) (built with v2.29.1)

------
**Authors:** Ander Díaz-Navarro, Pablo Bousquets-Muñoz & Xose S. Puente -- Universidad de Oviedo