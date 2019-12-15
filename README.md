# WES (Whole-Exome Sequencing) pipeline

## Introduction

Here, we provide a series of scripts for the analysis of Whole-Exome sequencing data of cancer samples.
The analysis is performed by standard tools as BWA, Picard, GATK and MuTect.
We also provide the script to perform the LODN filtering for the GATK results.
The LODN is based on the algorithm of MuTect and it evaluates the somatic status of single
nucleotide variants found in cancer samples.

See [wiki](https://bitbucket.org/BBDA-UNIBO/wes-pipeline/wiki/Home) for further details.

## Requirements

* Python

gatk_pipe.py:

* [BWA](http://bio-bwa.sourceforge.net/)
* [Picard](http://broadinstitute.github.io/picard/) (Version >= 2.1.1)
* [GATK](https://www.broadinstitute.org/gatk/download/)
* [Annovar](http://annovar.openbioinformatics.org/en/latest/)


run_lodn.py:

* [pysam](http://pysam.readthedocs.io/en/latest/api.html)


mutect_pipe.py:

* [MuTect](https://www.broadinstitute.org/cancer/cga/mutect)
* [Annovar](http://annovar.openbioinformatics.org/en/latest/)


## Installation

    git clone https://github.com/italodovalle/wes-pipeline.git
    cd wes-pipeline/
    python setup.py install

## Usage

### gatk_pipe.py

    usage: gatk_pipe.py [-h] -i [FASTQ [FASTQ ...]] -s NAME [-c CORES] -f CONFIG
                        [-d OUTDIR]

    Whole Exome Analysis Pipeline

    optional arguments:
    -h, --help            show this help message and exit
    -i [FASTQ [FASTQ ...]]
                            fastq files
    -s NAME               Sample name
    -c CORES              number of cores to use, default = 2
    -f CONFIG             config file yaml
    -d OUTDIR             output directory

### run_lodn.py

    usage: run_lodn.py [-h] -i INFILE -f FMT -b BAM -o OUTFILE
                    [-min_n_cov MIN_N_COV] [-min_t_cov MIN_T_COV]

    Run LODN filtering

    optional arguments:
    -h, --help            show this help message and exit
    -i INFILE             infile
    -f FMT                table or vcf
    -b BAM                normal bam file
    -o OUTFILE            outfile
    -min_n_cov MIN_N_COV  min normal coverage
    -min_t_cov MIN_T_COV  min tumor coverage


### mutect_pipe.py


    usage: mutect_pipe.py [-h] -n NORMAL -t TUMOR -s NAME -c CODE -f CONFIG
                        [-d OUTDIR]

    MuTect Analysis Pipeline

    optional arguments:
    -h, --help  show this help message and exit
    -n NORMAL   Normal bam file
    -t TUMOR    Tumor bam file
    -s NAME     Sample name
    -c CODE     Normal sample name
    -f CONFIG   config file yaml
    -d OUTDIR   output directory

## Citation

do Valle, Í. F. et al. Optimized pipeline of MuTect and GATK tools to improve the detection of somatic single nucleotide polymorphisms in whole-exome sequencing data. BMC Bioinformatics 17, 341 (2016).

[See paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5123378/)

## Contact

Ítalo Faria do Valle: italo.fariadovalle2@unibo.it

Daniel Remondini: daniel.remondini@unibo.it

Gastone Castellani: gastone.castellani@unibo.it

## License

Copyright (c) 2016 - Ítalo Faria do Valle, Daniel Remondini & Gastone Castellani

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
