ngsCAT2: a tool to assess the efficiency of targeted enrichment sequencing
=======================================


# Requirements
This instructions are for debian based linux distributions such as: Ubuntu, Linux mint, Lubuntu...  


## Python Packages Requirements

All standard python packages required are specified in
the `requirements.txt` and in `setup.py` files.

Requirements
------------

Python 3.5 or later to run ngsCAT2.

In Ubuntu, Mint and Debian you can install Python 3 like this:

    $ sudo apt-get install python3 python3-pip
    
- **Samtools**
samtools 1.7 using htslib 1.7.2
Samtools webpage http://www.htslib.org/download/
```
sudo apt-get install samtools=1.7-1

```
- **Bedops version 2.4.25 or greater**
Bedops webpage https://bedops.readthedocs.io/en/latest/content/revision-history.html#v2-4-35

 
## Installation
For the installation of the tool just only run:

```
pip3 install git+https://github.com/alegarsan2/ngsCAT2@master
```
# Usage instruction

```
Usage: 	
       	****************************************************************************************************************
       	Task: Assesses capture performance in terms of sensibility, specificity and uniformity of the coverage.
       	Output: An html report will be created at the path indicated with the --out option.
       	*****************************************************************************************************************
       	usage: ngscat2 --bams <filename> --bed <filename> --out <path> --annotation <filename> --reference <filename>  --tmp <path> --threads <integer>

Options:
  -h, --help            show this help message and exit
  --bams=BAMS           Required. Comma separated list of bam files (2
                        maximum). E.g.: --bams
                        /home/user/bam1.sorted.bam,/home/user/bam2.sorted.bam
  --bed=BED             Required. Full path to the bed file containing the
                        target regions.
  --out=OUTDIR          Required. Full path to the directory where results
                        will be saved.
  --reference=REFERENCE
                        Optional. String indicating the path to a .fasta file
                        containing the reference chromosomes. Default=None.
  --annotation=ANNOTATION
                        Optional. String indicating the path to a .bed file
                        containing annotated regions . Default=None.
  --coveragethrs=COVERAGETHRESHOLDS
                        Optional. Comma separated list of real numbers (do not
                        leave spaces between) indicating coverage thresholds
                        to be used when calculating percentages of covered
                        bases (first graph in the report).
                        Default=1,5,10,20,30.
  --tmp=TMP             Optional. String indicating the full path to a
                        temporary directory where temporary files will be
                        created. Default=/tmp/.
  --threads=NTHREADS    Optional. Integer indicating the number of concurrent
                        threads to launch. Default=cpu_count() - 1.
```

# Input data
Here an example of input data can he downloaded.

* Hg0097 alignment data BAM http://ngscat.clinbioinfosspa.es/_media/ngscat/download/hg00097.bam  
* Region of interest file BED ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/exome_pull_down_targets/20130108.exome.targets.bed
* Reference genome (for --reference option) http://ngscat.clinbioinfosspa.es/_media/ngscat/download/hg19.wochr.tar.gz
* Annotated bed (for --annotation option)  https://mega.nz/#!f2oRxI7Z!QUSRErtvwNGs2tfpwRWjiPZQkOEhRQhaWtXhkmXXpEM
#Output example

Exome example https://mega.nz/#!eywTGCBL!6HTDl7J9eLY4VhX2aFfe57Cp9mvwLw4JTcDcHr-zu1A
