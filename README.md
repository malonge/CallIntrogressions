# CallIntrogressions

## Description
The code in this repository was used to find putative introgressions in modern tomato genomes. This work is a part of a larger
panSV-genome analysis that is currently unpublished. This code is written expressly for this analysis, and is will not be
updated or maintaiened for future use.

## Dependencies
- python3
- numpy (developed, and run on v1.17.2)

## Installation and Usage
There is no installation required. Just run `python3 get_distances.py` to see the following usage message. 

```
usage: get_distances.py [-h] [-m 5]
                        <SVs.vcf> <chr_name> <group.txt> <SP>
                        <reference.fasta.fai> <100000>

Get Jaccard similarity between SLL and a comparison group.

positional arguments:
  <SVs.vcf>             SV vcf file with support vectors. Only one chromosome
                        at a time allowed.
  <chr_name>            Name of reference chromosome.
  <group.txt>           First column is the phylogenetic group (SLC, SP, GAL,
                        CHE, or SLL), second column is the accession
  <SP>                  Group to compare to SLL (SLC, SP, GAL, or CHE)
  <reference.fasta.fai>
                        Fasta index file for the reference genome
  <100000>              Introgression window size.

optional arguments:
  -h, --help            show this help message and exit
  -m 5                  minimum number of SVs needed to calculate Jaccard
```
