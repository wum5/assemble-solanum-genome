# assemble-solanum-genome

## Table of Contents

* [Overview](#overview)
* [Contributors](#contributors)
* [De novo Asseemble Genome](#de-novo-assemble-genome)
* [Genome Annotation](#genome-annotation)

## Overview
* This repository is to record the scripts to assemble and annotate some plant genomes in the genus Solanum. 
* The scripts are aimed to be reproducible that can be used in other plant genome projects. 
* You use -h to see the detailed information of how to run each bash script.
* We are still updating the scripts and this pipeline! 

## Contributors 
* Meng Wu
* Rafael Guerrero

## De novo Assemble Genome
##### Set the environment to add required programs to PATH
```
./set_path.sh
```
##### Assemble genome using MaSuRCA approach
```
./masurca_assembly.sh -d <working dir> -c <masurca configfile> -P <PacBio reads> \
-p <"pair-ended reads information"> -n 16 -K 2 -J 3000000000
```
##### Build contaminant databases (skip this step if you have deconseq databases being set)
```
./contaminants_dbbuild.sh -d <working dir> -n 16
```
##### Remove contaminant scaffolds
```
./contaminants_remove.sh -d <working dir> -g <genome file> -r <retention database> \
-R human,bacteria,viral -n 16 -i 80 -c 50
```
## Genome Annotation
##### Prepare species-specific repeats library
```
./repeats_annotation.sh -d <working dir> -g <genome file> -T <Tpases020812DNA> \
-P <alluniRefprexp070416> -n 16 -S 1234 
```

