# assemble-solanum-genome

## Table of Contents

* [Overview](#overview)
* [Contributors](#contributors)
* [De novo Asseemble Genome](#de-novo-assemble-genome)

## Overview
* This repository is to record the scripts to assemble and annotate some plant genomes in the genus Solanum. 
* The scripts are aimed to be reproducible that can be used in other plant genome projects. 
* detailed information of how to run this step are recorded in the corresponding bash script
* We are still updating the scripts and this pipeline! 

## Contributors 
* Meng Wu
* Rafael Guerrero

## De novo Assemble Genome
##### Assemble genome using MaSuRCA approach
```
qsub masurca_assembly.sh
```
##### Remove contaminant scaffolds
```
qsub contaminants_remove.sh
```