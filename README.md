# assemble-solanum-genome

## Table of Contents

* [Overview](#overview)
* [Contributors](#contributors)
* [De novo Asseemble Genome](#de-novo-assemble-genome)

## Overview
* This repository is to record the scripts to assemble and annotate some plant genomes in the genus Solanum. 
* The scripts are aimed to be reproducible that can be used in other plant genome projects. 
* Detailed information of how to run each step are recorded in the corresponding bash script (use -h to see).
* You need to edit the parameters and program PATH in each script before use.
* We are still updating the scripts and this pipeline! 

## Contributors 
* Meng Wu
* Rafael Guerrero

## De novo Assemble Genome
##### Assemble genome using MaSuRCA approach
```
sh masurca_assembly.sh -h
qsub masurca_assembly.sh
```
##### Remove contaminant scaffolds
```
sh contaminants_remove.sh -h
qsub contaminants_remove.sh
```