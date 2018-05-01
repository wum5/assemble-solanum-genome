# assemble-solanum-genome


## Table of Contents

* [Overview](#overview)
* [Contributors](#contributors)
* [De novo Asseemble Genome](#de-novo-assemble-genome)
* [Genome Annotation](#genome-annotation)


## Overview
* This repository is to record the scripts to assemble and annotate some plant genomes in the genus Solanum. 
* The scripts are aimed to be reproducible that can be used in other plant genome projects. 
* You could use -h to see the detailed information of how to run each bash script.
* We are still updating the scripts and this pipeline! 


## Contributors 
* Meng Wu
* Rafael Guerrero


## De novo Assemble Genome
##### Set the environment to add required programs to PATH
```
set_path.sh
```

##### Assemble genome using MaSuRCA approach
```
masurca -g masurca_config.txt   # run this command in the working directory
masurca masurca_config.txt   # need to manually edit "masurca_config" before to set required parameters

assemble.sh
```

##### Build contaminant databases (skip this step if you have deconseq databases being set)
```
contaminants_dbbuild.sh -d <working dir> -n 16   # need to manually edit "DeconSeqConfig.pm" then
```

##### Remove contaminant scaffolds
```
contaminants_remove.sh -d <working dir> -g <genome file> -r <retention database> \
-R human,bacteria,viral -n 16 -i 80 -c 50
```


## Genome Annotation
##### Prepare species-specific repeats library
```
repeats_annotation.sh -d <working dir> -g <genome file> -T <Tpases020812DNA> \
-P <alluniRefprexp070416> -n 16 -S 1234 
```

##### Annotate genome using MAKER2 pipeline
```
maker -EXE   # need to manually edit the file to get access to PATH of required programs then

gene_annotation.sh -d <working dir> -g <genome file> -f <forward RNA-seq> -r <reverse RNA-seq> \
-p <homologous protein sequences> -t <TE proteins sequences> -R <Repeats library> \
-B <BUSCO embryophyta directory> -M maker_exe.ctl -n 16 -S 1234
```

##### Check annotation statistics and pull out coding sequence with AED score <0.5
```
annotation_by_aed.py --genome <genome file> --gff <GFF output from MAKER2> --AED 0.5
```

