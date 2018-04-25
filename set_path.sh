#!/bin/bash
# a bash script to add the required programs to the PATH 


### Genome Assembly
PATH=$PATH:/N/u/wum5/Carbonate/softwares/masucra-3.2.2/bin  # PATH to masucra binary files


### Remove Contaminants by DeconSeq
PATH=$PATH:/N/u/wum5/Carbonate/softwares/deconseq-0.4.3    # PATH to deconseq
PATH=$PATH:/N/u/wum5/Carbonate/softwares/prinseq-lite-0.20.4  # PATH to prinseq


### required programs for building species-specific repeats library
module load cd-hit/4.6.8   # PATH to cd-hit
module load repeatmasker/4.0.7  # PATH to repeatmasker
PATH=$PATH:/N/u/wum5/Carbonate/softwares/hmmer-3.1b2/binaries  # PATH to hmmer
PATH=$PATH:/N/u/wum5/Carbonate/softwares/mite-hunter  # PATH to mite-hunter 
PATH=$PATH:/N/u/wum5/Carbonate/softwares/ltr_retriever  # PATH to ltr-retriever
PATH=$PATH:/N/u/wum5/Mason/bin/gt/bin   # PATH to ltr-harvest
PATH=$PATH:/N/u/wum5/Carbonate/softwares/deconseq-0.4.3  # PATH to deconseq
PATH=$PATH:/N/u/wum5/Carbonate/softwares/crl-scripts-1.0  # PATH to CRL scripts used for repeat annotation
PATH=$PATH:/N/u/wum5/Carbonate/softwares/repeat-modeler-1.0.11  # PATH to repeatmodeler
PATH=$PATH:/N/u/wum5/Carbonate/softwares/ncbi-blast-2.7.1+/bin  # PATH to blast 
PATH=$PATH:/N/u/wum5/Carbonate/softwares/ProtExcluder1.2  # PATH to ProtExcluder


### Predict gene models by MAKER2
module load python
module load bowtie2/intel/2.3.2 
module load tophat/2.1.1 
module load trinityrnaseq/2.4.0
module load maker/2.31.9

PATH=$PATH:/N/u/wum5/Carbonate/softwares/busco/scripts
PATH=$PATH:/gpfs/home/w/u/wum5/Carbonate/softwares/augustus.2.5.5/bin
PATH=$PATH:/N/u/wum5/Carbonate/softwares/gm_et_linux_64/gmes_petap
export AUGUSTUS_CONFIG_PATH=/gpfs/home/w/u/wum5/Carbonate/softwares/augustus.2.5.5/config/
AUGUSTUS_DB=/gpfs/home/w/u/wum5/Carbonate/softwares/augustus.2.5.5/config/species/