#!/bin/bash

#PBS -N Contaminant_RM
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=128gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


usage="\nThis script is to remove the potential contaminant scaffolds from the genome assembly, \
by searching the databases from human, bacteria and virus but retaining the scaffolds \
that have high similarity to the closely related species. \e[31mEdit the script before use!!!\e[39m \
You must ensure the Deconseq program has been correctly installed and the corresponding \
databases have been downloaded before running the script.

$(basename "$0") [-h] 
where:
    -h  show this help text\n\n"


###################################### USER DEFINED AREA ###########################################
## PATH of required programs, need to make sure the binary program and scripts used is executable
PATH=$PATH:/N/u/wum5/Carbonate/softwares/deconseq-0.4.3    # PATH to prinseq
WORKDIR=/N/dc2/projects/solanumgenome/sapp/drafts  # the directory for searching contaminants
GENOME=/N/dc2/projects/solanumgenome/sapp/drafts/sapp_1.0/sapp.fa  # input genome file

## parameters need to be set by the users
dbs_remove=human,bacteria,viral  # remove contaminants from human,bacteria and viral
dbs_retain=solanum  # retain scaffolds potentially from solanum species
ThreadN=8  # number of CPUs for parallel computing
alignidt=80  # alignment identity threshold in percentage	
aligncov=50  # alignment coverage threshold in percent
####################################################################################################


while getopts ':h' option; do
  case "$option" in
    h) printf "$usage"
       exit
       ;;
  esac
done


########## run DeconSeq to remove contaminants ########## 
cd $WORKDIR
if [ ! -d contaminants ]; then mkdir contaminants; fi
cd contaminants

genomeID=$(basename $GENOME ".fa") 
cp $GENOME $genomeID
splitFasta.pl -i $genomeID -n 32

function run_deconseq(){
dir=$(basename ${1} .fasta)  # {1}=input file
mkdir $dir
mv ${1} $dir
cd $dir  # ${2}=dbs_remove, ${3}=dbs_retain
deconseq.pl -f ${1} -dbs ${2} -dbs_retain ${3} -i ${alignidt} -c ${aligncov} -out_dir .
cd ..
}

export -f run_deconseq
find *fasta | xargs -n 1 | xargs -P $ThreadN -I{} bash -c 'run_deconseq {}' ${dbs_remove} ${dbs_retain}

cat ${genomeID}_c*/*_clean.fa > ${genomeID}"_clean.fa"
cat ${genomeID}_c*/*_cont.fa > ${genomeID}"_cont.fa"
cat ${genomeID}_c*/*_both.fa > ${genomeID}"_both.fa" 


