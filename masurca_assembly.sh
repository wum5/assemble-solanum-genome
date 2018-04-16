#!/bin/bash

#PBS -N MASURCA
#PBS -l nodes=1:ppn=24,walltime=320:00:00,vmem=320gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

usage="\nThis script is to build the initial genome assembly with PacBio long reads and \
and Illumina short reads by using MaSuRCA program. \e[31mEdit the script before use!!!\e[39m \
You must ensure the MaSuRCA program has been correctly installed and added to the PATH \
before running the script. All the parameters need to be set before running the script. 

$(basename "$0") [-h] 
where:
    -h  show this help text\n\n"


###################################### USER DEFINED AREA ###########################################
######## PATH of required programs (the binary program and scripts used should be executable) ########
WORKDIR=/N/dc2/projects/solanumgenome/spol  # the working directory for assembly
CONFIGFILE=config_spol.txt  # masurca config file (it should be in the working directory)
PATH=$PATH:/N/u/wum5/Carbonate/softwares/masucra-3.2.2/bin  # PATH to masucra binary files

######## Parameters setting (leave it blank if you don't have) ########
PACBIO=/N/dc2/projects/solanumgenome/spol/reads/pac/all_subreads.fa  # PATH to pacbio reads
JUMP=  # PATH to mate paired (jump) Illumia reads
PE="pe 180 35 /N/dc2/projects/solanumgenome/spol/reads/illumina/raw/malfempool_R1.fq /N/dc2/projects/solanumgenome/spol/reads/illumina/raw/malfempool_R2.fq"  # PE = two_letter_prefix mean stdev /PATH/fwd_reads.fastq /PATH/rev_reads.fastq
OTHER=  # PATH to the file with other information
NUM_THREADS=24  # number of cpus for parallel processing 
KMER_COUNT=2  # default=1, increase to 2 if illumina coverage >100
JF_SIZE=3000000000  # a safe value to be estimated_genome_size*estimated_coverage
####################################################################################################


while getopts ':h' option; do
  case "$option" in
    h) printf "$usage"
       exit
       ;;
  esac
done


######### A function to edit CONFIG file ######### 
function set_config(){
CONFIGFILE=$1
TARGET_KEY=$2
REPLACEMENT_VALUE=$3
sed -i "s|^\($TARGET_KEY\s*=\s*\).*\$|\1$REPLACEMENT_VALUE|" $CONFIGFILE
sed -i 's/\[s\]/ /g' $CONFIGFILE
}


######### Edit the config file of MaSuRCA #########
cd $WORKDIR

set_config ${CONFIGFILE} "PACBIO" "$PACBIO"
set_config ${CONFIGFILE} "JUMP" "$JUMP"
set_config ${CONFIGFILE} "OTHER" "$OTHER"
set_config ${CONFIGFILE} "PE" "$PE"

set_config ${CONFIGFILE} "NUM_THREADS" "$NUM_THREADS"
set_config ${CONFIGFILE} "KMER_COUNT_THRESHOLD" "$KMER_COUNT"
set_config ${CONFIGFILE} "JF_SIZE" "$JF_SIZE"

sed -i '/= $/ d' ${CONFIGFILE}
sed -i '/=$/ d' ${CONFIGFILE}


######### Run MasuRCA #########
#masurca config_spol.txt
#./assemble.sh

