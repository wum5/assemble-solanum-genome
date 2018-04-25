#!/bin/bash
# a bash script to assemble genome using MaSuRCA 


function usage(){
printf "\nThis script is to build the initial genome assembly with PacBio long reads and \
and Illumina short reads by using MaSuRCA program. You must ensure the MaSuRCA program has \
been correctly installed and added to the PATH before running the script. 

$(basename "$0") [-h] [-OPTIONS] -d <working dir> -c <masurca configfile> -P <PacBio reads> \
-p <'pair-ended reads information'>

where:
    -h  show this help text
    -d  the working directory for assembly
    -c  PATH to masurca config file 
    -P  PATH to pacbio reads
    -p  'two_letter_prefix mean stdev /PATH/fwd_reads.fastq /PATH/rev_reads.fastq' (paired-end reads)
    -j  'two_letter_prefix mean stdev /PATH/fwd_reads.fastq /PATH/rev_reads.fastq' (mate-paired reads; optional)
    -O  PATH to the file with other information (optional)
    -n  number of cpus for parallel processing (default=16)
    -K  increase to 2 if illumina coverage >100 (default=1)
    -J  a safe value to be estimated_genome_size*estimated_coverage (default=1000000000)\n\n"
}


######### Default parameters #########
WORKDIR=''  
CONFIGFILE=''  
PACBIO=''
PE=''
JUMP='' 
OTHER=''
NUM_THREADS=16
KMER_COUNT=1
JF_SIZE=1000000000 


######### Parse input #########
while getopts ':h:d:c:P:p:j:O:n:K:J' option; do
  case "$option" in
    h) usage
       exit
       ;;
    d) WORKDIR=${OPTARG}
       ;;
    c) CONFIGFILE=${OPTARG}
       ;;
    P) PACBIO=${OPTARG}
       ;;
    p) PE=${OPTARG}
       ;;
    j) JUMP=${OPTARG}
       ;;
    O) OTHER=${OPTARG}
       ;;
    n) NUM_THREADS=${OPTARG}
       ;;
    K) KMER_COUNT=${OPTARG}
       ;;
    J) JF_SIZE=${OPTARG}
       ;;
    *) usage
       exit
       ;;
  esac
done
shift $((OPTIND-1))


if [ -z "${WORKDIR}" ] || [ -z "${CONFIGFILE}" ] || [ -z "${PACBIO}" ] || [ -z "${PE}" ] ; then
    usage
    exit
fi

if [ $(echo "${PE}" | wc -w) != 5 ]; then 
    usage
    exit
fi



######### A function to edit CONFIG file ######### 
function set_config(){
CONFIGFILE=$1
TARGET_KEY=$2
REPLACEMENT_VALUE=$3
sed -i "s|^\($TARGET_KEY\s*=\s*\).*\$|\1$REPLACEMENT_VALUE|" $CONFIGFILE
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
masurca config_spol.txt
./assemble.sh

