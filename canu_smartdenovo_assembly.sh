#!/bin/bash
set -e
set -u
set -o pipefail


function usage(){
printf "\nDescription: this script is to assemble the genome using a hybrid-assembly approach. \
This pipeline can be splited into 6 major steps: 1: Canu correction, 2: SMARTdenovo assembly, \
3-5: pilon polish (3 iterations), and 6. scaffolding using SSPACE. Required software must be \
installed and added into the PATH, including Canu, Smart-denovo, BWA, Bowtie, Pilon, and SSPACE.  

$(basename "$0") [-h] [-OPTIONS] -d <working dir> -s <genome size> -P <PacBio reads> \
-f <foward reads> -r <reverse reads> -I <assembly prefix> 

where:
   -h  show this help text
   -d  the working directory for gene annotation
   -s  estimated genome size (e.g. 1.2g)
   -P  input PacBio reads 
   -I  the library file listing the illumina reads information
   -p  assembly prefix (no space in prefix)
   -n  number of threads for parallel processing (default=16)
   -S  the steps to processed (default=123456): 1: Canu; 2: SMARTdenovo; 3: Pilon1; \
   4: Pilon2; 5: Pilon3; 6: SSPACE \n\n"
}


######### Default parameters #########
WORKDIR=''
GenomeSZ=''
PACBIO=''
ILLUMINA=''
Prefix=''
ThreadN=16
ProcStep=123456


######### Parse input #########
while getopts 'h:d:s:P:I:p:n:S:' option; do
  case "$option" in
    h) usage
       exit
       ;;
    d) WORKDIR=${OPTARG}
       ;;
    s) GenomeSZ=${OPTARG}
       ;;
    P) PACBIO=${OPTARG}
       ;;
    I) ILLUMINA=${OPTARG}
       ;;
    p) Prefix=${OPTARG}
       ;;
    n) ThreadN=${OPTARG}
       ;;
    S) ProcStep=${OPTARG}
       ;;
    *) usage
       exit
       ;;
  esac
done
shift $((OPTIND-1))



if [ -z "${WORKDIR}" ] || [ -z "${GenomeSZ}" ] || [ -z "${PACBIO}" ] \
|| [ -z "${ILLUMINA}" ] || [ -z "${Prefix}" ] ; 
then
    usage
    exit
fi

cd ${WORKDIR}



############ Canu Correction ############
if [[ $ProcStep = *"1"* ]]; 
then
	if [ ! -d canu_correction ]; then mkdir canu_correction; fi
	cd canu_correction
	canu -correct -p ${Prefix} -d . genomeSize=${GenomeSZ} \
	merylMemory=48g merylThreads=6 batMemory=48g batThreads=4 \
	gridOptions="-l vmem=64gb,walltime=4:00:00" -pacbio-raw ${PACBIO}
	gunzip ${Prefix}.correctedReads.fasta.gz
	cd ..
fi



############ SMARTdenovo Assembly ############
if [[ $ProcStep = *"2"* ]]; 
then
	if [ ! -d smartdenvo_assembly ]; then mkdir smartdenvo_assembly; fi
	cp canu_correction/${Prefix}.correctedReads.fasta > smartdenvo_assembly/
	
	cd smartdenvo_assembly	
	smartdenovo.pl –c 1 –t ${ThreadN} –k 17 –p ${Prefix} \
	${Prefix}.correctedReads.fasta > ${Prefix}.make
	
	make –f ${Prefix}.make
	cd ..
fi



############ Pilon polish: 1st iteration ############
if [[ $ProcStep = *"3"* ]]; 
then
	if [ ! -d pilon_polish1 ]; then mkdir pilon_polish1; fi	
	cp smartdenvo_assembly/${Prefix}.dmo.cns pilon_polish/draft0.fasta
	cd pilon_polish1
	
	bwa index draft0.fasta
	while read -r line; 
	do
		arr=(echo $line)
		FORWARD=${arr[2]}
		REVERSE=${arr[3]}
		ReadPre=$(basename ${FORWARD} _R1.fq)
		
		bwa mem draft0.fasta ${FORWARD} ${REVERSE} -t ${ThreadN} | \
		samtools sort -@ ${ThreadN} -O BAM -m 10G -o ${ReadPre}.sort.bam - 
		samtools flagstat ${ReadPre}.sort.bam > ${ReadPre}_mapstats.txt    	
	done < ${ILLUMINA}
	
	samtools merge -@ ${ThreadN} pilon1.bam *.sort.bam 
	samtools index pilon1.bam

	mkdir pilon_data
	cp draft0.fasta pilon_data/draft
	cd pilon_data
	splitFasta.pl -i draft -s 50
	rm draft
		
	for file in *.fasta
	do
		filename=$(basename $file .fasta) 
		grep '>' ${file} | cut -f1 -d " " | tr --delete ">" | tr "\n" "," | \
		head -c -1 > ${filename}.list
	done
	
	rm *.fasta
		
	for file in *.list
	do
		prefix=$(basename $file .list)
		java -Xmx100g -jar ${PILON}/pilon-1.22.jar --genome ../draft0.fasta --jumps ../pilon1.bam \
		--output ${prefix}_pilon --targets ${prefix}.list --threads ${ThreadN}
	done
	
	cat *_pilon.fasta > ../pilon1.fasta
	cd ..
fi


############ Pilon polish: 2nd iteration ############
if [[ $ProcStep = *"4"* ]]; 
then
	if [ ! -d pilon_polish2 ]; then mkdir pilon_polish2; fi
	cp pilon_polish1/pilon1.fasta pilon_polish2/draft1.fasta
	
	cd pilon_polish2	
	bwa index draft1.fasta
	
	while read -r line; 
	do
		arr=(echo $line)
		FORWARD=${arr[2]}
		REVERSE=${arr[3]}
		ReadPre=$(basename ${FORWARD} _R1.fq)
		
		bwa mem draft1.fasta ${FORWARD} ${REVERSE} -t ${ThreadN} | \
		samtools sort -@ ${ThreadN} -O BAM -m 10G -o ${ReadPre}.sort.bam - 2> ${ReadPre}.error 
		samtools flagstat ${ReadPre}.sort.bam > ${ReadPre}_mapstats.txt   	
	done < ${ILLUMINA}
	
	samtools merge -@ ${ThreadN} pilon2.bam *.sort.bam 
	samtools index pilon2.bam

	mkdir pilon_data
	cp draft1.fasta pilon_data/draft
	cd pilon_data
	splitFasta.pl -i draft -s 50
	rm draft
		
	for file in *.fasta
	do
		filename=$(basename $file .fasta) 
		grep '>' ${file} | cut -f1 -d " " | tr --delete ">" | tr "\n" "," | \
		head -c -1 > ${filename}.list
	done
	
	rm *.fasta
		
	for file in *.list
	do
		prefix=$(basename $file .list)
		java -Xmx100g -jar ${PILON}/pilon-1.22.jar --genome ../draft1.fasta --jumps ../pilon2.bam \
		--output ${prefix}_pilon --targets ${prefix}.list --threads ${ThreadN}
	done
	
	cat *_pilon.fasta > ../pilon2.fasta
	cd ..
fi



############ Pilon polish: 3rd iteration ############
if [[ $ProcStep = *"5"* ]]; 
then
	if [ ! -d pilon_polish3 ]; then mkdir pilon_polish3; fi
	cp pilon_polish2/pilon2.fasta pilon_polish3/draft2.fasta
	
	cd pilon_polish3	
	bwa index draft2.fasta
		
	while read -r line; 
	do
		arr=(echo $line)
		FORWARD=${arr[2]}
		REVERSE=${arr[3]}
		ReadPre=$(basename ${FORWARD} _R1.fq)
		
		bwa mem draft2.fasta ${FORWARD} ${REVERSE} -t ${ThreadN} | \
		samtools sort -@ ${ThreadN} -O BAM -m 10G -o ${ReadPre}.sort.bam - 2> ${ReadPre}.error 	
		samtools flagstat ${ReadPre}.sort.bam > ${ReadPre}_mapstats.txt  	
	done < ${ILLUMINA}
	
	samtools merge -@ ${ThreadN} pilon3.bam *.sort.bam 
	samtools index pilon3.bam
	
	mkdir pilon_data
	cp draft2.fasta pilon_data/draft
	cd pilon_data
	splitFasta.pl -i draft -s 50
	rm draft
		
	for file in *.fasta
	do
		filename=$(basename $file .fasta) 
		grep '>' ${file} | cut -f1 -d " " | tr --delete ">" | tr "\n" "," | \
		head -c -1 > ${filename}.list
	done
	
	rm *.fasta
		
	for file in *.list
	do
		prefix=$(basename $file .list)
		java -Xmx100g -jar ${PILON}/pilon-1.22.jar --genome ../draft2.fasta --jumps ../pilon3.bam \
		--output ${prefix}_pilon --targets ${prefix}.list --threads ${ThreadN}
	done

	cat *_pilon.fasta > ../pilon3.fasta
	cd ..
fi



############ Scaffolding using SSPACE ############
if [[ $ProcStep = *"6"* ]]; 
then
	if [ ! -d sspace_scaffold ]; then mkdir sspace_scaffold; fi
	cp pilon_polish3/pilon3.fasta sspace_scaffold/${Prefix}_draft3.fasta
	cd sspace_scaffold
	SSPACE_Basic.pl -T ${ThreadN} -l ${ILLUMINA} -s ${Prefix}_draft3.fasta
	cd ..
fi
	
	
	
