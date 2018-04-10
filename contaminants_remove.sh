#!/bin/bash

#PBS -N Contaminant_RM
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=128gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

usage="\nThis script is to remove the potential contaminant scaffolds from the genome assembly, \
by searching the databases from human, bacteria and virus but retaining the scaffolds \
that have high similarity to the closely related species. \e[31mEdit the script before use!!!\e[39m \
You must ensure the Deconseq program has been correctly installed and the corresponding \
databases have been downloaded before running the script. Once the contaminant and retention \
databases have been set, set 'SetConfig=1'  and rerun the script.

$(basename "$0") [-h] 
where:
    -h  show this help text\n"

while getopts ':h' option; do
  case "$option" in
    h) printf "$usage"
       exit
       ;;
  esac
done

## PATH of required programs, need to make sure the binary program and scripts used is executable
PATH=$PATH:/N/u/wum5/Carbonate/softwares/prinseq-lite-0.20.4
PATH=$PATH:/N/u/wum5/Carbonate/softwares/deconseq-0.4.3
LIBDIR=/N/dc2/projects/solanumgenome/library   ## the directory for preparing searching datasets
WORKDIR=/N/dc2/projects/solanumgenome/sapp/drafts  ## the directory for searching contaminants
GENOME=${WORKDIR}'/sapp_1.0/sapp.fa'  ## input genome file

## parameters need to be set by the users
dbs_remove=human,bacteria,viral  ## remove contaminants from human,bacteria and viral
dbs_retain=solanum  ## retain scaffolds potentially from solanum species
SetConfig=0  ## change this to YES if the DeconSeqConfig.pm has been set, otherwise keep it as NO
ThreadN=8  ## number of CPUs for parallel computing

cd $LIBDIR
genomeID=$(basename $GENOME ".fa") 
sample_names=()


########## Contaminant database ##########
if [ $SetConfig -eq 0 ]
then
	if [ ! -d deconseq_DB ]; then mkdir deconseq_DB; fi
	cd deconseq_DB

	# human genome
	wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
	gunzip GRCh38_latest_genomic.fna.gz
	mv GRCh38_latest_genomic.fna human.fasta
	sample_names+="human "

	## a function to download genomic files
	function download_genome(){
	filename=$(basename ${1})'_genomic.fna.gz'
	wget -P ${2} ${1}'/'$filename
	}
	export -f download_genome

	## viral genomes
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt
	awk -F '\t' '{if($12=="Complete Genome") print $20}' assembly_summary.txt > viral_genomes.txt
	rm assembly_summary.txt
	mkdir GbVir
	cat viral_genomes.txt | xargs -n 1 | xargs -P $ThreadN -I{} bash -c 'download_genome {} GbVir'
	find GbVir/*.gz | xargs -n 1 | xargs -P $ThreadN -I{} gunzip {}
	cat GbVir/*.fna > viral.fasta
	sample_names+="viral "

	## bacterial genomes
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
	awk -F '\t' '{if($12=="Complete Genome") print $20}' assembly_summary.txt > bacteria_genomes.txt
	rm assembly_summary.txt
	mkdir GbBac
	cat bacteria_genomes.txt | xargs -n 1 | xargs -P $ThreadN -I{} bash -c 'download_genome {} GbBac'
	find GbBac/*.gz | xargs -n 1 | xargs -P $ThreadN -I{} gunzip {}
	cat GbBac/*.fna > bacteria.fasta
	mv bacteria.fasta bacteria
	splitFasta.pl -verbose -i bacteria -n 10
	sample_names+="bacteria "
fi


########## Retention database (this is for Solanum species; change this to your system) ##########
if [ $SetConfig -eq 0 && $dbs_retain == 'solanum' ]
then
	## S.lycopersicum
	wget ftp://ftp.solgenomics.net/genomes/Solanum_lycopersicum/annotation/ITAG2.4_release/S_lycopersicum_chromosomes.2.50.fa.gz
	gunzip S_lycopersicum_chromosomes.2.50.fa.gz 
	mv S_lycopersicum_chromosomes.2.50.fa.gz lyco.fasta
	sample_names+="lyco "

	## S.pennellii
	wget ftp://ftp.solgenomics.net/genomes/Solanum_pennellii/penn.fasta
	sample_names+="penn "

	## S.tuberosum
	wget ftp://ftp.solgenomics.net/genomes/Solanum_tuberosum/assembly/PGSC_DM_v4.03/PGSC_DM_v4.03_pseudomolecules.fasta.zip
	wget ftp://ftp.solgenomics.net/genomes/Solanum_tuberosum/assembly/PGSC_DM_v4.03/PGSC_DM_v4.03_unanchored_regions_chr00.fasta.zip
	unzip PGSC_DM_v4.03_pseudomolecules.fasta.zip
	unzip PGSC_DM_v4.03_unanchored_regions_chr00.fasta.zip
	cat PGSC_DM_v4.03_unanchored_regions_chr00.fasta PGSC_DM_v4.03_pseudomolecules.fasta > tube.fasta
	sample_names+="tube "
fi


########## Prepare the searching databases ########## 
if [ $SetConfig -eq 0 ]
then
	echo $sample_names

	function seq_split(){
	cat ${1}'.fasta' | perl -p -e 's/N\n/N/' | \
	perl -p -e 's/^N+//;s/N+$//;s/N{200,}/\n>split\n/' > ${1}'_split.fa'
	}
	export -f seq_split

	if [ "$sample_names" != '' ]
	then
		echo $sample_names | xargs -n 1 | xargs -P $ThreadN -I{} bash -c 'seq_split {}'

		echo $sample_names | xargs -n 1 | xargs -P $ThreadN -I{} -I{} prinseq-lite.pl \
		-log -verbose -fasta {}'_split.fa' -min_len 200 -ns_max_p 10 -derep 12345 \
		-out_good {}'_split_prinseq' -seq_id {} -rm_header -out_bad null

		echo $sample_names | xargs -n 1 | xargs -P $ThreadN -I{} bwa64 index -p {} \
		-a bwtsw {}'_split_prinseq.fasta'
	fi
fi


########## Edit the DeconSeqConfig.pm in the "deconseq" software folder (example below) ########## 
########## After editing parameters, change "SetConfig" to "YES", and rerun the script ##########

#use constant DB_DIR => '/N/dc2/projects/solanumgenome/library/deconseq_DB/';
#use constant TMP_DIR => '/N/dc2/projects/solanumgenome/sapp/drafts/contaminants_removed/tmp/';
#use constant OUTPUT_DIR => '/N/dc2/projects/solanumgenome/sapp/drafts/contaminants_removed/output/';

#use constant PROG_NAME => 'bwa64';  # should be either bwa64 or bwaMAC (based on your system architecture)
#use constant PROG_DIR => '';      # should be the location of the PROG_NAME file (use './' if in the same location at the perl script)

#use constant DBS => {human => {name => 'Human Reference GRCh38',  #database name used for display and used as input for -dbs and -dbs_retain
#                               db => 'human'},            #database name as defined with -p for "bwa index -p ..." (specify multiple database chunks separated with commas without space; e.g. hs_ref_s1,hs_ref_s2,hs_ref_s3)
#                     bact => {name => 'Bacterial genomes',
#                              db => 'bacteria_c1,bacteria_c2,bacteria_c3,bacteria_c4,bacteria_c5,bacteria_c6,bacteria_c7,bacteria_c8,bacteria_c9,bacteria_c10'},
#                     vir => {name => 'Viral genomes',
#                             db => 'viral'},
#                     solanum => {name => 'Solanum genomes',
#                             db => 'lyco,tube,penn'}};
#use constant DB_DEFAULT => 'human';

## exit the script/job if the DeconSeqConfig.pm has not yet been set
if [ $SetConfig -eq 0 ]; then echo "DeconSeqConfig.pm has not yet been set" && exit; fi


########## run DeconSeq to remove contaminants ########## 
cd $WORKDIR
if [ ! -d contaminants ]; then mkdir contaminants; fi
cd contaminants

cp $GENOME $genomeID
splitFasta.pl -i $genomeID -n 32

function run_deconseq(){
dir=$(basename ${1} .fasta)  # {1}=input file
mkdir $dir
mv ${1} $dir
cd $dir  # ${2}=dbs_remove, ${3}=dbs_retain
deconseq.pl -f ${1} -dbs ${2} -dbs_retain ${3} -i 80 -c 50 -out_dir .
cd ..
}

export -f run_deconseq
find *fasta | xargs -n 1 | xargs -P $ThreadN -I{} bash -c 'run_deconseq {}' ${dbs_remove} ${dbs_retain}

cat ${genomeID}_c*/*_clean.fa > ${genomeID}"_clean.fa"
cat ${genomeID}_c*/*_cont.fa > ${genomeID}"_cont.fa"
cat ${genomeID}_c*/*_both.fa > ${genomeID}"_both.fa" 


