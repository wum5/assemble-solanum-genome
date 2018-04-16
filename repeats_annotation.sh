#!/bin/bash

#PBS -N REPEAT-ANNOT
#PBS -l nodes=1:ppn=16,walltime=128:00:00,vmem=128gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


usage="\nThis script is to prepare species-specific library before performing gene \
annotation. \e[31mEdit the script at USER DEFINED AREA before use!!!\e[39m \
For example, you can control the steps need to be processed by editing the 'ProcStep'. \
You must also ensure the following programs have been correctly installed and added \
into PATH: MITE-hunter, LTR-harvest, LTR-retriever, Repeat-Modeler, Repeat-Masker, \
Deconseq, CD-HIT, Hmmer, Blast and ProtExcluder. 

$(basename "$0") [-h] 
where:
    -h  show this help text\n\n"

while getopts ':h' option; do
  case "$option" in
    h) printf "$usage"
       exit
       ;;
  esac
done    


###################################### USER DEFINED AREA ###########################################    
## PATH of required programs, need to make sure the binary program and scripts used is executable
module load hmmer/3.1  # PATH to hmmer
module load cd-hit/4.6.8   # PATH to cd-hit
module load repeatmasker/4.0.7  # PATH to repeatmasker
PATH=$PATH:/N/u/wum5/Carbonate/softwares/mite-hunter  # PATH to mite-hunter 
PATH=$PATH:/N/u/wum5/Carbonate/softwares/ltr_retriever  # PATH to ltr-retriever
PATH=$PATH:/N/u/wum5/Mason/bin/gt/bin   # PATH to ltr-harvest
PATH=$PATH:/N/u/wum5/Carbonate/softwares/deconseq-0.4.3  # PATH to deconseq
PATH=$PATH:/N/u/wum5/Carbonate/softwares/crl-scripts-1.0  # PATH to CRL scripts used for repeat annotation
PATH=$PATH:/N/u/wum5/Carbonate/softwares/repeat-modeler-1.0.11  # PATH to repeatmodeler
PATH=$PATH:/N/u/wum5/Carbonate/softwares/ncbi-blast-2.7.1+/bin  # PATH to blast 
PATH=$PATH:/N/u/wum5/Carbonate/softwares/ProtExcluder1.2  # PATH to ProtExcluder

WORKDIR=/N/dc2/projects/solanumgenome/sapp/drafts/repeats  # the working directory for assembly
seqfile=sapp_clean.fa  # the genome fasta file
transposase=/N/dc2/projects/solanumgenome/library/TE_proteins/Tpases020812DNA  # DNA transposase protein
plantprot=/N/dc2/projects/solanumgenome/library/TE_proteins/alluniRefprexp070416  # Plant protein database
fileindex=sapp  # genome index name
ThreadN=16  # number of cpus for parallel processing 
ProcStep=1234  # steps being processed. S1: Mite search; S2: LTR search; S3: Repeat-Modeler; S4: RepeatMasker
####################################################################################################



cd $WORKDIR

### MITE Hunter: take ~300 CPU hours for a genome (~1Gb) and might need large memories in step 4
if [[ $ProcStep = *"1"* ]]; 
then
	if [ ! -d mite_annotation ]; then mkdir mite_annotation; fi
	cd mite_annotation
	MITE_Hunter_manager.pl -i ../$seqfile -g $fileindex -c $ThreadN -S 5678
	cd ..
fi


### use LTR_retriever to classify intact LTR elements (finished within 24 CPU hours)
if [[ $ProcStep = *"2"* ]];
then
	if [ ! -d ltr_annotation ]; then mkdir ltr_annotation; fi
	cd ltr_annotation
	gt suffixerator -db $seqfile -indexname $fileindex -tis -suf -lcp -des -ssp -dna

	gt ltrharvest -index $fileindex -out $fileindex.out -outinner $fileindex.outinner \
	-gff3 $fileindex.gff -minlenltr 100 -maxlenltr 3000 -mintsd 4 -maxtsd 20 > $fileindex.result

	LTR_retriever -genome ../$seqfile -inharvest $fileindex.result -threads $ThreadN
	cd ..
fi


### Collecting unmasked sequences from previous steps and run RepeatModeler
if [[ $ProcStep = *"3"* ]];
then
	if [ ! -d other_repeats ]; then mkdir other_repeats; fi
	cp ltr_annotation/${seqfile}.mod.LTRlib.fa LTR.lib
	cp mite_annotation/${fileindex}_Step8_singlet.fa MITE.lib	
	RepeatMasker -lib MITE.lib -pa $ThreadN -dir other_repeats ltr_annotation/${seqfile}.mod.masked
	
	cd other_repeats
	rmaskedpart.pl ${seqfile}.mod.masked.masked 50  > umseqfile
	perl -e '$/="\n>"; while (<>) { s/>//g; my ($id, $seq) = split /\n/; print ">$_" if length $seq; }' \
	< umseqfile > umseq.fasta
	
	BuildDatabase -name umseqdb -engine ncbi umseq.fasta
	RepeatModeler -pa $ThreadN -engine ncbi -database umseqdb >& umseq.out
	
	cp RM_*/consensi.fa.classified ../
	cd ..
fi



### final processing (classify & double-check repeats) and prepare library
if [[ $ProcStep = *"4"* ]];
then
	repeatmodeler_parse.pl --fastafile consensi.fa.classified --unknowns repeatmodeler_unknowns.fasta  \
	--identities repeatmodeler_identities.fasta 

	makeblastdb -in $transposase -dbtype prot
	blastx -query repeatmodeler_unknowns.fasta -db $transposase -evalue 1e-10 -num_descriptions 10 \
	-num_threads $ThreadN -out modelerunknown_blast_results.txt

	transposon_blast_parse.pl --blastx modelerunknown_blast_results.txt --modelerunknown \
	repeatmodeler_unknowns.fasta

	mv unknown_elements.txt ModelerUnknown.lib
	cat identified_elements.txt repeatmodeler_identities.fasta > ModelerID.lib
	cat MITE.lib LTR.lib > allMITE_LTR.lib
	cat allMITE_LTR.lib ModelerID.lib > KnownRepeats.lib

	makeblastdb -in $plantprot -dbtype prot
	
	blastx -query ModelerUnknown.lib -db $plantprot -evalue 1e-10 -num_descriptions 10 \
	-num_threads $ThreadN -out ModelerUnknown.lib_blast_results.txt
	ProtExcluder.pl ModelerUnknown.lib_blast_results.txt ModelerUnknown.lib

	blastx -query KnownRepeats.lib -db $plantprot -evalue 1e-10 -num_descriptions 10 \
	-num_threads $ThreadN -out KnownRepeats.lib_blast_results.txt	
	ProtExcluder.pl KnownRepeats.lib_blast_results.txt KnownRepeats.lib

	cat KnownRepeats.libnoProtFinal ModelerUnknown.libnoProtFinal > allRepeats.lib
fi


