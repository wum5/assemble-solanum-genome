#!/bin/bash
# a bash script to build species-specific repeats library


function usage(){
printf "\nThis script is to prepare species-specific library before performing gene \
annotation. You can control the steps need to be processed by editing the flag '-S'. \
You must also ensure the following programs have been correctly installed and added \
into PATH: MITE-hunter, LTR-harvest, LTR-retriever, Repeat-Modeler, Repeat-Masker, \
Deconseq, CD-HIT, Hmmer, Blast and ProtExcluder. The files 'Tpases020812DNA' and \
'alluniRefprexp070416' can be downloaded from the website here: \
http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Repeat_Library_Construction-Advanced.

$(basename "$0") [-h] 
where:
    -h  show this help text
    -d  PATH to the working directory for assembly
    -g  the input genome fasta file
    -T  the DNA transposase protein file (Tpases020812DNA)
    -P  the plant protein database (alluniRefprexp070416)
    -n  number of cpus for parallel processing (default=16)
    -S  the steps to processed (default=1234): 1: Mite search; 2: LTR search; 3: Repeat-Modeler; 4: Final Processing\n\n"
}


######### Default parameters #########  
WORKDIR=''
seqfile=''
fileindex=''
transposase=''
plantprot=''
ThreadN=16 
ProcStep=1234 


######### Parse input #########
while getopts 'h:d:g:i:T:P:n:S' option; do
  case "$option" in
    h) usage
       exit
       ;;
    d) WORKDIR=${OPTARG}
       ;;
    g) seqfile=${OPTARG}
       ;;
    i) fileindex=${OPTARG}
       ;;
    T) transposase=${OPTARG}
       ;;
    P) plantprot=${OPTARG}
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


if [ -z "${WORKDIR}" ] || [ -z "${seqfile}" ] || [ -z "${fileindex}" ] || [ -z "${transposase}" ] || [ -z "${plantprot}" ] ; then
    usage
    exit
fi



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


