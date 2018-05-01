#!/bin/bash
set -e
set -u
set -o pipefail


function usage(){
printf "\nDescription: this script is to predict gene models from assembled genome using \
MAKER2 pipeline. You must ensure the follwing program has been correctly installed and added \
to the PATH: Bowtie2, Tophat2, Trinity, Maker2, Busco, Augustus, and GeneMark, Blast+, Hmmer. 

$(basename "$0") [-h] [-OPTIONS] -d <working dir> -g <genome file> -f <forward rna seq> \
-r <reverse rna seq> -p <homologous protein sequences> -t <TE proteins> -R <repeat library> \
-I <species name> -B <BUSCO searching database> -M <maker_exe.ctl>

where:
   -h  show this help text
   -d  the working directory for gene annotation
   -g  the input genome fasta file 
   -f  forward reads of paired-end RNA-seq
   -r  reverse reads of paired-end RNA-seq
   -p  the protein sequences from related species
   -t  all TE protein sequences 
   -R  the repeat library
   -I  species name (e.g. solanum_lycopersicum; no space in species name)
   -B  the searching database for BUSCO (e.g. embryophyta_odb9)
   -M  the maker_exe.ctl file indicating the PATH to required programs 
   -n  number of cpus for parallel processing (default=16)
   -S  the steps to processed (default=1234): 1: train Augustus; 2: train GeneMark; 3: train SNAP; 4. run MAKER\n\n"
}


######### Default parameters #########
WORKDIR=''
GENOME=''
RNAseq1=''
RNAseq2=''
BUSCODB=''
MAKEREXE=''
Protseq=''
TEprot=''
Repeatlib=''
SpeciesID=''
ThreadN=16
ProcStep=1234


######### Parse input #########
while getopts 'h:d:g:f:r:p:t:R:I:B:M:n:S:' option; do
  case "$option" in
    h) usage
       exit
       ;;
    d) WORKDIR=${OPTARG}
       ;;
    g) GENOME=${OPTARG}
       ;;
    f) RNAseq1=${OPTARG}
       ;;
    r) RNAseq2=${OPTARG}
       ;;
    p) Protseq=${OPTARG}
       ;;
    t) TEprot=${OPTARG}
       ;;
    R) Repeatlib=${OPTARG}
       ;;
    r) rRNAsDB=${OPTARG}
       ;;
    I) SpeciesID=${OPTARG}
       ;;
    B) BUSCODB=${OPTARG}
       ;;
    M) MAKEREXE=${OPTARG}
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



if [ -z "${WORKDIR}" ] || [ -z "${GENOME}" ] || [ -z "${RNAseq1}" ] || [ -z "${RNAseq2}" ] \
|| [ -z "${Protseq}" ] || [ -z "${TEprot}" ] || [ -z "${Repeatlib}" ] || [ -z "${SpeciesID}" ] \
|| [ -z "${BUSCODB}" ] || [ -z "${MAKEREXE}" ] ; 
then
    usage
    exit
fi


cd $WORKDIR 
Gindex=$(basename $GENOME .fa)



######### A function to edit CONFIG file ######### 
function set_config(){
CONFIGFILE=$1
TARGET_KEY=$2
REPLACEMENT_VALUE=$3
sed -i "s|^\($TARGET_KEY\s*=\).*\$|\1$REPLACEMENT_VALUE|" $CONFIGFILE
}


############ Train Augustus ############
if [[ $ProcStep = *"1"* ]]; 
then
	run_BUSCO.py -i $GENOME -o my_busco -l $BUSCODB -m genome -c $ThreadN -f
	cd run_my_busco/augustus_output
	
	sed -i "s/$(basename *_parameters.pbl _parameters.pbl)/${SpeciesID}/g" *
	basename *_exon_probs.pbl _exon_probs.pbl | xargs -I{} mv {}_exon_probs.pbl "${SpeciesID}"_exon_probs.pbl
	basename *_intron_probs.pbl _intron_probs.pbl | xargs -I{} mv {}_intron_probs.pbl "${SpeciesID}"_intron_probs.pbl
	basename *_igenic_probs.pbl _igenic_probs.pbl | xargs -I{} mv {}_igenic_probs.pbl "${SpeciesID}"_igenic_probs.pbl
	basename *_metapars.cfg _metapars.cfg | xargs -I{} mv {}_metapars.cfg "${SpeciesID}"_metapars.cfg
	basename *_metapars.utr.cfg _metapars.utr.cfg | xargs -I{} mv {}_metapars.utr.cfg "${SpeciesID}"_metapars.utr.cfg
	basename *_parameters.cfg _parameters.cfg | xargs -I{} mv {}_parameters.cfg "${SpeciesID}"_parameters.cfg
	basename *_weightmatrix.txt _weightmatrix.txt | xargs -I{} mv {}_weightmatrix.txt "${SpeciesID}"_weightmatrix.txt
			
	cd ..
	cp retraining_parameters $SpeciesID
	cp $SpeciesID "${AUGUSTUS_CONFIG_PATH}/species/${SpeciesID}"
	cd ..
fi


############ Train GeneMark ############
if [[ $ProcStep = *"2"* ]]; 
then
	if [ ! -d run_my_genemark ]; then mkdir run_my_genemark; fi
	cd run_my_genemark
	
	bowtie2-build $GENOME $Gindex
	tophat2 -p $ThreadN $Gindex $RNAseq1 $RNAseq2

	bet_to_gff.pl --bed tophat_out/junctions.bed --gff introns.gff --label Tophat2
	gmes_petap.pl --sequence $GENOME --ET introns.gff --cores $ThreadN 
	
	cp output/gmhmm.mod ../gmhmm.mod
	cd ..
fi


############ 1st round of training SNAP in MAKER2 ############
if [[ $ProcStep = *"3"* ]]; 
then
	if [ ! -d run_my_snap ]; then mkdir run_my_snap; fi
	cd run_my_snap
	
	Trinity --seqType fq --max_memory 64G --CPU $ThreadN --full_cleanup \
	--output $Gindex --left $RNAseq1 --right $RNAseq2
	
	genome_stat.py -i ${GENOME} 
	
	maker -OPTS
	maker -BOPTS
	cp ${MAKEREXE} "$( pwd )"/maker_exe.ctl
	
	set_config maker_opts.ctl "genome" "$( pwd )"/large_scaffolds.fa
	set_config maker_opts.ctl "est" "$( pwd )"/${Gindex}".Trinity.fasta"
	set_config maker_opts.ctl "protein" ${Protseq}
	set_config maker_opts.ctl "rmlib" ${Repeatlib}
	set_config maker_opts.ctl "repeat_protein" ${TEprot}
	set_config maker_opts.ctl "est2genome" "1"
	set_config maker_opts.ctl "protein2genome" "1"
	set_config maker_opts.ctl "keep_preds" "1"
	set_config maker_opts.ctl "single_exon" "1"
	set_config maker_opts.ctl "model_org" ""	
	
	mpiexec -n $ThreadN maker -base snap1
	
	cd snap1.maker.output
 	gff3_merge -d snap1_master_datastore_index.log
 	maker2zff snap1.all.gff
	fathom -categorize 1000 genome.ann genome.dna
	fathom -export 1000 -plus uni.ann uni.dna
	forge export.ann export.dna
	hmm-assembler.pl snap . > ../snap1.hmm	
	
	cd ..
	set_config maker_opts.ctl "snaphmm" "$( pwd )"/snap1.hmm
	set_config maker_opts.ctl "est2genome" "0"
	set_config maker_opts.ctl "protein2genome" "0" 		
 	
 	mpiexec -n $ThreadN maker -base snap2
 	
 	cd snap2.maker.output
 	gff3_merge -d snap2_master_datastore_index.log
 	maker2zff snap2.all.gff
	fathom -categorize 1000 genome.ann genome.dna
	fathom -export 1000 -plus uni.ann uni.dna
	forge export.ann export.dna
	hmm-assembler.pl snap . > ../snap2.hmm	
		
	cd ..
	set_config maker_opts.ctl "snaphmm" "$( pwd )"/snap2.hmm	
 	mpiexec -n $ThreadN maker -base snap3
 	
 	cd snap3.maker.output
 	gff3_merge -d snap3_master_datastore_index.log
 	maker2zff snap3.all.gff
	fathom -categorize 1000 genome.ann genome.dna
	fathom -export 1000 -plus uni.ann uni.dna
	forge export.ann export.dna
	hmm-assembler.pl snap . > ../snap3.hmm	
	
	cd ..
	cp snap3.hmm ../snap3.hmm
	cp maker_opts.ctl ../maker_opts.ctl 
	cp maker_bopts.ctl ../maker_bopts.ctl 
	cp maker_exe.ctl ../maker_exe.ctl 
	cd ..
fi


############ Predict gene models by MAKER2 ############	
if [[ $ProcStep = *"4"* ]]; 
then
	set_config maker_opts.ctl "genome" ${GENOME}
	set_config maker_opts.ctl "augustus_species" ${SpeciesID}
	set_config maker_opts.ctl "gmhmm" "$( pwd )"/gmhmm.mod
	set_config maker_opts.ctl "snaphmm" "$( pwd )"/snap3.hmm
	set_config maker_opts.ctl "pred_stats" "1"
	set_config maker_opts.ctl "min_protein" "30"
	set_config maker_opts.ctl "alt_splice" "1"	
 	
 	mpiexec -n $ThreadN maker 
 	
 	cd ${Gindex}".maker.output"
 	gff3_merge -d ${Gindex}"_master_datastore_index.log"
 	fasta_merge -d ${Gindex}"_master_datastore_index.log"
	cd ..
fi


