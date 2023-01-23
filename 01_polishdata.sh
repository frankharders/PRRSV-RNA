#!/bin/bash

##  activate the environment for this downstream analysis
eval "$(conda shell.bash hook)";
conda activate 02-POPPUNK


cnt=$(cat samples.txt | wc -l);

##  list of directories used with this pipeline

OUTPUT=$PWD/01_polished/;
LOG=$PWD/LOGS/;
TEMP=$PWD/TEMP/;

## create directories for general use of the script

mkdir -p "$LOG" "$TEMP" "$OUTPUT";

##  remove all files from previous analysis


#rmdir "01_polished";
#rmdir "dbMAPPING"_"$ToDay";
#rmdir "$LOG";
#rmdir "$TEMP";





## used variables read trimming

KTRIM=r;
TRIMQ=20;
QTRIM='rl';
k=29;

while read SAMPLE; do 

echo $SAMPLE;


SEQDIR=./RAWREADS/;

OUTDIR=./01_polished/;

echo "$SAMPLE";

## LOG files for debugging
LOG0="$LOG"/"$SAMPLE".filterbytile.log;
LOG1="$LOG"/"$SAMPLE".adapterFind.log;
LOG3="$LOG"/"$SAMPLE".adapterclip.log;
LOG4="$LOG"/"$SAMPLE".qualitytrim.log;
LOG5="$LOG"/"$SAMPLE".pairedinfo.log;
LOG6="$LOG"/"$SAMPLE".insertSize.log;
LOG7="$LOG"/"$SAMPLE"."$REF".readmapping.log;

#R1=$PWD/RAWREADS/R1.1.fastq.gz;
#R2=$PWD/RAWREADS/R2.1.fastq.gz;

R1=$PWD/RAWREADS/"$SAMPLE"_R1.fastq.gz;
R2=$PWD/RAWREADS/"$SAMPLE"_R2.fastq.gz;


	

## intermediate files stored in a TEMp directory
## files will be deleted when re-starting the script 
FILTERED1="$TEMP"/"$SAMPLE"_R1.filterbytile.fq.gz;
FILTERED2="$TEMP"/"$SAMPLE"_R2.filterbytile.fq.gz;
ADAPTER="$TEMP"/"$SAMPLE"_adapters.fa;
ADAPTERout1="$TEMP"/"$SAMPLE"_R1.adapter.fq.gz;
ADAPTERout2="$TEMP"/"$SAMPLE"_R2.adapter.fq.gz;
NEXTERA=/home/harde004/pipelines/resources/nextera.fa.gz;
	
OUTPUT1="$OUTDIR"/"$SAMPLE"_R1.QTR.adapter.nextera.k"$k".fq.gz;
OUTPUT2="$OUTDIR"/"$SAMPLE"_R2.QTR.adapter.nextera.k"$k".fq.gz;

MERGED="$TEMP"/"$SAMPLE"-merged.fq.gz;
UNMERGED="$TEMP"/"$SAMPLE"-unmerged.fq.gz;
INSERTout="$LOG"/"$SAMPLE"-insertSize.tab;


echo -e "R1=$R1";
echo -e "R2=$R2";



		# 1
		##  filter blurry parts from flowcell
#		filterbytile.sh -Xmx24g in1="$R1" in2="$R2" out1="$FILTERED1" out2="$FILTERED2" ow=t > "$LOG0" 2>&1;

		# 2
		##  Discover adapter sequence for this library based on read overlap.
		##  You can examine the adapters output file afterward if desired;
		##  If there were too few short-insert pairs this step will fail (and you can just use the default Illumina adapters).
#		bbmerge.sh in1="$FILTERED1" in2="$FILTERED2" outa="$ADAPTER" ow reads=1m > "$LOG1" 2>&1;

		# 3
		##  Perform adapter-trimming on the reads.
		##  Also do quality trimming and filtering.
		##  If desired, also do primer-trimming here by adding, e.g., 'ftl=20' to to trim the leftmost 20 bases.
		##  If the prior adapter-detection step failed, use "ref=adapters"
#     	bbduk.sh -Xmx24g in1="$FILTERED1" in2="$FILTERED2" out1="$ADAPTERout1" out2="$ADAPTERout2" ktrim="$KTRIM" ref="$NEXTERA" k="$k" mink=21 ow=t > "$LOG3" 2>&1;

     	bbduk.sh -Xmx24g in1="$R1" in2="$R2" out1="$ADAPTERout1" out2="$ADAPTERout2" ktrim="$KTRIM" ref="$NEXTERA" k="$k" mink=21 ow=t > "$LOG3" 2>&1;
		bbduk.sh -Xmx48g in1="$ADAPTERout1" in2="$ADAPTERout2" out1="$OUTPUT1" out2="$OUTPUT2" qtrim="$QTRIM" trimq="$TRIMQ" ow=t > "$LOG4" 2>&1;

		# 4
		##  verify read pairs
		reformat.sh verifypaired=t in1="$OUTPUT1" in2="$OUTPUT2" ow=t > "$LOG5" 2>&1;

		# 5
		##  verify insertSize of NGS libs
#		bbmerge-auto.sh in1="$OUTPUT1" in2="$OUTPUT2" out="$MERGED" outu="$UNMERGED" outinsert="$INSERTout" iterations=10 k="$k" ecct qtrim2=r trimq="$TRIMQ" strict minlength=150 ow=t > "$LOG6" 2>&1;

		# 6
		##  map the polished reads against the clustered database file
#		bbmap.sh -Xmx48g ref="$reference" build=1 in1="$OUTPUT1" in2="$OUTPUT2" outm="$TEMP"/"$SAMPLE"."$REF".sam untrim=t showprogress=0 ambig=best local=t sam=1.3 fast=t ow=t > "$LOG7" 2>&1 ;
		
#		samtools view -@ "$NODES" -S "$TEMP"/"$SAMPLE"."$REF".sam -b -o "$TEMP"/"$SAMPLE"."$REF".bam;
#		samtools sort -@ "$NODES"  -o "$dbMAPPING"/"$SAMPLE"."$REF"_sorted.bam "$TEMP"/"$SAMPLE"."$REF".bam;
#		samtools index "$dbMAPPING"/"$SAMPLE"."$REF"_sorted.bam "$dbMAPPING"/"$SAMPLE"."$REF"_sorted.bai;

	let cnt--;	
	echo -e "$cnt samples to go!";
	echo "NEXT";





#else 

#echo "reference check if the file is missing from the original location";


#fi 












	
		done < samples.txt

echo "read polishing is done for all samples";
echo -e "output files for downstream processing can be found in the directory $PROCESSED";

		
exit 1

