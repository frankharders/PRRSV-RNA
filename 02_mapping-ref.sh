#!/bin/bash

##  activate the environment for this downstream analysis
eval "$(conda shell.bash hook)";
conda activate 02-POPPUNK

## no minimal coverage is selected because this part of the pipeline is only used for testing the sample prep what extraction will be used (Shield or Trizol)

cnt=$(cat samples.txt | wc -l);

mkdir -p ./02_mapping-ref/;
mkdir -p ./LOGS/;

LOG=./LOGS/;

MAPPING=./02_mapping-ref/;


reference=./references/ref.fa;
ref=$(basename $reference | cut -f1 -d'.');


count0=1;
countS=$(cat samples.txt | wc -l);


while [ $count0 -le $countS ];do

		SAMPLE=$(cat samples.txt | awk 'NR=='$count0 );

#R1=./01_polished/$SAMPLE'_R1.QTR.adapter.nextera.k29.fq.gz';
#R2=./01_polished/$SAMPLE'_R2.QTR.adapter.nextera.k29.fq.gz';




LOG1=$LOG/$SAMPLE.mapping.$ref.log;
LOG2=$LOG/$SAMPLE.dedupedbymapping.$ref.log;
LOG3=$LOG/$SAMPLE.rename.fasta.$ref.log;
LOG4=$LOG/$SAMPLE.mapping.consensus.log;


OUT1=$MAPPING/$SAMPLE.$ref.sam.gz;
OUT1b=$MAPPING/$SAMPLE.$ref.deduped.sam.gz;

OUTM=$MAPPING/$SAMPLE.consensus.sam.gz;
OUTd1=$MAPPING/$SAMPLE.consensus.trimclip.sam.gz;


OUT2=$MAPPING/$SAMPLE.$ref.sorted.bam;
OUTc2=$MAPPING/$SAMPLE.consensus.sorted.bam;

CONSout=$MAPPING/$SAMPLE.consensus.fa;
RENAME=$MAPPING/$SAMPLE.consensus.SL.fa;

VCF=$MAPPING/$SAMPLE.consensus.vcf;
VCFfilter=$MAPPING/$SAMPLE.consensus.filter.vcf;

##  align reads to selected reference, standard genotype 3c refernece sequence. MT362711.1 Hepeviridae isolate HEVgt3c_NL_serum_906010, partial genome

#		reformat.sh in=./RAWREADS/"$SAMPLE"_R#.fastq.gz out=./RAWREADS/"$SAMPLE".fq.gz ow=t ;#> "$LOG1" 2>&1;
		echo 'clumpify';
		# 2
#		clumpify.sh in=./RAWREADS/"$SAMPLE".fq.gz out=./01_polished/"$SAMPLE"_clumped.fq.gz zl=9 dedupe s=2 passes=4 -Xmx31g ;#> "$LOG2" 2>&1;
		echo 'bbduk_1';
		# 3
#		bbduk.sh in=./01_polished/"$SAMPLE"_clumped.fq.gz out=./01_polished/"$SAMPLE"_trimmed.fq.gz minlen=60 ktrim=r k=21 mink=9 hdist=2 hdist2=1 ref=./adapter.fa altref=adapters maq=14 qtrim=r trimq=10 maxns=0 tbo tpe ow -Xmx1g ftm=5 ;#> "$LOG3" 2>&1;
		echo 'bbmap'
		# 4
#		bbmap.sh ref="$reference" in=./01_polished/"$SAMPLE"_trimmed.fq.gz outm="$MAPPING"/"$SAMPLE"_mapped.sam.gz nodisk local maxindel=500 -Xmx4g ow k=12 ;#> "$LOG4" 2>&1;




##  reduce data by dedupe reads for consensus calling, this can't be used for quasi species determination
#     2
echo -e "dedupe reads by mapping to reduce data for draft consensus sequence";
#dedupebymapping.sh -Xmx384G in="$MAPPING"/"$SAMPLE"_mapped.sam.gz out="$OUT1b"  ow ;#> "$LOG2" 2>&1;

##  use samtools to create a consensus sequence
#     3
echo -e "samtools create consensus"; 
samtools view -bShu "$OUT1b" | samtools sort -m 2G -@ 3 - -o "$OUT2"; 
samtools index "$OUT2";
#     4
samtools consensus -@ 48 -a -q -f fasta --show-ins yes --show-del yes "$OUT2"  -o "$CONSout";

echo -e "rename consensus file header and make in single line fasta file";
#     5
#bbrename.sh in="$CONSout" out="$RENAME" fastawrap=1000000 prefixonly="$SAMPLE" ignorejunk ow > "$LOG3" 2>&1;


## align reads to consensus sequence for SNP calling and selecting read for assembly (optional)

echo -e "align reads to the consensus sequence and create vcf file for manual inspection";
#     6
#bbmap.sh -Xmx4g ref="$RENAME" in1="$R1" in2="$R2" outm="$OUTM" nodisk local maxindel=500 k=12 ow > "$LOG4" 2>&1;

#     7
#samtools view -bShu "$OUTM" | samtools sort -m 2G -@ 3 - -o "$OUTc2"; 
#samtools index "$OUTc2";




#freebayes -f "$RENAME" --haplotype-length 0 --min-alternate-count 1 --min-alternate-fraction 0 --pooled-continuous --report-monomorphic "$OUTc2" > "$VCF";
#vcffilter  -f "SRP > 20" -f "SAP > 20" -f "EPP > 20" -f "QUAL > 20" -f "DP > 20" "$VCF" > "$VCFfilter";

MINCOV=1;

MINEDISTMAX=30;
MINEDIST=16
MINALLELFRAC=0.05;
MINQUALMAX=15;
MINSCORE=15;

BASECOV="$MAPPIN"/"$SAMPLE".consensus.basecov.txt;

## consensus logs

LOG8=./LOGS/"$SAMPLE".trimclip.log;
LOG9=./LOGS/"$SAMPLE".callvariants.log;
LOG10=./LOGS/"$SAMPLE".pileup.log;



#     8
#		bbduk.sh -Xmx2g in="$OUTM" trimclip out="$OUTd1" ow > "$LOG8" 2>&1;
#		echo 'callvariants'
		# 9
#		callvariants.sh -Xmx4g in="$OUTd1" ref="$RENAME" out="$VCF" strandedcov usebias=f minstrandratio=0 maf=0.05 minreads="$MINCOV" mincov="$MINCOV" minedistmax="$MINEDISTMAX" minedist="$MINEDIST" flagnearby ow > "$LOG9" 2>&1;
#		echo 'pileup'
		# 10
#		pileup.sh -Xmx4g in="$OUTd1" basecov="$BASECOV" border=5 ow > "$LOG10" 2>&1;
#		echo 'applyvariants'
		# 11
#		applyvariants.sh in=$RENAME out="$OUT"/"$SAMPLE"_genome.fa vcf="$OUT"/"$SAMPLE"_vars.vcf basecov="$OUT"/"$SAMPLE"_basecov_border5.txt mindepth="$MINCOV" ow;#> "$LOG11" 2>&1;
#		echo 'mapping original reads to "consensus sequence"'
		# 12
#		bbrename.sh in="$OUT"/"$SAMPLE"_genome.fa out="$OUT"/"$SAMPLE"_genome.consensus.fa prefix="$SAMPLE";


count0=$((count0+1));

	let cnt--;	
	echo -e "$cnt samples to go!";
	echo "NEXT";

done


exit 1


#Usage: samtools consensus [options] <in.bam>

#Options:
#  -r, --region REG      Limit query to REG. Requires an index
#  -f, --format FMT      Output in format FASTA, FASTQ or PILEUP [FASTA]
#  -l, --line-len INT    Wrap FASTA/Q at line length INT [70]
#  -o, --output FILE     Output consensus to FILE
#  -m, --mode STR        Switch consensus mode to "simple"/"bayesian" [bayesian]
#  -a                    Output all bases (start/end of reference)
#  --rf, --incl-flags STR|INT
#                        Only include reads with any flag bit set [0]
#  --ff, --excl-flags STR|INT
#                        Exclude reads with any flag bit set
#                        [UNMAP,SECONDARY,QCFAIL,DUP]
#  --min-MQ INT          Exclude reads with mapping quality below INT [0]
#  --show-del yes/no     Whether to show deletion as "*" [no]
#  --show-ins yes/no     Whether to show insertions [yes]
#  -A, --ambig           Enable IUPAC ambiguity codes [off]
#
#For simple consensus mode:
#  -q, --(no-)use-qual   Use quality values in calculation [off]
#  -c, --call-fract INT  At least INT portion of bases must agree [0.75]
#  -d, --min-depth INT   Minimum depth of INT [1]
#  -H, --het-fract INT   Minimum fraction of 2nd-most to most common base [0.5]
#
#For default "Bayesian" consensus mode:
#  -C, --cutoff C        Consensus cutoff quality C [10]
#      --(no-)adj-qual   Modify quality with local minima [on]
#      --(no-)use-MQ     Use mapping quality in calculation [on]
#      --(no-)adj-MQ     Modify mapping quality by local NM [on]
#      --NM-halo INT     Size of window for NM count in --adj-MQ [50]
#      --scale-MQ FLOAT  Scale mapping quality by FLOAT [1.00]
#      --low-MQ  INT     Cap minimum mapping quality [1]
#      --high-MQ INT     Cap maximum mapping quality [60]
#      --P-het FLOAT     Probability of heterozygous site[1.0e-04]
#
#Global options:
#      --input-fmt-option OPT[=VAL]
#               Specify a single input file format option in the form
#               of OPTION or OPTION=VALUE
#  -@, --threads INT
#               Number of additional threads to use [0]
#      --verbosity INT
#               Set level of verbosity





