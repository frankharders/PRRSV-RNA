#!/bin/bash

##  activate the environment for this downstream analysis
eval "$(conda shell.bash hook)";
conda activate 02-POPPUNK

mkdir -p ./02_polished;
mkdir -p ./dbMAPPING;
mkdir -p ./LOGS;


MINCOV=30;

MINEDISTMAX=30;
MINEDIST=16
MINALLELFRAC=0.05;
MINQUALMAX=15;
MINSCORE=15;
PLOIDY=10;
FRACTION=0.3; # 1 is alles
MINLEN=100;
NORM=10000;
FRAC=0.6;

count0=1;
countS=$(cat samples.txt | wc -l);

REF=ref.fa;

mkdir -p $PWD/dbMAPPING/"$REF"/;

OUT=$PWD/dbMAPPING/"$REF"/;
OUT2=$PWD/dbMAPPING;


		while [ $count0 -le $countS ];do

SAMPLE=();

		SAMPLE=$(cat samples.txt | awk 'NR=='$count0 );

		LOG1=./LOGS/"$SAMPLE".reformat.log;
		LOG2=./LOGS/"$SAMPLE".clumpify.log;
		LOG3=./LOGS/"$SAMPLE".bbduk_1.log;
		LOG4=./LOGS/"$SAMPLE".bbmap.log;
		LOG5=./LOGS/"$SAMPLE".dedupe.log;
		LOG6=./LOGS/"$SAMPLE".filtersam_1.log;
		LOG7=./LOGS/"$SAMPLE".filtersam_2.log;
		LOG8=./LOGS/"$SAMPLE".bbduk_2.log;
		LOG9=./LOGS/"$SAMPLE".callvariants.log;
		LOG10=./LOGS/"$SAMPLE".pileup.log;
		LOG11=./LOGS/"$SAMPLE".applyvariants.log;
		LOG12=./LOGS/"$SAMPLE".bbmap_2.log;
		LOG13=./LOGS/"$SAMPLE".callvariants.consensus.log;

		LOG=./LOGS/"$SAMPLE".general.log;
		echo "general log for sample $SAMPLE" > "$LOG";

		echo -e "$SAMPLE";
		echo -e "used reference $REF";

		echo 'reformat';
		# 1 combine R1 & R2 into 1 file
		echo "reformat interleave R1 and R2 files" >> "$LOG";
		reformat.sh in=./RAWREADS/"$SAMPLE"_R#.fastq.gz out=./02_polished/"$SAMPLE".fq.gz ow=t -Xmx384g >> "$LOG" 2>&1;
		# 2 normalize the reads to $NORM
		echo "normalize reads to $NORM depth" >> "$LOG";
		bbnorm.sh in=./02_polished/"$SAMPLE".fq.gz out=./02_polished/"$SAMPLE".norm"$NORM".fq.gz -Xmx384g ow passes=2 target="$NORM" fixspikes=t >> "$LOG" 2>&1;
		echo 'clumpify';
		# 2 get rid of duplicates
		echo "clumpify reads" >> "$LOG";
		clumpify.sh in=./02_polished/"$SAMPLE".norm"$NORM".fq.gz out=./02_polished/"$SAMPLE"_clumped.fq.gz zl=9 dedupe s=2 passes=4 -Xmx384g ow >> "$LOG" 2>&1;
		echo 'bbduk_1';
		# 3 trim adapters from file adapters.fa
		echo "trim adapter from reads" >> "$LOG";
		bbduk.sh in=./02_polished/"$SAMPLE"_clumped.fq.gz out=./02_polished/"$SAMPLE"_trimmed.fq.gz minlen=60 ktrim=r k=21 mink=9 hdist=2 hdist2=1 ref=./adapter.fa altref=adapters maq=14 qtrim=r trimq=10 maxns=0 tbo tpe ow -Xmx384g ftm=5 >> "$LOG" 2>&1;
		echo 'bbmap'
		# 4 align reads against the reference
		echo "align reads to $reference" >> "$LOG";
		bbmap.sh ref="$PWD"/references/"$REF" in=./02_polished/"$SAMPLE"_trimmed.fq.gz outm="$OUT"/"$SAMPLE"_mapped.sam.gz nodisk local maxindel=500 -Xmx384g ow k=12 >> "$LOG" 2>&1;
		echo 'dedupe'
		# 5 get rid of duplicates based on aligning the sequences. Do'nt use if minorities must be scored, consensus only
		echo "dedupe reads, get rid of duplicate reads" >> "$LOG";
		dedupebymapping.sh in="$OUT"/"$SAMPLE"_mapped.sam.gz out="$OUT"/"$SAMPLE"_deduped.sam.gz -Xmx384g ow >> "$LOG" 2>&1;
		echo 'filtersam_1'
		# 6 filter sam file from uniq bad reads conatining uniq kmers
		echo "filter sam file for artefacs" >> "$LOG";
		filtersam.sh ref="$PWD"/references/"$REF" ow in="$OUT"/"$SAMPLE"_deduped.sam.gz out="$OUT"/"$SAMPLE"_filtered.sam.gz mbad=1 del sub=f mbv=0 -Xmx384g >> "$LOG" 2>&1;
		echo 'filtersam_2'
		# 7
		echo "filter sam file for artefacs" >> "$LOG";
		filtersam.sh ref="$PWD"/references/"$REF" ow in="$OUT"/"$SAMPLE"_filtered.sam.gz out="$OUT"/"$SAMPLE"_filtered2.sam.gz mbad=1 sub mbv=2 -Xmx384g >> "$LOG" 2>&1;
		echo 'bbduk_2'
		# 8 softclip the reads
		echo "softclip reads in sam file" >> "$LOG";
		bbduk.sh in="$OUT"/"$SAMPLE"_filtered2.sam.gz trimclip out="$OUT"/"$SAMPLE"_trimclip.sam.gz -Xmx384g ow >> "$LOG" 2>&1;
		echo 'samtools view & index'
		# 9 from sam to sorted/indexed bam using samtools
		samtools view -bShu "$OUT"/"$SAMPLE"_trimclip.sam.gz | samtools sort -m 2G -@ 3 - -o "$OUT"/"$SAMPLE"_sorted.ref.bam;
		samtools index "$OUT"/"$SAMPLE"_sorted.ref.bam;
		echo 'samtools consesnsus'
		# 10 create a consensus sequence using samtools 
		samtools consensus -f fasta --call-fract "$FRAC" --ambig -a -d 3 -@ 4 "$OUT"/"$SAMPLE"_sorted.ref.bam -o "$OUT"/"$SAMPLE".samtools.fasta;
		echo 'rename header consensus to sample'
		# 11 rename header consensus
		echo "rename header consensus to sample without an additional number" >> "$LOG";
		bbrename.sh in="$OUT"/"$SAMPLE".samtools.fasta out="$OUT"/"$SAMPLE"_genome.consensus.fa prefix="$SAMPLE" prefixonly=t;


## align reads against consensus sequence
		echo bbmap to consensus
		# 12 align reads against the reference
		echo "align reads to consensus" >> "$LOG";
		bbmap.sh ref="$OUT"/"$SAMPLE"_genome.consensus.fa in=./02_polished/"$SAMPLE"_trimmed.fq.gz outm="$OUT"/"$SAMPLE"_mapped.consensus.sam.gz nodisk local mintrimlength="$MINLEN" samplerate="$FRACTION" maxindel=500 -Xmx384g ow k=12 >> "$LOG" 2>&1;
		echo 'dedupe'
		# 13 get rid of duplicates based on aligning the sequences. Don't use if minorities must be scored, consensus only
		echo "dedupe reads, get rid of duplicate reads, consensus" >> "$LOG";
		dedupebymapping.sh in="$OUT"/"$SAMPLE"_mapped.consensus.sam.gz out="$OUT"/"$SAMPLE"_deduped.consensus.sam.gz -Xmx384g ow >> "$LOG" 2>&1;
		echo 'filtersam_1'
		# 14 filter sam file from uniq bad reads containing uniq kmers 
		echo "filter sam file for artefacs, consensus" >> "$LOG";
		filtersam.sh ref="$OUT"/"$SAMPLE"_genome.consensus.fa ow in="$OUT"/"$SAMPLE"_deduped.consensus.sam.gz out="$OUT"/"$SAMPLE"_filtered.consensus.sam.gz mbad=1 del sub=f mbv=0 -Xmx384g >> "$LOG" 2>&1;
		echo 'filtersam_2'
		# 15 filter sam file from uniq bad reads containing uniq kmers
		echo "filter sam file for artefacs, consensus" >> "$LOG";
		filtersam.sh ref="$OUT"/"$SAMPLE"_genome.consensus.fa ow in="$OUT"/"$SAMPLE"_filtered.consensus.sam.gz out="$OUT"/"$SAMPLE"_filtered2.consensus.sam.gz mbad=1 sub mbv=2 -Xmx384g >> "$LOG" 2>&1;
		echo 'bbduk_2'
		# 16 softclip the reads 
		echo "softclip reads in sam file, consensus" >> "$LOG";
		bbduk.sh in="$OUT"/"$SAMPLE"_filtered2.consensus.sam.gz trimclip out="$OUT"/"$SAMPLE"_trimclip.consensus.sam.gz -Xmx384g ow >> "$LOG" 2>&1;
		echo 'samtools view & index'
		# 17 from sam to sorted/indexed bam using samtools
		samtools view -bShu "$OUT"/"$SAMPLE"_trimclip.consensus.sam.gz | samtools sort -m 2G -@ 3 - -o "$OUT"/"$SAMPLE"_sorted.consensus.bam;
		samtools index "$OUT"/"$SAMPLE"_sorted.consensus.bam;
		echo 'callvariants consensus'
		# 18 call variants from consensus
		echo "call variants from mapping to consensus" >> "$LOG";		
		callvariants.sh in="$OUT"/"$SAMPLE"_trimclip.consensus.sam.gz ref="$OUT"/"$SAMPLE"_genome.consensus.fa out="$OUT"/"$SAMPLE"_vars.consensus.vcf -Xmx4g ow strandedcov usebias=f minstrandratio=0 maf=0.05 minreads="$MINCOV" mincov="$MINCOV" minedistmax="$MINEDISTMAX" minedist="$MINEDIST" flagnearby >> "$LOG" 2>&1;


cp "$PWD"/references/"$REF" "$OUT";



rm ./02_polished/*.fq.gz;
rm "$OUT"/*.sam.gz;









#		# 9
#		callvariants.sh in="$OUT"/"$SAMPLE"_mapped.sam.gz ref="$PWD"/references/"$REF" out="$OUT"/"$SAMPLE"_vars.vcf -Xmx4g ow usebias=f callsub=t calldel=t callins=t minreads="$MINCOV" mincov="$MINCOV" maf=0.6 minedistmax="$MINEDISTMAX" minedist="$MINEDIST" ;#flagnearby > "$LOG9" 2>&1;
#		echo 'pileup'
#		# 10
#		pileup.sh in="$OUT"/"$SAMPLE"_mapped.sam.gz basecov="$OUT"/"$SAMPLE"_basecov_border5.txt -Xmx4g ow border=1 > "$LOG10" 2>&1;
#		echo 'applyvariants'
#		# 11
#		applyvariants.sh in=$PWD/references/"$REF" out="$OUT"/"$SAMPLE"_genome.fa vcf="$OUT"/"$SAMPLE"_vars.vcf basecov="$OUT"/"$SAMPLE"_basecov_border5.txt ow noframeshifts=f mindepth="$MINCOV" > "$LOG11" 2>&1;
#		echo 'mapping original reads to "consensus sequence"'
#		# 12
#		bbrename.sh in="$OUT"/"$SAMPLE"_genome.fa out="$OUT"/"$ID"_genome.consensus.fa prefix="$ID" prefixonly=t;



#CONS1="$OUT"/"$ID"_genome.consensus.fa;
#cp $CONS1 $OUT2;

		#13
#		bbmap.sh ref="$CONS1" in=./02_polished/"$SAMPLE"_trimmed.fq.gz outm="$OUT"/"$ID"_mapped.consensus2.sam.gz nodisk local maxindel=500 -Xmx4g ow k=12;# > "$LOG12" 2>&1;
		#14
#		filtersam.sh ref="$CONS1" ow in="$OUT"/"$ID"_mapped.consensus2.sam.gz out="$OUT"/"$ID"_filtered.cons1.sam.gz mbad=1 del sub=f mbv=0 -Xmx4g ;#> "$LOG6" 2>&1;
#		echo 'filtersam_2'
		# 7
#		filtersam.sh ref="$CONS1" ow in="$OUT"/"$ID"_filtered.cons1.sam.gz out="$OUT"/"$ID"_filtered2.cons1.sam.gz mbad=1 sub mbv=2 -Xmx4g ;#> "$LOG7" 2>&1;
#		echo 'bbduk_2'
		# 8
#		bbduk.sh in="$OUT"/"$ID"_filtered2.cons1.sam.gz trimclip out="$OUT"/"$ID"_trimclip.cons1.sam.gz -Xmx1g ow ;#> "$LOG8" 2>&1;
#		echo 'callvariants'
		# 9
#		callvariants.sh in="$OUT"/"$ID"_trimclip.cons1.sam.gz ref="$CONS1" out="$OUT"/"$ID"_vars.cons1.vcf -Xmx4g ow strandedcov usebias=f minstrandratio=0 maf=0.6 minreads="$MINCOV" mincov="$MINCOV" minedistmax="$MINEDISTMAX" minedist="$MINEDIST" flagnearby ;#> "$LOG9" 2>&1;
#		echo 'pileup'
		# 10
#		pileup.sh in="$OUT"/"$ID"_trimclip.cons1.sam.gz basecov="$OUT"/"$ID"_basecov_border5.cons1.txt -Xmx4g ow border=5 ;#> "$LOG10" 2>&1;
#		echo 'applyvariants'
		# 11
#		applyvariants.sh in="$CONS1" out="$OUT"/"$ID"_genome.cons1.fa vcf="$OUT"/"$ID"_vars.cons1.vcf basecov="$OUT"/"$ID"_basecov_border5.cons1.txt ow mindepth="$MINCOV" ;#> "$LOG11" 2>&1;
#		echo 'mapping original reads to "consensus sequence"'

#		bbrename.sh in="$OUT"/"$ID"_genome.cons1.fa out="$OUT"/"$ID"_genome.cons2.fa prefix="$ID" prefixonly=t;

#CONS2="$OUT"/"$ID"_genome.cons2.fa;
#cp $CONS2 $OUT2;




		count0=$((count0+1));

	done


exit 1


#Variant-Calling Cutoffs:
#minreads=2              (minad) Ignore variants seen in fewer reads.
#maxreads=BIG            (maxad) Ignore variants seen in more reads.
#mincov=0                Ignore variants in lower-coverage locations.
#maxcov=BIG              Ignore variants in higher-coverage locations.
#minqualitymax=15        Ignore variants with lower max base quality.
#minedistmax=20          Ignore variants with lower max distance from read ends.
#minmapqmax=0            Ignore variants with lower max mapq.
#minidmax=0              Ignore variants with lower max read identity.
#minpairingrate=0.1      Ignore variants with lower pairing rate.
#minstrandratio=0.1      Ignore variants with lower plus/minus strand ratio.
#minquality=12.0         Ignore variants with lower average base quality.
#minedist=10.0           Ignore variants with lower average distance from ends.
#minavgmapq=0.0          Ignore variants with lower average mapq.
#minallelefraction=0.1   Ignore variants with lower allele fraction.  This
#                        should be adjusted for high ploidies.
#minid=0                 Ignore variants with lower average read identity.
#minscore=20.0           Ignore variants with lower Phred-scaled score.
#clearfilters            Clear all filters.  Filter flags placed after
#                        the clearfilters flag will still be applied.






