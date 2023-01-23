#!/bin/bash

cnt=$(cat samples.txt | wc -l);

POLISHED=./01_polished/;

# (temp) variables

DEPTH=100;
MINlen=800;
MINcov=100;

SHOVILL=./02_assembly/;

mkdir -p $SHOVILL;

TARGET=60;


export PATH=$PATH:/home/harde004/miniconda3/pkgs/gsl-2.7-he838d99_0/lib/;



while read SAMPLE; do 

R1=$POLISHED/$SAMPLE'_R1.QTR.adapter.nextera.k29.fq.gz';
R2=$POLISHED/$SAMPLE'_R2.QTR.adapter.nextera.k29.fq.gz';
R1n=$POLISHED/$SAMPLE'_R1.QTR.adapter.nextera.norm.k29.fq.gz';
R2n=$POLISHED/$SAMPLE'_R2.QTR.adapter.nextera.norm.k29.fq.gz';

bbnorm.sh -Xmx12g in1=$R1 in2=$R2 out1=$R1n out2=$R2n target=$TARGET min=5 ow=t fixspikes=t;#> $LOG6 2>&1;




	echo $SAMPLE;
	
OUTPUTdir=$SHOVILL/$SAMPLE/;

shovill --outdir $OUTPUTdir --depth $DEPTH --minlen $MINlen --mincov $MINcov --keepfiles --namefmt $SAMPLE'_contig%05d' --force --R1 $R1n --R2 $R2n; 

FASTAin=$OUTPUTdir/contigs.fa;
FASTAout=$OUTPUTdir/$SAMPLE'_contigs.fa';

#reformat.sh in=$FASTAin out=$FASTAout fastawrap=10000000;
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $FASTAin > $FASTAout;
	
	let cnt--;	
	echo -e "$cnt samples to go!";
	echo "NEXT";



				
done < samples.txt

echo "genome assembly is done using $ASSEMBLER for all samples";
echo -e "output files for downstream processing can be found in the directory $SHOVILL";


exit 1

#####
#SYNOPSIS
#  De novo assembly pipeline for Illumina paired reads
#USAGE
#  shovill [options] --outdir DIR --R1 R1.fq.gz --R2 R2.fq.gz
#GENERAL
#  --help          This help
#  --version       Print version and exit
#  --check         Check dependencies are installed
#INPUT
#  --R1 XXX        Read 1 FASTQ (default: '')
#  --R2 XXX        Read 2 FASTQ (default: '')
#  --depth N       Sub-sample --R1/--R2 to this depth. Disable with --depth 0 (default: 100)
#  --gsize XXX     Estimated genome size eg. 3.2M <blank=AUTODETECT> (default: '')
#OUTPUT
#  --outdir XXX    Output folder (default: '')
# --force         Force overwite of existing output folder (default: OFF)
#  --minlen N      Minimum contig length <0=AUTO> (default: 0)
#  --mincov n.nn   Minimum contig coverage <0=AUTO> (default: 2)
#  --namefmt XXX   Format of contig FASTA IDs in 'printf' style (default: 'contig%05d')
#  --keepfiles     Keep intermediate files (default: OFF)
#RESOURCES
#  --tmpdir XXX    Fast temporary directory (default: '')
#  --cpus N        Number of CPUs to use (0=ALL) (default: 16)
#  --ram n.nn      Try to keep RAM usage below this many GB (default: 32)
#ASSEMBLER
#  --assembler XXX Assembler: velvet megahit skesa spades (default: 'spades')
#  --opts XXX      Extra assembler options in quotes eg. spades: "--untrusted-contigs locus.fna" ... (default: '')
#  --kmers XXX     K-mers to use <blank=AUTO> (default: '')
#MODULES
#  --trim          Enable adaptor trimming (default: OFF)
#  --noreadcorr    Disable read error correction (default: OFF)
#  --nostitch      Disable read stitching (default: OFF)
#  --nocorr        Disable post-assembly correction (default: OFF)
#HOMEPAGE
#  https://github.com/tseemann/shovill - Torsten Seemann

