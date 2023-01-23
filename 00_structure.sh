#!/bin/bash

#create several directories
while getopts "i:o:" opt; do
  case $opt in
    w)
      echo "-w was triggered! $OPTARG"
      WORKDIR="`echo $(readlink -m $OPTARG)`"
      echo $WORKDIR
      ;;
     a)
      echo "-a was triggered! $OPTARG"
      RAW_FASTQ="`echo $(readlink -m $OPTARG)`"
      echo $RAW_FASTQ
      ;;
     b)
      echo "-b was triggered! $OPTARG"
      RAWSTATS="`echo $(readlink -m $OPTARG)`"
      echo $RAWSTATS
      ;;
	 c)
      echo "-c was triggered! $OPTARG"
      POLISHED="`echo $(readlink -m $OPTARG)`"
      echo $POLISHED
      ;;
	 d)
      echo "-d was triggered! $OPTARG"
      TRIMMEDSTATS="`echo $(readlink -m $OPTARG)`"
      echo $TRIMMEDSTATS
      ;;
	 e)
      echo "-e was triggered! $OPTARG"
      SHOVILL="`echo $(readlink -m $OPTARG)`"
      echo $SHOVILL
      ;;
	 f)
      echo "-f was triggered! $OPTARG"
      QUAST="`echo $(readlink -m $OPTARG)`"
      echo $QUAST
      ;;
	 g)
      echo "-g was triggered! $OPTARG"
      QUASTparse="`echo $(readlink -m $OPTARG)`"
      echo $QUASTparse
      ;;
     h)
      echo "-h was triggered! $OPTARG"
      MLST="`echo $(readlink -m $OPTARG)`"
      echo $MLST
      ;;
     i)
      echo "-i was triggered! $OPTARG"
      MLSTparse="`echo $(readlink -m $OPTARG)`"
      echo $MLSTparse
      ;;	  
	 j)
      echo "-j was triggered! $OPTARG"
      ABRICATE="`echo $(readlink -m $OPTARG)`"
      echo $ABRICATE
      ;;
     k)
      echo "-k was triggered! $OPTARG"
      ABRICATEparse="`echo $(readlink -m $OPTARG)`"
      echo $ABRICATEparse
      ;;	  
	 l)
      echo "-l was triggered! $OPTARG"
      TMP="`echo $(readlink -m $OPTARG)`"
      echo $TMP
      ;;	  
	 n)
      echo "-n was triggered! $OPTARG"
      LOG="`echo $(readlink -m $OPTARG)`"
      echo $LOG
      ;;
     r)
      echo "-r was triggered! $OPTARG"
      REPORTING="`echo $(readlink -m $OPTARG)`"
      echo $REPORTING
      ;;	  
    \?)
      echo " Let's start some analysis" >&2
      ;;

  esac
done

if [ "x" == "x$WORKDIR" ] || [ "x" == "x$RAW_FASTQ" ] || [ "x" == "x$RAWSTATS" ]  || [ "x" == "x$POLISHED" ]  || [ "x" == "x$TRIMMEDSTATS" ] ; then
    echo "-w $WORKDIR -a $RAW_FASTQ -b $RAWSTATS -c $POLISHED -d $TRIMMEDSTATS -e $SHOVILL -f QUAST -g $QUASTparse -h $MLST -i $MLSTparse -j $ABRICATE -k $ABRICATEparse -l $TMP -n $LOG -r $REPORTING"
    echo "-w, -a, -b, -c, -d, -e, -f, -g, -h, -i, -j, -k, -l, -n, -r  [options] are required"
    exit 1
fi


### standard project directories

mkdir -p $RAWSTATS $POLISHED $TRIMMEDSTATS $SHOVILL $TMP $LOG $REPORTING ;

while read SAMPLE;do


mkdir -p $REPORTING/$SAMPLE;

done < samples.txt


exit 1
