#!/bin/bash

##  loop over all used conda envs for this pipeline so all verions can be recorded for the used software and databases

## central file for reporting all PROGVERSIONS of software ans databases that will be used for the several steps in the pipeline

mkdir -p "$PWD"/REPORTING/;

PROGVERSIONS="$PWD"/REPORTING/versions.txt;

rm -f "$PROGVERSIONS";		
		

##  00_fastqc.sh
##  fastqc was used

	eval "$(conda shell.bash hook)";
	conda activate NGS-QC;

		ENV=amr-QC;
		PROG=$(fastqc -v);
				
			echo -e "script: 01_fastqc.sh" >> "$PROGVERSIONS";
			echo -e "conda environment=$ENV" >> "$PROGVERSIONS";			
			echo -e "programm\tversion" >> "$PROGVERSIONS";
			echo -e "fastQC\t$PROG" >> "$PROGVERSIONS";
			echo -e "\n" >> "$PROGVERSIONS";

##  PRRSV-consensus.sh
##  BBMap suite is used

	eval "$(conda shell.bash hook)";
	conda activate 02-POPPUNK;

		ENV=();
		ENV=02-POPPUNK;
		PROG=();
		PROG=$(conda list | grep '^bbmap' | sed 's/  */;/g' | cut -f2 -d';');
		PROG1=$(conda list | grep '^samtools' | sed 's/  */;/g' | cut -f2 -d';');

			echo -e "script: PRRSV-consensus.sh" >> "$PROGVERSIONS";
			echo -e "conda environment=$ENV" >> "$PROGVERSIONS";			
			echo -e "programm\tversion" >> "$PROGVERSIONS";			
			echo -e "bbmap-Suite\t$PROG" >> "$PROGVERSIONS";
			echo -e "samtools\t$PROG1" >> "$PROGVERSIONS";
			echo -e "\n" >> "$PROGVERSIONS";









exit 1
