#!/bin/sh
#$ -S /bin/bash
#$ -t 1-48
#$ -cwd
#qsub -l medium -l s_vmem=10G -l mem_req=10G -pe def_slot 10 j1_fastp.sh

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

DATADIR=${PROJECTDIR}/Data/rawdata/
ANALYSISDIR=${PROJECTDIR}/Analysis/

cd $DATADIR
FILENAME=`ls -1 *_1.fq.gz | awk -v line=$SGE_TASK_ID '{if (NR == line) print $0}'`
SAMPLE=$(echo $FILENAME | rev | cut -c 9- | rev | uniq)

FASTQ1=${DATADIR}/${SAMPLE}_1.fq.gz
FASTQ2=${DATADIR}/${SAMPLE}_2.fq.gz

cd $ANALYSISDIR

fastp -i ${FASTQ1} -I ${FASTQ2} -3\
 -o ${ANALYSISDIR}/1_post_fastp/${SAMPLE}_1.fq.gz -O ${ANALYSISDIR}/1_post_fastp/${SAMPLE}_2.fq.gz\
 -h ${ANALYSISDIR}/1_post_fastp/report_${SAMPLE}.html -q 30 -u 30 -n 10 -l 20 -w 10 -f 1 -F 1

echo ending at
date
