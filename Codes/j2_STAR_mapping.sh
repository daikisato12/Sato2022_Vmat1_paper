#!/bin/sh
#$ -S /bin/bash
#$ -t 1-48
#$ -cwd
#qsub -l medium -l s_vmem=30G -l mem_req=30G -pe def_slot 16 j2_STAR_mapping.sh

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

DATADIR=${PROJECTDIR}/Data/rawdata/
ANALYSISDIR=${PROJECTDIR}/Analysis/
GENOMEDIR=${PROJECTDIR}/Data/Index_STAR

cd $DATADIR
FILENAME=`ls -1 *_1.fq.gz | awk -v line=$SGE_TASK_ID '{if (NR == line) print $0}'`
SAMPLE=$(echo $FILENAME | rev | cut -c 9- | rev | uniq)

RESDIR=${ANALYSISDIR}/2_post_STAR
cd ${RESDIR}

FASTQ1=${ANALYSISDIR}/1_post_fastp/${SAMPLE}_1.fq.gz
FASTQ2=${ANALYSISDIR}/1_post_fastp/${SAMPLE}_2.fq.gz
STAR --genomeDir ${GENOMEDIR} --readFilesIn $FASTQ1 $FASTQ2 --readFilesCommand zcat --runThreadN 16 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $SAMPLE

echo ending at
date
