#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#qsub -l medium -l s_vmem=30G -l mem_req=30G -pe def_slot 16 j3_featureCounts.sh

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

DATADIR=${PROJECTDIR}/Data/Ensembl_v100/GTF
ANALYSISDIR=${PROJECTDIR}/Analysis

LIST=""
for SAMPLE in `ls -1 ${ANALYSISDIR}/2_post_STAR/*Aligned.sortedByCoord.out.bam`; do
        LIST="$LIST ${SAMPLE}"
done

RESDIR=${ANALYSISDIR}/3_featureCounts
cd ${RESDIR}

featureCounts -Q 20 -p -B -T 16 -t exon -g gene_id -a ${DATADIR}/Mus_musculus.GRCm38.100.gtf -o Counts.txt $LIST

echo ending at
date
