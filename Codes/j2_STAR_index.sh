#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#qsub -l medium -l s_vmem=100G -l mem_req=100G -pe def_slot 16 j2_STAR_index.sh

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

REFDIR=${PROJECTDIR}/Data/Ensembl_v100
ANALYSISDIR=${PROJECTDIR}/Analysis/2_post_STAR
GENOMEDIR=${PROJECTDIR}/Data/Index_STAR

STAR --runMode genomeGenerate --genomeDir ${GENOMEDIR} --genomeFastaFiles ${REFDIR}/fasta/Mus_musculus.GRCm38.dna.toplevel.fa --sjdbGTFfile ${REFDIR}/GTF/Mus_musculus.GRCm38.100.gtf --sjdbOverhang 100 --runThreadN 16 --limitGenomeGenerateRAM 100000000000

echo ending at
date
