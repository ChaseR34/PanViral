#!/bin/env bash

#SBATCH --job-name=Pdenoise
#SBATCH --mail-type=All
#SBATCH --mail-user=clr96@nau.edu
#SBATCH --output=/scratch/clr96/PanFlavi/slurmout/%x_%A.out
#SBATCH --error=/scratch/clr96/PanFlavi/slurmout/%x_%A.err
#SBATCH --time=00:45:00
#SBATCH --mem=5000
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

module purge
module load qiime2/2020.11

echo "start time:  `date`"

export data=$1
bname=$(basename $data)
outname=${bname%%-*}

echo $data
echo $outname


qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $data \
  --p-trim-left-f 0 \
  --p-trunc-len-f 215 \
  --p-trim-left-r 0  \
  --p-trunc-len-r 210 \
  --o-representative-sequences ${outname}-rep-seqs-dada2.qza \
  --o-table ${outname}-table-dada2.qza \
  --o-denoising-stats ${outname}-stats-dada2.qza

echo "end time:  `date`"
