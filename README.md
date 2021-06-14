# PanViral
# PanFlavi Classification Tutorial

## Setup Project Folder

```
cd /scratch/<nauid> # this can be any folder you want
mkdir qiime2_tutorial
cd qiime2_tutorial
```

Adding folders for sequences

```
mkdir linked_sequences trimmed_sequences scripts

```

Copy over useful bash scripts to scripts folder

```
cp /scratch/clr96/qiime_tutuorial/scripts/make_manifest_file.sh scripts
cp /scratch/clr96/qiime_tutuorial/scripts/trimReads.sh scripts

#check scripts are in folder

ls scripts

```

## Preparing Sequences for qiime2

### Create symbolic-link of the original sequencess into the *linked_sequences* folder

```
ln -s /scratch/clr96/qiime2_tutorial/sequences/*fastq* linked_sequences

#check sequences are in the folder

ls linked_sequences

```
The reason you link the sequences instead of copying or using the direct path to where the sequences are stored is safety and saves space. If you accidentally delete or attempt to modify the linked files, the original files will not be changed.


### Trim Sequences using bbduk

bbduk will trim off the primers by matching the supplied primer sequences to the first 30bp of the sequencing reads. The primer sequences are supplied to it via a reference file. In our case it will be the following file:
```
ln -s /scratch/clr96/qiime2_tutorial/primers.fasta .
```

The trimReads.sh bash script calls bbduk for you. The script takes three arguments, <path-to-sequences-dir>, <path-to-primer-file>, and <path-to-output-dir>

```
bash ./scripts/trimReads.sh linked_sequences primers.fasta trimmed_sequences

#check trimmed_sequences
ls trimmed_sequences
```


## Importing Sequences in qiime2

### Create manifest file

This file is used by qiime2 to import the seqeunces into a qiime2 specific .qza file. Documentation on this import can be found at: [link]( https://docs.qiime2.org/2020.6/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq)
scroll down to the section titled *“Fastq manifest” formats*

The make_manifest_file.sh bash script will be used to generate your manifest:

```
bash ./scripts/make_manifest_file.sh -d trimmed_sequences -o AA_manifest

#check the generated manifest file

cat AA_manifest | column -t

```


### qiime2 import command

```
module load qiime2
 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path AA_manifest \
  --output-path paired-end-demux-AA.qza \
  --input-format PairedEndFastqManifestPhred33V2

```

Notes on inputs: 

--type: tells qiime2 what type of sequence data is being imported, e.g., single or paired end

--input-path: path to manifest file

--output-path: name of output file

--input-format tells qiime2 what format the sequencing data is. See link above for more information

Now you should have a file titled paired-end-demux-AA.qza.

###  Visualization of sequence quality

This step produces plots of the average quality per base for the forward and reverse reads. This information is used to trim off the bases of low quality within the reads.

To produce the qiime2 visualization file, .qzv, use the following command:

```
 qiime demux summarize \
  --i-data paired-end-demux-AA.qza \
  --o-visualization demux_AA.qzv
```



To download the file from monsoon to your computer using a terminal:

1. open a new local terminaluse the following command:
2. run the following command

```
# scp <nauid>@monsoon.hpc.nau.edu:/scratch/<nauid>/qiime2_tutorial/demux_AA.qzv <path-to-local-directory>

scp clr96@monsoon.hpc.nau.edu:/scratch/clr96/qiime2_tutorial/demux_AA.qzv /Users/Curly/Desktop
```
The demux_AA.qzv can be visualized at the following website [link](https://view.qiime2.org/)

import or drag the file into the open vew.qiime2.org webpage.

## Clustering sequences using DADA2

This section is where we cluster the reads into 100% OTUs using the algorithm implemented by DADA2

Using the demux_AA.qzv we will decide where to trim and truncate our forward and reverse reads. 

```
qiime dada2 denoise-paired \
 --i-demultiplexed-seqs paired-end-demux-AA.qza \
 --p-trim-left-f 0 \
 --p-trim-left-r 0 \
 --p-trunc-len-f 245 \
 --p-trunc-len-r 100 \
 --o-representative-sequences rep-seqs-dada2-AA.qza \
 --o-table table-dada2-AA.qza \
 --o-denoising-stats stats-dada2-AA.qza

```

Notes on inputs

*  --p-trim-left-f Trims all of the bases upto the base given for the forward reads. e.g., --p-trim-left-f 25 will trim base position 1-25 off all the reads.
*  --p-trim-left-r same as above but for reverse reads.
*  --p-trunc-len-f Truncates all of the bases after the base given for the forward reads. e.g., --p-trunc-len-f 200 will remove bases 201-end of read.
*  --p-trunc-len-r same as above but for reverse reads.
  

### Generate Metadata file using Keemei  

### Visualizing feature table 

```
qiime feature-table summarize \
  --i-table table-dada2-AA.qza \
  --o-visualization table-AA.qzv \
  --m-sample-metadata-file AA_meta.tsv
  
qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2-AA.qza \
  --o-visualization rep-seqs-AA.qzv
  
```

## Classification

Using a naive bayes classifier, the features are classified against all publically available virus genomes.

```

qiime feature-classifier classify-sklearn   \
--i-classifier /scratch/clr96/classifier.qza  \
 --i-reads rep-seqs-dada2-AA.qza   \
 --o-classification taxonomy-AA.qza
 
```

## Species Bar plot

```
 
 qiime taxa barplot \
 --i-table paired-table-dada2.qza \ 
 --i-taxonomy trainingClassifier/Taxonomy/taxonomy-AA.qza \
 --m-metadata-file AA_meta.tsv \
 --o-visualization taxa-bar-plots.qzv

``` 

## Create Combined Table 

This section we will produce a table that contains the feature id, sequences, and feature counts per sample.

1. we have to reorient the feature table so the columns and rows match the taxonomy_AA.qza and paired-rep-seqs-dada2.qza. This is accomplished using the 
feature-table transpose command:

```
qiime feature-table transpose 
  --i-table paired-table-dada2.qza \
  --o-transposed-feature-table tansposed_table.qza \

```

2. Now we generate the table:

```
qiime metadata tabulate 
  --m-input-file taxonomy_AA.qza \
  --m-input-file paired-rep-seqs-dada2.qza \
  --m-input-file transposed_table.qza \
  --o-visualization feat-tax-rep.qzv

```


## TODO Creating your own classifier
## TODO Diversity Statistics

