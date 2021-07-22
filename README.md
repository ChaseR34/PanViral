# 

# PanViral

## Purpose

Using the Qiime2 framework to develop a pipeline to classify pan-viral samples against a database of known viruses. The database will be generated by the user. 



## Pipeline Steps Overview

1. Create Custom Classifier
   1. Following steps found on Qiime2 website: https://docs.qiime2.org/2021.4/tutorials/feature-classifier/
      1. import reference sequences
      2. Extract reference reads
      3. Train Classifier
      4. Test Classifier
2. Classifying Samples (Modified Procedure from Qiime2 Website)
   1. Preparing Sequences
      1. Trim sequences
      2. Import sequences
      3. Check quality of reads
   2. Clustering Sequences (Dada2)
   3. Classify Clusters
   4. Generate Visualizations

Below is a tutorial classifying flaviviruses using pan-flavi primers. 



# PanFlavi Classification Tutorial



Note: This tutorial expects that you have Qiime2 already installed. If you do not, please follow the directions on their website: 

https://docs.qiime2.org/2021.4/install/



## Clone git repository

```
git clone https://github.com/ChaseR34/PanViral

cd PanViral/PanFlaviTutorial
```

The reference sequences to used are based on the flaviviruses described in Moureau et al. 2015 (https://doi.org/10.1371/journal.pone.0117849). The sequences are found in the directory **example_viral_db** in the file ***moureau_2015_ref_sequences.fasta***

The example samples were produced from mosquito samples using the pan-flavi primer set described in Vina-Rodriguez et al. 2017 (https://doi.org/10.1155/2017/4248756). The primer name and sequence found in the table below:

| Primer Name | Sequence                              | Location in Genome (AF196835) |
| ----------- | ------------------------------------- | ----------------------------- |
| PFlav-fAAR  | TACAACATGATGGGAAAG**A**GAGAGAA**R**AA | 9040 -- 9068                  |
| PFlav-rKR   | GTGTCCCA**K**CC**R**GC**T**GTGTCATC   | 9305 -- 9283                  |



## Part 1: Create Custom Classifier

1. ### Import Reference Sequences

   1. ```bash
      qiime tools import \
        --type 'FeatureData[Sequence]' \
        --input-path example_viral_db/moureau_2015_ref_sequences.fasta \
        --output-path classifier/panflavi_ref_sequences.qza
      
      qiime tools import \
        --type 'FeatureData[Taxonomy]' \
        --input-format HeaderlessTSVTaxonomyFormat \
        --input-path example_viral_db/moureau_2015_taxonomy.tsv \
        --output-path classifier/panflavi_ref_taxonomy.qza
      ```

2. ### Extract Reference Reads

   1. ```bash
      qiime feature-classifier extract-reads \
        --i-sequences classifier/panflavi_ref_sequences.qza \
        --p-f-primer TACAACATGATGGGAAAGAGAGAGAARAA \
        --p-r-primer GTGTCCCAKCCRGCTGTGTCATC \
        --p-min-length 200 \
        --p-max-length 300 \
        --o-reads classifier/panflavi_ref_seqs_redueced.qza
      ```

      

3. ### Train Classifier

   1. ```bash
      qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads classifier/panflavi_ref_seqs_redueced.qza \
        --i-reference-taxonomy classifier/panflavi_ref_taxonomy.qza \
        --o-classifier classifier/panflavi_classifier.qza
      ```





## Part 2: Classifying Samples

### objective:

We will be using the classifier generated in part1 to classify the 5 example samples found in the folder **example_sequences**



## Preparing Sequences

1. ####  Trim Sequences using bbduk

   bbduk will trim off the primers by matching the supplied primer sequences to the first 30bp of the sequencing reads. The primer sequences are supplied to it via a reference file. The file can be found in the **primers** directory and the file name is ***primers.fasta***.

   The trimReads.sh bash script calls bbduk for you. The script takes three arguments, <path-to-sequences-dir>, <path-to-primer-file>, and <path-to-output-dir>

   ```bash
   
   mkdir example_sequences_trimmed
   
   bash ./scripts/trimReads.sh linked_sequences primers/primers.fasta example_sequences_trimmed
   
   #check trimmed_sequences
   ls trimmed_sequences
   ```

2. ### Importing Sequences into Qiime2

   1. #### Create manifest file

   This file is used by qiime2 to import the sequences into a qiime2 specific .qza file. Documentation on this import can be found at: [link]( https://docs.qiime2.org/2020.6/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq)
   scroll down to the section titled *“Fastq manifest” formats*

   The make_manifest_file.sh bash script will be used to generate your manifest:

   ```bash
   bash ./scripts/make_manifest_file.sh -d trimmed_sequences -o panflavi_manifest
   
   #check the generated manifest file
   
   cat panflavi_manifest | column -t
   ```

   

3. #### 	qiime2 import command

   ```bash
   module load qiime2
    
   qiime tools import \
     --type 'SampleData[PairedEndSequencesWithQuality]' \
     --input-path panflavi_manifest \
     --output-path paired_end_demux_panflavi.qza \
     --input-format PairedEndFastqManifestPhred33V2
   ```

   

   1. #### Notes on inputs: 

      - --type: tells qiime2 what type of sequence data is being imported, e.g., single or paired end

      - --input-path: path to manifest file

      - --output-path: name of output file

      - --input-format tells qiime2 what format the sequencing data is. See link above for more information

      

4. ### Visualization of sequence quality

   This step produces plots of the average quality per base for the forward and reverse reads. This information is used to trim off the bases of low quality within the reads.

   To produce the qiime2 visualization file, .qzv, use the following command:

   ```bash
    qiime demux summarize \
     --i-data paired-end-demux-panflavi.qza \
     --o-visualization visual_output/demux_panflavi.qzv
   ```

   To visualize the .qzv file you have to visit https://view.qiime2.org/

   

5. ## Clustering sequences using DADA2

   1. ​	This section is where we cluster the reads into 100% OTUs using the algorithm implemented by DADA2

      Using the demux_AA.qzv we will decide where to trim and truncate our forward and reverse reads. 

      ```
      qiime dada2 denoise-paired \
           --i-demultiplexed-seqs paired_end_demux_panflavi.qza \
           --p-trim-left-f 0 \
           --p-trim-left-r 0 \
           --p-trunc-len-f 245 \
           --p-trunc-len-r 100 \
           --o-representative-sequences rep_seqs_dada2_panflavi.qza \
           --o-table table_dada2_panflavi.qza \
           --o-denoising-stats stats_dada2_panflavi.qza
      
      ```

      #### Notes on inputs:

   - --p-trim-left-f Trims all of the bases upto the base given for the forward reads. e.g., --p-trim-left-f 25 will trim base position 1-25 off all the reads.

   *  --p-trim-left-r same as above but for reverse reads.
   *  --p-trunc-len-f Truncates all of the bases after the base given for the forward reads. e.g., --p-trunc-len-f 200 will remove bases 201-end of read.
   *  --p-trunc-len-r same as above but for reverse reads.

6. ### Generate Metadata file using Keemei  

   1. See instructions on the qiime2 website: https://docs.qiime2.org/2021.4/tutorials/metadata/
   2. save file as ***panflavi_meta.tsv***

   

7. ### Visualizing feature table 

   1. ```bash
      qiime feature-table summarize \
        --i-table table_dada2_panflavi.qza \
        --o-visualization visual_output/table_panflavi.qzv \
        --m-sample-metadata-file panflavi_meta.tsv
      
      qiime feature-table tabulate-seqs \
        --i-data rrep_seqs_dada2_panflavi.qza \
        --o-visualization visual_output/rep_seqs_panflavi.qzv
        
      ```



## Classification

1. Using a naïve Bayes classifier, the features are classified against all publicly available virus genomes.

   1. ```bash
      qiime feature-classifier classify-sklearn   \
      --i-classifier classifier/panflavi_classifier.qza  \
       --i-reads rep_seqs_dada2_panflavi.qza   \
       --o-classification taxonomy_panflavi.qza
      ```

      

2. ## Visualizing Results

   1. #### Create Taxa Bar plot

   2. ```bash
       qiime taxa barplot \
       --i-table table_dada2_panflavi.qza \ 
       --i-taxonomy taxonomy_panflavi.qza \
       --m-metadata-file panflavi_meta.tsv \
       --o-visualization visual_ouput/taxa_bar_plots_panflavi.qzv
      
      ```

3. ## Create Combined Table 

   This section we will produce a table that contains the feature id, sequences, and feature counts per sample.

   1. we have to reorient the feature table so the columns and rows match the taxonomy_AA.qza and paired-rep-seqs-dada2.qza. This is accomplished using the 
     feature-table transpose command:

   ```bash
   qiime feature-table transpose 
     --i-table table_dada2_panflavi.qza \
     --o-transposed-feature-table transposed_table_dada2_panflavi.qza \
   
   ```

   Now we generate the table:

   ```bash
   qiime metadata tabulate 
     --m-input-file taxonomy_panflavi.qza \
     --m-input-file rep_seqs_dada2_panflavi.qza \
     --m-input-file transposed_table_dada2_panflavi.qza \
     --o-visualization visual_output/feat-tax-rep.qzv
   ```

   ​	

   
