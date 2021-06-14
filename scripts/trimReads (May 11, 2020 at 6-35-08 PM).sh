#!/bin/env bash

 for i in /scratch/clr96/PanFlavi/KomarSet/*R1*.fastq.gz; do
       # i=/scratch/clr96/PanFlavi/KomarSet/UNKNOWN-PanFlavxxA022xx43xxPE18xx225-xx-xx-xxx-xxxx-143-CH_S83_L001_R1_001.fastq.gz;
       input1=$i;
       input2=${input1/R1_001.fastq.gz/R2_001.fastq.gz};
       output1=$(basename ${input1/R1_001.fastq.gz/R1_output.fastq});
       output2=$(basename ${input2/R2_001.fastq.gz/R2_output.fastq});
       mismatch1=$(basename ${input1/R2_001.fastq.gz/R1_mismatch.fastq});
       mismatch2=$(basename ${input2/R2_001.fastq.gz/R2_mismatch.fastq});
       stats=$(basename ${input1/R1_001.fastq.gz/R1_stats.txt});
       echo;
       echo "input1: $input1";
       echo "input2: $input2";
       echo "output1: $output1";
       echo "output2: $output2";
       echo "stats: $stats";
       echo;

    /packages/bbmap/38.61b/bbduk.sh bbduk.sh -Xmx6g -da \
                                            in=$input1 \
                                            in2=$input2 \
                                            out=$output1 \
                                            out2=$output2 \
                                            outm=$mismatch1 \
                                            outm2=$mismatch2 \
                                            ref=../primers.fastq \
                                            stats=$stats \
                                            statscolumns=5 \
                                            ottm=t \
                                            restrictleft=30 \
                                            ordered=t k=7 \
                                            minlen=80 \
                                            ktrim=l \
                                            edist=3 \
                                            edist2=1 \
                                            ftm=5 \
                                            mink=5 \
                                            copyundefined \
                                            fixjunk;
 done
