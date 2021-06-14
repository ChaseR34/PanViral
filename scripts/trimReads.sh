#!/bin/env bash
module load bbmap 

  seq_folder=$1
  primer_file=$2
  output=$3

  [ ! -d "${seq_folder}" ] && echo "This ${seq_folder} does not " && exit 1
  [ ! -f "${primer_file}" ] && echo "Primer file: ${seq_folder} does not " && exit 1

 for i in "${seq_folder}"/*R1*.fastq.gz; do

       input1="$i";
       input2="${input1/R1_001.fastq.gz/R2_001.fastq.gz}";
       output1="$(basename ${input1/R1_001.fastq.gz/R1_output.fastq})";
       output2="$(basename ${input2/R2_001.fastq.gz/R2_output.fastq})";
       mismatch1="$(basename ${input1/R2_001.fastq.gz/R1_mismatch.fastq})";
       mismatch2="$(basename ${input2/R2_001.fastq.gz/R2_mismatch.fastq})";
       stats="$(basename ${input1/R1_001.fastq.gz/R1_stats.txt})";
       echo;
       echo "input1: $input1";
       echo "input2: $input2";
       echo "output1: $output1";
       echo "output2: $output2";
       echo "stats: $stats";
       echo;

   if  [[ ! -f "${input1}" ]]; then

   echo "${seq_folder}"
   else

    /scratch/clr96/packages/src/bbmap/bbduk.sh -Xmx6g -da \
                                            in="$input1" \
                                            in2="$input2" \
                                            out="${output}/$output1" \
                                            out2="${output}/$output2" \
                                            outm="${output}/$mismatch1" \
                                            outm2="${output}/$mismatch2" \
                                            ref="${primer_file}" \
                                            stats="${output}/$stats" \
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
  fi

 done
