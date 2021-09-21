import qiime2
import os
import numpy as np
from biom.table import Table
import glob
import tempfile
from Bio import SeqIO

def filter_fastq_by_size(fastq_file_path: str, filter_size: int) -> list:
    out_records = list()
    with open(fastq_file_path) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if len(record) >= filter_size:
                record.letter_annotations = {}
                record.seq = record.seq[0:16000]
                out_records.append(record)

    return out_records


def _concatenate_fastq_(barcodes_directory: str, metadata_file: str = None) -> str:
    output_sequence_directory = os.path.join(barcodes_directory, "output_sequences")
    bar_dirs = glob.glob(os.path.join(barcodes_directory, "barcode*"))
    combined_file_path = os.path.join(output_sequence_directory, "combined_all.fasta")
    with open(combined_file_path, 'w') as combined_file:
        combined_seqs = list()
        for bar in bar_dirs:
            barcode_name = os.path.basename(bar)
            seq_file_name = os.path.join(output_sequence_directory, barcode_name + ".fasta")
            print(seq_file_name)

            with open(seq_file_name, 'w') as outfile:

                out_seqs = list()
                fastq_files = glob.glob(os.path.join(bar, "*.fastq"))
                for fastq in fastq_files:

                    out_seqs += filter_fastq_by_size(fastq, 16000)

                combined_seqs += out_seqs

                SeqIO.write(out_seqs, outfile, "fasta")
        SeqIO.write(combined_seqs, combined_file, "fasta")

    return output_sequence_directory

 def _generate_feature_table_(sequences_directory: str, path_to_combined_fasta:str) -> Table:
    sample_files = glob.glob(os.path.join(sequences_directory, "barcode*"))

    sample_names = [(os.path.basename(i)).replace(".fasta","") for i in sample_files]
    sequence_data = dict()

    with open(path_to_combined_fasta, "r") as combined_file:

        for record in SeqIO.parse(combined_file, "fasta"):
            sequence_data.update({record.id: 0})





def _minion_generate_manifest_file_(barcode_directory: str, manifest_file_dir: str) -> None:
    HEADER = ['sample-id', 'absolute-filepath']
    sample_directory = _concatenate_fastq_(barcodes_directory=barcode_directory)

    sample_paths = glob.glob(os.path.join(sample_directory, '*.fastq'))
    sample_names = [(os.path.basename(x)).replace(".fastq", "").strip() for x in sample_paths]

    manifest_file = os.path.join(manifest_file_dir, 'MANIFEST')

    with open(manifest_file, 'w') as man_file:
        man_file.write('\t'.join(HEADER) + '\n')
        for out in zip(sample_names, sample_paths):
            man_file.write('\t'.join(out) + '\n')

def minion_sequence_import(barcode_directory: str):

    pass




if __name__ == '__main__':
    barcode_directory = '/home/chase/DissertationProjects/Qiime2SummerCamp/PanViral/CulexMitoGenome/sequence_reads'
    output_sequence_directory = '/home/chase/DissertationProjects/Qiime2SummerCamp/PanViral/CulexMitoGenome/sequence_reads/output_sequences/'
    #_concatenate_fastq_(barcode_directory)
    _minion_generate_manifest_file_(barcode_directory=barcode_directory, manifest_file_dir=barcode_directory)
