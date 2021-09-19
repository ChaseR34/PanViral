import qiime2
import os
import glob
import tempfile
from Bio import SeqIO

def filter_fastq_by_size(fastq_file_path: str, filter_size: int) -> list:
    out_records = list()
    with open(fastq_file_path) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if len(record) >= filter_size:
                out_records.append(record)

    return out_records


def _concatenate_fastq_(barcodes_directory: str, metadata_file: str = None) -> str:
    output_sequence_directory = os.path.join(barcodes_directory, "output_sequences")
    bar_dirs = glob.glob(os.path.join(barcodes_directory, "barcode*"))

    for bar in bar_dirs:
        barcode_name = os.path.basename(bar)
        seq_file_name = os.path.join(output_sequence_directory, barcode_name + ".fastq")
        print(seq_file_name)
        with open(seq_file_name, 'w') as outfile:
            out_seqs = list()
            fastq_files = glob.glob(os.path.join(bar, "*.fastq"))
            for fastq in fastq_files:

                out_seqs += filter_fastq_by_size(fastq, 10000)

            SeqIO.write(out_seqs, outfile, "fastq")

    return output_sequence_directory


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






if __name__ == '__main__':
    barcode_directory = '/home/chase/DissertationProjects/Qiime2SummerCamp/PanViral/CulexMitoGenome/sequence_reads'
    output_sequence_directory = '/home/chase/DissertationProjects/Qiime2SummerCamp/PanViral/CulexMitoGenome/sequence_reads/output_sequences/'
    _minion_generate_manifest_file_(barcode_directory=barcode_directory, manifest_file_dir=barcode_directory)
