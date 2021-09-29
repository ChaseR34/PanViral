import biom
import qiime2
import os
import numpy as np
import pandas as pd
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


def _get_sample_names_(barcode_names: list, metadata_file:str) -> list:
    metadata = pd.read_csv(metadata_file, sep='\t')
    sample_names = metadata["sampleid"]
    barcode_id = metadata["Barcode_ID"]


    bar_name_tmp = ['BC' + x.replace('barcode', '') for x in barcode_names]

    sample_out = list()

    sample_out += [sample_names[barcode_id.str.contains(bar)].astype(str).values.tolist()[0] for bar in bar_name_tmp]

    return sample_out



def _concatenate_fastq_(barcodes_directory: str, metadata_file: str = None) -> str:
    output_sequence_directory: str = os.path.join(barcodes_directory, "output_sequences")
    bar_dirs = glob.glob(os.path.join(barcodes_directory, "barcode*"))
    bar_names = [os.path.basename(x) for x in bar_dirs]
    if metadata_file:
        bar_tmp = [os.path.basename(x) for x in bar_dirs]
        sample_names = _get_sample_names_(barcode_names=bar_tmp, metadata_file=metadata_file)
    else:
        sample_names = [os.path.basename(x) for x in bar_dirs]
    combined_file_path = os.path.join(output_sequence_directory, "combined_all.fq")
    with open(combined_file_path, 'w') as combined_file:
        combined_seqs = list()
        for bar_name, sample_name in zip(bar_names, sample_names):
            seq_file_name = os.path.join(output_sequence_directory, sample_name + ".fastq")

            print(seq_file_name)

            with open(seq_file_name, 'w') as outfile:

                out_seqs = list()


                fastq_files =glob.glob(os.path.join(barcodes_directory, bar_name, "*.fastq"))

                for fastq in fastq_files:
                    out_seqs += filter_fastq_by_size(fastq, 16000)

                combined_seqs += out_seqs

                SeqIO.write(out_seqs, outfile, "fastq")
        SeqIO.write(combined_seqs, combined_file, "fastq")

    return output_sequence_directory


def _generate_feature_table_(sequences_directory: str, path_to_combined_fasta: str) -> Table:
    sample_files = glob.glob(os.path.join(sequences_directory, "*.fasta"))
    sample_names = [(os.path.basename(i)).replace(".fasta", "") for i in sample_files]
    sequence_data = dict()

    with open(path_to_combined_fasta, "r") as combined_file:

        for record in SeqIO.parse(combined_file, "fasta"):
            sequence_data.update({record.id: 0})

    sequence_names = [i for i in sequence_data.keys()]

    table_data = []

    for indx, sample_file in enumerate(sample_files):
        sd_tmp_dict = sequence_data.copy()
        with open(sample_file, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                sd_tmp_dict[record.id] = 1

            sd_tmp_values = [i for i in sd_tmp_dict.values()]
            table_data.append(sd_tmp_values)

    seq_names_array = np.array(sequence_names)
    samp_names_array = np.array(sample_names)
    td_array = np.transpose(np.array(table_data))



    table_out = Table(data=td_array, observation_ids=seq_names_array, sample_ids=samp_names_array, table_id="minion_frequency_table")

    return table_out


def _minion_generate_manifest_file_(barcode_directory: str, manifest_file_dir: str, metadata_file:str = None) -> None:
    HEADER = ['sample-id', 'absolute-filepath']
    sample_directory = _concatenate_fastq_(barcodes_directory=barcode_directory, metadata_file=metadata_file)

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

    import biom
    barcode_directory = '/home/chase/DissertationProjects/Qiime2SummerCamp/PanViral/CulexMitoGenome/sequence_reads/testing/'
    output_sequence_directory = '/home/chase/DissertationProjects/Qiime2SummerCamp/PanViral/CulexMitoGenome/sequence_reads/testing/output_sequences/'
    metadata_file = '/home/chase/DissertationProjects/Qiime2SummerCamp/PanViral/CulexMitoGenome/Qiime_MetaData_MIP01_09162021_ZAB.tsv'

    # _concatenate_fastq_(barcode_directory)
    _minion_generate_manifest_file_(barcode_directory=barcode_directory, manifest_file_dir=barcode_directory, metadata_file=metadata_file)

    combined_fasta = os.path.join(output_sequence_directory, 'combined_all.fas')
    t_out = _generate_feature_table_(output_sequence_directory, combined_fasta)

    with biom.util.biom_open(os.path.join(output_sequence_directory, "minion_table.biom"), 'w') as ff:
        t_out.to_hdf5(ff, "minion frequency table")


