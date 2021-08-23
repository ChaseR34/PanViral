import qiime2
import os
import tempfile

def _return_names_(file_path_names: list) -> list:

    names_raw = [os.path.basename(x) for x in file_path_names]

    names_out = list()

    for name in names_raw:
        if '-' in name:
            name_split = name.split('-')
            name_1_len = len(name_split[0])

            if name_1_len < 8:
                names_out.append(name_split[1])
            else:
                names_out.append(name_split[0])
        else:
            names_out = name.split('_')[0]
    return names_out

def _generate_manifest_file_(sample_dir_path: str, manifest_file_dir: str) -> tuple:
    """Generates a manifest file in a temporary directory to be used in sequence reads upload"""

    HEADER = ['sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath']

    sequence_paths = os.listdir(sample_dir_path)
    forward_paths = [os.path.join(sample_dir_path, x) for x in sequence_paths if "_R1" in x]
    reverse_paths_tmp = [os.path.join(sample_dir_path, x) for x in sequence_paths if "_R2" in x]

    names = _return_names_(forward_paths)

    reverse_paths = []

    for name in names:
        for rev in reverse_paths_tmp:
            if name in rev:
                reverse_paths.append(rev)

    with open(os.path.join(manifest_file_dir, "manifest"), "w") as manifest_file:

        manifest_file.write('\t'.join(HEADER) + '\n')

        for out in zip(names, forward_paths, reverse_paths):

            manifest_file.write('\t'.join(out) + '\n')

    return tuple([os.path.join(manifest_file_dir, "manifest"), names])


def _write_metadata_template_(sample_names: list, output_dir: str ) -> None:
    HEADER = "sample-id"

    if output_dir and not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    elif not output_dir:
        output_dir = os.getcwd()

    with open(os.path.join(output_dir, "metadata_template.tsv"), "w") as metadata_file:
        metadata_file.write(HEADER + '\n')
        metadata_file.write('\n'.join(sample_names))


def prep_sequence_reads(ctx, sequences_directory, metadata_template_dir=None, primer_f=None, primer_r=None):
    results = []
    sequences_directory = os.path.abspath(sequences_directory)

    # importing external plugins to be used later
    cut_adapt = ctx.get_action('cutadapt', 'trim_paired')
    create_table_viz = ctx.get_action('demux', 'summarize')

    #importing sequences

    temp_dir = tempfile.TemporaryDirectory()

    manifest_file_path, names = _generate_manifest_file_(sequences_directory, temp_dir.name)

    read_seqs = qiime2.Artifact.import_data(type='SampleData[PairedEndSequencesWithQuality]',
                                            view=manifest_file_path,
                                            view_type='PairedEndFastqManifestPhred33V2')
    #write metadata template
    _write_metadata_template_(sample_names=names, output_dir=metadata_template_dir)


    # using plugins to trim reads and create reads visualization

    if primer_f and primer_r:
        trimmed_reads = cut_adapt(demultiplexed_sequences=read_seqs, front_f=[primer_f], front_r=[primer_r])
        table_viz = create_table_viz(data=trimmed_reads.trimmed_sequences)

        results += trimmed_reads
    else:
        table_viz = create_table_viz(data=read_seqs.trimmed_sequences)
        results += [read_seqs]

    results += table_viz

    return tuple(results)
