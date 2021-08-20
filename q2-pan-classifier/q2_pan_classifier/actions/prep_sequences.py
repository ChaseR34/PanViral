import qiime2
import glob
import os
import tempfile



def generate_manifest_file(sample_dir_path: str) -> str:
    """Generates a manifest file in a temporary directory to be used in sequence reads upload"""

    HEADER = ['sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath']
    temp_dir = tempfile.TemporaryDirectory()

    forward_paths = [os.path.abspath(x) for x in glob.glob(os.path.join(sample_dir_path,"*R1*fastq*"))]
    reverse_paths = [os.path.abspath(x) for x in glob.glob(os.path.join(sample_dir_path,"*R1*fastq*"))]

    # sample_names = [x.split]

    return os.path.join(temp_dir.name, "manifest")



def prep_sequence_reads(ctx, manifest_file_path, primer_f=None, primer_r=None):
    results = []

    # importing external plugins to be used later
    cut_adapt = ctx.get_action('cutadapt', 'trim_paired')
    create_table_viz = ctx.get_action('demux', 'summarize')

    #importing sequences
    read_seqs = qiime2.Artifact.import_data(type='SampleData[PairedEndSequencesWithQuality]',
                                            view=manifest_file_path,
                                            view_type='PairedEndFastqManifestPhred33V2')

    # using plugins to trim reads and create reads visualization
    trimmed_reads = cut_adapt(demultiplexed_sequences=read_seqs, front_f=[primer_f], front_r=[primer_r])
    table_viz = create_table_viz(data=trimmed_reads.trimmed_sequences)

    results += trimmed_reads
    results += table_viz

    return tuple(results)