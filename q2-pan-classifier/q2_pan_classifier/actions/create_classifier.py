#   Copyright 2021 Chase Ridenour
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

from os import path

import pandas as pd
import pkg_resources
from jinja2 import Environment, FileSystemLoader

import qiime2
from q2_types.feature_data import FeatureData, Sequence


def generate_taxonomy(ref_seqs: pd.Series) -> list:

    seq_names = [name.metadata['id'] for name in ref_seqs]

    return seq_names

def create_classifier(ctx,
                      ref_seqs_file,
                      ref_tax=None,
                      ref_tax_file=None,
                      f_primer=None,
                      r_primer=None,
                      min_len=None,
                      max_len=None):

    #importing external plugins to be used later
    extract_refs = ctx.get_action('feature_classifier', 'extract_reads')
    train_classifier = ctx.get_action('feature_classifier', 'fit_classifier_naive_bayes')

    results = []

    #importing reference sequences and taxonomy
    ref_seqs = qiime2.Artifact.import_data(type='FeatureData[Sequence]',
                                           view=ref_seqs_file,
                                           view_type=None)

    if ref_tax and ref_tax_file:
        raise ValueError("Please only provide ref_tax OR ref_tax_file not both")

    if ref_tax:
        ref_tax_out = ref_tax
    elif ref_tax_file:
        ref_tax_out = qiime2.Artifact.import_data(type='FeatureData[Taxonomy]',
                                                  view=ref_tax_file,
                                                  view_type=None)
    else:
        ref_tax_out = qiime2.Artifact.import_data(type='FeatureData[Taxonomy]',
                                                  view=ref_seqs_file,
                                                  view_type='DNAFastaNCBIFormat')

    # using imported plugins to extract reference and train classifier
    if f_primer and r_primer:
        trimmed_refs = extract_refs(sequences=ref_seqs,
                                    f_primer=f_primer,
                                    r_primer=r_primer,
                                    min_length=min_len,
                                    max_length=max_len)

        trained_class = train_classifier(reference_reads=trimmed_refs.reads,
                                         reference_taxonomy=ref_tax_out)
        results += trimmed_refs
    else:
        trained_class = train_classifier(reference_reads=ref_seqs,
                                         reference_taxonomy=ref_tax_out)
        results += [ref_seqs]


    results += [ref_tax_out]
    results += trained_class

    return tuple(results)


def prep_sequence_reads(ctx, manifest_file_path, primer_f, primer_r):
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


def classify_reads(ctx, samp_reads, trunc_len_f, trunc_len_r, trained_classifier):
    results = []

    # action importing
    dada2 = ctx.get_action('dada2', 'denoise_paired')
    classify_sklearn = ctx.get_action('feature_classifier', 'classify_sklearn')
    barplot = ctx.get_action('taxa', 'barplot')
    transpose = ctx.get_action('feature_table', 'transpose')
    tabulate = ctx.get_action('metadata', 'tabulate')
    # vis_test = ctx.get_action('pan_classifier', 'visualization_final')

    # getting some output
    dada2_table, dada2_rep_seqs, dada2_stats = dada2(demultiplexed_seqs=samp_reads,
                                                     trunc_len_f=trunc_len_f,
                                                     trunc_len_r=trunc_len_r
                                                     )
    classified, = classify_sklearn(classifier=trained_classifier, reads=dada2_rep_seqs)
    barplot_taxonomy = barplot(table=dada2_table, taxonomy=classified)

    # transposing table and getting metadata
    tt, = transpose(table=dada2_table)

    tt_m = tt.view(view_type=qiime2.Metadata)
    dr_m = dada2_rep_seqs.view(view_type=qiime2.Metadata)
    c_m = classified.view(view_type=qiime2.Metadata)

    # merging table
    merge_table = tabulate(c_m.merge(dr_m, tt_m))

    results += [dada2_table, dada2_rep_seqs, dada2_stats, classified]
    results += barplot_taxonomy
    results += merge_table

    return tuple(results)

def visualization_final(output_dir: str) -> None:

    # temp_dir = tempfile.TemporaryDirectory()
    template_data = pkg_resources.resource_filename('q2_pan_classifier', 'templates')
    jin_env = Environment(loader=FileSystemLoader(template_data), auto_reload=True)
    # shutil.copy2(path.join(template_data, 'base.html'), temp_dir.name)

    # jin_env = Environment(loader=FileSystemLoader("templates"))
    jin_out = jin_env.get_template('base.html').render(title="This is my title")

    with open(path.join(output_dir, 'index.html'), 'w') as f:
        f.write(jin_out)
