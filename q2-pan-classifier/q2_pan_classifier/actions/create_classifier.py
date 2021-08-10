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
import shutil
import tempfile

from jinja2 import Environment, FileSystemLoader

import qiime2
from qiime2.sdk import Results
from qiime2.core.type.primitive import Str
from q2_types.sample_data import SampleData


from q2_pan_classifier.format_types import MyStringFormat

from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import PairedEndFastqManifestPhred33V2, PairedEndSequencesWithQuality, \
    SingleLanePerSamplePairedEndFastqDirFmt
from q2_types.feature_data import Sequence

from qiime2.plugin import model


def test_function() -> str:
    # print("testing one, two, three")
    # cool = MyStringFormat()
    # # cool.__name__ = "Chase"
    # with open(cool.path, 'w') as ff:
    #     ff.write("Chase's Cool Project")

    return "Chase's Cool Project"


def import_sequences() -> PairedEndFastqManifestPhred33V2:
    pass
    # directory = SingleLanePerSamplePairedEndFastqDirFmt()
    #
    # out.export_data(directory.path)
    # out.save()
    #
    #
    # return out


# , : str, forward_primer: str, reverse_primer: str

def create_classifier(ctx,
                      ref_seqs_file,
                      ref_tax_file,
                      f_primer,
                      r_primer,
                      min_len,
                      max_len):

    extract_refs = ctx.get_action('feature_classifier', 'extract_reads')
    train_classifier = ctx.get_action('feature_classifier', 'fit_classifier_naive_bayes')

    results = []

    ref_seqs = qiime2.Artifact.import_data(type='FeatureData[Sequence]',
                                           view=ref_seqs_file,
                                           view_type=None)

    ref_tax = qiime2.Artifact.import_data(type='FeatureData[Taxonomy]',
                                          view=ref_tax_file,
                                          view_type='HeaderlessTSVTaxonomyFormat')

    trimmed_refs = extract_refs(sequences=ref_seqs,
                                f_primer=f_primer,
                                r_primer=r_primer,
                                min_length=min_len,
                                max_length=max_len)

    trained_class = train_classifier(reference_reads=trimmed_refs.reads,
                                     reference_taxonomy=ref_tax)

    results += trimmed_refs
    results += trained_class

    return tuple(results)


def prep_sequence_reads(ctx, manifest_file_path, primer_f, primer_r):
    results = []

    cut_adapt = ctx.get_action('cutadapt', 'trim_paired')
    create_table_viz = ctx.get_action('demux', 'summarize')

    read_seqs = qiime2.Artifact.import_data(type='SampleData[PairedEndSequencesWithQuality]',
                                            view=manifest_file_path,
                                            view_type='PairedEndFastqManifestPhred33V2')

    trimmed_reads = cut_adapt(demultiplexed_sequences=read_seqs, front_f=[primer_f], front_r=[primer_r])
    table_viz = create_table_viz(data=trimmed_reads.trimmed_sequences)

    results += trimmed_reads
    results += table_viz

    return tuple(results)


def classify_reads(ctx, samp_reads, trunc_len_f, trunc_len_r, trained_classifier):
    results = []

    # action importing
    dada2 = ctx.get_action('dada2', 'denoise_paired')
    classify = ctx.get_action('feature_classifier', 'classify_sklearn')
    taxa_barplot = ctx.get_action('taxa', 'barplot')
    trans_table = ctx.get_action('feature_table', 'transpose')
    merge_tables = ctx.get_action('metadata', 'tabulate')
    vis_test = ctx.get_action('pan_classifier', 'visualization_final')

    # getting some output
    dada2_table, dada2_rep_seqs, dada2_stats = dada2(demultiplexed_seqs=samp_reads,
                                                     trunc_len_f=trunc_len_f,
                                                     trunc_len_r=trunc_len_r
                                                     )
    classified, = classify(classifier=trained_classifier, reads=dada2_rep_seqs)
    barplot_taxonomy = taxa_barplot(table=dada2_table, taxonomy=classified)

    # transposing table and getting metadata
    tt, = trans_table(table=dada2_table)

    tt_m = tt.view(view_type=qiime2.Metadata)
    dr_m = dada2_rep_seqs.view(view_type=qiime2.Metadata)
    c_m = classified.view(view_type=qiime2.Metadata)

    # merging table
    merge_table_out = merge_tables(c_m.merge(tt_m, dr_m))

    # compiling results
    results += Results(['table'], [dada2_table])
    results += Results(['representative_sequences'], [dada2_rep_seqs])
    results += Results(['denoising_stats'], [dada2_stats])
    results += Results(['classification'], [classified])
    results += barplot_taxonomy
    results += merge_table_out

    return tuple(results)

def visualization_final(output_dir: str) -> None:

    temp_dir = tempfile.TemporaryDirectory()
    template_data = pkg_resources.resource_filename('q2_pan_classifier', 'templates')
    jin_env = Environment(loader=FileSystemLoader(temp_dir.name), auto_reload=True)
    shutil.copy2(path.join(template_data, 'base.html'), temp_dir.name)

    # jin_env = Environment(loader=FileSystemLoader("templates"))
    jin_out = jin_env.get_template('base.html').render(title="This is my title")

    with open(path.join(output_dir, 'index.html'), 'w') as f:
        f.write(jin_out)
