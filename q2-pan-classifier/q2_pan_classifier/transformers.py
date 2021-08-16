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
import pandas as pd
import skbio

from q2_pan_classifier.format_types import DNAFastaNCBIFormat
from q2_pan_classifier.plugin_setup import plugin

from q2_types.feature_data import TSVTaxonomyFormat, FeatureData, Sequence, DNAFASTAFormat



@plugin.register_transformer
def _2(ref_seqs: list) -> TSVTaxonomyFormat:

    tax_out = TSVTaxonomyFormat()

    with open(tax_out.path, 'w') as ff:
        ff.write('\t'.join(tax_out.HEADER) + '\n')
        for name in ref_seqs:
            ff.write('\t'.join([name, 'virus']) + '\n')

    return tax_out


@plugin.register_transformer
def _3(ref_seqs: DNAFASTAFormat) -> TSVTaxonomyFormat:

    seq_names = [name.metadata['id'] for name in ref_seqs.view(pd.Series)]

    tax_out = TSVTaxonomyFormat()

    with open(tax_out.path, 'w') as ff:
        ff.write('\t'.join(tax_out.HEADER) + '\n')
        for name in seq_names:
            ff.write('\t'.join([name, 'virus']) + '\n')

    return tax_out



# def _read_dna_fasta(path):
#     return skbio.read(path, format='fasta', constructor=skbio.DNA)
# #
#
# @plugin.register_transformer
# def _4(ff: DNAFastaNCBIFormat) -> pd.Series:
#     data = {}
#     for sequence in _read_dna_fasta(str(ff)):
#         data[sequence.metadata['id']] = sequence
#     return pd.Series(data)


@plugin.register_transformer
def _4(ff: DNAFastaNCBIFormat) -> list:
    return ff.get_accession_numbers()

# @plugin.register_transformer
# def _5(ref_seqs: DNAFastaNCBIFormat) -> TSVTaxonomyFormat:
#
#     seq_names = [name.metadata['id'] for name in ref_seqs.view(pd.Series)]
#
#     tax_out = TSVTaxonomyFormat()
#
#     with open(tax_out.path, 'w') as ff:
#         ff.write('\t'.join(tax_out.HEADER) + '\n')
#         for name in seq_names:
#             ff.write('\t'.join([name, 'virus']) + '\n')
#
#     return tax_out

@plugin.register_transformer
def _5(ref_seqs: DNAFastaNCBIFormat) -> TSVTaxonomyFormat:
    ref_seqs.get_accession_numbers()
    seq_names = ref_seqs.accession_numbers

    tax_out = TSVTaxonomyFormat()

    with open(tax_out.path, 'w') as ff:
        ff.write('\t'.join(tax_out.HEADER) + '\n')
        for name in seq_names:
            ff.write('\t'.join([name, 'virus_list2']) + '\n')

    return tax_out