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

from q2_pan_classifier.format_types import MyStringFormat
from q2_pan_classifier.plugin_setup import plugin

from q2_types.feature_data import TSVTaxonomyFormat, FeatureData, Sequence, DNAFASTAFormat

@plugin.register_transformer
def _1(data: str) -> MyStringFormat:
    print("testing one, two, three")
    cool = MyStringFormat()
    # cool.__name__ = "Chase"
    with open(cool.path, 'w') as ff:
        ff.write(data)

    return cool

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
