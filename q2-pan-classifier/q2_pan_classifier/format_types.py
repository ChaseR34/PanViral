#   Copyright 2021 Chase Ridenour

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
import time

import skbio
from Bio import Entrez
import re
import qiime2.plugin.model as model
from q2_types.feature_data import DNAFASTAFormat
from qiime2.plugin import SemanticType

DNAFastaNCBI = SemanticType('DNAFastaNCBI')


class DNAFastaNCBIFormat(DNAFASTAFormat):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.PIPE = '|'
        self.accession_numbers = []
        self.names = []
        self.taxonomy = []
        self.email = "clr96@nau.edu"
        Entrez.email = self.email

    def get_accession_numbers(self) -> None:
        samps = skbio.read(str(self), format="fasta")

        for samp in samps:
            name = samp.metadata['id']
            self.names.append(name)
            if "gb:" and self.PIPE in name:
                a_n_tmp = name.split(":")[1].split(self.PIPE)[0]
                self.accession_numbers.append(a_n_tmp)
            elif self.PIPE in name:
                a_n_tmp = name.split(self.PIPE)[1]
                self.accession_numbers.append(a_n_tmp)
            else:
                a_n_tmp = name.split(" ")[0]
                self.accession_numbers.append(a_n_tmp)

    def get_taxonomy(self):

        if not self.accession_numbers:
            raise ValueError("Missing accession numbers")

        tax_set = set()

        def _accession_number_split(sub_list_size: int = 100) -> list:

            out_list = list()
            ac_len = len(self.accession_numbers)

            total_count = 0
            sublist = list()

            while total_count < ac_len:

                index = 0
                while index < sub_list_size and total_count < ac_len:
                    sublist.append(self.accession_numbers[total_count])
                    index += 1
                    total_count += 1

                out_list.append(sublist.copy())
                sublist.clear()

            return out_list

        def _check_set(tax_set_input: set, taxon: str) -> str:

            for tax in tax_set_input:
                if taxon in tax:
                    return tax
            return taxon

        def _get_taxonomy(ac_numbers_subset: list):

            handle = Entrez.efetch(db="nuccore",
                                   id=ac_numbers_subset,
                                   rettype="gb",
                                   retmode="xml")
            records = Entrez.read(handle)

            for rec in records:
                tax_tmp = rec['GBSeq_taxonomy']
                tax_split = tax_tmp.split(';')

                if len(tax_split) > 8:
                    tax_tmp = ';'.join(tax_split[0:8])
                elif len(tax_split) < 8:
                    tax_tmp = prev_taxonomy
                if tax_tmp in tax_set:
                    taxonomy = tax_tmp
                else:
                    taxonomy = _check_set(tax_set, tax_tmp)

                prev_taxonomy = taxonomy
                tax_set.add(taxonomy)

                scientific_name = rec['GBSeq_organism']
                self.taxonomy.append(taxonomy + "; " + scientific_name)

            handle.close()

        accession_number_subsets = _accession_number_split()

        for subset in accession_number_subsets:
            _get_taxonomy(subset)
            time.sleep(1)

    def _validate_(self, level):
        super()._validate_(level=level)

        # fasta file check
        # get it from DNAFASTA
        # Check if sequence names have pipes
        # Check if all accession numbers are valid
        # if not raise error
        # use model.ValidationError
        # TODO: make validatoin function to check if ncbi taxonomic names are there


DNAFastaNCBIDirFormat = model.SingleFileDirectoryFormat('DNAFastaNCBIDirFormat', 'taxonomy.tsv', DNAFastaNCBIFormat)
