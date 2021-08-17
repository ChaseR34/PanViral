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

import skbio
from Bio import Entrez
import re
import qiime2.plugin.model as model
from qiime2.plugin import SemanticType
import qiime2.core.path as qpath

DNAFastaNCBI = SemanticType('DNAFastaNCBI')


class DNAFastaNCBIFormatError(ValueError):
	"""Custom error for invalid format """
	#TODO: Finish making this error class



class DNAFastaNCBIFormat(model.TextFileFormat):

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
			if self.PIPE in name:
				a_n_tmp = name.split(self.PIPE)[1]
				self.accession_numbers.append(a_n_tmp)
			else:
				self.accession_numbers.append(name)

	def get_taxonomy(self, accession_numbers):

		if not isinstance(accession_numbers, list):
			accession_numbers = [accession_numbers]

		handle = Entrez.efetch(db="nuccore",
							   id=accession_numbers,
							   rettype="gb",
							   retmode="xml")
		records = Entrez.read(handle)

		for rec in records:
			taxonomy = rec['GBSeq_taxonomy']
			scientific_name = rec['GBSeq_organism']
			self.taxonomy.append(taxonomy + "; " + scientific_name)

		handle.close()

	def _validate_(self, level):
		#TODO: make validatoin function to check if ncbi taxonomic names are there
		pass



DNAFastaNCBIDirFormat = model.SingleFileDirectoryFormat('DNAFastaNCBIDirFormat', 'taxonomy.tsv', DNAFastaNCBIFormat)

