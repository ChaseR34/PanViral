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
import qiime2.plugin.model as model
from qiime2.plugin import SemanticType
import qiime2.core.path as qpath

DNAFastaNCBI = SemanticType('DNAFastaNCBI')


class DNAFastaNCBIFormatError(ValueError):
	"""Custom error for invalid format """



class DNAFastaNCBIFormat(model.TextFileFormat):

	def __init__(self, path: str = None, mode: str = 'w') -> None:
		self.PIPE = '|'
		self.accession_numbers = []
		self.names = []

		if path is None:
			if mode != 'w':
				raise ValueError("A path must be provided when reading.")
		else:
			if mode != 'r':
				raise ValueError("A path must be omitted when writing.")

		if mode == 'w':
			self.path = qpath.OutPath(
				# TODO: parents shouldn't know about their children
				dir=isinstance(self, model.DirectoryFormat),
				prefix='q2-%s-' % self.__class__.__name__)
		else:
			self.path = qpath.InPath(path)

		self._mode = mode

	def get_accession_numbers(self) -> None:
		samps = skbio.read(str(self), format="fasta")

		for samp in samps:
			name = samp.metadata['id']
			if self.PIPE in name:
				a_n_tmp = name.split(self.PIPE)[1]
				self.accession_numbers.append(a_n_tmp)
			else:
				self.accession_numbers.append(name)


	def _validate_(self, level):
		pass



DNAFastaNCBIDirFormat = model.SingleFileDirectoryFormat('DNAFastaNCBIDirFormat', 'taxonomy.tsv', DNAFastaNCBIFormat)

