from Bio import Entrez
import re
import skbio
import threading
import time
PIPE = '|'
class NCBI_Tax:
    def __init__(self, email):
        self.email = email
        self.taxid = []
        self.taxonomy = []

        Entrez.email = self.email

    def get_taxid(self, accession_numbers):

        if not isinstance(accession_numbers, list):
            accession_numbers = [accession_numbers]

        handle = Entrez.efetch(db="nuccore",
                               id=accession_numbers,
                               mode="gb",
                               rettype="text")
        records = Entrez.read(handle)

        for rec in records['Bioseq-set_seq-set']:
            # self.taxid = records['Bioseq-set_seq-set'][0]['Seq-entry_set']['Bioseq-set']['Bioseq-set_descr']["Seq-descr"][0]["Seqdesc_source"]["BioSource"]["BioSource_org"]["Org-ref"]["Org-ref_db"][0]["Dbtag_tag"]["Object-id"]["Object-id_id"]
            tax_match = re.search('Object-id_id.*?([0-9]+)', str(rec))

            self.taxid.append(int(tax_match.group(1)))

        handle.close()

        return

    def get_taxonomy(self, taxid):

        handle = Entrez.efetch(db="taxonomy", id=int(taxid), mode="text", rettype="xml")
        records = Entrez.read(handle)
        scientific_name = records[0]["ScientificName"]
        taxonomy = records[0]["Lineage"]
        self.taxonomy.append(taxonomy + "; " + scientific_name)
        handle.close()

        return

def validate_ncbi_names(path: str=None) -> list:
    accession_numbers = []
    samps = skbio.read(path, format="fasta")

    for samp in samps:
        name = samp.metadata['id']
        if PIPE in name:
            a_n_tmp = name.split(PIPE)[1]
            accession_numbers.append(a_n_tmp)
        else:
            accession_numbers.append(name)

    return accession_numbers


if __name__ == "__main__":
    # a_n = ['NZ_JAEUAI010000014', 'MW183195.1', 'MW183195.1',
    #        'NZ_JAEUAI010000014', 'MW183195.1','NZ_JAEUAI010000014',
    #        'MW183195.1','NZ_JAEUAI010000014', 'MW183195.1',
    #        'NZ_JAEUAI010000014','MW183195.1'
    #        ]
    # ncbi_tax = NCBI_Tax("clr96@nau.edu")
    #
    # ncbi_tax.get_taxid(a_n)
    #
    # threads = []
    #
    # out_list = []
    #
    # for index, i in enumerate(ncbi_tax.taxid):
    #
    #     if index % 3 == 0:
    #         time.sleep(1)
    #
    #     x = threading.Thread(target=ncbi_tax.get_taxonomy, args=(i,))
    #     threads.append(x)
    #     x.start()
    #
    # for t in threads:
    #     t.join()
    # [print(x) for x in ncbi_tax.taxonomy]

    names = validate_ncbi_names("../PanFlaviTutorial/example_viral_db/moureau_2015_ref_sequences.fasta")
    print(names)

