"""
Python script for c-Myc
"""

import sys
import subprocess # for running r through later
import re
import math
import csv

psm_file = sys.argv[1]


def protein_get_name(protein):
    return protein['name']


class Protein:
    def __init__(self, name, description, gene_name, hours):
        self.name = name # Protein accession
        self.description = description
        self.gene_name = gene_name
        self.hours = hours

    # For double checking: will display the prot accession and description if
    # called as a string object
    def __str__(self):
        return "{0} {1}".format(self.name, self.description)

    def get_name(self):
        return self.name

    def get_hr(self, hr):
        return self.hours[hr]

    # Class method to directly access FC log 2
    def get_fclog2(self, hr):
        return math.log(self.hours[hr] / self.hours[0], 2)


prot_dict = {}


def get_first_protein_accession(raw_accessions):
    return re.match("(^.*?)(?:;|\Z)", raw_accessions).group(0)


with open(psm_file) as psm:
    psmreader = csv.reader(psm)

    next(psmreader)  # Skip header
    for row in psmreader:
        # Parse protein accession
        raw_prot_accession = row[12]
        if raw_prot_accession == "sp":
            continue
        mast_prot_accession = get_first_protein_accession(raw_prot_accession)

        # Parse description
        full_desc = row[14]
        # Get everything up to the first optional ;
        prot_description = re.match("(^.*?)(?:;|\Z)", full_desc).group(0)
        if len(prot_description) == 0:
            continue

        # Parse the protein description
        results = re.findall("(^.+)OS=.*GN=(\S*)", prot_description)
        if len(results) == 0:
            continue
        result = results[0]

        desc = result[0]
        gene = result[1]

        # Storing all timepoints in a dictionary, to be contained w/in object
        hours_all = {0: row[51], 4: row[52], 8: row[53], 12: row[54],
                     24: row[55], 48: row[56]}

        prot = Protein(mast_prot_accession, desc, gene, hours_all)

        prot_dict[mast_prot_accession] = prot

for value in prot_dict:
    print(prot_dict[value])
