"""
Python script for c-Myc

Master script designed to be used in conjunction with R_pull.py and create_
plots.py, and located in the same directory as these files.

Generates pathway images, regression plots,

"""

import sys
import re
import math
import csv
import R_pull
import create_plots


# windows is python3 path_analysis.py
## C:/Program Files/R/R-3.3.2/bin/Rscript.exe 90mfull_PSM.csv
# Linux is python3 path_analysis.py Rscript 90full_PSM.csv


class Protein:
    def __init__(self, name: str, description: str, gene_name: str, hours):
        self.name = name  # Protein accession
        self.description = description
        self.gene_name = gene_name
        self.hours = hours

    def get_name(self):
        return self.name

    def get_gene_name(self):
        return self.gene_name

    def get_hr(self, hr):
        return float(self.hours[hr])

    # Class method to directly access FC log 2
    def get_fclog2(self, hr):
        return math.log(self.get_hr(hr) / self.get_hr(0), 2)

    def get_all_fclog2(self):
        return list(
            map(lambda hrlog2: self.get_fclog2(hrlog2), [0, 4, 8, 12, 24, 48]))


def get_first_protein_accession(raw_accessions: str):
    """
    Parses the raw accession, which may include ; or whitespace artifacts.
    :param raw_accessions: Provide the raw raw protein accession as a string
    :return: Cleaned accession without artifacts
    """
    return re.match("(^.*?)(?=(?:;|\s|$))", raw_accessions).group(0)


def plot_pathways(pathway_set, is_upreg: bool, foldchange_hr: int):
    """
    Works w/ create_plots.py to plot all proteins provided to the program

    :param foldchange_hr: Foldchange on which to determine the significantly
     changed protein set
    :param pathway_set: genes present in the pathway of interest
    :param is_upreg: Is this an upregulated or downregulated pathway set?
    :return: Saves 1 image for each pathway to file of the 9 differentially
    expressed proteins within that set.

    """
    for pathway, genes in pathway_set.items():
        geneset = set(genes)
        intersect_genes = create_plots.search_current_data(available_genes,
                                                           geneset)
        # Find top genes in intersect_genes by fold change
        # Top 9 upregulated
        top_9 = sorted(list(intersect_genes),
                       key=lambda x: prot_dict[x].get_fclog2(foldchange_hr),
                       reverse=is_upreg)[:9]

        create_plots.plot(prot_dict, top_9, pathway, is_upreg)


psm_file = sys.argv[2]
prot_dict = {}
timepoints = [0, 4, 8, 12, 24, 48]
cleaned_csv_name = "clean_data.csv"

with open(psm_file) as psm:
    psmreader = csv.reader(psm)

    next(psmreader)  # Skip header
    for row in psmreader:

        # Parse protein accession
        raw_prot_accession = row[12]
        if raw_prot_accession == "sp":
            continue
        if raw_prot_accession == "":
            continue
        mast_prot_accession = get_first_protein_accession(raw_prot_accession)

        # Parse description to pull out gene name and protein description
        full_desc = row[14]
        results = re.findall("(^.+?)(?:OS.+?GN=)(.+?)(?:;|\s)", full_desc)
        if len(results) == 0:
            continue
        result = results[0]

        desc = result[0]
        gene = result[1]

        # Storing all timepoints in a dictionary, to be contained w/in object
        hours_all = {}
        data_valid = True
        for i, time in enumerate(timepoints):
            hour = row[51 + i]
            if hour == "":
                data_valid = False
                break
            hours_all[time] = hour
        if data_valid is False:
            continue

        prot_dict[gene] = Protein(mast_prot_accession, desc, gene, hours_all)

with open(cleaned_csv_name, 'w') as newfile:
    writefile = csv.writer(newfile)
    writefile.writerow(
        ["Prot Accession", "Gene Name", "0h", "4h", "8h", "12h", "24h", "48h"])

    for key, value in prot_dict.items():
        spread_fc = []
        for time in timepoints:
            spread_fc.append(str(value.get_fclog2(time)))
        line_to_write = [value.get_name(), value.get_gene_name()]
        line_to_write += spread_fc
        writefile.writerow(line_to_write)

# Streams the input through R and generates pathways significantly impacted
# after 48 h time points.
upreg, downreg = R_pull.pathway_analysis(sys.argv[1], 'pathwayscript.R')

# Generates a list of genes available within the dataset
available_genes = set(prot_dict.keys())

plot_pathways(upreg, True, 48)
plot_pathways(downreg, False, 48)
