import subprocess
from typing import Tuple, List, Dict

Pathway = Dict[str, List[str]]


def __pathway_parse(path_to_rscript: str, path_to_script: str):
    """
    Run the R script to develop pathway analysis.
    :param path_to_script: String of the path to script.
    :return: Printed output from R: Upregulated genes and downregulate gene
    sets. Each gene is on a new line, pathway sets are separated by '#' headers
    """
    cmd = [path_to_rscript, path_to_script, 'clean_data.csv']

    x = subprocess.check_output(cmd, universal_newlines=True)
    return x


def __split_r_output(output_lines: List[str]) -> Tuple[List[str], List[str]]:
    """
    Splits R output by new line.

    :param output_lines: Full output printed from R, where upregulated and
    downregulated gene sets are separated by a blank line
    :return: Two lists: one for upregulated and downregulated genes
    """
    upreg = []
    downreg = []
    is_upreg = True

    for line in output_lines:
        if line == "":
            is_upreg = False  # Find where the separation between upreg and
            # downreg are.
        elif is_upreg is True:
            upreg.append(line)
        else:
            downreg.append(line)
    return upreg, downreg


def __parse_gene_list(gene_list: List[str]) -> Pathway:
    """
    Parses the given gene list giving a list of Pathway objects (dict
    where key is the name of the pathway and includes a list of gene names as
    strings)
    :param gene_list: Gene list, new line separated
    :return: Pathway object.
    """
    pathways = {}
    current_pathway = ""

    for line in gene_list:
        if line.startswith("# "):
            line = line.lstrip("# ")
            pathways[line] = []
            current_pathway = line
        else:
            pathways[current_pathway] += [line.rstrip("\n\t ")]
    return pathways


def pathway_analysis(path_to_rscript: str, path_to_script: str) -> \
        Tuple[Pathway, Pathway]:
    """
    Master function to complete all pathway analysis and parse the list of
    upregulated and downregulated genes printed from R

    :param path_to_rscript: String of pathway to R scripts
    :param path_to_script: String of pathway to R script (differs if on windows
    or Linux)
    :return: List of upregulated and downregulate gene pathways. R pathway
    analysis results are saved as CSV, and images saved in the directory file.
    """
    output = __pathway_parse(path_to_rscript, path_to_script)
    output_lines = output.split("\n")

    # Remove artifacts of the R output (error messages and such)
    output_lines = list(filter(lambda line: line.startswith("[1]") is False,
                               output_lines))

    upreg, downreg = __split_r_output(output_lines)
    upreg_pathways = __parse_gene_list(upreg)
    downreg_pathways = __parse_gene_list(downreg)
    return upreg_pathways, downreg_pathways
