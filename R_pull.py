import subprocess
from typing import Tuple, List, Dict

Pathway = Dict[str, List[str]]


def __pathway_parse(path_to_rscript: str, path_to_script: str):
    """

    :param path_to_script:
    :return:
    """
    cmd = [path_to_rscript, path_to_script, 'clean_data.csv']

    x = subprocess.check_output(cmd, universal_newlines=True)
    return x


def __split_r_output(output_lines: List[str]) -> Tuple[List[str], List[str]]:
    """

    :param output_lines:
    :return:
    """
    upreg = []
    downreg = []
    is_upreg = True
    for line in output_lines:
        if line == "":
            is_upreg = False
        elif is_upreg is True:
            upreg.append(line)
        else:
            downreg.append(line)
    return upreg, downreg


def __parse_gene_list(gene_list: List[str]) -> Pathway:
    """

    :param gene_list:
    :return:
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

    :param path_to_rscript:
    :param path_to_script:
    :return:
    """
    output = __pathway_parse(path_to_rscript, path_to_script)
    output_lines = output.split("\n")
    output_lines = list(filter(lambda line: line.startswith("[1]") is False,
                               output_lines))
    upreg, downreg = __split_r_output(output_lines)
    upreg_pathways = __parse_gene_list(upreg)
    downreg_pathways = __parse_gene_list(downreg)
    return upreg_pathways, downreg_pathways
