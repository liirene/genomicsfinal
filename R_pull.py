import subprocess
from typing import Tuple, List, Dict
def pathway_parse(path_to_script):
    command = 'C:/Program Files/R/R-3.3.2/bin/Rscript.exe'
    args = ['clean_data.csv']
    cmd = [command, path_to_script] + args

    x = subprocess.check_output(cmd, universal_newlines=True)
    return x

# Upreg is before downreg
def split_r_output(output_lines: List[str]) -> Tuple[List[str], List[str]]:
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
    return (upreg, downreg)

def parse_gene_list(gene_list: List[str]) -> Dict[str, List[str]]:
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


output = pathway_parse('pathwayscript.R')
#print(output)
output_lines = output.split("\n")
output_lines = list(filter(lambda line: line.startswith("[1]") is False, output_lines))
upreg, downreg = split_r_output(output_lines)
upreg_pathways = parse_gene_list(upreg)
downreg_pathways = parse_gene_list(downreg)

# Do this for each group of pathways that you're interested in (1 figure for
# each group of pathways)
#pathwayname = 'pathhere'
#title = str(pathwayname)

# get each protein name into the tuple form


#array = np.ndarray()
#f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3,
                                                                      #squeeze=False,
                                                                      #sharey='row')

#plt.suptitle(title, fontsize=12)

# include however many you want w/in the tuple

# Include subplot (and making into multiple plots) here

# Add pathway title name at top of the entire figure

#path = 'path'
# name this as the name of the pathway itself
#plt.savefig('path' + title + '.png', bbox_inches='tight')
# Make multiple graphs at once?
# Make 8 at once to be default

#plt.close('all')  # closes all open figures
