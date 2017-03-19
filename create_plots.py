"""
create_plots.py
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Set, Dict, List

def search_current_data(available_genes: Set[str], genes_in_pathway: Set[str]) -> Set[str]:
    """
    Finds intersection between two datasets
    """
    return available_genes.intersection(genes_in_pathway)

def plot_protein(protein, subplot): # Gene name is a protein object from cleandata.py
    # Will return one subplot object at a time

    make_subplot(subplot, protein.get_gene_name(), protein.get_hours_fclog2())

def plot(full_protdict: Dict, prots_to_plot: List[str], pathway_name: str):

    # Create subplot array here
    fig, subplots = create_master(pathway_name)
    for i in range(0,9):
        plot_protein(full_protdict[prots_to_plot[i]], subplots[int(i / 3)][i % 3])

    fig.savefig(pathway_name + '.png', dpi=200)


def make_subplot(subplot, gene_name, y, num=3, x=np.array([0,4,8,12,24,48]),
                 ylabel='Log2 Fold Change'):
    """

    :param prot_name:
    :param y:
    :param num:
    :param x:
    :param ylabel:
    :return:
    """

    subplot.scatter(x= x, y = y, c= 'red') # change color here if using color map
    subplot.set_title(gene_name)
    subplot.set_ylabel(ylabel)
    subplot.set_ylim(0, max(y) + 1)
    calc_regression(x, y, subplot)


def create_master(title:str): # Subplot is a list of plotting objects
    """

    :param subplots:
    :return:
    """
    fig, ax = plt.subplots(3,3, squeeze=False, sharex='row')  # initializes figure and axes object
    plt.suptitle(title, fontsize=12)
    plt.xlim(0, 50.0)
    plt.xlabel('Hours')
    return (fig, ax)

def calc_regression(x, y, subplot):
    """

    :param x:
    :param y:
    :param degree:
    :return:
    """
    fit = np.polyfit(x=x, y=y, deg = 2)
    # figure this out here

    subplot.plot(x, fit[0] * x ** 2 + fit[1] * x ** 1 + fit[2], color='blue')

    #subplot.plot(x, fit[0] * x ** 3 + fit[1] * x ** 2 + fit[2] * x + fit[3],
            #color='blue')  # add regression line to plot

# calculate coefficients and store them: array outputs are:
# polynomial coeffs, highest power first (3,2,1,0)

# save shit to disk

#plt.close('fig')
