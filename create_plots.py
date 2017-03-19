"""
create_plots.py
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Set, Dict, List


def search_current_data(available_genes: Set[str],
                        genes_in_pathway: Set[str]) -> Set[str]:
    """
    Finds intersection between two datasets
    """
    return available_genes.intersection(genes_in_pathway)


# Gene name is a protein object from cleandata.py
def plot_protein(protein, subplot, is_upreg):
    # Will return one subplot object at a time

    make_subplot(subplot, protein.get_gene_name(), protein.get_hours_fclog2(),
                 is_upreg)


def plot(full_protdict: Dict, prots_to_plot: List[str], pathway_name: str,
         is_upreg: bool):
    file_prefix = "upreg_"
    if is_upreg is False:
        file_prefix = "downreg_"
    # Create subplot array here
    fig, subplots = __create_master(pathway_name)
    for i in range(0, 9):
        plot_protein(full_protdict[prots_to_plot[i]],
                     subplots[int(i / 3)][i % 3], is_upreg)

    plt.subplots_adjust(hspace=0.4, wspace=0.40)

    fig.savefig(file_prefix + pathway_name + '.png', dpi=300)

    # TODO: Add a legend


def make_subplot(subplot, gene_name: str, y, is_upreg: bool,
                 x=np.array([0, 4, 8, 12, 24, 48])):
    """

    :param subplot:
    :param gene_name:
    :param y:
    :param x:
    :return:
    """
    scatter_color = 'green'
    if is_upreg is False:
        scatter_color = 'red'

    y_ticks = np.arange(-2, 10, 0.5)
    subplot.scatter(x=x, y=y, s=3, c=scatter_color)
    subplot.set_title(gene_name, fontsize=7, weight='demibold')


    subplot.tick_params(which='both', direction='in')
    subplot.set_xticks(x)
    subplot.set_xticklabels(x, fontsize=4, color="grey")
    subplot.set_yticks(y_ticks)
    subplot.set_yticklabels(y_ticks, fontsize=4, color="grey")
    subplot.set_ylim(min(y) - 0.5, max(y) + 1)

    subplot.axhline(y=0, c='grey', lw=0.5, ls='dashed')

    __calc_regression(x, y, subplot)


def __create_master(title: str):  # Subplot is a list of plotting objects
    """

    :param title:
    :return:
    """
    # initializes figure and axes object
    fig, ax = plt.subplots(3, 3, squeeze=False)
    # Specifying tight layout for multiple subplots

    plt.suptitle(title, fontsize=15)
    plt.xlim(-5, 50)

    # Set X and Y Axes Labels

    fig.text(x=.45, y=.04, s='Timepoints',  size=10, weight='semibold')
    fig.text(x=.04, y=.6, s='Log2 Fold Change', size=10, rotation='vertical',
             weight='semibold')

    return fig, ax


def __calc_regression(x, y, subplot):
    """

    :param x:
    :param y:
    :param subplot:
    :return:
    """
    fit = np.polyfit(x=x, y=y, deg=2)

    # TODO: Calculate R values and other parameters

    subplot.plot(x, fit[0] * x ** 2 + fit[1] * x ** 1 + fit[2], lw=0.5, c='b')

    # subplot.plot(x, fit[0] * x ** 3 + fit[1] * x ** 2 + fit[2] * x + fit[3],
    # color='blue')  # add regression line to plot

# calculate coefficients and store them: array outputs are:
# polynomial coeffs, highest power first (3,2,1,0)


# plt.close('fig')
