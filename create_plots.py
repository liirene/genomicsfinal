"""
create_plots.py

Contains methods for plotting, statistics, and calculation of gene sets
provided to the script.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from typing import Set, Dict, List


def search_current_data(available_genes: Set[str],
                        genes_in_pathway: Set[str]) -> Set[str]:
    """
    Finds intersection of genes between two datasets, allowing us to pull out
    which genes are present in the proteomic dataset (the dataset unfortunately
    does not contain all genes that we are interested in).
    :param available_genes: Set of gene names present in the proteomic dataset
    :param genes_in_pathway: Set of genes associated with the pathway of
    interest.
    :return: Set of gene names common to both the genes that we have proteomic
    data on and are implicated in the pathway of interest
    """
    return available_genes.intersection(genes_in_pathway)


def plot_protein(protein, subplot, is_upreg: bool):
    """
    Plot one protein's fold change values on the provided subplot.
    :param protein: Protein object in the proteomic dataset
    :param subplot: Axes object specified by the its location in the master
    plot array
    :param is_upreg: Is the protein part of the upregulated or downregulated
    dataset?
    :return: Plotted object on specified subplot
    """
    make_subplot(subplot, protein.get_gene_name(), protein.get_all_fclog2(),
                 is_upreg)


def plot(full_protdict: Dict, prots_to_plot: List[str], pathway_name: str,
         is_upreg: bool):
    """
    Master plotting function, takes in other module wrappers to generate a file
    containing 9 of the most hihgly upregulated proteins within the pathway
    of interest

    :param full_protdict: Dictionary of protein classes
    :param prots_to_plot: List of 9 proteins to be plotted
    :param pathway_name: Name of the pathway, to be included as title
    :param is_upreg: Is the protein part of the upregulated or downregulated
    dataset?
    :return: Saves figure to file.
    """

    # Allows for customization of the file name.
    file_prefix = "upreg_"
    if is_upreg is False:
        file_prefix = "downreg_"

    fig, subplots = __create_master(pathway_name)
    for i in range(0, 9):
        plot_protein(full_protdict[prots_to_plot[i]],
                     subplots[int(i / 3)][i % 3], is_upreg)

    # Add legend lines
    line1 = mlines.Line2D([], [], linewidth=2, color='b')
    plt.figlegend(handles=[line1], labels=['2nd Deg \nRegression'],
                  loc='lower right', fontsize=8)

    # Adds padding to make all labels visible
    plt.subplots_adjust(hspace=0.4, wspace=0.40)

    # Saves figure as .png file. DPI can be adjusted
    fig.savefig(file_prefix + pathway_name + '.png', dpi=300)


def make_subplot(subplot, gene_name: str, y: List[int], is_upreg: bool,
                 x=np.array([0, 4, 8, 12, 24, 48])):
    """
    Creates a subplot for the specified gene name, calculates and plots the
    associated regression line with the dataset

    :param is_upreg: Is the protein in the upregulated or downregulated array?
    :param subplot: Subplot to be plotted (by location in array)
    :param gene_name: Gene name (string)
    :param y: Log2 FC values.
    :param x: Timecourse values. Defaults to the provided values, but others
    can be specified for future analyses
    :return: Plotted values within the master figure object.
    """

    # Changes the color of scatter objects if the protein is part of the upreg
    # ulated or downregulated datasets
    scatter_color = 'green'
    if is_upreg is False:
        scatter_color = 'red'

    # Style parameters and plotting
    y_ticks = np.arange(-2, 10, 0.5)
    subplot.scatter(x=x, y=y, s=3, c=scatter_color)
    subplot.set_title(gene_name, fontsize=7, weight='demibold')

    subplot.tick_params(which='both', direction='in')
    subplot.set_xticks(x)
    subplot.set_xticklabels(x, fontsize=4, color="grey")
    subplot.set_yticks(y_ticks)
    subplot.set_yticklabels(y_ticks, fontsize=4, color="grey")

    # Expand subplot boundaries past the min and max w/in the dataset
    subplot.set_ylim(min(y) - 0.5, max(y) + 1)

    # Sets horizontal line at 0 (for easier visual comparison)
    subplot.axhline(y=0, c='grey', lw=0.5, ls='dashed')

    __calc_regression(x, y, subplot)


def __create_master(title: str):  # Subplot is a list of plotting objects
    """
    Creates a master plot with associated and subplots in a 3x3 array. Contains
    9 plots total.

    :param title: Title for overall plot (master)
    :return: Figure object for the master figure and numpy array containing a
    list of 3x3 axes subplots.
    """
    # initializes figure and axes object
    fig, ax = plt.subplots(3, 3, squeeze=False)

    # Set label and horizontal axis limits (are the same for all subplots)
    plt.suptitle(title, fontsize=15)
    plt.xlim(-5, 50)

    # Set X and Y Axes Labels and style
    fig.text(x=.45, y=.04, s='Timepoints',  size=10, weight='semibold')
    fig.text(x=.04, y=.6, s='Log2 Fold Change', size=10, rotation='vertical',
             weight='semibold')

    return fig, ax


def __calc_regression(x, y, subplot):
    """
    Calculate second degree regression given a list of values and plots them
    directly on to the subplot specified.

    :param x: Independent values for regression calculation
    :param y: Dependent values (mapping directly to the independent values
    given in x)
    :param subplot: Subplot associated with values to be analyzed
    :return: Plots regression line for the specific subplot
    """
    fit = np.polyfit(x=x, y=y, deg=2)

    # FOrmula for second degree regression line plotting
    subplot.plot(x, fit[0] * x ** 2 + fit[1] * x ** 1 + fit[2], lw=0.5, c='b')
