#!/usr/bin/env python3

"""
Last update: April 2021
Author: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzoluigi@gmail.com
GitHub project repository: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Library for the generation of statistical analysis plots
"""

import os, sys
import pandas as pd
import stats_tools as stats
import matplotlib.pyplot as plt
import matplotlib as mpl
import ntpath

# Set publication-quality settings for the statistical plots
plt.style.use('ggplot')
plt.rc('font', family='serif')
mpl.rcParams['axes.titlesize'] = 7
mpl.rcParams['axes.titleweight'] = 'bold'
mpl.rcParams['axes.labelsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['font.size'] = 5
mpl.rcParams['axes.facecolor'] = 'whitesmoke'
pd.options.mode.chained_assignment = None


class Stats:
    def __init__(self, targets, decoys, match, isotope, FDR, Sp_cutoff, sequence_lengths,
                 targets_without_decoys, box,
                 scatter_Sp_vs_dSp, scatterbox_lengths_z):
        self.targets_csv, self.decoys_csv, self.match_file, self.light_heavy = targets, decoys, match, isotope
        self.FDR_input, self.lengths = FDR, sequence_lengths
        self.Sp_cutoff = Sp_cutoff
        self.targets_without_decoys, self.hide_info_box = targets_without_decoys, box
        self.scatter_Sp_dSp, self.scatterbox_lengths_z = scatter_Sp_vs_dSp, scatterbox_lengths_z

    def graphs_plot(self):
        """
        Statistical plot generation from the csv-format data
        """
        csv_files = (read_csv_input(self.targets_csv, self.light_heavy),
                     read_csv_input(self.decoys_csv, self.light_heavy))

        if csv_files[0].empty:
            print("ERROR! CSV input file {} is empty. Execution terminated without generating any output".format(
                self.targets_csv))
            sys.exit(1)

        name_input = filename_from_path(self.targets_csv)[8:-4]

        # Some dirty tricks to incorporate options without changing the stats_tools from the command line version
        if self.targets_without_decoys:
            targets_decoys = 'n'
        else:
            targets_decoys = 'y'

        if self.hide_info_box:
            bbox = 'n'
        else:
            bbox = 'y'

        if self.lengths == 'all':
            lengths = 'all'
        else:
            lengths = [int(x) for x in self.lengths.split(',')]

        # Histogram with the distribution of ppm error offsets for MS1 and MS2 matches
        plt.subplot(2, 1, 1)
        stats.ppm_errors_histogram(csv_files, self.Sp_cutoff, 'MS1', self.match_file, bbox)

        plt.subplot(2, 1, 2)
        stats.ppm_errors_histogram(csv_files, self.Sp_cutoff, 'MS2', self.match_file, bbox)
        save_graph("Hist_ppm_offset_{}.png".format(name_input))
        plt.close()

        # Histogram with the distribution of the Sp scores for targets/decoys and top target/decoys
        plt.subplot(2, 1, 1)
        stats.hist_Sp(csv_files, lengths, self.Sp_cutoff, targets_decoys, bbox)

        plt.subplot(2, 1, 2)
        stats.hist_top_Sp(csv_files, lengths, self.Sp_cutoff, targets_decoys, bbox)
        plt.subplots_adjust(hspace=0.25)
        save_graph("Hist_Sp_{}.png".format(name_input))
        plt.close()

        # Histogram with the distribution of dSp2 for top targets with a competing target
        stats.hist_second_dSp(csv_files, lengths, self.Sp_cutoff, targets_decoys, bbox)
        save_graph("Hist_dSp2_{}.png".format(name_input))
        plt.close()

        # # Scatter + boxplot of the Sp scores vs the sequences length
        plt.subplot(2, 1, 1)
        stats.scatter_nts_vs_score(csv_files, 'Sp', 0, 0, lengths, self.Sp_cutoff, bbox, targets_decoys)
        plt.subplot(2, 1, 2)
        stats.box_nts_vs_score(csv_files, 'Sp', 0, 0, lengths, self.Sp_cutoff,
                               targets_decoys)

        save_graph("Scatterbox_Sp_length_{}.png".format(name_input))
        plt.close()

        # Scatter + boxplot of the Sp scores vs the sequence length separated by charge state
        if self.scatterbox_lengths_z:
            try:
                plt.subplot(2, 1, 1)
                stats.scatter_nts_z_vs_score(csv_files, 'Sp', 0, 0, lengths, self.Sp_cutoff, targets_decoys, bbox)
                plt.subplot(2, 1, 2)
                stats.box_nts_z_vs_score(csv_files, 'Sp', 0, 0, lengths, self.Sp_cutoff, targets_decoys)
                save_graph("Scatterbox_Sp_length_charge_{}.png".format(name_input))
                plt.close()

            except:
                print('Scatter + boxplot of Sp scores of matches vs length grouped by charge could '
                      'not be produced. This is probably due to the parameters chosen')

        # Additional graphs are plotted only if the decoys csv file is not empty
        if not csv_files[1].empty:
            # Scatter plot of Sp of targets vs Sp of their highest-score competing decoys
            stats.scatter_Sp_vs_Spdecoy(csv_files, lengths, self.Sp_cutoff, bbox)
            save_graph("Scatter_Sp_Spdecoy_{}.png".format(name_input))
            plt.close()

            # Scatter plot of Sp of targets vs Sp of decoys of the given targets, divided by sequence length
            stats.scatter_Sp_vs_Spdecoy(csv_files, ['analysis'], self.Sp_cutoff, bbox)
            save_graph("Scatter_Sp_Spdecoy_bylength_{}.png".format(name_input))
            plt.close()

            # Scatter plot with the dSp scores vs Sp scores for targets and decoys.
            # FDR and percentile lines are optional
            if self.scatter_Sp_dSp:
                stats.scatter_dSp_vs_Sp(csv_files, lengths, self.Sp_cutoff, bbox, targets_decoys)
                save_graph("Scatter_dSp_Sp_{}.png".format(name_input))
                plt.close()

            # Edit the FDR table if additional targets (without decoys) are used in the plots
            if self.targets_without_decoys:
                stats.FDR_update(csv_files, lengths, self.Sp_cutoff, name_input)

    def final_output(self):
        self.graphs_plot()


def read_csv_input(csv_file, isotopes):
    """
    Read the csv input file creating a pandas dataframe
    Use the matches from the light, heavy channel or both based on the user preferences
    """
    df = pd.read_csv(csv_file)

    # Applies the filter of light/heavy ions specified by the user
    if isotopes == 'light':
        df = df[df['isotope'] == 'light']
    elif isotopes == 'heavy':
        df = df[df['isotope'] == 'heavy']
    else:
        pass

    return df


def save_graph(namefile):
    """
    Save the matplotlib graphs and creates copies of already existing files
    """
    if os.path.exists(namefile):
        numbers = list(range(100))
        for n in numbers:
            if not os.path.exists("_".join(namefile.split("_")[:-1]) + "_{}.png".format(n)):
                plt.savefig("_".join(namefile.split("_")[:-1]) + "_{}.png".format(n), bbox_inches='tight', dpi=300)
                break

    else:
        plt.savefig(namefile, bbox_inches='tight', dpi=300)


def filename_from_path(path):
    """
    Extract the name of the output from the given input file
    """
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)
