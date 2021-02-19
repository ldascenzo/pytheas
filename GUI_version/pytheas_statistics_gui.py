#!/usr/bin/env python3

"""
Last update: December 2020
Credits to: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzo@scripps.edu
"""

import os, sys
import pandas as pd
import stats_tools as stats
import matplotlib.pyplot as plt
import matplotlib as mpl
import ntpath

# Publication-quality settings for the graphs
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
    def __init__(self, targets, decoys, match, isotope, FDR, y_min, y_max, lengths, percentile, FDR_line, Sp_cutoff,
                 targets_without_decoys, box, MS1_offset_cutoff_min, MS1_offset_cutoff_max):
        self.targets_csv, self.decoys_csv, self.match_file, self.light_heavy = targets, decoys, match, isotope
        self.FDR_input, self.y_min, self.y_max, self.lengths = FDR, y_min, y_max, lengths
        self.percentiles, self.FDR_line_values, self.Sp_cutoff = percentile, FDR_line, Sp_cutoff
        self.targets_without_decoys, self.hide_info_box = targets_without_decoys, box
        self.MS1_offset_min, self.MS1_offset_max = MS1_offset_cutoff_min, MS1_offset_cutoff_max

    def graphs_plot(self):
        """
        Plot the graphs from the input csv data
        """
        csv_files = (read_csv_input(self.targets_csv, self.light_heavy),
                     read_csv_input(self.decoys_csv, self.light_heavy))

        if csv_files[0].empty:
            print("ERROR! CSV input file {} is empty. Execution terminated without generating any output".format(
                self.targets_csv))
            sys.exit(1)

        if csv_files[1].empty:
            print("ERROR! CSV input file {} is empty. Execution terminated without generating any output".format(
                self.decoys_csv))
            sys.exit(1)

        if self.FDR_input:
            name_input = filename_from_path(self.FDR_input)[:-4]

        else:
            name_input = "output"

        # Some dirty tricks to incorporate options without changing the stats_tools from the command line version
        if self.percentiles:
            percentile_option = 'y'
        else:
            percentile_option = 'n'

        if self.targets_without_decoys:
            targets_decoys = 'n'
        else:
            targets_decoys = 'y'

        if self.hide_info_box:
            bbox = 'n'
        else:
            bbox = 'y'

        if 'all' in self.lengths:
            lengths = self.lengths
        else:
            lengths = self.lengths[0].split(',')

        # Scatter plot of Sp of targets vs Sp of decoys of the given targets
        stats.scatter_Sp_vs_Spdecoy(csv_files, lengths, self.Sp_cutoff, self.MS1_offset_min, self.MS1_offset_max,
                                    bbox)
        save_graph("Scatter_Sp_Spdecoy_{}.png".format(name_input))
        plt.close()

        # Scatter plot of Sp of targets vs Sp of decoys of the given targets, divided by nucleotide lengths of
        # precursors
        stats.scatter_Sp_vs_Spdecoy(csv_files, ['analysis'], self.Sp_cutoff, self.MS1_offset_min, self.MS1_offset_max,
                                    bbox)
        save_graph("Scatter_Sp_Spdecoy_bylength_{}.png".format(name_input))
        plt.close()

        # TEST
        stats.barchart_global_normalized_scores(csv_files, lengths, self.Sp_cutoff, self.MS1_offset_min,
                                                self.MS1_offset_max, bbox)
        save_graph("Barchart_NormalizedSp_targets&decoys_{}.png".format(name_input))
        plt.close()

        # Histogram with the distribution of the Sp scores for targets/decoys and top target/decoys
        plt.subplot(2, 1, 1)
        stats.hist_Sp(csv_files, lengths, self.Sp_cutoff, targets_decoys, self.MS1_offset_min,
                      self.MS1_offset_max, bbox)

        plt.subplot(2, 1, 2)
        stats.hist_top_Sp(csv_files, lengths, self.Sp_cutoff, targets_decoys, self.MS1_offset_min,
                          self.MS1_offset_max, bbox)
        plt.subplots_adjust(hspace=0.25)
        save_graph("Hist_Sp_{}.png".format(name_input))
        plt.close()

        # Histogram with the distribution of the dSp values for the second best targets
        stats.hist_second_dSp(csv_files, lengths, self.Sp_cutoff, targets_decoys,
                              self.MS1_offset_min, self.MS1_offset_max, bbox)
        save_graph("Hist_dSp_2ndbest_{}.png".format(name_input))
        plt.close()

        # Histogram with the distribution of ppm error offsets for MS1 and MS2 matches
        plt.subplot(2, 1, 1)
        stats.ppm_errors_histogram(csv_files, lengths, self.Sp_cutoff, 'MS1', self.match_file, self.MS1_offset_min,
                                   self.MS1_offset_max, bbox)

        plt.subplot(2, 1, 2)
        stats.ppm_errors_histogram(csv_files, lengths, self.Sp_cutoff, 'MS2', self.match_file, self.MS1_offset_min,
                                   self.MS1_offset_max, bbox)
        save_graph("Hist_ppm_offset_{}.png".format(name_input))
        plt.close()

        # Scatter with the dSp scores vs Sp scores for targets and decoys. FDR and percentile lines are optional
        stats.scatter_dSp_vs_Sp(csv_files, lengths, percentile_option, self.FDR_line_values, self.Sp_cutoff,
                                name_input, self.FDR_input, bbox, targets_decoys,
                                self.MS1_offset_min, self.MS1_offset_max)
        save_graph("Scatter_dSp_Sp_{}.png".format(name_input))
        plt.close()

        # Scatter + boxplot of the Sp scores vs the matches ions length
        plt.subplot(2, 1, 1)
        stats.scatter_nts_vs_score(csv_files, 'Sp', self.y_min, self.y_max, lengths, percentile_option,
                                   self.FDR_line_values, self.Sp_cutoff, name_input, self.FDR_input, bbox,
                                   targets_decoys, self.MS1_offset_min, self.MS1_offset_max)
        plt.subplot(2, 1, 2)
        stats.box_nts_vs_score(csv_files, 'Sp', self.y_min, self.y_max, lengths, self.Sp_cutoff,
                               targets_decoys, self.MS1_offset_min, self.MS1_offset_max)

        save_graph("Scatterbox_Sp_nts_{}.png".format(name_input))
        plt.close()

        # Scatter + boxplot of the Sp scores vs the matches ions length grouped by charge
        try:
            plt.subplot(2, 1, 1)
            stats.scatter_nts_z_vs_score(csv_files, 'Sp', self.y_min, self.y_max, lengths, self.Sp_cutoff,
                                         self.MS1_offset_min, self.MS1_offset_max, targets_decoys)
            plt.subplot(2, 1, 2)
            stats.box_nts_z_vs_score(csv_files, 'Sp', self.y_min, self.y_max, lengths, self.Sp_cutoff,
                                     self.MS1_offset_min, self.MS1_offset_max, targets_decoys)
            save_graph("Scatterbox_Sp_nts_z_{}.png".format(name_input))
            plt.close()

        except:
            print('Scatter + boxplots of Sp scores of matches vs length grouped by charge could not be produced. '
                  'This is probably due to the parameters chosen')

        # Edit the FDR table with the option used in the statistics
        stats.FDR_update(csv_files, lengths, self.Sp_cutoff, name_input, self.MS1_offset_min, self.MS1_offset_max,
                         targets_decoys)

    def final_output(self):
        if self.FDR_line_values and not self.FDR_input:
            print(
                "ERROR! Please specify an FDR input file to output the FDR lines in the graphs. "
                "Execution terminated without generating any output")
            sys.exit(1)

        if self.FDR_input and not self.FDR_line_values:
            print(
                "WARNING! Input FDR file was selected but no output lines values "
                "(via option --FDR_line_values) was specified")

        self.graphs_plot()


def read_csv_input(csv_file, isotopes):
    """
    Reads the csv input file creating a pandas dataframe
    Uses the matches from the light, heavy channel or both based on the user preferences
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
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)
