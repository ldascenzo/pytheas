#!/usr/bin/env python3

"""
Last update: August 2020
Credits to: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzo@scripps.edu

DESCRIPTION
The script creates a series of plots based on the data from the matching, useful for statistical evaluation and
assessment of the general data quality

USAGE
python pytheas_statistics.py --OPTIONS

OPTIONS
--targets_csv (Required) -> Input file with targets information, in csv format
--decoys_csv  (Required) -> Input file with decoys information, in csv format
--match_file (Required) -> Input file with the matching/scoring output, in .txt format
--light_heavy (Optional) -> Choose if only light matches, heavy matches or all have to be used for statistics
--FDR_input (Optional, if FDR is requested on scatterplots) -> Input file with information on the FDR for the input
  dataset, in csv format
--y_min (Optional, default = auto) -> Minimum value for the y axis for visualization on scatter plots with Sp vs nts
--y_max (Optional, default = auto) -> Maximum value for the y axis for visualization on scatter plots with Sp vs nts
--lengths (Optional, default = all - NOTE: use without '=' and with no comma between values) -> Input the nucleotide
  lenght for the fragments to be used for statistical analysis. Option is not applied for the ppm offset histogram
--percentiles (Optional, default = n) -> Select (y/n) if percentiles have to be calculated, output as csv file and
  shown in the scatter plot
--FDR_line_values (Optional, default = None - NOTE: use without '=' and with no comma between values) -> Choose the
  FDR (%) lines to be shown in the scatterplots. Lines are shown only if at least 100 target matches are identified
--Sp_cutoff (Optional, default = 0) -> Select an Sp score cutoff for which only target matches above the threshold
  will be in the graphs output and to edit the FDR table
--only_targets_with_decoys (Optional, default = y) -> Select (y/n) if only targets with at least one decoy will be
  used for the scatterplots and boxplots
--info_box (Optional, default = y) -> Select (y/n) if boxes with info on the plots (statistics, parameters...) have
  to be shown
--MS1_offset_min (Optional, default = None) -> Minimum value for the MS1 offset as all the matches below this value are
  excluded for statistics - NOTE: not an absolute value, usually a negative value
--MS1_offset_max (Optional, default = None) -> Maximum value for the MS1 offset as all the matches above this value are
  excluded for statistics - usually a positive value


OUTPUT & HELP
Following graphs in .png format:
- Scatter_Sp_Spdecoys -> Scatterplot with the distribution of Sp scores for the top scoring (per each unique match)
  targets and decoys
- Scatter_multiplot
- Barchart highest scoring target/decoy couples per unique precursor ion matching
- Scatter_dSp_Sp -> Scatterplot with the distribution of dSp scores vs Sp for top scoring (per each unique match)
  targets and decoys
- Hist_Sp -> Double histogram with the distribution of Sp for all targets and decoys (top) and only targets and
  decoys with dSp = 0 (bottom)
- Hist_dSp_2nd_best -> Histogram with the dSp values of the second best target matches for targets with dSp = 0
- Hist_ppm_offset -> Double histogram with the distribution of ppm offset errors for MS1 matches (top) and MS2 matches
  (bottom)
- Scatterbox_Sp_nts -> Scatterplot (top) and boxplot (bottom) with the distribution of Sp scores vs the lengths of
  matches for targets with their respective decoys
- Scatterbox_Sp_nts_z -> Scatterplot (top) and boxplot (bottom) with the distribution of Sp scores vs the lengths of
  matches (divided by different charges) for targets with their respective decoys
"""

import os, sys
import pandas as pd
import argparse
import stats_tools as stats
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

# Initialize and define launch options
parser = argparse.ArgumentParser(description='List of available options')
parser.add_argument('--targets_csv', required=True,
                    help='Input file with targets information, in csv format (Required)')
parser.add_argument('--decoys_csv', required=True, help='Input file with decoys information, in csv format (Required)')
parser.add_argument('--match_file', required=True, help='Input file obtained from matching/scoring script (Required)')
parser.add_argument('--light_heavy', default='all', choices=['all', 'light', 'heavy'],
                    help='Choose if only light, only heavy or all matches will be used for statistics (Optional)')
parser.add_argument('--FDR_input', default=None,
                    help='Input file with information on the FDR for the input dataset, in csv format (Optional)')
parser.add_argument('--y_min', default=0, type=np.float64,
                    help='Minimum value for the y axis for visualization on scatter plots with Sp vs nts '
                         '(Optional, default = auto)')
parser.add_argument('--y_max', default=0, type=np.float64,
                    help='Maximum value for the y axis for visualization on scatter plots with Sp vs nts '
                         '(Optional, default = auto)')
parser.add_argument('--lengths', nargs='*', default=["all"], type=int,
                    help='Input the nucleotide lenght for the fragments to be used for statistical analysis. '
                         'Option is not used for ppm offset histogram (Optional, default = all)')
parser.add_argument('--percentiles', default='n', choices=['y', 'n'],
                    help='Select if percentiles have to be calculated, output as csv file and shown in the '
                         'scatter plot (Optional, default = n)')
parser.add_argument('--FDR_line_values', nargs='*', default=0, type=np.float64,
                    help='Choose the FDR (%) lines to be shown in the scatterplots. Lines are shown only if at '
                         'least 100 target matches are identified (Optional, default = none)')
parser.add_argument('--Sp_cutoff', default=0, type=np.float64,
                    help='Select an Sp score cutoff for which only target matches above the threshold will be '
                         'in the graphs output and to edit the FDR table (Optional, default = 0')
parser.add_argument('--only_targets_with_decoys', default='y', choices=['y', 'n'],
                    help='Select (y/n) if only targets with at least one decoy will be used for the scatterplots '
                         'and boxplots (Optional, default = y)')
parser.add_argument('--info_box', default='y', choices=['y', 'n'],
                    help='Select (y/n) if boxes with info on the plots (statistics, parameters...) have '
                         'to be shown (Optional, default = y)')
parser.add_argument('--MS1_offset_min', default=None, type=int,
                    help='Minimum value for the MS1 offset as all the matches below this value are excluded '
                         'for statistics - NOTE: not an absolute value, usually a negative value  '
                         '(Optional, default = None')
parser.add_argument('--MS1_offset_max', default=None, type=int,
                    help='Maximum value for the MS1 offset as all the matches above this value are excluded for '
                         'statistics - usually a positive value (Optional, default = None')

args = parser.parse_args()

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


def graphs_plot():
    """
    Plot the graphs from the input csv data
    """
    csv_files = (read_csv_input(args.targets_csv, args.light_heavy), read_csv_input(args.decoys_csv, args.light_heavy))

    if csv_files[0].empty:
        print("ERROR! CSV input file {} is empty. Execution terminated without generating any output".format(
            args.targets_csv))
        sys.exit(1)

    if csv_files[1].empty:
        print("ERROR! CSV input file {} is empty. Execution terminated without generating any output".format(
            args.decoys_csv))
        sys.exit(1)

    if args.FDR_input:
        name_input = args.FDR_input[4:-4]

    else:
        name_input = "output"

    # Scatter plot of Sp of targets vs Sp of decoys of the given targets 
    stats.scatter_Sp_vs_Spdecoy(csv_files, args.lengths, args.Sp_cutoff, args.MS1_offset_min, args.MS1_offset_max,
                                args.info_box)
    save_graph("Scatter_Sp_Spdecoy_{}.png".format(name_input))
    plt.close()

    # Scatter plot of Sp of targets vs Sp of decoys of the given targets, divided by nucleotide lengths of precursors
    stats.scatter_Sp_vs_Spdecoy(csv_files, ['analysis'], args.Sp_cutoff, args.MS1_offset_min, args.MS1_offset_max,
                                args.info_box)
    save_graph("Scatter_Sp_Spdecoy_bylength_{}.png".format(name_input))
    plt.close()

    # TEST
    stats.barchart_global_normalized_scores(csv_files, args.lengths, args.Sp_cutoff, args.MS1_offset_min,
                                            args.MS1_offset_max, args.info_box)
    save_graph("Barchart_NormalizedSp_targets&decoys_{}.png".format(name_input))
    plt.close()

    # Histogram with the distribution of the Sp scores for targets/decoys and top target/decoys
    plt.subplot(2, 1, 1)
    stats.hist_Sp(csv_files, args.lengths, args.Sp_cutoff, args.only_targets_with_decoys, args.MS1_offset_min,
                  args.MS1_offset_max, args.info_box)

    plt.subplot(2, 1, 2)
    stats.hist_top_Sp(csv_files, args.lengths, args.Sp_cutoff, args.only_targets_with_decoys, args.MS1_offset_min,
                      args.MS1_offset_max, args.info_box)
    plt.subplots_adjust(hspace=0.25)
    save_graph("Hist_Sp_{}.png".format(name_input))
    plt.close()

    # Histogram with the distribution of the dSp values for the second best targets
    stats.hist_second_dSp(csv_files, args.lengths, args.Sp_cutoff, args.only_targets_with_decoys, args.MS1_offset_min,
                          args.MS1_offset_max, args.info_box)
    save_graph("Hist_dSp_2ndbest_{}.png".format(name_input))
    plt.close()

    # Histogram with the distribution of ppm error offsets for MS1 and MS2 matches
    plt.subplot(2, 1, 1)
    stats.ppm_errors_histogram(csv_files, args.lengths, args.Sp_cutoff, 'MS1', args.match_file, args.MS1_offset_min,
                               args.MS1_offset_max, args.info_box)

    plt.subplot(2, 1, 2)
    stats.ppm_errors_histogram(csv_files, args.lengths, args.Sp_cutoff, 'MS2', args.match_file, args.MS1_offset_min,
                               args.MS1_offset_max, args.info_box)
    save_graph("Hist_ppm_offset_{}.png".format(name_input))
    plt.close()

    # Scatter with the dSp scores vs Sp scores for targets and decoys. FDR and percentile lines are optional
    stats.scatter_dSp_vs_Sp(csv_files, args.lengths, args.percentiles, args.FDR_line_values, args.Sp_cutoff, name_input,
                            args.FDR_input, args.info_box, args.only_targets_with_decoys, args.MS1_offset_min,
                            args.MS1_offset_max)
    save_graph("Scatter_dSp_Sp_{}.png".format(name_input))
    plt.close()

    # Scatter + boxplot of the Sp scores vs the matches ions length    
    plt.subplot(2, 1, 1)
    stats.scatter_nts_vs_score(csv_files, 'Sp', args.y_min, args.y_max, args.lengths, args.percentiles,
                               args.FDR_line_values, args.Sp_cutoff, name_input, args.FDR_input, args.info_box,
                               args.only_targets_with_decoys, args.MS1_offset_min, args.MS1_offset_max)
    plt.subplot(2, 1, 2)
    stats.box_nts_vs_score(csv_files, 'Sp', args.y_min, args.y_max, args.lengths, args.Sp_cutoff,
                           args.only_targets_with_decoys, args.MS1_offset_min, args.MS1_offset_max)

    save_graph("Scatterbox_Sp_nts_{}.png".format(name_input))
    plt.close()

    # Scatter + boxplot of the Sp scores vs the matches ions length grouped by charge
    try:
        plt.subplot(2, 1, 1)
        stats.scatter_nts_z_vs_score(csv_files, 'Sp', args.y_min, args.y_max, args.lengths, args.Sp_cutoff,
                                     args.MS1_offset_min, args.MS1_offset_max, args.only_targets_with_decoys)
        plt.subplot(2, 1, 2)
        stats.box_nts_z_vs_score(csv_files, 'Sp', args.y_min, args.y_max, args.lengths, args.Sp_cutoff, args.MS1_offset_min,
                                 args.MS1_offset_max, args.only_targets_with_decoys)
        save_graph("Scatterbox_Sp_nts_z_{}.png".format(name_input))
        plt.close()

    except:
        print('Scatter + boxplots of Sp scores of matches vs length grouped by charge could not be produced. '
              'This is probably due to the parameters chosen')

    # Edit the FDR table with the option used in the statistics
    stats.FDR_update(csv_files, args.lengths, args.Sp_cutoff, name_input, args.MS1_offset_min, args.MS1_offset_max,
                     args.only_targets_with_decoys)


if __name__ == "__main__":
    if args.FDR_line_values and not args.FDR_input:
        print(
            "ERROR! Please specify an FDR input file to output the FDR lines in the graphs. "
            "Execution terminated without generating any output")
        sys.exit(1)

    if args.FDR_input and not args.FDR_line_values:
        print(
            "WARNING! Input FDR file was selected but no output lines values "
            "(via option --FDR_line_values) was specified")

    graphs_plot()
