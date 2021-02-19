#!/usr/bin/python3

"""
Last update: December 2020
Credits to: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzo@scripps.edu
GitHub project: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Pytheas statistics, GUI version

DESCRIPTION
The script creates a series of plots based on the data from the matching, useful for statistical evaluation and
assessment of the general data quality

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

from gooey import Gooey, GooeyParser
from pytheas_statistics_gui import Stats
import numpy as np


@Gooey(dump_build_config=True, program_name='Pytheas descriptive statistical plots', default_size=(1920, 1080))
def statistical_plots():
    description = 'Elaboration of Pytheas descriptive statistical plots'
    parser = GooeyParser(description=description)

    # Required arguments
    parser.add_argument('Targets', widget="FileChooser",
                        help='Targets .csv file obtained after matching')
    parser.add_argument('Decoys', widget="FileChooser",
                        help='Decoys .csv file obtained after matching')
    parser.add_argument('Match_output', widget="FileChooser",
                        help="Match output file (match_output_xxx.txt")

    # Optional Arguments
    parser.add_argument('--light_heavy', default='all', choices=['all', 'light', 'heavy'],
                        help='Use light, heavy or all isotope matches for statistics (default=all)')
    parser.add_argument('--Sp_cutoff', default=0, type=float,
                        help='Sp minimum score cutoff for targets to be considered for plots')
    parser.add_argument('--FDR_input', default=None, widget="FileChooser",
                        help='FDR .csv file obtained after matching')
    parser.add_argument('--y_min', default=0, type=float,
                        help='Minimum value for the y axis for visualization on scatter plots with Sp vs nts '
                             '(default=auto)')
    parser.add_argument('--y_max', default=0, type=float,
                        help='Maximum value for the y axis for visualization on scatter plots with Sp vs nts '
                             '(default=auto)')
    parser.add_argument('--lengths', nargs='*', default="all",
                        help='Nucleotide lengths for matches to be used for statistical analysis, separated by commas '
                             '(default = all)')
    parser.add_argument('--FDR_line_value', default=0, type=float,
                        help='FDR (%) threshold to show as line in the scatter plots')
    parser.add_argument('--MS1_offset_min', default=None, type=int,
                        help='Exclude all the matches with MS1 offset below the given value (usually <0)')
    parser.add_argument('--MS1_offset_max', default=None, type=int,
                        help='Exclude all the matches with MS1 offset above the given value')
    parser.add_argument('--targets_without_decoys', action='store_true', default=False,
                        help='Use targets without any competing decoy for scatter plots and box plots')
    parser.add_argument('--hide_box_graphs', action='store_true', default=False,
                        help='Hide boxes with info on the plots (statistics, parameters...)')
    parser.add_argument('--percentiles', action='store_true', default=False,
                        help='Add percentiles in the output and in the scatter plots')

    ####################################################
    args = parser.parse_args()

    statistics = Stats(args.Targets, args.Decoys, args.Match_output, args.light_heavy, args.FDR_input, args.y_min,
                       args.y_max, args.lengths, args.percentiles, args.FDR_line_value, args.Sp_cutoff,
                       args.targets_without_decoys, args.hide_box_graphs, args.MS1_offset_min, args.MS1_offset_max)

    statistics.final_output()


if __name__ == '__main__':
    statistical_plots()
