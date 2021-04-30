#!/usr/bin/python3

"""
Last update: April 2021
Author: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzoluigi@gmail.com
GitHub project repository: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Pytheas statistical analysis tools. Given information on target and decoy sequences (obtained from the matching and
scoring step), relevant statistical plots are generated. More information on the parameters and the output plots
can be found in the Pytheas manual.
"""

from gooey import Gooey, GooeyParser
from statistics_library import Stats


@Gooey(dump_build_config=True, program_name='Pytheas statistical analysis', default_size=(1920, 1080))
def statistical_plots():
    description = 'Generate statistical plots for quality assessment of the database matching process'
    parser = GooeyParser(description=description)

    # Required arguments
    parser.add_argument('Targets', widget="FileChooser",
                        help='Pytheas matching and scoring output file with a list of top scoring '
                             'target sequences (targets_[dataset].csv)')
    parser.add_argument('Decoys', widget="FileChooser",
                        help='Pytheas matching and scoring output file with a list of top scoring '
                             'decoy sequences (decoys_[dataset].csv)')
    parser.add_argument('Match_output', widget="FileChooser",
                        help="Pytheas matching and scoring .txt output file")

    # Optional Arguments
    parser.add_argument('--isotopic_species', default='all', choices=['light', 'heavy', 'all'],
                        help='Unlabeled (light) or isotopically labeled target/decoy sequences to include in the '
                             'statistical plots. By default, both are included')
    parser.add_argument('--Sp_cutoff', default=0, type=float,
                        help='Minimum Sp score cutoff for target sequences to use for statistical analysis')
    parser.add_argument('--sequence_lengths', default='all',
                        help='Sequence length values to use for statistical analysis. '
                             'Input multiple values separated by commas (default = all)')
    parser.add_argument('--targets_without_decoys', action='store_true', default=False,
                        help='Use all the available targets for statistical analysis and FDR re-calculation. '
                             'By default, only targets with competing decoys are used')
    parser.add_argument('--hide_box_graphs', action='store_true', default=False,
                        help='Hide descriptive statistics and parameters summary box')
    parser.add_argument('--Sp_vs_dSp_plot', action='store_true', default=False,
                        help='Generate an additional Sp vs dSp scatter for target and decoy sequences')
    parser.add_argument('--Sp_length_charge_plot', action='store_true', default=False,
                        help='Generate an additional scatter and box plots of Sp vs sequence_length_charge for'
                             'target and decoy sequences')

    ####################################################
    args = parser.parse_args()

    statistics = Stats(args.Targets, args.Decoys, args.Match_output, args.isotopic_species,
                       args.Sp_cutoff, args.sequence_lengths, args.targets_without_decoys,
                       args.hide_box_graphs, args.Sp_vs_dSp_plot, args.Sp_length_charge_plot)

    statistics.final_output()


if __name__ == '__main__':
    statistical_plots()
