#!/usr/bin/python3

"""
Last update: December 2020
Credits to: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzo@scripps.edu
GitHub project: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Pytheas visualization, GUI version

DESCRIPTION
Visualization script
"""

from gooey import Gooey, GooeyParser
from pytheas_visualization_gui import Visualize


@Gooey(dump_build_config=True, program_name='Pytheas spectra visualization of annotated spectra',
       default_size=(1920, 1080))
def visualizer():
    description = 'Elaboration annotated spectra plots'
    parser = GooeyParser(description=description)

    # Required arguments
    parser.add_argument('Digest_input', widget="FileChooser", help='Digest file used for matching')
    parser.add_argument('Input_mgf', widget="FileChooser", help='.mgf input file to visualize')
    parser.add_argument('Match_file', widget="FileChooser", help='Output file of the matching routine (match_output_'
                                                                 'xx.txt)')

    # Optional arguments
    parser.add_argument('--visualize_decoys', action='store_true', default=False,
                        help='Visualize decoys spectra')
    parser.add_argument('--modified_only', action='store_true', default=False,
                        help='Visualize only spectra for matches with nucleotide modifications')
    parser.add_argument('--only_unique_positions', action='store_true', default=False,
                        help='Add only matches to unique positions in the final output')
    parser.add_argument('--remove_redundant_sequences_with_X', action='store_true', default=False,
                        help='Remove redundant matches containing X, keeping only the highest ranking')
    parser.add_argument('--rank_max', default=99, type=int,
                        help='Maximum rank of spectra to visualize (default=99[all])')
    parser.add_argument('--Sp_minimum_cutoff', default=0, type=float,
                        help='Value cutoff of the Sp score, all matched spectra Sp >= cutoff will be included in '
                             'the output (DEFAULT = 0)')
    parser.add_argument('--dSp_maximum_cutoff', default=1, type=float,
                        help='Value cutoff of the dSp parameter, all matched spectra dSp <= cutoff will be included in '
                             'the output (DEFAULT = 1)')
    parser.add_argument('--highest_peaks', default='all',
                        help='Number of most intense MS2+ peaks to be shown in the output spectrum (default=all)')
    parser.add_argument('--x_min', default=-1, type=int, help='Minimum value to use for x axis on the spectra plotting'
                                                              '(default=-1[auto])')
    parser.add_argument('--x_max', default=-1, type=int, help='Maximum value to use for x axis on the spectra plotting'
                                                              '(default=-1[auto])')
    parser.add_argument('--digest_spectra', action='store_true', default=False,
                        help='Visualize all the digest spectra')


    ####################################################
    args = parser.parse_args()

    visualization = Visualize(args.Digest_input, args.Input_mgf, args.Match_file, args.highest_peaks,
                              args.x_min, args.x_max, args.digest_spectra, args.Sp_minimum_cutoff,
                              args.dSp_maximum_cutoff, args.visualize_decoys, args.modified_only, args.rank_max,
                              args.only_unique_positions, args.remove_redundant_sequences_with_X)

    visualization.final_output()


if __name__ == '__main__':
    visualizer()
