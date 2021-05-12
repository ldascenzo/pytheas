#!/usr/bin/python3

"""
Last update: April 2021
Author: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzoluigi@gmail.com
GitHub project repository: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Pytheas filtering algorithm to generate an easy-to-use and organized output file.
Additional information on the output file and the parameters can be found in the Pytheas manual
"""

from gooey import Gooey, GooeyParser
from final_report_library import Pytheas_Output


@Gooey(dump_build_config=True, program_name='Pytheas final report', default_size=(1920, 1080))
def pytheas_output():
    description = 'Generate a filtered list of identified sequences'
    parser = GooeyParser(description=description)

    # Required arguments
    parser.add_argument('Match_output', widget="FileChooser",
                        help="Pytheas matching and scoring .txt output file")

    # Optional Arguments
    parser.add_argument('--rank_maximum', default=99, type=int,
                        help='Maximum rank value for matched sequences (default=99[all])')
    parser.add_argument('--Sp_minimum', default=0, type=float,
                        help='Minimum Sp score cutoff for matched sequences (default=0)')
    parser.add_argument('--dSp_maximum', default=1, type=float,
                        help='Maximum dSp score cutoff for matched sequences (default=1)')
    parser.add_argument('--dSp2_minimum', default=0, type=float,
                        help='Minimum dSp2 score cutoff for matched sequences (default=0)')
    parser.add_argument('--MS1_ppm_cutoff', default=None, type=float,
                        help='Precursor ion mass tolerance cutoff. Only matches with precursor ion m/z within '
                             'specified +- ppm window are included (default=no cutoff)')
    parser.add_argument('--output_decoys', action='store_true', default=False,
                        help='Output decoy matches. By default, only targets are included')
    parser.add_argument('--modified_only', action='store_true', default=False,
                        help='Only include matches containing nucleotide modifications')
    parser.add_argument('--unique_positions_only', action='store_true', default=False,
                        help='Only include matches mapping to unique positions within input RNA sequence(s)')
    parser.add_argument('--remove_redundant_SeqX_matches', action='store_true', default=False,
                        help='When SeqX is used for sequence consolidation, output only the highest ranking matches '
                             'with X')

    ####################################################
    args = parser.parse_args()

    output = Pytheas_Output(args.Match_output, args.Sp_minimum, args.dSp_maximum, args.output_decoys,
                            args.modified_only, args.unique_positions_only, args.rank_maximum,
                            args.remove_redundant_SeqX_matches, args.dSp2_minimum, args.MS1_ppm_cutoff)

    output.parse_match_file()


if __name__ == '__main__':
    pytheas_output()
