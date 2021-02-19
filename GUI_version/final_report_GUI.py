"""

XXX

"""

from gooey import Gooey, GooeyParser
from pytheas_final_report_gui import Pytheas_Output


@Gooey(dump_build_config=True, program_name='Pytheas final report generator', default_size=(1920, 1080))
def pytheas_output():
    description = 'Pytheas final report generator'
    parser = GooeyParser(description=description)

    # Required arguments
    parser.add_argument('Match_file', widget="FileChooser",
                        help="Match file (match_output_xxx.txt)")

    # Optional Arguments
    parser.add_argument('--visualize_decoys', action='store_true', default=False,
                        help='Add decoys in the final output')
    parser.add_argument('--modified_only', action='store_true', default=False,
                        help='Add only matches with nucleotide modifications in the final output')
    parser.add_argument('--only_unique_positions', action='store_true', default=False,
                        help='Add only matches to unique positions in the final output')
    parser.add_argument('--remove_redundant_sequences_with_X', action='store_true', default=False,
                        help='Remove redundant matches containing X, keeping only the highest ranking')
    parser.add_argument('--Sp_minimum_cutoff', default=0, type=float,
                        help='Value cutoff of the Sp score, all matched spectra Sp >= cutoff will be included in '
                             'the output (DEFAULT = 0)')
    parser.add_argument('--dSp_maximum_cutoff', default=1, type=float,
                        help='Value cutoff of the dSp parameter, all matched spectra dSp <= cutoff will be included in '
                             'the output (DEFAULT = 1)')
    parser.add_argument('--rank_max', default=99, type=int,
                        help='Maximum rank of matches to include in the output (default=99[all])')


    ####################################################
    args = parser.parse_args()

    output = Pytheas_Output(args.Match_file, args.Sp_minimum_cutoff, args.dSp_maximum_cutoff, args.visualize_decoys,
                            args.modified_only, args.only_unique_positions, args.rank_max,
                            args.remove_redundant_sequences_with_X)

    output.parse_match_file()


if __name__ == '__main__':
    pytheas_output()
