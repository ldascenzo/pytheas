#!/usr/bin/python3

"""
Last update: December 2020
Credits to: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzo@scripps.edu
GitHub project: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Pytheas matching and scoring, GUI version
"""

from gooey import Gooey, GooeyParser
from pytheas_match_gui import Match
from datetime import datetime


@Gooey(dump_build_config=True, program_name='Pytheas matching and scoring', default_size=(1920, 1080))
def matching_scoring():
    description = 'Pytheas matching and scoring routine'
    parser = GooeyParser(description=description)

    # Required Arguments
    parser.add_argument('Input_Digest', help='MS2 theoretical digest file (from Pytheas digestion workflow)',
                        widget="FileChooser")
    parser.add_argument('Input_mgf', help='experimental measured peaks in .mgf format file', widget="FileChooser")

    # Optional Arguments
    parser.add_argument('--light_heavy', default='all', choices=['light', 'heavy', 'all'],
                        help='Isotopic species to use for matching/scoring (default=all)')
    parser.add_argument('--MS2_ions_minimum_absolute_intensity', default='None',
                        help='Minimum absolute intensity (in the mgf file) of MS2 ions to be selected '
                             'for matching (default=None)')
    parser.add_argument('--MS2_most_intense_ions', default='all',
                        help='Number of the most intense MS2 ions of the same precursor (from the mgf file)'
                             'to be selected for matching (default=all)')
    parser.add_argument('--precursor_window_removal', default=2.0, type=float,
                        help='Exclusion MS2 matching m/z window (+-Da) around the precursor ion (default=2.0)')
    parser.add_argument('--losses_window_removal', default=1.5, type=float,
                        help='Exclusion MS2 matching m/z window (+-Da) around M-xx losses and free bases (default=1.5)')
    parser.add_argument('--beta_increment', default=0.075, type=float,
                        help='Increment value for consecutive sequence ions matches in the scoring function '
                             '(default=0.075)')
    parser.add_argument('--alpha', default=0, type=float,
                        help='Incremental step value for each consecutive sequence ion match within a series '
                             '(default=0)')
    parser.add_argument('--MS2_normalized_intensity_cutoff', default=5, type=int,
                        help='Cutoff for MS2 ions intensity NORMALIZED to top sequence ion match (default=5)')
    parser.add_argument('--MS1_mz_minimum', default=400, type=int,
                        help='Minimum m/z for MS1 matching window (default=400)')
    parser.add_argument('--MS1_mz_maximum', default=2000, type=int,
                        help='Maximum m/z for MS1 matching window (default=2000)')
    parser.add_argument('--MS2_mz_minimum', default=300, type=int,
                        help='Minimum m/z for MS2 matching window (default=300)')
    parser.add_argument('--MS2_mz_maximum', default=2000, type=int,
                        help='Maximum m/z for  MS2 matching window (default = 2000)')
    parser.add_argument('--MS1_ppm', default=30, type=int,
                        help='MS1 ions matching ppm window (default=+-30)')
    parser.add_argument('--MS2_ppm', default=50, type=int,
                        help='MS1 ions matching ppm window (default = +-50)')
    parser.add_argument('--MS1_ppm_offset', default=0, type=int,
                        help='ppm offset to "center" all the MS1 ions measured m/z (default=0)')
    parser.add_argument('--MS2_ppm_offset', default=0, type=int,
                        help='ppm offset to "center" all the MS2 ions measured m/z (default = 0)')
    parser.add_argument('--FDR_light_heavy', default='all', choices=['all', 'light', 'heavy'],
                        help='Select if only light (light), only heavy (heavy) isotopes or all the matching ions '
                             'have to be used to calculate FDR values and in the targets/decoys csv output '
                             '(Optional, default=all)')
    parser.add_argument('--precursor_isotopologues', action='store_true', default=False,
                        help='Include the +-1 isotopologue peaks of precursor ions for matching')
    parser.add_argument('--ignore_charges_mgf', action='store_true', default=False,
                        help='Ignore the charges in the mgf input file to select precursor ions for matching')
    parser.add_argument('--targets_without_decoys', action='store_true', default=False,
                        help='Use also targets without any competing decoys for the FDR table')

    ####################################################
    args = parser.parse_args()

    match = Match(args.Input_Digest, args.Input_mgf, args.light_heavy, args.MS1_mz_minimum, args.MS1_mz_maximum,
                  args.MS2_mz_minimum, args.MS2_mz_maximum, args.MS1_ppm, args.MS2_ppm, args.MS1_ppm_offset,
                  args.MS2_ppm_offset, args.MS2_ions_minimum_absolute_intensity, args.MS2_most_intense_ions,
                  args.precursor_window_removal, args.losses_window_removal, args.beta_increment, args.alpha,
                  args.MS2_normalized_intensity_cutoff, args.precursor_isotopologues, args.ignore_charges_mgf,
                  args.FDR_light_heavy, args.targets_without_decoys)

    match.final_output()


if __name__ == '__main__':
    matching_scoring()
