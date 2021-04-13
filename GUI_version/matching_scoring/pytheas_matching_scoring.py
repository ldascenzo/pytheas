#!/usr/bin/python3

"""
Last update: April 2021
Author: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzoluigi@gmail.com
GitHub project repository: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Pytheas matching and scoring algorithm. Alongside the matching output, FDR estimation file and lists of target/decoys
are generated.
Additional information on the output files and the parameters can be found in the Matching&Scoring section of the
Pytheas manual
"""

from gooey import Gooey, GooeyParser
from match_library import Match


@Gooey(dump_build_config=True, program_name='Pytheas matching and scoring', default_size=(1920, 1080))
def matching_scoring():
    description = 'Matching and scoring workflow'
    parser = GooeyParser(description=description)

    # Required Arguments
    parser.add_argument('Theoretical_digest', help='Input file obtained from the Pytheas in silico '
                                                   'digestion workflow', widget="FileChooser")
    parser.add_argument('Input_mgf', help='Experimental measured peaks in mgf file format', widget="FileChooser")

    # Optional Arguments
    parser.add_argument('--isotopic_species', default='all', choices=['light', 'heavy', 'all'],
                        help='Isotopically labeled (heavy) or unlabeled (light) sequences to include in the matching. '
                             'By default, both are included')
    parser.add_argument('--MS2_ions_minimum_intensity', default='None',
                        help='Minimum absolute intensity threshold of MS2 ions from the mgf file to be included in '
                             'the matching')
    parser.add_argument('--beta', default=0.075, type=float,
                        help='beta parameter value in the scoring function')
    parser.add_argument('--alpha', default=0, type=float,
                        help='alpha parameter value in the scoring function')
    parser.add_argument('--MS2_most_intense_ions', default='all',
                        help='Number of the most intense MS2 ions of the same precursor from the mgf file'
                             'to be included in the matching')
    parser.add_argument('--MS2_normalized_intensity_cutoff', default=5, type=int,
                        help='Minimum intensity threshold normalized to the most intense sequence-defining ion to be'
                             'included in the matching')
    parser.add_argument('--precursor_window_removal', default=2.0, type=float,
                        help='Exclusion window in Da centered around the precursor ion for matching of MS2 ions')
    parser.add_argument('--losses_window_removal', default=1.5, type=float,
                        help='Exclusion window in Da centered around losses (M-xx and free bases ions) for matching of '
                             'MS2 ions')
    parser.add_argument('--MS1_ppm', default=30, type=float,
                        help='Matching tolerance window (ppm) for precursor ions matching')
    parser.add_argument('--MS2_ppm', default=50, type=float,
                        help='Matching tolerance window (ppm) for MS2 ions matching')
    parser.add_argument('--MS1_offset', default=0, type=float,
                        help='Offset (ppm) to apply to all instances of precursor ions matching')
    parser.add_argument('--MS2_offset', default=0, type=float,
                        help='Offset (ppm) to apply to all instances of MS2 ions matching')
    parser.add_argument('--sequence_lengths_FDR', default='all',
                        help='Sequence length values to use for FDR estimation. Insert the values separated by a comma'
                             '(default = all)')
    parser.add_argument('--FDR_isotopic_species', default='all', choices=['all', 'light', 'heavy'],
                        help='Isotopically labeled (heavy) or unlabeled (light) sequences to consider for FDR '
                             'estimation. By default, both are included')
    parser.add_argument('--MS1_mz_minimum', default=400, type=int,
                        help='Lower end of the matching window for precursor ions in m/z')
    parser.add_argument('--MS1_mz_maximum', default=2000, type=int,
                        help='Higher end of the matching window for precursor ions in m/z')
    parser.add_argument('--MS2_mz_minimum', default=300, type=int,
                        help='Lower end of the matching window (m/z) for MS2 ions')
    parser.add_argument('--MS2_mz_maximum', default=2000, type=int,
                        help='Higher end of the matching window (m/z) for MS2 ions')
    parser.add_argument('--precursor_isotopologues', action='store_true', default=False,
                        help='Include the +-1 isotopologue peaks of precursor ions for matching')
    parser.add_argument('--ignore_charges_mgf', action='store_true', default=False,
                        help='Do not use charge information from the mgf to select precursor ions for matching')
    parser.add_argument('--targets_without_decoys', action='store_true', default=False,
                        help='Use targets without competing decoys to estimate FDR')

    ####################################################
    args = parser.parse_args()

    match = Match(args.Theoretical_digest, args.Input_mgf, args.isotopic_species, args.MS1_mz_minimum,
                  args.MS1_mz_maximum, args.MS2_mz_minimum, args.MS2_mz_maximum, args.MS1_ppm, args.MS2_ppm,
                  args.MS1_offset, args.MS2_offset, args.sequence_lengths_FDR,
                  args.MS2_ions_minimum_intensity, args.MS2_most_intense_ions,
                  args.precursor_window_removal, args.losses_window_removal, args.beta, args.alpha,
                  args.MS2_normalized_intensity_cutoff, args.precursor_isotopologues, args.ignore_charges_mgf,
                  args.FDR_isotopic_species, args.targets_without_decoys)

    match.final_output()


if __name__ == '__main__':
    matching_scoring()
