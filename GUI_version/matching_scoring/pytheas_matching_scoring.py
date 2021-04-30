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


@Gooey(dump_build_config=True, program_name="Pytheas matching and scoring", default_size=(1920, 1080))
def matching_scoring():
    description = 'Generate a list of oligonucleotide spectra matches against the theoretical digest library'
    parser = GooeyParser(description=description)

    # Required Arguments
    parser.add_argument('Theoretical_digest', help='Pytheas in silico digestion output file', widget="FileChooser")
    parser.add_argument('MS_data', help='Experimental data peak list in mgf format', widget="FileChooser")

    # Optional Arguments
    parser.add_argument('--isotopic_species', default='all', choices=['light', 'heavy', 'all'],
                        help='Unlabeled (light) or isotopically labeled (heavy) sequences to include in the search. '
                             'By default, both are included')
    parser.add_argument('--MS2_absolute_peak_intensity', default='None',
                        help='Minimum absolute intensity threshold for MS2 peaks')
    parser.add_argument('--MS2_most_intense_peaks', default='all',
                        help='Number of the most intense MS2 peaks to be used for matching')
    parser.add_argument('--MS2_relative_peak_intensity', default=5, type=int,
                        help='Minimum relative intensity threshold for MS2 peaks. Normalization is done '
                             'to the most intense non-precursor-derived peak')
    parser.add_argument('--precursor_exclusion_window', default=2.0, type=float,
                        help='Mass range centered around the precursor ion m/z for MS2 peaks exclusion')
    parser.add_argument('--precursor_losses_exclusion_window', default=1.5, type=float,
                        help='Mass range centered around the precursor ion losses m/z for MS2 peaks exclusion')
    parser.add_argument('--MS1_ppm', default=30, type=float,
                        help='Precursor ion matching tolerance in ppm')
    parser.add_argument('--MS2_ppm', default=50, type=float,
                        help='MS2 ions matching tolerance in ppm')
    parser.add_argument('--MS1_offset', default=0, type=float,
                        help='Precursor ion mass offset correction in ppm')
    parser.add_argument('--MS2_offset', default=0, type=float,
                        help='MS2 ions mass offset correction in ppm')
    parser.add_argument('--beta', default=0.075, type=float,
                        help='beta parameter value in the scoring function')
    parser.add_argument('--alpha', default=0, type=float,
                        help='alpha parameter value in the scoring function')
    parser.add_argument('--FDR_isotopic_species', default='all', choices=['all', 'light', 'heavy'],
                        help='Unlabeled (light) or isotopically labeled (heavy) sequences to use for FDR. '
                             'By default, both are included')
    parser.add_argument('--MS1_mz_minimum', default=400, type=int,
                        help='Minimum m/z value for precursor ion matching')
    parser.add_argument('--MS1_mz_maximum', default=2000, type=int,
                        help='Maximum m/z value for precursor ion matching')
    parser.add_argument('--MS2_mz_minimum', default=300, type=int,
                        help='Minimum m/z value for MS2 ions matching')
    parser.add_argument('--MS2_mz_maximum', default=2000, type=int,
                        help='Maximum m/z value for MS2 ions matching')
    parser.add_argument('--isotopologue_precursor_matching', action='store_true', default=False,
                        help='Include the mass of +-1 isotopologue for precursor ion matching')
    parser.add_argument('--exclude_precursor_charge_matching', action='store_true', default=False,
                        help='Precursor ion charges in mgf are not used for matching')
    parser.add_argument('--targets_without_decoys', action='store_true', default=False,
                        help='Use all available target sequences for FDR estimation. By default, only targets with'
                             ' competing decoys are used')

    ####################################################
    args = parser.parse_args()

    match = Match(args.Theoretical_digest, args.MS_data, args.isotopic_species, args.MS1_mz_minimum,
                  args.MS1_mz_maximum, args.MS2_mz_minimum, args.MS2_mz_maximum, args.MS1_ppm, args.MS2_ppm,
                  args.MS1_offset, args.MS2_offset,
                  args.MS2_absolute_peak_intensity, args.MS2_most_intense_peaks,
                  args.precursor_exclusion_window, args.precursor_losses_exclusion_window, args.beta, args.alpha,
                  args.MS2_relative_peak_intensity, args.isotopologue_precursor_matching,
                  args.exclude_precursor_charge_matching,
                  args.FDR_isotopic_species, args.targets_without_decoys)

    match.final_output()


if __name__ == '__main__':
    matching_scoring()
