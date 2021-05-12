#!/usr/bin/python3

"""
Last update: April 2021
Author: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzoluigi@gmail.com
GitHub project repository: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Pytheas visualization algorithm.
Additional information on the output files and the parameters can be found in the Matching&Scoring section of the
Pytheas manual

***OUTPUT***
1) visualization_output -> html file containing a table with all the visualized targets/decoys and the spectra
"""

from gooey import Gooey, GooeyParser
from visualization_library import Visualize


@Gooey(dump_build_config=True, program_name='Pytheas visualization of annotated spectra',
       default_size=(1920, 1080))
def visualizer():
    description = 'Generate a summary table and graphical output of spectra matched to candidate sequences'
    parser = GooeyParser(description=description)

    # Required arguments
    parser.add_argument('Theoretical_digest', widget="FileChooser", help='Pytheas in silico digestion output file')
    parser.add_argument('MS_data', widget="FileChooser", help='Experimental data peak list in mgf format')
    parser.add_argument('Match_output', widget="FileChooser", help='Pytheas matching and scoring .txt output file')

    # Optional arguments
    parser.add_argument('--visualize_decoys', action='store_true', default=False,
                        help='Output decoy matches. By default, only targets are shown')
    parser.add_argument('--modified_only', action='store_true', default=False,
                        help='Only include matches containing nucleotide modifications')
    parser.add_argument('--unique_positions_only', action='store_true', default=False,
                        help='Only include matches mapping to unique positions within input RNA sequence(s)')
    parser.add_argument('--remove_redundant_SeqX_matches', action='store_true', default=False,
                        help='When SeqX is used for sequence consolidation, output only the highest ranking matches '
                             'with X')
    parser.add_argument('--rank_maximum', default=99, type=int,
                        help='Maximum rank value for matched sequences (default=99)')
    parser.add_argument('--Sp_minimum', default=0, type=float,
                        help='Minimum Sp score cutoff for matched sequences (default=0)')
    parser.add_argument('--dSp_maximum', default=1, type=float,
                        help='Maximum dSp score cutoff for matched sequences (default=1)')
    parser.add_argument('--dSp2_minimum', default=0, type=float,
                        help='Minimum dSp2 score cutoff for matched sequences (default=0)')
    parser.add_argument('--MS1_ppm_cutoff', default=None, type=float,
                        help='Precursor ion mass tolerance cutoff. Only matches with precursor ion m/z within '
                             'specified +- ppm window are included (default=no cutoff)')
    parser.add_argument('--MS2_peaks_number', default='all',
                        help='Number of most intense MS2 peaks in the spectra. By default all peaks are shown')
    parser.add_argument('--spectra_mz_minimum', default=None, type=int, help='Minimum m/z value for the spectra '
                                                                             '(default=based on mgf)')
    parser.add_argument('--spectra_mz_maximum', default=None, type=int, help='Maximum m/z value for the spectra '
                                                                             '(default=based on mgf)')

    ####################################################
    args = parser.parse_args()

    visualization = Visualize(args.Theoretical_digest, args.MS_data, args.Match_output, args.MS2_peaks_number,
                              args.spectra_mz_minimum, args.spectra_mz_maximum, args.Sp_minimum,
                              args.dSp_maximum, args.visualize_decoys, args.modified_only, args.rank_maximum,
                              args.unique_positions_only, args.remove_redundant_SeqX_matches,
                              args.dSp2_minimum, args.MS1_ppm_cutoff)

    visualization.final_output()


if __name__ == '__main__':
    visualizer()
