#!/usr/bin/python3

"""
Last update: April 2021
Author: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzoluigi@gmail.com
GitHub project repository: https://github.com/ldascenzo/pytheas

DESCRIPTION
Mapping routine of the Pytheas workflow

OUTPUT
1) mapping_output -> html file with information on the RNA sequences and the mapped sequences
"""

from gooey import Gooey, GooeyParser
from mapping_library import Mapping


@Gooey(dump_build_config=True, program_name='Pytheas mapping of matching targets on the input RNA sequence(s)',
       default_size=(1920, 1080))
def mapper():
    description = 'Target sequences are mapped over the original RNA sequence using an informative color code scheme'
    parser = GooeyParser(description=description)

    # Required arguments
    parser.add_argument('Final_report',  help='Final Report file in csv format', widget="FileChooser")
    parser.add_argument('Alphabet', widget="FileChooser", help='Excel spreadsheet file with nucleotides and '
                                                               'modification alphabet used for the Digest')
    parser.add_argument('RNA_sequence', widget="FileChooser", help='RNA sequence(s) (fasta format) used for the in '
                                                                   'silico digestion')

    # Optional arguments
    parser.add_argument('--minimum_length', default=3, type=int, help='Minimum target nucleotide length to be '
                                                                      'mapped on the RNA sequence')
    parser.add_argument('--Sp_score_cutoff', default=0, type=float, help='Minimum Sp score cutoff for targets'
                                                                         'to be mapped on the RNA sequence')

    ####################################################
    args = parser.parse_args()

    map_sequences = Mapping(args.Alphabet, args.Final_report, args.RNA_sequence, args.minimum_length, args.Sp_score_cutoff)

    map_sequences.final_output()


if __name__ == '__main__':
    mapper()
