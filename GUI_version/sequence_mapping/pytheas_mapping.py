#!/usr/bin/python3

"""
Last update: January 2022
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


@Gooey(dump_build_config=True, program_name='Pytheas sequence mapping', default_size=(1920, 1080))
def mapper():
    description = 'Generate a graphical output of identified targets mapped on the input RNA sequence(s)'
    parser = GooeyParser(description=description)

    # Required arguments
    parser.add_argument('Final_report',  help='Pytheas final report output file', widget="FileChooser")
    parser.add_argument('Nucleotides_list', widget="FileChooser", help='Elemental composition file for standard and '
                                                                       'modified nucleotides (Excel spreadsheet)')
    parser.add_argument('RNA_sequence', widget="FileChooser", help='Input RNA sequence in fasta format')

    # Optional arguments
    parser.add_argument('--minimum_targets_length', default=3, type=int, help='Minimum targets length to be mapped on '
                                                                              'the RNA sequence')
    parser.add_argument('--Sp_cutoff', default=0, type=float, help='Minimum Sp score cutoff for targets'
                                                                         'to be mapped on the RNA sequence')

    ####################################################
    args = parser.parse_args()

    map_sequences = Mapping(args.Nucleotides_list, args.Final_report, args.RNA_sequence, args.minimum_targets_length,
                            args.Sp_cutoff)

    map_sequences.final_output()


if __name__ == '__main__':
    mapper()
