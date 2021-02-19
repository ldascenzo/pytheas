#!/usr/bin/python3

"""
Last update: December 2020
Credits to: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzo@scripps.edu
GitHub project: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Pytheas mapping, GUI version

DESCRIPTION
Mapping routine
"""

from gooey import Gooey, GooeyParser
#from pytheas_filter_match_gui import Match_Filtering
from pytheas_mapping_gui import Mapping



@Gooey(dump_build_config=True, program_name='Pytheas mapping of matching targets on the input RNA sequence(s)',
       default_size=(1920, 1080))
def mapper():
    description = 'Mapping of the matching targets on the input RNA sequence(s)'
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

    #filtering = Match_Filtering(args.Match_output)

    #filtering.final_output()

    # processed_xx file obtained from filtering the matching output
    #if platform.system() == 'Windows':
    #    naming_output = args.Match_output.split(';')

    #else:
    #    naming_output = args.Match_output.split(':')

    #processed_file = "processed_{}.csv".format("_".join(naming_output[0].split('/')[-1].split('.')[0].split("_")[2:]))

    map_sequences = Mapping(args.Alphabet, args.Final_report, args.RNA_sequence, args.minimum_length, args.Sp_score_cutoff)

    map_sequences.final_output()


if __name__ == '__main__':
    mapper()
