#!/usr/bin/python3

"""
Last update: March 2021
Author: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzoluigi@gmail.com
GitHub project repository: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Pytheas in silico digestion generation pipeline. Additional information on the steps and parameters/input file can be
found in the Digest section of the Pytheas manual
"""

from gooey import Gooey, GooeyParser
from enzyme_gui import Enzyme_cleavage
from modify_gui import Modifications
from consolidate_gui import Consolidation
from decoy_gui import Decoys
from calc_mass_gui import Masses


@Gooey(dump_build_config=True, program_name='Pytheas in silico digestion', default_size=(1920, 1080))
def in_silico_digest():
    description = 'Theoretical digest library generation workflow'
    parser = GooeyParser(description=description)

    # Required Arguments
    parser.add_argument("RNA_sequences", help="Input RNA sequence(s) (fasta)",
                        widget="MultiFileChooser")
    parser.add_argument("Enzyme", choices=['A', 'T1', 'U2', 'Cus', 'none', 'nonspecific'],
                        help='RNA endonuclease for the in silico digestion', gooey_options={'validator': {
            'test': "user_input != 'Select Option'", 'message': 'Enzyme is required'}})
    parser.add_argument('Nucleotides_alphabet',
                        help='Alphabet file for standard and modified nucleotides (Excel spreadsheet)',
                        widget="FileChooser", default='nts_alphabet_light_std.xlsx')
    parser.add_argument('Ion_mode', choices=['+', '-'], default='-',
                        help='Negative (-) or Positive (+)')
    parser.add_argument('MS1_charges_table', default='charges_MS1.txt', widget='FileChooser',
                        help='Charge table for precursor ions (MS1)')
    parser.add_argument('MS2_charges_table', default='charges_MS2.txt', widget='FileChooser',
                        help='Charge table for fragment ions (MS2)')

    # Optional Arguments
    parser.add_argument('--decoys', action='store_true', default=False, help='Add decoy sequences')
    parser.add_argument('--modification_profile', help='Input file with nucleotide modifications info',
                        widget="FileChooser")
    parser.add_argument("--missed_cleavages", choices=['0', '1', '2', '3', '4'], default='0',
                        help='Allowed consecutive missed cleavages for the selected nuclease')
    parser.add_argument('--isotopic_labeled_alphabet', default=None, widget='FileChooser',
                        help='Alphabet file for isotopic labeled nucleotides and modifications (heavy)')
    parser.add_argument("--cleaved_fragments_5end_chem", nargs='*', default=['OH'], choices=['OH', 'P'],
                        help="5' chemistry of the nucleolytic RNA fragments (Hold Ctrl for multiple selection)",
                        widget="Listbox")
    parser.add_argument("--cleaved_fragments_3end_chem", nargs='*', default=['P'], choices=['OH', 'P', 'cP'],
                        help="3' chemistry of the nucleolytic RNA fragments (Hold Ctrl for multiple selection)",
                        widget="Listbox")
    parser.add_argument("--RNA_5end_chem", nargs='*', default=['P'], choices=['OH', 'P'],
                        help="5' chemistry of the whole input RNA molecule(s) "
                             "(Hold Ctrl for multiple selection)", widget="Listbox")
    parser.add_argument("--RNA_3end_chem", nargs='*', default=['OH'], choices=['OH', 'P', 'cP'],
                        help="3' chemistry of the whole input RNA molecule(s) "
                             "(Hold Ctrl for multiple selection)", widget="Listbox")
    parser.add_argument('--CID_series', nargs='*', widget='Listbox',
                        default=["c", "y", "a", "a-B", "w", "b", "x", "d", "z", "y-P", "z-P"],
                        choices=["c", "y", "a", "a-B", "w", "b", "x", "d", "z", "y-P", "z-P"],
                        help="CID fragment ion series for MS2 ions m/z calculation (Hold Ctrl for multiple "
                             "selection)")
    parser.add_argument('--mz_consolidation', action='store_true', default=False,
                        help="Consolidate sequences differing only by nucleotides below a given ppm threshold")
    parser.add_argument('--MS1_ppm_consolidation', default=0, type=float,
                        help='ppm threshold for precursor ion masses for mz_consolidation')
    parser.add_argument('--MS2_ppm_consolidation', default=0, type=float,
                        help='ppm threshold for MS2 masses for mz_consolidation')
    parser.add_argument('--minimum_length_consolidation', default=3, type=int,
                        help='Minimum length of redundant sequences to consolidate')
    parser.add_argument('--MS', choices=['1', '2', 'MS1', 'MS2'], default='MS2',
                        help='Choose the MS level for the final output')
    parser.add_argument('--MS1_mzmin', default=400, type=int, help='Minimum value of precursor ions m/z window')
    parser.add_argument('--MS1_mzmax', default=2000, type=int,
                        help='Maximum value of the precursor ions m/z window')
    parser.add_argument('--MS2_mzmin', default=300, type=int,
                        help='Minimum value of the fragment ions m/z window')
    parser.add_argument('--MS2_mzmax', default=2000, type=int,
                        help='Maximum value of the fragment ions m/z window')
    parser.add_argument('--nonspecific_min_length', type=int, default=3,
                        help='Minimum sequence length if nonspecific cleavage is selected')
    parser.add_argument('--nonspecific_max_length', type=int, default=10,
                        help='Maximum sequence length if nonspecific cleavage is selected')



    ####################################################
    args = parser.parse_args()

    # Step 1/4 - Enzymatic Cleavage
    enzyme = Enzyme_cleavage(args.RNA_sequences, args.Enzyme, args.missed_cleavages, args.nonspecific_min_length,
                             args.nonspecific_max_length, args.cleaved_fragments_5end_chem,
                             args.cleaved_fragments_3end_chem,
                             args.RNA_5end_chem, args.RNA_3end_chem)
    enzyme.final_output()

    # Step 2/4 - Add modifications
    mods = Modifications(args.modification_profile, args.Nucleotides_alphabet)
    mods.final_output()

    # Step 3/4 - Consolidate at MS1/MS2 levels
    consolidate = Consolidation(args.MS, args.minimum_length_consolidation)
    consolidate.final_output()

    # Step 3.5/4 - Add Decoys (if selected)
    if args.decoys:
        decoy = Decoys()
        decoy.final_output()

    # Step 4/4 - Calculate final masses
    calc_mass = Masses(args.Ion_mode, args.Nucleotides_alphabet, args.isotopic_labeled_alphabet, args.MS1_charges_table,
                       args.MS2_charges_table, args.MS1_mzmin, args.MS1_mzmax, args.MS2_mzmin, args.MS2_mzmax,
                       args.CID_series, args.mz_consolidation, args.MS1_ppm_consolidation, args.MS2_ppm_consolidation)
    calc_mass.final_output()


if __name__ == '__main__':
    in_silico_digest()
