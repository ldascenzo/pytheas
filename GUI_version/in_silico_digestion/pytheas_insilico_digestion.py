8#!/usr/bin/python3

"""
Last update: May 2021
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
    description = 'Generate theoretical digest library'
    parser = GooeyParser(description=description)

    # Required Arguments
    parser.add_argument("RNA_sequences", help="Input RNA sequence(s) in fasta format",
                        widget="MultiFileChooser")
    parser.add_argument("Enzyme", choices=['A', 'T1', 'U2', 'Cus', 'MC1', 'MAZ', 'none', 'nonspecific'],
                        help='RNA endonuclease for the in silico digestion', gooey_options={'validator': {
            'test': "user_input != 'Select Option'", 'message': 'Enzyme is required'}})
    parser.add_argument('Nucleotides_light',
                        help='Elemental composition file for standard and modified nucleotides (Excel spreadsheet)',
                        widget="FileChooser", default='nts_light.xlsx')
    parser.add_argument('Ion_mode', choices=['+', '-'], default='-',
                        help='Negative (-) or Positive (+)')
    parser.add_argument('MS1_charges_table', default='charges_MS1.txt', widget='FileChooser',
                        help='Charge table for precursor ions (MS1)')
    parser.add_argument('MS2_charges_table', default='charges_MS2.txt', widget='FileChooser',
                        help='Charge table for fragment ions (MS2)')

    # Optional Arguments
    parser.add_argument('--decoys', action='store_true', default=False, help='Add decoy sequences')
    parser.add_argument('--list_of_known_RNA_modifications', help='List of nucleotide modifications info',
                        widget="FileChooser")
    parser.add_argument("--enzyme_missed_cleavages", choices=['0', '1', '2', '3', '4'], default='0',
                        help='Allowed consecutive missed cleavages for the selected nuclease')
    parser.add_argument('--Nucleotides_heavy', default=None, widget='FileChooser',
                        help='Elemental composition file for isotopically labeled nucleotides and modifications '
                             '(Excel spreadsheet))')
    parser.add_argument("--cleaved_RNA_5end_chemistry", nargs='*', default=['OH'], choices=['OH', 'P'],
                        help="5' chemistry of the nucleolytic RNA fragments (Hold Ctrl for multiple selection)",
                        widget="Listbox")
    parser.add_argument("--cleaved_RNA_3end_chemistry", nargs='*', default=['P'], choices=['OH', 'P', 'cP'],
                        help="3' chemistry of the nucleolytic RNA fragments (Hold Ctrl for multiple selection)",
                        widget="Listbox")
    parser.add_argument("--RNA_5end_chemistry", nargs='*', default=['P'], choices=['OH', 'P'],
                        help="5' chemistry of the input RNA molecule(s) "
                             "(Hold Ctrl for multiple selection)", widget="Listbox")
    parser.add_argument("--RNA_3end_chemistry", nargs='*', default=['OH'], choices=['OH', 'P', 'cP'],
                        help="3' chemistry of the input RNA molecule(s) "
                             "(Hold Ctrl for multiple selection)", widget="Listbox")
    parser.add_argument('--CID_HCD_series', nargs='*', widget='Listbox',
                        default=["c", "y", "a", "a-B", "w", "b", "x", "d", "z", "y-P", "z-P"],
                        choices=["c", "y", "a", "a-B", "w", "b", "x", "d", "z", "y-P", "z-P"],
                        help="MS2 ion series for m/z calculation (Hold Ctrl for multiple "
                             "selection)")
    parser.add_argument('--SeqX_consolidation', action='store_true', default=False,
                        help="Consolidate sequences that cannot be uniquely identified below specified "
                             "MS1_SeqX and MS2_SeqX thresholds")
    parser.add_argument('--MS1_SeqX', default=0, type=float,
                        help='Precursor ion mass threshold for SeqX consolidation (in ppm) ')
    parser.add_argument('--MS2_SeqX', default=0, type=float,
                        help='MS2 ions mass threshold for SeqX consolidation (in ppm)')
    parser.add_argument('--MS1_mz_min', default=400, type=int, help='Minimum value of the precursor ions m/z')
    parser.add_argument('--MS1_mz_max', default=2000, type=int,
                        help='Maximum value of the precursor ions m/z')
    parser.add_argument('--MS2_mz_min', default=300, type=int,
                        help='Minimum value of the fragment ions m/z')
    parser.add_argument('--MS2_mz_max', default=2000, type=int,
                        help='Maximum value of the fragment ions m/z')
    parser.add_argument('--nonspecific_cleavage_min_length', type=int, default=3,
                        help='Minimum sequence length if nonspecific cleavage is selected')
    parser.add_argument('--nonspecific_cleavage_max_length', type=int, default=10,
                        help='Maximum sequence length if nonspecific cleavage is selected')


    ####################################################
    args = parser.parse_args()

    # Step 1/4 - Enzymatic Cleavage
    enzyme = Enzyme_cleavage(args.RNA_sequences, args.Enzyme, args.enzyme_missed_cleavages,
                             args.nonspecific_cleavage_min_length,
                             args.nonspecific_cleavage_max_length, args.cleaved_RNA_5end_chemistry,
                             args.cleaved_RNA_3end_chemistry,
                             args.RNA_5end_chemistry, args.RNA_3end_chemistry)
    enzyme.final_output()

    # Step 2/4 - Add modifications
    mods = Modifications(args.list_of_known_RNA_modifications, args.Nucleotides_light)
    mods.final_output()

    # Step 3/4 - Consolidate at MS1/MS2 levels
    consolidate = Consolidation()
    consolidate.final_output()

    # Step 3.5/4 - Add Decoys (if selected)
    if args.decoys:
        decoy = Decoys()
        decoy.final_output()

    # Step 4/4 - Calculate final masses
    calc_mass = Masses(args.Ion_mode, args.Nucleotides_light,
                       args.Nucleotides_heavy,
                       args.MS1_charges_table,
                       args.MS2_charges_table, args.MS1_mz_min, args.MS1_mz_max, args.MS2_mz_min, args.MS2_mz_max,
                       args.CID_HCD_series, args.SeqX_consolidation, args.MS1_SeqX,
                       args.MS2_SeqX, args.RNA_sequences)
    calc_mass.final_output()


if __name__ == '__main__':
    in_silico_digest()
