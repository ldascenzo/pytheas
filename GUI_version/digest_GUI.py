#!/usr/bin/python3

"""
Last update: December 2020
Credits to: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzo@scripps.edu
GitHub project: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Pytheas in silico digestion, GUI version
"""

from gooey import Gooey, GooeyParser
from enzyme_gui import Enzyme_cleavage
from modify_gui import Modifications
from consolidate_gui import Consolidation
from decoy_gui import Decoys
from calc_mass_gui import Masses


@Gooey(dump_build_config=True, program_name='Pytheas in silico digestion', default_size=(1920, 1080))
def in_silico_digest():
    description = 'Pytheas in silico digestion workflow'
    parser = GooeyParser(description=description)

    # Required Arguments
    parser.add_argument("RNA_sequences", help="Input RNA sequence(s) (fasta format)",
                        widget="MultiFileChooser")
    parser.add_argument("Enzyme", choices=['A', 'T1', 'U2', 'Cus', 'none', 'nonspecific'],
                        help='Nuclease to perform in-silico digestion', gooey_options={'validator': {
            'test': "user_input != 'Select Option'", 'message': 'Enzyme is required'}})
    parser.add_argument('Alphabet_file',
                        help='Excel spreadsheet with the standard & modified nucleotides '
                             'alphabet (default = nts_alphabet_light_std.xlsx)',
                        widget="FileChooser", default='nts_alphabet_light_std.xlsx')
    parser.add_argument('Ion_mode', choices=['+', '-'], default='-',
                        help='Positive (+) or negative (-) ion mode')
    parser.add_argument('MS1_charges_table', default='charges_MS1.txt', widget='FileChooser',
                        help='Charge table for MS1 ions (default = charges_MS1.txt)')
    parser.add_argument('MS2_charges_table', default='charges_MS2.txt', widget='FileChooser',
                        help='Charge table for MS2 ions (default = charges_MS2.txt)')

    # Optional Arguments
    parser.add_argument('--add_decoys', action='store_true', default=False, help='Add decoy sequences to the digest')
    parser.add_argument('--mz_consolidation', action='store_true', default=False,
                        help="Consolidate nucleotides within a given ppm threshold (MSX_ppm_consolidation options)")
    parser.add_argument('--modification_profile', help='File with nucleotide-specific RNA modifications',
                        widget="FileChooser")
    parser.add_argument("--missed_cleavages", choices=['0', '1', '2', '3', '4'], default='0',
                        help='Consecutive missed cleavages to include, up to the selected value')
    parser.add_argument('--alphabet_heavy', default=None, widget='FileChooser',
                        help='Alphabet for isotopic labeled nucleotides and modifications')
    parser.add_argument("--cleaved_fragments_5end_chem", nargs='*', default=['OH'], choices=['OH', 'P'],
                        help="5' chemistry of the cleaved RNA fragments (Ctrl for multiple selection)",
                        widget="Listbox")
    parser.add_argument("--cleaved_fragments_3end_chem", nargs='*', default=['P'], choices=['OH', 'P', 'cP'],
                        help="3' chemistry of the cleaved RNA fragments (Ctrl for multiple selection)",
                        widget="Listbox")
    parser.add_argument("--RNA_5end_chem", nargs='*', default=['P'], choices=['OH', 'P'],
                        help="5' chemistry of the whole input RNA molecule(s) "
                             "(Ctrl for multiple selection)", widget="Listbox")
    parser.add_argument("--RNA_3end_chem", nargs='*', default=['OH'], choices=['OH', 'P', 'cP'],
                        help="3' chemistry of the whole input RNA molecule(s) "
                             "(Ctrl for multiple selection)", widget="Listbox")
    parser.add_argument('--CID_series', nargs='*', widget='Listbox',
                        default=["c", "y", "a", "a-B", "w", "b", "x", "d", "z", "y-P", "z-P"],
                        choices=["c", "y", "a", "a-B", "w", "b", "x", "d", "z", "y-P", "z-P"],
                        help="CID fragments to include when calculating MS2 ions (Ctrl for multiple selection)"
                             "(default=c,y,a,a-B,w,b,x,d,z,y-P,z-P)")
    parser.add_argument('--MS', choices=['1', '2', 'MS1', 'MS2'], default='MS2',
                        help='Choose to generate a digest at the MS1 or MS2 level (default = MS2)')
    parser.add_argument('--min_length', default=3, type=int,
                        help='Minimum length of redundant sequences to be consolidated together (default = 3)')
    parser.add_argument('--MS1_ppm_consolidation', default=0, type=int,
                        help='ppm threshold for MS1 masses for mz_consolidation (default = 0)')
    parser.add_argument('--MS2_ppm_consolidation', default=0, type=int,
                        help='ppm threshold for MS2 masses for mz_consolidation (default = 0)')
    parser.add_argument('--MS1_mzmin', default=400, type=int,
                        help='Minimum value of the MS1 ions m/z window (default=400)')
    parser.add_argument('--MS1_mzmax', default=2000, type=int,
                        help='Maximum value of the MS1 ions m/z window (default=2000)')
    parser.add_argument('--MS2_mzmin', default=300, type=int,
                        help='Minimum value of the MS2 ions m/z window (default=300)')
    parser.add_argument('--MS2_mzmax', default=2000, type=int,
                        help='Maximum value of the MS2 ions m/z window (default=2000)')
    parser.add_argument('--nonspecific_min_length', type=int, default=3,
                        help='Minimum nts length from nonspecific cleavage (default = 3)')
    parser.add_argument('--nonspecific_max_length', type=int, default=10,
                        help='Maximum nts length from nonspecific cleavage (default = 10)')

    ####################################################
    args = parser.parse_args()

    # Step 1/4 - Enzymatic Cleavage
    enzyme = Enzyme_cleavage(args.RNA_sequences, args.Enzyme, args.missed_cleavages, args.nonspecific_min_length,
                             args.nonspecific_max_length, args.cleaved_fragments_5end_chem,
                             args.cleaved_fragments_3end_chem,
                             args.RNA_5end_chem, args.RNA_3end_chem)
    enzyme.final_output()

    # Step 2/4 - Add modifications
    mods = Modifications(args.modification_profile, args.Alphabet_file)
    mods.final_output()

    # Step 3/4 - Consolidate at MS1/MS2 levels
    consolidate = Consolidation(args.MS, args.min_length)
    consolidate.final_output()

    # Step 3.5/4 - Add Decoys (if selected)
    if args.add_decoys:
        decoy = Decoys()
        decoy.final_output()

    # Step 4/4 - Calculate final masses
    calc_mass = Masses(args.Ion_mode, args.Alphabet_file, args.alphabet_heavy, args.MS1_charges_table,
                       args.MS2_charges_table, args.MS1_mzmin, args.MS1_mzmax, args.MS2_mzmin, args.MS2_mzmax,
                       args.CID_series, args.mz_consolidation, args.MS1_ppm_consolidation, args.MS2_ppm_consolidation)
    calc_mass.final_output()


if __name__ == '__main__':
    in_silico_digest()
