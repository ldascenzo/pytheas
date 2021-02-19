#!/usr/bin/env python3

"""
Last update: January 2021
Credits to: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzoluigi@gmail.com
GitHub project: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Script to elaborate the final text output in the Pytheas workflow, GUI version.

"""

import os, sys
import pandas as pd
import ntpath


class Pytheas_Output:
    def __init__(self, match_file, Sp_cutoff, dSp_cutoff, visualize_decoys, only_modified, only_unique_positions, rank,
                 redundant_x):
        self.match_file, self.Sp_cutoff, self.dSp_cutoff = match_file, Sp_cutoff, dSp_cutoff
        self.visualize_decoys, self.only_modified = visualize_decoys, only_modified
        self.only_unique_positions = only_unique_positions
        self.rank, self.redundant_x = rank, redundant_x

    def parse_match_file(self):
        """
        Parses the match_output file creating a dataframe with all the info
        """
        d = {'m/z': [], 'RT': [], 'MS1_offset(ppm)': [], 'isotopologue match': [],'length': [], "5'-end": [],
             'sequence': [], "3'-end": [], 'charge': [], 'isotope': [], 'Rank': [], 'Score (Sp)': [], 'dSp': [],
             'n/L': [], 'sequence_mods': [], 'molecule ID': [], 'sequence_location': [],
             'residue_start(unique)': [], 'residue_end(unique)': [], 'modifications_positions': []}

        sequences_with_x = []

        if not os.path.exists(self.match_file):
            print("File {} does not exist. Execution terminated without output".format(self.match_file))
            sys.exit(1)

        # Prepare the variables for the test on the options about decoys, modified matches and unique positions
        if self.visualize_decoys:
            decoys_string = '@#$%'
        else:
            decoys_string = 'decoy'

        if self.only_modified:
            mod_string = '['
        else:
            mod_string = ''

        if self.only_unique_positions:
            unique_length = 1
        else:
            unique_length = 999

        with open(self.match_file, 'r') as input_file:
            for line in input_file:

                if 'PRECURSOR_ION' in line:
                    unique_precursor = line.split('=')[1]

                if line[0].isdigit():
                    sp = line.split()

                    score, dSp, molecule, mod = sp[8], sp[6], sp[17], sp[14]
                    unique_list = sp[17].split(';')
                    rank = int(sp[-2])

                    # Only include data within requested cutoffs for Sp and dSp
                    if float(score) >= self.Sp_cutoff and float(dSp) <= self.dSp_cutoff \
                            and decoys_string not in molecule and mod_string in mod and \
                            len(unique_list) <= unique_length and rank <= self.rank:
                        mz, rt, seq = sp[0], sp[1].split('=')[1], sp[13]

                        # Make sure that only the top scoring/lowest MS1 offset sequence containing 'X' is
                        # used for output
                        flag_x = False

                        if self.redundant_x:
                            if 'X' in seq:
                                if '{}_{}_{}'.format(seq, unique_precursor, rt) not in sequences_with_x:
                                    sequences_with_x.append('{}_{}_{}'.format(seq, unique_precursor, rt))
                                else:
                                    flag_x = True

                        if not flag_x:

                            # Add the info about isotopologues in the appropriate column, removing the * from the
                            # m/z
                            if '*' in mz:
                                mz_corrected = mz[:-1]
                                d['isotopologue match'].append('y')
                            else:
                                mz_corrected = mz
                                d['isotopologue match'].append('')

                            d['m/z'].append(mz_corrected)

                            # Insert all the other info for all ions in preparation for the table output
                            d['Rank'].append(sp[-2])
                            d['RT'].append(rt), d['Score (Sp)'].append(score), d['dSp'].append(dSp), d['charge'].append(
                                sp[12]),
                            d['isotope'].append(sp[11]), d['sequence_location'].append(sp[17]), \
                            d['sequence'].append(seq), d['sequence_mods'].append(sp[14]),
                            d["5'-end"].append(sp[16]), d["3'-end"].append(sp[15]), d['MS1_offset(ppm)'].append(
                                sp[3].split('p')[0]), d['length'].append(len(seq)), \
                            d['n/L'].append(round(int(sp[-1].split(';')[1].split('=')[1]) /
                                                  int(sp[-1].split(';')[-2].split('=')[1]), 3))
                            d['molecule ID'].append(extract_molecule(molecule))

                            # Add the information about residue start and end for unique matches
                            if len(unique_list) == 1 and 'decoy' not in sp[17]:
                                d['residue_start(unique)'].append(sp[17].split(',')[1])
                                d['residue_end(unique)'].append(sp[17].split(',')[2])

                            else:
                                d['residue_start(unique)'].append('')
                                d['residue_end(unique)'].append('')

                            # Add the modifications positions if modified nucleotides are present in the match
                            if sp[14] != '-' and 'X' not in sp[13] and 'decoy' not in sp[17]:
                                d['modifications_positions'].append(modification_position(sp[17], sp[13], sp[14]))
                            else:
                                d['modifications_positions'].append('')

        df = pd.DataFrame(data=d)
        df.index += 1

        return df.to_csv('final_report_{}.csv'.format(
            ("".join(self.match_file.split('match_output_')[1]).split('.')[0])))


def extract_molecule(string):
    """
    Extract the molecule(s) from an input string in the format moleculeID,startres,endres;moleculeID,startres,endres
    """
    unique_molecule_id = []

    for molecule_info in string.split(';'):
        molecule_id = molecule_info.split(',')[0]
        if molecule_id not in unique_molecule_id:
            unique_molecule_id.append(molecule_id)

    return ";".join(unique_molecule_id)


def modification_position(molecule, modification_short, modification_ext):
    """
    Extract a list of the positions for the modified nucleotides inside a given match
    """
    standard_bases = ['A', 'G', 'U', 'C']
    modification_positions, output_strings = [], []

    for index, nt in enumerate(modification_short):
        if nt not in standard_bases:
            modification_positions.append(index)

    modification_order = []
    for mod in modification_ext.split(']'):
        if '[' in mod:
            modification_order.append(mod.split('[')[-1])

    for position, mod in enumerate(modification_order):
        for match in molecule.split(';'):
            split = match.split(',')
            output_strings.append("{}:{}[{}]".format(split[0], int(split[1]) + modification_positions[position], mod))

    return "|".join(output_strings)