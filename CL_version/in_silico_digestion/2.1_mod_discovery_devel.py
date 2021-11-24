#!/usr/bin/python3

"""
Last update: March 2021
Author: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzoluigi@gmail.com
GitHub project repository: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Preliminary work-work-in-progress step towards the Pytheas support of discovery mode. At the present state, the user can
specify one or more modifications that will be added in any possible combination to the available sequences extracted
from the output.2 file.
If for instance a methylated C is requested, all the instances of C with all possible combination of modified/unmodified
are computed. Note that the process can become computationally intensive if many "discovery" modifications
are added together for long sequences (it is advised to use only one for this version)

***OPTIONS***
--nts_light (REQUIRED): Excel spreadsheet with the standard nucleotides
                                                                      and optional modification alphabet.
--mod_discovery (REQUIRED):  text file with RNA modifications, must be in the same directory of the script. Details
                             on format and usage can be found in the Digest section of the Pytheas manual
--max_mods_per_fragment (OPTIONAL): Specify the maximum number of variable modifications allowed per RNA fragment
                                    (Default = 2)

***OUTPUT***
1) output.2 -> same format of output.1 with added info on modified nucleotides.

"""

import argparse
import os
import sys
from itertools import permutations
from itertools import product

import pandas as pd

# Initialize and define launch options
parser = argparse.ArgumentParser(description='List of available options')
parser.add_argument('--mod_discovery', required=True,
                    help='File with molecule-specific RNA modifications '
                         '(Required, be sure to fill the Originating_base column!!!!)')
parser.add_argument('--nts_light', required=True,
                    help='Excel spreadsheet with the standard nucleotides and optional modification alphabet '
                         '(Required)')
parser.add_argument('--max_mods_per_fragment', default=2, type=int,
                    help='Specify the maximum number of variable modifications allowed per RNA fragment (Default = 2)')
args = parser.parse_args()


def read_excel_input(nts_alph=args.nts_light):
    """
     Create three dictionaries mod_alphabet, mod_origin and mod_partial
     * mod_alphabet contains all ID : ID_ext couples, one letter and extended IDs or all modified base
     * mod_origin contains all the unmodified nucleotides and is used for modification positions validation
     * mod_partial contains all the partial modified nucleotides for each modification (as from the alphabet file)
     in case the option for partial modifications is selected
     """
    mod_origin, mod_alphabet = {}, {}

    if nts_alph:
        # Checking that the nts_light file given in argument exists
        if not os.path.exists(nts_alph):
            print("ERROR! File {} does not exist. Execution terminated without output".format(nts_alph))
            sys.exit(1)

        # Create a dataframe with info from Excel spreadsheet
        df = pd.read_excel(nts_alph, header=12)

        # Drop rows with the 4 standard nucleobases
        df = df[df.ID != 'C']
        df = df[df.ID != 'G']
        df = df[df.ID != 'A']
        df = df[df.ID != 'U']

        # Drop rows with NaN values
        df = df[pd.notnull(df['ID'])]

        # Transform all ID values in string (so numbers can be used as one letter code for bases)
        df = df.astype({"ID": str})

        # Create a dictionary with the pair ID : ID_ext as modification alphabet
        mod_alphabet = dict(zip(df.ID, df.ID_ext))

        # Create a dictionary with the ID : Originating_base for quality control purposes
        mod_origin = dict(zip(df.ID, df.Originating_base))

    return mod_origin, mod_alphabet


def parse_input(modfile):
    """
    Create a list with all the modifications to consider for discovery mode
    """
    # Checking that the nts_light file given in argument exists
    if not os.path.exists(modfile):
        print("ERROR! File {} does not exist. Execution terminated without output".format(modfile))
        sys.exit(1)

    with open(modfile, 'r') as f:
        for line in f:
            yield line.strip().split()


def mod_input(mods, header):
    """
    Create a dataframe with all the input modification to be included into discovery mode
    """
    df = pd.DataFrame(mods)

    new_header = df.iloc[header]
    df = df[(header + 1):]
    df.columns = new_header

    return df

def make_patterns(seq, mods):
    """
    Output all the combinations of base modifications in a given sequence
    """
    seq_list, new_seqs = list(seq), []

    # Finds indices of key letters in seq_list
    indices = [i for i, c in enumerate(seq_list) if c in mods]

    # Set the maximum number of modifications per fragment
    if len(indices) > args.max_mods_per_fragment:

        mod_sites = permutations(indices, args.max_mods_per_fragment)

        for sites in mod_sites:

            for i in range(len(sites) + 1):

                out_seq = list(seq)

                for position in sites[:i]:
                    out_seq[position] = mods[0]

                new_seqs.append("".join(out_seq))
                i -= 1

    # If the fragment contains less than the specified maximum number of variable modifications, apply all    
    else:
        for t in product(mods, repeat=len(indices)):
            for i, c in zip(indices, t):
                seq_list[i] = c

            new_seqs.append(''.join(seq_list))

    return list(set(new_seqs))


def new_fragments(mod_origin, mod_alph, mods=mod_input(parse_input(args.mod_discovery), 0),
                  input_lines=mod_input(parse_input("output.2"), 7)):
    """
    Generate all the possible combination of base modifications for the input fragments
    
    Combinations are generated via cartesian product between the original sequence and the
    modified/unmodified nucleobase
    """
    out_lines = []

    # Creates a list with all the input fragment sequences
    input_nts = list(input_lines['Seq'])
    new_seqs = list(input_lines['Seq'])

    # Creates a list with all the input modifications for discovery mode
    input_mods = list(mods['ID'])

    # Adds all the new modification combinations
    for mod in input_mods:
        for nts in input_nts:

            mol, nstart, nend, miss, chem3, chem5 = (input_lines.loc[input_lines['Seq'] == nts, 'Molecule'].iloc[0],
                                       input_lines.loc[input_lines['Seq'] == nts, 'Nstart'].iloc[0],
                                       input_lines.loc[input_lines['Seq'] == nts, 'Nend'].iloc[0],
                                       input_lines.loc[input_lines['Seq'] == nts, 'Miss'].iloc[0],
                                       input_lines.loc[input_lines['Seq'] == nts, "3'chem"].iloc[0],
                                       input_lines.loc[input_lines['Seq'] == nts, "5'chem"].iloc[0])

            for new_seq in make_patterns(nts, mod + mod_origin[mod]):

                if new_seq not in new_seqs:

                    # Determine the sequence written in human readable format
                    modification_input, mod_flag = '', 0
                    for s in new_seq:
                        if s == 'A' or s == 'U' or s == 'G' or s == 'C':
                            modification_input += s

                        else:
                            modification_input += mod_alph[s]
                            mod_flag = 1

                    if mod_flag == 0:
                        modification_input = '-'

                    # Prepare the final line for the output    
                    line = "{} {} {} {} {} {} {} {}\n".format(mol, new_seq, nstart, nend, miss, chem5, chem3,
                                                              modification_input)

                    new_seqs.append(new_seq)
                    out_lines.append(line)

    return (out_lines)


if __name__ == "__main__":

    if not os.path.exists(os.getcwd() + "/output.2"):
        print("ERROR! File output.2 is missing. Execution terminated.")
        sys.exit(1)

    mod_origin, mod_alph = read_excel_input()

    with open("output.2", 'a') as infile:
        infile.writelines(new_fragments(mod_origin, mod_alph))
