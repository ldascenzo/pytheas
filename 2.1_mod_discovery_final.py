#!/usr/bin/python3

"""
Last update: May 2020
Credits to: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzo@scripps.edu

***DESCRIPTION***
This script is a preliminary step towards the support for discovery mode. It allows to specify one or more
modifications that will be added in any possible combination to the available sequences from the output.2 file.
If for instance a methylated C is requested, all the instances of C with all possible combination
will be added as modified. Note that the process can become computationally intensive if many "discovery" modifications
are added together for long sequences (it is advised to use only one for this test version)

***USAGE***
python 2.1_mod_discovery.py --OPTIONS

***OPTIONS***
--nts_alphabet_light (REQUIRED): Excel spreadsheet with the standard nucleotides and optional modification alphabet,
                                 the standard version is nts_alphabet_light.xslx
                                 (Be sure to fill the Originating_base column!!!!)
--mod_discovery (REQUIRED):  text file with RNA modifications, must be in the same directory of the script
--max_mods_per_fragment (OPTIONAL): Specify the maximum number of variable modifications allowed per RNA fragment
                                    (Default = 2)

Entries in the modification profile file are in the format: "Molecule Position Modification Include" |One entry per line

Molecule -> The same of fasta header (special care for uppercase/lowercase characters)
Modification -> One of the modified bases within modification list


***INPUT***
mod_discovery file with the following format:

Mol ID ID_ext
xxx x xxx


***RNA MODIFICATIONS LIST***
The nts_alphabet_light file must have the structure present in the standard nts_alphabet_light.xlsx. Namely:
1) One letter name (ID) for modified nt -> uppercase for ribose modifications and lowercase for base mods/both base and
   backbone
2) Longer ID (ID_ext) introduced in the output for readability
3) Extended chemical name of the modified nt
4) Originating base for the given modification
5) Possible alternative IDs in case of partial modification
6) Atomic composition of the nucleobase and backbone (use convention in attached scheme "Std_atomic_composition")

Std backbone composition has PO3H- 3' end and -CH2 5' end (see attached scheme "Std_atomic_composition)

RNA modification database: http://modomics.genesilico.pl/

***OUTPUT***
output.2 - in addition to the info from output.1, info about modification are added in the fragments
(both one letter and human-readable formats)

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
parser.add_argument('--nts_alphabet_light', required=True,
                    help='Excel spreadsheet with the standard nucleotides and optional modification alphabet '
                         '(Required)')
parser.add_argument('--max_mods_per_fragment', default=2, type=int,
                    help='Specify the maximum number of variable modifications allowed per RNA fragment (Default = 2)')
args = parser.parse_args()


def read_excel_input(nts_alph=args.nts_alphabet_light):
    """
     Produces three dictionaries mod_alphabet, mod_origin and mod_partial to be used by the rest of the program
     mod_alphabet contains all ID : ID_ext couples, thus the one letter and extended codes for each modified base
     mod_origin contains all the unmodified originating species for each modification and is used by the program
     to check if invalid modification positions are given by the user
      mod_partial collects all the partial modified nucleotides for each modification,
      in case the option for partial modifications is selected
     """
    mod_origin, mod_alphabet = {}, {}

    if nts_alph:
        # Checking that the nts_alphabet_light file given in argument exists
        if not os.path.exists(nts_alph):
            print("ERROR! File {} does not exist. Execution terminated without output".format(nts_alph))
            sys.exit(1)

        # Creates a dataframe with info from Excel spreadsheet
        df = pd.read_excel(nts_alph, header=12)

        # Drops rows with the 4 standard nucleobases
        df = df[df.ID != 'C']
        df = df[df.ID != 'G']
        df = df[df.ID != 'A']
        df = df[df.ID != 'U']

        # Drops rows with NaN values
        df = df[pd.notnull(df['ID'])]

        # Transform all ID values in string (so numbers can be used as one letter code for bases)
        df = df.astype({"ID": str})

        # Creates a dictionary with the pair ID : ID_ext as modification alphabet
        mod_alphabet = dict(zip(df.ID, df.ID_ext))

        # Creates a dictionary with the ID : Originating_base for quality control purposes
        mod_origin = dict(zip(df.ID, df.Originating_base))

    return mod_origin, mod_alphabet


def parse_input(modfile):
    """
    Creates a list with all the modifications to consider for discovery mode
    """
    # Checking that the nts_alphabet_light file given in argument exists
    if not os.path.exists(modfile):
        print("ERROR! File {} does not exist. Execution terminated without output".format(modfile))
        sys.exit(1)

    with open(modfile, 'r') as f:
        for line in f:
            yield line.strip().split()


def mod_input(mods, header):
    """
    Creates a dataframe with all the input modification to be included into discovery mode
    """
    df = pd.DataFrame(mods)

    new_header = df.iloc[header]
    df = df[(header + 1):]
    df.columns = new_header

    return df


def make_patterns(seq, mods):
    """
    Outputs all the combinations of base modifications in a given sequence
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
    Generates all the possible combination of base modifications for the input fragments
    
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
