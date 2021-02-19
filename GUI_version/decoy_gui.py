#!/usr/bin/python3

"""
Last update: May 2020
Credits to: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzo@scripps.edu

***DESCRIPTION***
This script generates decoy sequences for the theoretical digest RNA fragments obtained after running the
script 3_consolidate.
Decoys are obtained shuffling the nts except for the last of the original sequence,
with up to one per original sequence (zero if all combinations already present).
NOTE that running the script is optional and will only modify the input/output file output.3.MS2 adding the new decoys
"""

import os
import sys
from random import shuffle
import pandas as pd
from shutil import move


class Decoys:
    def __init__(self):
        self.input_file = "output.3.MS2"

    def final_output(self):
        open('./output.3.5', 'w').writelines(output_lines(self.input_file, decoy_lines(self.input_file,
                                                                                       list(sequence_list(
                                                                                           self.input_file)),
                                                                                       read_csv())))
        os.remove('./' + self.input_file)
        move('./output.3.5', './' + self.input_file)


def read_csv(input_csv='./nts_light.csv'):
    """
     Produces a dictionary nts_alphabet
     mod_alphabet contains all ID : ID_ext couples, thus the one letter and extended codes for each nucleobase
     """
    if not os.path.isfile(input_csv):
        print("ERROR! File {} with info on nucleotides from script 2_modify.py is missing."
              "Execution terminated without generating decoys".format(input_csv))
        sys.exit(1)

    else:
        # Read the csv file with the nucleotides dictionary
        df = pd.read_csv(input_csv, usecols=['ID', 'ID_ext'])

        # Drops rows with NaN values
        df = df[pd.notnull(df['ID'])]

        return dict(zip(df.ID, df.ID_ext))


def sequence_list(input_file):
    """
    Creates a list with all the sequences available in the digest (without duplicates)
    """
    # Checks if the input file output.3.MS2 is present in the directory
    if not os.path.isfile("./" + input_file):
        print("ERROR! Input file output.3.MS2 is missing. Execution terminated without generating decoys")
        sys.exit(1)

    # Extracts all the sequences present in the input file
    sequence_list = []
    with open("./" + input_file, 'r') as infile:
        for line in infile:
            if line[0] != '#' and line.split()[0] != 'Seq':
                split = line.split()
                sequence_list.append((split[0], split[3], split[4]))

    return set(sequence_list)


def decoy_lines(input_file, seq_list, nts):
    """
    Creates a list with all the decoy lines to be appended on the original file
    Each decoy line is obtained taking line per line the sequences in the input file
    and randomizing the residues of the positions excluded the last (that remains conserved
    since it is the specific enzymatic cleavage site). The first shuffled sequence + last residue
    that is not present in the input file is appended as decoy, then the following line is considered.
    If no unique sequences are found the line is skipped and no decoy is created
    """
    outlines = []

    # Cycle among the lines of the input file
    with open("./output.3.MS2", 'r') as infile:
        for line in infile:

            # Lines of the header are excluded
            if line[0] != '#' and line.split()[0] != 'Seq' and 'decoy' not in line:
                seq, chem3, chem5 = line.split()[0], line.split()[3], line.split()[4]

                # Only sequences at least 3 nts long are considered to create decoys
                if len(seq) > 2:
                    entry = list(seq[:-1])

                    # Sets the maximum number of shuffled sequences to be tried, being the length of the entry ^ 4
                    for i in range(1, len(entry) ** 4):
                        shuffle(entry)
                        new_seq = ''.join(entry)

                        # If the shuffled sequence is not in the input it is added in a new line
                        # and moves to next input line
                        if (new_seq + seq[-1], chem3, chem5) not in seq_list:

                            seq_list.append((new_seq + seq[-1], chem3, chem5))
                            if line.split()[1] == '-':
                                outlines.append(
                                    "{} {} 1 decoy decoy - -\n".format(new_seq + seq[-1], ' '.join(line.split()[1:5])))
                                break

                            else:
                                final_seq = ""

                                # Add the extended output for modified nucleotides
                                for nt in new_seq + seq[-1]:
                                    final_seq += nts[nt]

                                outlines.append("{} {} {} 1 decoy decoy - -\n".format(new_seq + seq[-1], final_seq,
                                                                                      ' '.join(line.split()[2:5])))

                                break

    return outlines


def output_lines(input_file, decoy_lines):
    """
    Merges the lines of the input and the decoy, ready to be written in the output file
    """
    # Checks if the input file output.3.MS2 is present in the directory
    if not os.path.isfile("./" + input_file):
        print("ERROR! Input file output.3.MS2 is missing. Execution terminated without generating decoys")
        sys.exit(1)

    with open("./output.3.MS2", 'r') as infile:
        input_lines = infile.readlines()

    return input_lines + decoy_lines
