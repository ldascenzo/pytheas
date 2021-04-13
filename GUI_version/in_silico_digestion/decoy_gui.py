#!/usr/bin/python3

"""
Last update: March 2021
Author: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzoluigi@gmail.com
GitHub project repository: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Optional step of the Pytheas in silico digest library generation. Decoy sequences are generated (if possible) for every
competing target sequence. Decoy sequences are obtained shuffling the nucleotides except for the last.

***OUTPUT***
-> output.3.MS2 -> the same file used for input, with appended decoy sequences
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
        """
        Generate the output file and replace the output.3.MS2
        """
        open('./output.3.5', 'w').writelines(output_lines(self.input_file, decoy_lines(self.input_file,
                                                                                       list(sequence_list(
                                                                                           self.input_file)),
                                                                                       read_csv())))
        os.remove('./' + self.input_file)
        move('./output.3.5', './' + self.input_file)


def read_csv(input_csv='./nts_light.csv'):
    """
    Create a dictionary containing all ID : ID_ext couples, one letter and extended IDs for all nucleotides
    """
    if not os.path.isfile(input_csv):
        print("ERROR! File {} with info on nucleotides from script 2_modify.py is missing."
              "Execution terminated without generating decoys".format(input_csv))
        sys.exit(1)

    else:
        # Read the csv file with the nucleotides dictionary
        df = pd.read_csv(input_csv, usecols=['ID', 'ID_ext'])

        # Drop rows with NaN values
        df = df[pd.notnull(df['ID'])]

        return dict(zip(df.ID, df.ID_ext))


def sequence_list(input_file):
    """
    Create a list with all the sequences available in the digest (without duplicates)
    """
    # Check if the input file output.3.MS2 is present in the directory
    if not os.path.isfile("./" + input_file):
        print("ERROR! Input file output.3.MS2 is missing. Execution terminated without generating decoys")
        sys.exit(1)

    # Extract all the sequences present in the input file
    sequence_list = []
    with open("./" + input_file, 'r') as infile:
        for line in infile:
            if line[0] != '#' and line.split()[0] != 'sequence':
                split = line.split()
                sequence_list.append((split[0], split[3], split[4]))

    return set(sequence_list)


def decoy_lines(input_file, seq_list, nts):
    """
    Generate a list with all the decoy lines to be appended on the original file
    Each decoy line is obtained taking line per line the sequences in the input file
    and shuffling the residues of the positions excluded the last. The first shuffled sequence + last residue
    that is not present in the input file is appended as decoy, then the following line is considered.
    If no unique sequences are found the line is skipped and no decoy is created
    """
    outlines = []

    # Cycle among the lines of the input file
    with open(input_file, 'r') as infile:
        for line in infile:

            # Lines of the header are excluded
            if line[0] != '#' and line.split()[0] != 'sequence' and 'decoy' not in line:
                seq, chem3, chem5 = line.split()[0], line.split()[3], line.split()[4]

                # Only sequences at least 3 nts long are considered to create decoys
                if len(seq) > 2:
                    entry = list(seq[:-1])

                    # Set the maximum number of shuffled sequences to be tried, being the length of the entry ^ 4
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
    Merge the lines of the input and the decoy, ready to be written in the output file
    """
    # Check if the input file output.3.MS2 is present in the directory
    if not os.path.isfile("./" + input_file):
        print("ERROR! Input file output.3.MS2 is missing. Execution terminated without generating decoys")
        sys.exit(1)

    with open("./output.3.MS2", 'r') as infile:
        input_lines = infile.readlines()

    return input_lines + decoy_lines
