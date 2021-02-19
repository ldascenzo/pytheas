#!/usr/bin/python3

"""
Last update: May 2020
Credits to: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzo@scripps.edu

***DESCRIPTION***
This script performs an in silico digest of given RNA sequence(s), using an enzyme given as input and
allowing a number of consecutive missed cleavages up to 4. Information on the 3'-end chemistry of the cleaved RNA
fragment is also required. Information on the 5'-end of cleaved RNA fragments and of the 3' and 5'ends of the whole RNA
molecules are optional.
The output "output.1" file contains all the RNA fragments with info on their sequence,
starting/ending nucleotides, number of missed cleavages and 3'/5'-end chemistry of the cleaved fragments.

***USAGE***
python 1_enzyme.py --OPTIONS

***OPTIONS***
--RNA_sequences (REQUIRED) -> Input RNA sequence file(s). Please use fasta format and input the file names without "="
                              after the option. First string of the header for each sequence will be used as id'
--enzyme (REQUIRED) -> type of nuclease utilized:  A (cuts after C or U),T1 (cuts after G), U2 (cuts after A,G), 'none'
                       (no RNA digestion), 'nonspecific' and 'Cus'
--miss (OPTIONAL, DEFAULT = 0) -> number of possible consecutive missed cleavages

--nonspecific_min_length (OPTIONAL, DEFAULT = 3) -> Maximum length for the oligos obtained from nonspecific cleavage.
                                                Default = 3
--nonspecific_max_length (OPTIONAL, DEFAULT = 10) -> Maximum length for the oligos obtained from nonspecific cleavage.
                                                Default = 10

--cleaved_fragments_5end_chem (OPTIONAL, DEFAULT='OH') -> Set the 5' chemistry of the RNA fragments cleaved from the
                      chosen endonuclease. Possible choices=['OH', 'P']. Input of values spaced
                      without '=' nor commas
--cleaved_fragments_3end_chem (REQUIRED, DEFAULT='P') -> Set the 3' chemistry of the RNA fragments cleaved from the
                      chosen endonuclease. Ppossible choices=['OH', 'P', 'cP', 'CP']. Input of values spaced
                      without '=' nor commas

--RNA_5end_chem (OPTIONAL, DEFAULT='P') -> Set the 5' chemistry of the input RNA molecule(s) [options are 'OH' or 'P'].
                      Note that this option refers to the whole RNA molecule and not to the fragments
                      generated by endonuclease cleavage. Please Input values spaced without '=' nor commas
--RNA_3end_chem (OPTIONAL, DEFAULT='OH') -> Set the 3' chemistry of the input RNA molecule(s) [options are 'OH',
                      'P' or 'CP' ]. Note that this option refers to the whole RNA molecule and not to the fragments
                      generated by endonuclease cleavage. Please Input values spaced without '=' nor commas

***INPUT***
RNA sequences files in fasta format

***OUTPUT***
output.1 - text file in the script directory where fragments are listed in lines with header: Molecule Seq Nstart Nend Miss
seq_output - accessory output file with the fasta sequences of the given RNAs, used to validate modification sites

"""

import argparse
import os
import re
import sys
from Bio import SeqIO

# Initialize and define launch options
parser = argparse.ArgumentParser(description='List of available options')
parser.add_argument('--RNA_sequences', nargs='*', required=True,
                    help='Input RNA sequence file(s). Please use fasta format and input the file names without '
                         '"=" after the option. First string of the header for each sequence will be used as id')
parser.add_argument('--enzyme', choices=['A', 'T1', 'U2', 'none', 'nonspecific', 'Cus'], required=True,
                    help='Nuclease enzyme used for digestion')
parser.add_argument('--miss', type=int, choices=[0, 1, 2, 3, 4], default=0,
                    help='Number of possible consecutive miss cleavages to consider (Max 4), up to given value')
parser.add_argument('--nonspecific_min_length', type=int, default=3,
                    help='Minimum length for the oligos obtained from nonspecific cleavage. Default = 3')
parser.add_argument('--nonspecific_max_length', type=int, default=10,
                    help='Maximum length for the oligos obtained from nonspecific cleavage. Default = 10')

parser.add_argument("--cleaved_fragments_5end_chem", nargs='*', default=['OH'], choices=['OH', 'P'],
                    help="Set the 5' chemistry of the RNA fragments cleaved from the chosen endonuclease to be 'OH' or "
                         "'P' . Input of values spaced without '=' nor commas (default = OH)")
parser.add_argument("--cleaved_fragments_3end_chem", nargs='*', required=True, default=['P'], choices=['OH', 'P', 'cP',
                                                                                                       'CP'],
                    help="Set the 3' chemistry of the RNA fragments cleaved from the chosen endonuclease to be 'OH', "
                         "'P' or 'CP'. Input of values spaced without '=' nor commas (default = P)")

parser.add_argument("--RNA_5end_chem", nargs='*', default=['P'], choices=['OH', 'P'],
                    help="Set the 5' chemistry of the input RNA molecule(s) [options are 'OH' or 'P' ]. "
                         "Note that this option refers to the whole RNA molecule and not to oligomers after digestion "
                         "Input of values spaced without '=' nor commas (default = P)")
parser.add_argument("--RNA_3end_chem", nargs='*', default=['OH'], choices=['OH', 'P', 'cP', 'CP'],
                    help="Set the 3' chemistry of the input RNA molecule(s) [options are 'OH', 'P' or 'CP']. "
                         "Note that this option refers to the whole RNA molecule and not to oligomers after digestion "
                         "Input of values spaced without '=' nor commas (default = OH)")
args = parser.parse_args()


def enzyme_cut(enzyme):
    """
    Defines endonuclease cutting sites (cleavage on 3' of the given nucleotide)
    """
    if enzyme == 'A':
        cut_nts = ['C', 'U']

    elif enzyme == 'T1':
        cut_nts = ['G']

    elif enzyme == 'U2':
        cut_nts = ['A', 'G']

    elif enzyme == 'Cus':
        cut_nts = ['C']

    elif enzyme == 'none':
        cut_nts = []

    return cut_nts


def print_ReSites(id, sequence, cut_nts):
    """
    Extracts all RNA cleavages fragments for the nuclease. Fragments are given with the format:
    Molecule Sequence Nstart Nend 5'chem 3'chem
    """
    output_lines = []

    if len(cut_nts) == 1:
        # Regular expression matching all occurrences of the given base except if at the end of the sequence
        pattern = r"{0}(?!$)".format(cut_nts[0])

    elif len(cut_nts) == 2:
        # Regular expression matching all occurrences of the given base except if at the end of the sequence
        pattern = r"({0}|{1})(?!$)".format(cut_nts[0], cut_nts[1])

    elif len(cut_nts) == 0:
        # Regular expression matching no nucleotide in case of no enzymatic digestion
        pattern = 'X1X'

    sites = [str(m.start()) for m in re.finditer(pattern, sequence)]

    if sites:
        # Manually adding the first fragment to the fragments list
        output_lines.append("{} {} {} {} {}\n".format(id, sequence[:int(sites[0]) + 1], str(1), str(int(sites[0]) + 1),
                                                      str(0)))

    else:
        # Adding only one fragment in the case of no enzymatic digestion
        output_lines.append("{} {} {} {} {}\n".format(id, sequence, str(1), str(len(sequence)), str(
            0)))

    sites.append(str(len(sequence)))

    # Cycle to add all the remaining fragments
    for start, end in zip(sites, sites[1:]):

        if int(end) < len(sequence):
            output_lines.append(
                "{} {} {} {} {}\n".format(id, sequence[int(start) + 1:int(end) + 1], str(int(start) + 2),
                                          str(int(end) + 1), str(0)))

        elif int(end) == len(sequence):
            output_lines.append(
                "{} {} {} {} {}\n".format(id, sequence[int(start) + 1:int(end) + 1], str(int(start) + 2), str(int(end)),
                                          str(0)))

    return output_lines


def miss_1(input_list):
    """
    Adds the fragments with 1 missed cleavage
    """
    output_list = []

    for x, y in zip(input_list, input_list[1:]):
        output_list.append(
            "{} {}{} {} {} {}\n".format(x.split()[0], x.split()[1], y.split()[1], x.split()[2], y.split()[3], str(1)))

    return input_list + output_list


def miss_2(input_list):
    """
    Adds the fragments with 2 missed cleavages
    """
    output_list = []

    for x, y, z in zip(input_list, input_list[1:], input_list[2:]):
        if z.split()[-1][0] == "0":
            output_list.append(
                "{} {}{}{} {} {} {}\n".format(x.split()[0], x.split()[1], y.split()[1], z.split()[1], x.split()[2],
                                              z.split()[3], str(2)))

    return input_list + output_list


def miss_3(input_list):
    """
    Adds the fragments with 3 missed cleavages
    """
    output_list = []

    for x, y, z, w in zip(input_list, input_list[1:], input_list[2:], input_list[3:]):
        if w.split()[-1][0] == "0":
            output_list.append(
                "{} {}{}{}{} {} {} {}\n".format(x.split()[0], x.split()[1], y.split()[1], z.split()[1], w.split()[1],
                                                x.split()[2], w.split()[3], str(3)))

    return input_list + output_list


def miss_4(input_list):
    """
    Adds the fragments with 4 missed cleavages
    """
    output_list = []

    for x, y, z, w, v in zip(input_list, input_list[1:], input_list[2:], input_list[3:], input_list[4:]):
        if v.split()[-1][0] == "0":
            output_list.append(
                "{} {}{}{}{}{} {} {} {}\n".format(x.split()[0], x.split()[1], y.split()[1], z.split()[1], w.split()[1],
                                                  v.split()[1], x.split()[2], v.split()[3], str(4)))

    return input_list + output_list


def nonspecific(rna_id, sequence, min_length, max_length):
    """
    Obtains all the fragment sequences in case of nonspecific cleavage
    """
    output_sequences, seq_list = [], list(sequence)

    for i in range(min_length, max_length + 1):

        if i <= len(sequence):
            for position in range(0, len(sequence) - i + 1):
                seq_to_add = ''.join(seq_list[position:position + i])
                output_sequences.append(
                    "{} {} {} {} {}\n".format(rna_id, seq_to_add, position + 1, position + i, 0))

    return output_sequences


def generate_output():
    """
    Generates the final output file with lines in the format: Molecule Sequence Nstart Nend Miss 5'end 3'end
    """
    final_lines, seq_output, unique_ids = [], [], []
    sequences_dictionary = {}

    # Cycles through the fasta files given as input
    for fasta_file in args.RNA_sequences:

        if fasta_file[0].isdigit() or fasta_file[0].isalpha():
            with open(os.getcwd() + "/" + fasta_file.rstrip(), 'r') as handle:

                # Extract and process the sequences within fasta files
                for seq in SeqIO.parse(handle, "fasta"):

                    if seq.id in unique_ids:
                        print(
                            "ERROR!! The molecule id {} is used to identify multiple molecules. Please edit the "
                            "molecules id in the fasta file to be uniques. Execution terminated without output".format(
                                seq.id))
                        sys.exit(1)

                    else:
                        unique_ids.append(seq.id)
                        sequence = str(seq.seq.ungap("-"))
                        sequences_dictionary[str(seq.id)] = sequence
                        seq_output.append(str(seq.id) + " " + sequence + "\n")

                        # Adds the nonspecific cleavage lines or the missed cleavages lines
                        if args.enzyme == 'nonspecific':
                            final_lines = final_lines + nonspecific(str(seq.id), sequence, args.nonspec_min_length,
                                                                    args.nonspec_max_length)

                        else:

                            if args.miss == 0:
                                final_lines = final_lines + print_ReSites(str(seq.id), sequence,
                                                                          enzyme_cut(args.enzyme))

                            elif args.miss == 1:
                                final_lines = final_lines + miss_1(
                                    print_ReSites(str(seq.id), sequence, enzyme_cut(args.enzyme)))

                            elif args.miss == 2:
                                final_lines = final_lines + miss_2(
                                    miss_1(print_ReSites(str(seq.id), sequence, enzyme_cut(args.enzyme))))

                            elif args.miss == 3:
                                final_lines = final_lines + miss_3(
                                    miss_2(miss_1(print_ReSites(str(seq.id), sequence, enzyme_cut(args.enzyme)))))

                            elif args.miss == 4:
                                final_lines = final_lines + miss_4(miss_3(
                                    miss_2(miss_1(print_ReSites(str(seq.id), sequence, enzyme_cut(args.enzyme))))))

    open(os.getcwd() + "/seq_output", 'w').writelines(seq_output)

    # Add the information about 3' and 5' chemistry to the cleaved fragments
    output_lines = []
    for line in final_lines:
        output_lines.extend(assign_chemistry(line, sequences_dictionary))

    return output_lines


def assign_chemistry(fragment, sequences):
    """
    Assigns the 3' and 5' chemistry to the cleaved fragments
    """
    values = fragment.split()
    outlines = []

    if values[2] == '1':
        for end5 in args.RNA_5end_chem:
            if int(values[3]) == len(sequences[values[0]]):
                for end3 in args.RNA_3end_chem:
                    outlines.append("{} {} {}\n".format(fragment.rstrip(), end5, end3))
            else:
                for end3 in args.cleaved_fragments_3end_chem:
                    outlines.append("{} {} {}\n".format(fragment.rstrip(), end5, end3))
    else:
        for end5 in args.cleaved_fragments_5end_chem:
            if int(values[3]) == len(sequences[values[0]]):
                for end3 in args.RNA_3end_chem:
                    outlines.append("{} {} {}\n".format(fragment.rstrip(), end5, end3))
            else:
                for end3 in args.cleaved_fragments_3end_chem:
                    outlines.append("{} {} {}\n".format(fragment.rstrip(), end5, end3))

    return outlines


if __name__ == "__main__":
    """
    Writes the output file output.1
    """

    # The header depends on the choice of nonspecific cleavage or an endonuclease by the user
    if args.enzyme == 'nonspecific':
        starting_lines = ["#ENZYME {}\n#NONSPECIFIC_MIN_LENGTH {}\n#NONSPECIFIC_MAX_LENGTH {}"
                          "\n#CLEAVED_FRAGMENTS_5'CHEMISTRY"
                          " {}\n#CLEAVED_FRAGMENTS_3'CHEMISTRY {}\n#WHOLE_RNA_5'CHEMISTRY {}\n#WHOLE_RNA_3'CHEMISTRY "
                          "{}\nMolecule Seq Nstart Nend Miss 5'chem 3'chem\n".format(args.enzyme,
                                                                                     args.nonspec_min_length,
                                                                                     args.nonspec_max_length,
                                                                                     ','.join(
                                                                                         args.cleaved_fragments_5end_chem),
                                                                                     ','.join(
                                                                                         args.cleaved_fragments_3end_chem),
                                                                                     ','.join(args.RNA_5end_chem),
                                                                                     ','.join(args.RNA_3end_chem))]
    else:
        starting_lines = ["#ENZYME {}\n#MISSED_CLEAVAGES {}\n#CLEAVED_FRAGMENTS_5'CHEMISTRY {}\n"
                          "#CLEAVED_FRAGMENTS_3'CHEMISTRY {}\n#WHOLE_RNA_5'CHEMISTRY {}\n#WHOLE_RNA_3'CHEMISTRY {}\n"
                          "Molecule Seq Nstart Nend Miss 5'chem 3'chem\n".format(args.enzyme, args.miss,
                                                                                 ','.join(
                                                                                     args.cleaved_fragments_5end_chem),
                                                                                 ','.join(
                                                                                     args.cleaved_fragments_3end_chem),
                                                                                 ','.join(args.RNA_5end_chem),
                                                                                 ','.join(args.RNA_3end_chem))]

    open(os.getcwd() + "/output.1", 'w').writelines(starting_lines + generate_output())

    print("Done! Output file(s) -> output.1 seq_output")