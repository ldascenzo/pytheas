#!/usr/bin/python3

"""
Last update: March 2021
Author: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzoluigi@gmail.com
GitHub project repository: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
First step of the Pytheas in silico digest library generation. Given RNA sequences(s) in fasta format are cleaved
with a specific RNA endonuclease enzyme (or not if the user choses so). Additional parameters such as missed cleavages
(up to 4), 3' and 5' chemistry [P, OH or cP] for the nucleolytic fragments and the whole RNA molecule are optional.
If non-specific nucleoltic cleavage is requested, the values of minimum and maximum nucleolytic sequences length
has to be specified.

***OPTIONS***
--RNA_sequences (REQUIRED) -> Input RNA sequence file(s) in fasta format. First string of the header for each sequence
                              will be used as sequence id. NOTE: file names input without "=" nor comma.
--enzyme (REQUIRED) -> RNA endonuclease selected for the digestion: options are A (cleaves after C or U),T1 (cleaves
                        after G), U2 (cleaves after A or G), Cus' (cuts after C), 'none' (no RNA digestion) and
                        'nonspecific' (cleaves after every nucleotide, generating sequences of given length).
--miss (OPTIONAL, DEFAULT = 0) -> number of possible consecutive missed cleavages (up to).
--nonspecific_min_length (OPTIONAL, DEFAULT = 3) -> Maximum length for the nucleolytic fragments obtained if
                                                    nonspecific cleavage is selected.
--nonspecific_max_length (OPTIONAL, DEFAULT = 10) -> Maximum length for the nucleolytic fragments obtained if
                                                    nonspecific cleavage is selected.
--cleaved_fragments_5end_chem (OPTIONAL, DEFAULT='OH') -> 5' chemistry of the RNA nucleolytic fragments
                            Possible choices=['OH', 'P']. NOTE: parameters input without "=" nor comma.
--cleaved_fragments_3end_chem (OPTIONAL, DEFAULT='P') -> 3' chemistry of the RNA nucleolytic fragments
                            Possible choices=['OH', 'P', 'cP']. NOTE: parameters input without "=" nor comma.
--RNA_5end_chem (OPTIONAL, DEFAULT='P') -> 5' chemistry of the input RNA molecule(s) ['OH' or 'P'].
                      NOTE: this refers to the whole RNA molecule and not to the nucleolytic fragments. Parameters
                      input without "=" nor comma.
--RNA_3end_chem (OPTIONAL, DEFAULT='OH') -> 3' chemistry of the input RNA molecule(s) ['OH' or 'P'].
                      NOTE: this refers to the whole RNA molecule and not to the nucleolytic fragments. Parameters
                      input without "=" nor comma.

***OUTPUT***
1) output.1 file contains all the RNA nucleolytic fragments with info on their sequence, numbering referred to the input
RNA sequence, number of missed cleavages and 3'/5'-end chemistry.
2) seq_output containing a compact version of the input RNA sequences, used for modification sites validation in later
steps of the in silico digestion generation

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
parser.add_argument("--cleaved_fragments_3end_chem", nargs='*', default=['P'], choices=['OH', 'P', 'cP'],
                    help="Set the 3' chemistry of the RNA fragments cleaved from the chosen endonuclease to be 'OH', "
                         "'P' or 'cP'. Input of values spaced without '=' nor commas (default = P)")

parser.add_argument("--RNA_5end_chem", nargs='*', default=['P'], choices=['OH', 'P'],
                    help="Set the 5' chemistry of the input RNA molecule(s) [options are 'OH' or 'P' ]. "
                         "Note that this option refers to the whole RNA molecule and not to oligomers after digestion "
                         "Input of values spaced without '=' nor commas (default = P)")
parser.add_argument("--RNA_3end_chem", nargs='*', default=['OH'], choices=['OH', 'P', 'cP'],
                    help="Set the 3' chemistry of the input RNA molecule(s) [options are 'OH', 'P' or 'cP']. "
                         "Note that this option refers to the whole RNA molecule and not to oligomers after digestion "
                         "Input of values spaced without '=' nor commas (default = OH)")
args = parser.parse_args()


def enzyme_cut(enzyme):
    """
    Definition of the RNA endonuclease cleavages (on 3' of the given nucleotide)
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
    Computation of all the RNA cleavages fragments for the nuclease.
    Fragment sequences are given with the format:
    molecule sequence residue_start residue_end 5'end 3'end
    """
    output_lines = []

    if len(cut_nts) == 1:
        # All occurrences of the given base are matched except if at the end of the sequence
        pattern = r"{0}(?!$)".format(cut_nts[0])

    elif len(cut_nts) == 2:
        # All occurrences of the given base are matched except if at the end of the sequence
        pattern = r"({0}|{1})(?!$)".format(cut_nts[0], cut_nts[1])

    elif len(cut_nts) == 0:
        # No nucleotides are matched in case of no enzymatic digestion
        pattern = 'X1X'

    # List of all the cleavage sites for the given enzyme
    sites = [str(m.start()) for m in re.finditer(pattern, sequence)]

    if sites:
        # Dirty trick: manually adding the first fragment to the fragments list
        output_lines.append("{} {} {} {} {}\n".format(id, sequence[:int(sites[0]) + 1], str(1), str(int(sites[0]) + 1),
                                                      str(0)))

    else:
        # Adding only one fragment in the case of no enzymatic digestion
        output_lines.append("{} {} {} {} {}\n".format(id, sequence, str(1), str(len(sequence)), str(
            0)))

    sites.append(str(len(sequence)))

    # Loop to add all the remaining fragments to the output list
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
    Generate the fragments with 1 missed cleavage
    """
    output_list = []

    for x, y in zip(input_list, input_list[1:]):
        output_list.append(
            "{} {}{} {} {} {}\n".format(x.split()[0], x.split()[1], y.split()[1], x.split()[2], y.split()[3], str(1)))

    return input_list + output_list


def miss_2(input_list):
    """
    Generate the fragments with 2 missed cleavages
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
    Generate the fragments with 3 missed cleavages
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
    Generate the fragments with 4 missed cleavages
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
    Compute all the fragment sequences in case of nonspecific cleavage, based on the info selected by the user
    on minimum and maximum length for the sequences generated from nonspecific cleavage
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
    Final output file generation with lines in the format:
    molecule sequence residue_start residue_end miss 5'end 3'end
    """
    final_lines, seq_output, unique_ids = [], [], []
    sequences_dictionary = {}

    # Loop through the fasta files given as input extracting the sequences
    for fasta_file in args.RNA_sequences:

        if fasta_file[0].isdigit() or fasta_file[0].isalpha():
            with open(os.getcwd() + "/" + fasta_file.rstrip(), 'r') as handle:

                # Extract and process the input fasta sequences
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

                        # Append the nonspecific cleavage lines
                        if args.enzyme == 'nonspecific':
                            final_lines = final_lines + nonspecific(str(seq.id), sequence, args.nonspecific_min_length,
                                                                    args.nonspecific_max_length)

                        # Append the missed cleavages lines based on the selected values
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
    Assign the 3' and 5' chemistry to the cleaved fragments based on the user-specific options
    """
    values = fragment.split()
    outlines = []

    # Add the 5' and 3' chemistry for the starting and ending nucleotides of each RNA sequence
    if values[2] == '1':
        for end5 in args.RNA_5end_chem:
            if int(values[3]) == len(sequences[values[0]]):
                for end3 in args.RNA_3end_chem:
                    outlines.append("{} {} {}\n".format(fragment.rstrip(), end5, end3))
            else:
                for end3 in args.cleaved_fragments_3end_chem:
                    outlines.append("{} {} {}\n".format(fragment.rstrip(), end5, end3))
    # Add the chemistry info for all the other sequences
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
    Generate the output file output.1
    """

    # The header info differ based on the choice of nonspecific or specific cleavage
    if args.enzyme == 'nonspecific':
        starting_lines = ["#ENZYME {}\n#NONSPECIFIC_MIN_LENGTH {}\n#NONSPECIFIC_MAX_LENGTH {}"
                          "\n#CLEAVED_FRAGMENTS_5'CHEMISTRY"
                          " {}\n#CLEAVED_FRAGMENTS_3'CHEMISTRY {}\n#WHOLE_RNA_5'CHEMISTRY {}\n#WHOLE_RNA_3'CHEMISTRY "
                          "{}\nmolecule sequence residue_start residue_end miss 5'end 3'end\n".format(args.enzyme,
                                                                                     args.nonspecific_min_length,
                                                                                     args.nonspecific_max_length,
                                                                                     ','.join(
                                                                                         args.cleaved_fragments_5end_chem),
                                                                                     ','.join(
                                                                                         args.cleaved_fragments_3end_chem),
                                                                                     ','.join(args.RNA_5end_chem),
                                                                                     ','.join(args.RNA_3end_chem))]
    else:
        starting_lines = ["#ENZYME {}\n#MISSED_CLEAVAGES {}\n#CLEAVED_FRAGMENTS_5'CHEMISTRY {}\n"
                          "#CLEAVED_FRAGMENTS_3'CHEMISTRY {}\n#WHOLE_RNA_5'CHEMISTRY {}\n#WHOLE_RNA_3'CHEMISTRY {}\n"
                          "molecule sequence residue_start residue_end miss 5'end 3'end\n".format(args.enzyme, args.miss,
                                                                                 ','.join(
                                                                                     args.cleaved_fragments_5end_chem),
                                                                                 ','.join(
                                                                                     args.cleaved_fragments_3end_chem),
                                                                                 ','.join(args.RNA_5end_chem),
                                                                                 ','.join(args.RNA_3end_chem))]

    open(os.getcwd() + "/output.1", 'w').writelines(starting_lines + list(set(generate_output())))

    print("Done! Output file(s) -> output.1 seq_output")