#!/usr/bin/python3

"""
Last update: March 2021
Author: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzoluigi@gmail.com
GitHub project repository: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Third step of the Pytheas in silico digest library generation. Nucleolytic fragments with the same sequences above
a user-defined length threshold are consolidated together, preserving the RNA sequence and residue numbering info.

***OUTPUT***
1) output.3.MS2 (or MS1 if selected by the user) -> collection of all nucleolytic sequences with additional info on
                                multiple sequence information in case of consolidation, otherwise similar to output.2.
"""

import sys
from collections import defaultdict


class Consolidation:
    def __init__(self):
        self.MS_level, self.min_length = 'MS2', 3

    def fragment_list(self):
        """
        Create a list of the output.2 file lines
        """
        output_lines = []

        with open("./output.2", 'r') as infile:
            for line in infile:
                if line[0] != "#":
                    split = line.split()
                    if split[2].isdigit():

                        # Only truly identical fragment longer or equal to min_length are consolidated
                        if self.MS_level[-1] == "2" and len(split[1]) >= self.min_length:
                            output_lines.append(
                                "{} {} {} {} {} {} {},{},{} {} {} {}".format(split[1], split[7], split[4],
                                                                             split[6], split[5], 1,
                                                                             split[0], split[2], split[3],
                                                                             split[0], split[2], split[3]))

                        # Append the original sequence as info in case of MS1 fragments consolidated by masses
                        else:
                            output_lines.append(
                                "{} {} {} {} {} {} {},{},{},{} {} {} {}".format(split[1], split[7], split[4],
                                                                                split[6], split[5], 1,
                                                                                split[0], split[2], split[3],
                                                                                split[1], split[0], split[2],
                                                                                split[3]))

        return output_lines

    def redundant_dic(self):
        """
        Create a dictionary with all the lines with redundant fragments
        """
        output_lines = []
        with open("./output.2", 'r') as infile:
            for line in infile:
                if line[0] != "#":
                    split = line.split()
                    if split[2].isdigit():
                        output_lines.append('{} {} {} {} {} {} {}'.format(split[1], split[7], split[6], split[5],
                                                                          split[2], split[3], split[0]))

        d = defaultdict(list)

        for i, item in enumerate(output_lines):

            # In case of MS2 only fragments >= min_length and with identical sequences
            # are consolidated (cannot be discriminated by MS2)
            if self.MS_level[-1] == "2" and len(item.split()[0]) >= self.min_length:
                d[' '.join(item.split()[:4])].append(i)

            # In case MS1 is selected, consolidate all fragments with equal masses writing alphabetized sequences
            else:
                d["{} {} {}".format("".join(sorted(item.split()[0])), item.split()[2], item.split()[3])].append(i)

        # The dictionary contains info on the redundant fragment coupled with its line number
        d = {k: v for k, v in d.items() if len(v) > 1}

        return d

    def consol(self, fragments, min_length, redundant_dictionary):
        """
        Consolidate fragments and remove redundant lines
        """
        if min_length < 0:
            # Check that the input min_length value is a positive value
            print("ERROR!!! Min length value must be positive!!! ")
            sys.exit(1)

        indices = []
        # Loop through the dictionary containing the redundant fragments/lines
        # Lines with redundant fragments are consolidated together
        for key in redundant_dictionary:
            for i in range(1, len(redundant_dictionary[key])):
                fragments[redundant_dictionary[key][0]] = (
                        " ".join(fragments[redundant_dictionary[key][0]].split()[0:7])
                        + ";" + fragments[redundant_dictionary[key][i]].split()[6] + " "
                        + " ".join(fragments[redundant_dictionary[key][0]].split()[-3:]))
                indices.append(redundant_dictionary[key][i])

        # Delete the redundant occurrences from the fragments list
        final_frag = [i for j, i in enumerate(fragments) if j not in indices]

        # Add the correct number of redundant occurrences per line (when > 1)
        for i, line in enumerate(final_frag):
            if ";" in line:

                unique_mol = []

                # Modify the last three columns of the line in case of redundant fragments
                for x in line.split()[6].split(';'):
                    if x.split(',')[0] not in unique_mol:
                        unique_mol.append(x.split(',')[0])

                # When the redundance is between two fragments of different molecules the output is a - followed by 0 0
                if len(unique_mol) > 1:
                    line_end = str((final_frag[i].split()[6])) + " - 0 0"

                # When the redundance is only between fragments of the same molecule that molecule is specified before
                # the 0 0 closing the line
                else:
                    line_end = "{} {} 0 0".format((final_frag[i].split()[6]), unique_mol[0])

                # Assemble the final lines for MS1 and MS2 output
                if self.MS_level[-1] == "1":
                    final_frag[i] = "{} {} {} {}".format("".join(sorted(final_frag[i].split()[0])),
                                                         " ".join(final_frag[i].split()[1:5]),
                                                         line.count(';') + 1, line_end)

                else:
                    final_frag[i] = "{} {} {} {}".format("".join(final_frag[i].split()[0]),
                                                         " ".join(final_frag[i].split()[1:5]),
                                                         line.count(';') + 1, line_end)

        return final_frag

    def output_file(self, lines, outfile, min_length):
        """
        Prepare lines for the output
        """
        final_lines = []

        # Only in case of MS2, info about the maximum length of fragments specified
        # as parameter is added to the output file header
        if self.MS_level[-1] == "2":
            final_lines.append("#MIN_LENGTH_CONSOLIDATE " + str(min_length) + "\n")

        with open(outfile, 'r') as infile:
            for line in open(outfile, 'r'):
                if line[0] == "#":
                    final_lines.append(line)

                else:
                    if not line.split()[2].isdigit() and line.split()[2] != '-':
                        # Write the header line for the columns of the output
                        final_lines.append("sequence mod miss 3'end 5'end num_copy sequence_location"
                                           " molecule residue_start residue_end\n")

        for line in lines:

            # Output only the fragments longer than min length
            if self.MS_level[-1] == "2":
                if len(line.split()[0]) >= min_length:
                    final_lines.append(line + "\n")
            else:
                final_lines.append(line + "\n")

        return final_lines

    def final_output(self):
        # Writing the lines in the output.3.MSx file
        out_files = []

        if self.MS_level[-1] == "1":
            open("./output.3.MS1", 'w').writelines(
                self.output_file(self.consol(self.fragment_list(), 10000, self.redundant_dic()),
                                 "./output.2", self.min_length))
            out_files.append("output.3.MS1")
        else:
            open("./output.3.MS2", 'w').writelines(
                self.output_file(self.consol(self.fragment_list(), self.min_length, self.redundant_dic()),
                                 "./output.2", self.min_length))
            out_files.append("output.3.MS2")
