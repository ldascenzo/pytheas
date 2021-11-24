#!/usr/bin/python3

"""
Last update: March 2021
Author: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzoluigi@gmail.com
GitHub project repository: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Second step of the Pytheas in silico digest library generation. RNA nucleolytic fragments obtained from the in silico
digestion are edited with user-specified chemical modifications.
The nucleotide chemical dictionary is defined for modified and unmodified nucleotides from the input alphabet file.
The standard alphabet file, "nts_light.xslx" with standard nucleotides, common modifications and IDs is
provided. More info on the formatting and usage of the alphabet file found in the Digest section of the Pytheas manual.
A template of the modification file, "modfile_template.txt", is provided. More info on the formatting and usage of the
modfile can be found in the Digest section of the Pytheas manual.

***OUTPUT***
1) output.2 -> same format of output.1 with added info on modified nucleotides.
2) nts_light.csv -> nucleotide alphabet in csv format

"""

import itertools
import os
import sys
import pandas as pd


class Modifications:
    def __init__(self, mod_profile, nts_alphabet_light):
        self.modfile, self.alphabet = mod_profile, nts_alphabet_light

    def profile_mod(self):
        """
        Modifications from the mod_profile input file are extracted and put in a list
        """
        # Checking that the modification file given in argument exists
        if self.modfile:
            if not os.path.exists(self.modfile):
                print("ERROR! File {} does not exist. Execution terminated without output".format(self.modfile))
                sys.exit(1)

            profile_mod = []

            if self.modfile:
                with open(self.modfile, 'r') as infile:
                    for line in infile:
                        if line[0] != '\n' and line[0] != '#':
                            if line.split()[1].isdigit():
                                profile_mod.append(line.replace('\n', ''))

            return profile_mod

    def read_excel_input(self):
        """
        Create three dictionaries mod_alphabet, mod_origin and mod_partial
        * mod_alphabet contains all ID : ID_ext couples, one letter and extended IDs or all modified base
        * mod_origin contains all the unmodified nucleotides and is used for modification positions validation
        * mod_partial contains all the partial modified nucleotides for each modification (as from the alphabet file)
        in case the option for partial modifications is selected
        """
        # Checking that the nts_alphabet_light file given in argument exists
        if not os.path.exists(self.alphabet):
            print("ERROR! File {} does not exist. Execution terminated without output".format(self.alphabet))
            sys.exit(1)

        # Create a dataframe with info from the input alphabet file
        df = pd.read_excel(self.alphabet, header=12)
        out_alph = df

        # Drop the first 4 rows with the standard nucleotides
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

        red_df = df[pd.notnull(df['Partial_modifications'])]

        # Create a dictionary to add the partial modified nucleotides
        mod_partial = dict(zip(red_df.ID, red_df.Partial_modifications))

        # Assign multiple partial modifications for a base to a list
        for key, value in mod_partial.items():
            mod_partial[key] = list(value.split(','))

        return mod_alphabet, mod_origin, mod_partial, out_alph

    def modification_check(self, mod_origin):
        """
        Control check that the modification positions from the modification file are correct
        (e.g.: mG should be G in fasta and not A/U/C)
        """
        if not os.path.exists("./seq_output"):
            print("ERROR! Problem reading the file seq_output, quality check for modification will not be performed")

        else:
            seq_list = fasta_sequences(open("./seq_output", 'r'))
            for seq in seq_list:
                for mod in self.profile_mod():
                    if seq.split()[0] == mod.split()[0]:
                        if mod_origin[mod.split()[2]] != list(seq.split()[1])[int(mod.split()[1]) - 1]:
                            print(
                                "WARNING!! Modification '{}' in input modfile does not correspond to the right "
                                "unmodified nucleotide in the input fasta sequence (it is {} instead of {}). "
                                "Specified position will be modified anyway, please check the input "
                                "files.".format(mod, list(seq.split()[1])[int(mod.split()[1]) - 1],
                                                mod_origin[mod.split()[2]]))

    def mod_0_1_2_mode(self, input_lines, mod_alphabet):
        """
        Add modification lines for nucleotides.
        Copies of the unmodified and modified sequences are appended if requested
        """
        # Add modified bases with the one-letter and extended IDs
        with open(self.modfile, 'r') as infile:
            for modline in infile:

                if modline.split()[-1].isdigit():
                    # Check if the modification is "on" (1 or 2 from the modification file)
                    if int(modline.split()[-1]) > 0:
                        mod_lines, positions, del_lines = [], [], []

                        for i, line in enumerate(input_lines):
                            if int(line.split()[2]) <= int(modline.split()[1]) <= int(
                                    line.split()[3]) and modline.split()[0] == line.split()[0]:

                                # Implement only the modified line if modification is not requested, deleting the
                                # non-modified entry line
                                if modline.split()[-1] == "1":
                                    split = line.split()
                                    s = list(split[1])
                                    s[int(modline.split()[1]) - int(split[2])] = modline.split()[2]
                                    split[1] = "".join(s)

                                    # Tracking trick: add @s at the end of the line for each modified base in the
                                    # fragment
                                    if "@" in split[-1]:
                                        s = split[1]

                                        for a in mod_alphabet:
                                            s = s.replace(a, mod_alphabet[a])
                                        split[-1] = s + "@" * split[-1].count('@') + "@"

                                    else:
                                        s[int(modline.split()[1]) - int(split[2])] = mod_alphabet[modline.split()[2]]
                                        split[-1] = split[-1] + " " + "".join(s) + "@"

                                    del_lines.append(line)
                                    mod_lines.append(" ".join(split)), positions.append(i)

                                # Implement both modified and non-modified lines if requested
                                elif modline.split()[-1] == "2":
                                    split = line.split()
                                    s = list(split[1])
                                    s[int(modline.split()[1]) - int(split[2])] = modline.split()[2]
                                    split[1] = "".join(s)

                                    # Tracking trick: add @s at the end of the line for each modified base in the
                                    # fragment
                                    if "@" in split[-1]:
                                        s = list(split[1])
                                        for i, p in enumerate(s):
                                            if p in mod_alphabet.keys():
                                                s[i] = mod_alphabet[p]

                                        split[-1] = "".join(s) + "@" * split[-1].count('@') + "@"

                                    else:
                                        s[int(modline.split()[1]) - int(split[2])] = mod_alphabet[modline.split()[2]]
                                        split[-1] = split[-1] + " " + "".join(s) + "@"

                                    mod_lines.append(" ".join(split)), positions.append(i + 1)

                                # Raise an error if an invalid value is specified for the modification option
                                # (value different from 0,1,2)
                                else:
                                    print(
                                        "ERROR! Invalid line {} in the modification file. "
                                        "Execution terminated without output".format(modline))
                                    sys.exit(1)

                        assert (len(mod_lines) == len(positions))
                        acc = 0

                        # Add lines with modifications
                        for i in range(len(mod_lines)):
                            input_lines.insert(positions[i] + acc, mod_lines[i])
                            acc += 1

                        # Delete unmodified lines when option is 1
                        for d in del_lines:
                            input_lines.remove(d)

        return input_lines

    def output_alphabet(self, outfile):
        """
        Transform the input dictionary in csv format for output
        """
        self.read_excel_input()[3].to_csv(outfile)

    def output_file(self, lines, outfile):
        """
        Prepare the lines for the output file
        """
        final_lines = []

        if self.modfile:
            # Add info on the file header
            final_lines.append("#MODIFICATIONS_PROFILE {}\n".format(self.modfile))

        # Copy information in the header of the input file output.1
        with open(outfile, 'r') as infile:
            for line in infile:
                if line[0] == "#":
                    final_lines.append(line)

                else:
                    if not line.split()[2].isdigit():
                        final_lines.append(line[:-1] + " mod\n")

        # Prepare the final lines to be written, stripping the @ added as modification flags
        for line in lines:
            if "@" in line:
                final_lines.append(line.rstrip('@') + "\n")

            else:
                # Add a dash to the ID column when the fragment does not contain modifications
                final_lines.append(line + " -\n")

        return final_lines

    def final_output(self):
        # Quality check of the modifications with the given modification file/fasta sequences
        if self.modfile:
            self.modification_check(self.read_excel_input()[1])

        lines = []

        # Create a list with the lines from output.1 file
        with open("./output.1", 'r') as infile:
            for l in infile:
                if l[0] != "#":
                    if l.split()[2].isdigit():
                        lines.append(l.rstrip())

        if self.modfile:
            input_lines = self.mod_0_1_2_mode(lines, self.read_excel_input()[0])
        else:
            input_lines = lines

        # Write the lines in the output.1 file
        open("./output.2", 'w').writelines(
            unique_list(self.output_file(mod_nts_exceptions(input_lines, self.read_excel_input()[2],
                                                            self.read_excel_input()[0]), "./output.1")))

        # Write the nts light alphabet in nts_light.csv
        self.output_alphabet('nts_light.csv')


def mod_nts_exceptions(input_lines, mod_partial, mod_alphabet):
    """
    Append additional lines for multiple modification states (e.g mmA can be A, mA or mmA)
    """
    mod_lines, positions = [], []

    for i, line in enumerate(input_lines):
        # If only one residue is modified the combinations are only the number of possible modified states
        # of the residue
        if line.count("@") == 1:

            for nt in mod_partial.keys():
                if nt in line.split()[1]:
                    for x in mod_partial[nt]:
                        s, t = line.split()[1].replace(nt, x), line.split()[-1].replace(mod_alphabet[nt],
                                                                                        mod_alphabet[x])
                        mod_lines.append('{} {} {} {}'.format(line.split()[0], s, " ".join(line.split()[2:7]), t))
                        positions.append(i + 1)

        # If more than one residue has modifications the combination of possible fragment is more complex
        # and requires further operations
        elif line.count('@') > 1:
            modifs, pos = [], []

            # Identify the positions with partial modifications to be added
            for a in mod_partial:
                if a in line.split()[1]:
                    for n, ch in enumerate(line.split()[1]):
                        if ch == a:
                            modifs.append(a), pos.append(n)

            new_modifs = []

            # Create the list of the modifications for the cartesian product
            for mod in modifs:
                new_modifs.append(mod + "".join(mod_partial[mod]))

            # Cartesian product to compute all possible modification combinations
            permut_mod = itertools.product(*new_modifs)

            # Add the info on modifications with extended ID for human readability
            for x in permut_mod:

                string = list(line.split()[1])

                for n, t in enumerate(x):
                    string[pos[n]] = t
                string2, string1 = string, "".join(string)

                for y, s in enumerate(string):
                    if s in mod_alphabet:
                        string2[y] = mod_alphabet[s]

                mod_lines.append(
                    line.split()[0] + " " + string1 + " " + " ".join(line.split()[2:7]) + " " + "".join(
                        string2) + "@" *
                    line.split()[-1].count('@')), positions.append(i + 1)

    assert (len(mod_lines) == len(positions))
    acc = 0

    # Add lines with extended modification possibilities to the list
    for i in range(len(mod_lines)):
        input_lines.insert(positions[i] + acc, mod_lines[i])
        acc += 1

    return input_lines


def unique_list(seq):
    """
    Delete redundant entrances within a list, keeping the original order
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def fasta_sequences(fastafile):
    """
    Creates a list with the fasta sequences for all the RNA molecules
    """
    seq_lines = []

    for line in fastafile:
        seq_lines.append(line)

    return seq_lines
