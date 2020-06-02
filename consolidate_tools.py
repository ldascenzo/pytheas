#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 13:35:54 2019

@author: luigi
"""

import os
import sys
from copy import deepcopy
import numpy as np
import pandas as pd

##########
# ELEMENT MASS DICTIONARY
##########
"""
Note that for each element with no specified isotope the mass of the most abundant isotope is used. 
In the case of O, C, H and N the most abundant isotope is >=99% abundant. 
S32 has 95% abundance
Se80 has only 49% abundance, so you could consider adding multiple isotopes if working with Se in case of bad data 
matching
"""
ele_mass = {"C": 12.0000000, "H": 1.007825032, "N": 14.0030740044, "O": 15.9949146196, "S": 31.9720711744,
            "P": 30.9737619984, "Se": 79.9165218, "H2": 2.0141017781, "C13": 13.0033548351, "N15": 15.0001088989,
            "O18": 17.9991596128}

ion_series = ['a', 'b', 'c', 'd', 'w', 'x', 'y', 'z']


def ppm_range(value, difference):
    """
    Calculates incertitude on MS1/MS2 masses equivalent to given ppm value
    """
    return difference * 1000000 / value


def read_excel_input(nts_file):
    """
    Produces a dataframe with all the info on the nucleobases from the input file nts_alphabet_light
    """
    # Checking that the nts_alphabet_light file given in argument exists
    if not os.path.exists(nts_file):
        print("ERROR! File " + nts_file + " does not exist. Execution terminated without generating any output")
        sys.exit(1)

    # Creates a dataframe with info from Excel spreadsheet
    df = pd.read_excel(nts_file, header=12)

    # Drops rows with NaN values
    df = df[pd.notnull(df['ID'])]

    # Transform all ID values in string (so numbers can be used as one letter code for bases)
    df = df.astype({"ID": str})

    return df


def nts_mass(df):
    """
    Creates a dictionary with the masses of all nts based on their atomic composition
    """
    mass_nts, mass_base, mass_backbone = {}, {}, {}

    for index, row in df.iterrows():
        # Calculates mass of the base only
        mass_b = (np.float64(row['C']) * ele_mass["C"] + np.float64(row['O']) * ele_mass["O"] + np.float64(row['H']) *
                  ele_mass["H"] + np.float64(row['N']) * ele_mass["N"] +
                  np.float64(row['P']) * ele_mass["P"] + np.float64(row['S']) * ele_mass["S"] + np.float64(row['Se']) *
                  ele_mass["Se"] + np.float64(row['C13']) * ele_mass["C13"]
                  + np.float64(row['O18']) * ele_mass["O18"] + np.float64(row['N15']) * ele_mass["N15"] + np.float64(
                    row['H2']) * ele_mass["H2"])

        mass_base[row['ID']] = mass_b

        # Calculate mass of the backbone only
        mass_back = (np.float64(row['C.1']) * ele_mass["C"] + np.float64(row['O.1']) * ele_mass["O"] + np.float64(
            row['H.1']) * ele_mass["H"] + np.float64(row['N.1']) * ele_mass["N"] +
                     np.float64(row['P.1']) * ele_mass["P"] + np.float64(row['S.1']) * ele_mass["S"] + np.float64(
                    row['Se.1']) * ele_mass["Se"] + np.float64(row['C13.1']) * ele_mass["C13"]
                     + np.float64(row['O18.1']) * ele_mass["O18"] + np.float64(row['N15.1']) * ele_mass[
                         "N15"] + np.float64(row['H2.1']) * ele_mass["H2"])

        mass_backbone[row['ID']] = mass_back

        # Calculate mass of the whole nucleotide
        mass_nts[row['ID']] = mass_b + mass_back

    return mass_nts, mass_base, mass_backbone


def nucleotides_to_consolidate(alphabet_file, ppm_threshold):
    """
    Check if any of the nucleotides specified on the light or heavy channel are within a set amount of ppm
    and have to be highlighted for redundancy
    """
    dic = nts_mass(read_excel_input(alphabet_file))[0]

    nt_pairs = []

    for key1 in dic.keys():
        for key2 in dic.keys():

            # Reduces the amount of nucleotide pairs in the search space to only the ones that are within the ppm
            # specified range at 5000 m/z
            if ppm_range(2000, abs(dic[key1] - dic[key2])) <= ppm_threshold and key1 != key2 and (
                    key1, key2) not in nt_pairs and (key2, key1) not in nt_pairs:
                nt_pairs.append((key1, key2))

    return nt_pairs


def check_Da_nucleotides(alphabet_file):
    """
    Check if any of the nucleotides specified on the light or heavy channel are within a set amount of ppm
    and have to be highlighted for redundancy
    """
    dic = nts_mass(read_excel_input(alphabet_file))[0]

    nt_pairs_warning, nt_pairs = [], []

    for key1 in dic.keys():
        for key2 in dic.keys():
            if abs(dic[key1] - dic[key2]) < 0.5 and key1 != key2 and (key1, key2) not in nt_pairs and (
                    key2, key1) not in nt_pairs:
                nt_pairs_warning.append((key1, key2, str(round(abs(dic[key1] - dic[key2]), 3)) + " Da"))
                nt_pairs.append((key1, key2))

    print(
        "WARNING!!! The following nucleotide couples from the alphabet file {} have masses within 0.5 Da, "
        "which may be an issue for discrimination of MS2 fragments. "
        "Please run again the script with the --mz_consolidation option set as 'y' and choose "
        "a ppm threshold with --ppm_consolidation for "
        "isobaric/close fragments to be consolidated \n{}\n\n".format(alphabet_file, nt_pairs_warning))


def digest_lines(digest_input):
    """
    Creates a list with all the lines of the Digest file output
    """
    outlines = []
    with open(digest_input, 'r') as infile:
        for line in infile:
            outlines.append(line)
    return outlines


def replaceMultiple(mainString, toBeReplaces, newString):
    """
    Replace a set of multiple sub strings with a new string in main string.
    """
    # Iterate over the strings to be replaced
    for elem in toBeReplaces:
        # Check if string is in the main string
        if elem in mainString:
            # Replace the string
            mainString = mainString.replace(elem, newString)

    return mainString


def Average(lst):
    return sum(lst) / len(lst)


def avg_masses_lines(lines_list):
    """
    Averages the MS1 and all the MS2 masses of ions that will be consolidated
    """
    dic_masses = {}
    for line in lines_list:
        lineSplit = line.split()
        if dic_masses:
            dic_masses['prec'].append(np.float64(lineSplit[0]))
        else:
            dic_masses['prec'] = [np.float64(lineSplit[0])]

        for ms2_ion in lineSplit[13:]:
            ion_series = ms2_ion.split(':')[0]

            if ion_series in dic_masses.keys():
                dic_masses[ion_series].append(np.float64(ms2_ion.split(':')[1]))
            else:
                dic_masses[ion_series] = [np.float64(ms2_ion.split(':')[1])]

    list_output = []
    for key in dic_masses.keys():
        list_output.append(key + ":" + str(Average(dic_masses[key])))

    return list_output


def diff_list(first, second):
    second = set(second)
    return [item for item in first if item not in second]


def check_ppm_threshold(lines_list, MS1_ppm_threshold):
    """
    Keeps in the list of precursor lines to be consolidated only the ones with masses < input threshold
    """
    out_list = []
    for i in range(0, len(lines_list)):

        ref_mass = np.float64(lines_list[i].split()[0])
        for line in lines_list:

            mss = np.float64(line.split()[0])
            if ppm_range(mss, abs(ref_mass - mss)) <= MS1_ppm_threshold:
                out_list.append(line)

    return out_list


def mz_consolidate(alphabet_file, digest_file, channel, MS1_ppm_threshold, MS2_ppm_threshold):
    """
    Consolidates the lines with residues closer in m/z than the ppm window that will be used for MS2 matching. 
    It replaces the competing residues with a X and concatenate the location information of the consolidated fragments
    """
    nt_pairs, digest_list, unique_lines = nucleotides_to_consolidate(alphabet_file, MS1_ppm_threshold), [], []

    # Only if any pair of nucleotides within the mass threshold are found the other operations follow
    if nt_pairs:
        digest_list = digest_lines(digest_file)
        digest_list_out = deepcopy(digest_list)
        for pair in nt_pairs:
            lines_dict = {}

            for line_digest in digest_list:
                if line_digest[0].isdigit() and channel in line_digest:
                    sequence, charge = line_digest.split()[7], line_digest.split()[5]

                    if pair[0] in sequence or pair[1] in sequence:

                        key_sequence = replaceMultiple(sequence, [pair[0], pair[1]], "X")

                        if key_sequence + "_" + charge not in lines_dict.keys():
                            lines_dict["_".join([key_sequence, charge])] = [line_digest]
                        else:
                            lines_dict["_".join([key_sequence, charge])].append(line_digest)

            for key in lines_dict.keys():

                lines = deepcopy(lines_dict[key])

                # Check if the masses of the multiple precursors are within the input ppm threshold 
                to_consolidate = check_ppm_threshold(lines, MS1_ppm_threshold)
                if len(to_consolidate) > 1:
                    for line in lines_dict[key]:
                        if 'decoy' in line.split()[2]:
                            lines.remove(line), digest_list_out.remove(line)

                    for i in range(0, len(lines)):

                        position_molecules, ref_mass = [lines[i].split()[12]], np.float64(lines[i].split()[0])

                        if lines[i].split()[8] != '-':
                            seq_molecules = [lines[i].split()[8]]
                        else:
                            seq_molecules = [lines[i].split()[7]]

                        ms2_ions, newLine_1, newLine_2 = {}, '', ''

                        for x in lines[i].split()[13:]:
                            ms2_ions[x.split(':')[0]] = np.float64(x.split(':')[1])

                        # Creates the string to keep the information on the positions of the fragments
                        # on the sequence when consolidating
                        for line in lines:

                            consolidation_flag = 0
                            if ppm_range(np.float64(line.split()[0]),
                                         abs(ref_mass - np.float64(line.split()[0]))) <= MS1_ppm_threshold and line != \
                                    lines[i]:

                                consolidation_flag = 1
                                # Check MS2 level threshold for consolidation
                                for x in line.split()[13:]:
                                    ion, ms2_mass = x.split(':')[0], np.float64(x.split(':')[1])

                                    if ion in ms2_ions.keys() and ion[0] in ion_series:
                                        if ppm_range(ms2_mass, abs(ms2_mass - ms2_ions[ion])) > MS2_ppm_threshold:
                                            consolidation_flag = 0
                                            break

                            # Creating the human readable sequence with information
                            # on which nucleotides have been consolidated
                            if consolidation_flag == 1:
                                if line.split()[12] not in position_molecules:
                                    position_molecules.append(line.split()[12])

                                if line.split()[8] != '-':
                                    if line.split()[8] not in seq_molecules:
                                        seq_molecules.append(line.split()[8])

                                else:
                                    if line.split()[7] not in seq_molecules:
                                        seq_molecules.append(line.split()[7])

                                # Determine the consensus sequences with X replacing consolidated nts
                                consensus_sequence = find_sequence_with_x(key.split('_')[0], seq_molecules)

                                newLine_1 = " ".join(lines[i].split()[:3]) + " 0 0 " + " ".join(lines[i].split()[5:7]) \
                                            + " " + consensus_sequence

                                if newLine_1 not in unique_lines:
                                    unique_lines.append(newLine_1)

                                newLine_2 = " " + " ".join(lines[i].split()[9:12]) + " " + ";".join(
                                    position_molecules) + " " + " ".join(lines[i].split()[13:])

                                # Eliminate the lines at the end of the cycle through the nt pairs,
                                # in order to consider all the possible
                                if lines[i] in digest_list_out:
                                    digest_list_out.remove(lines[i])

                        if newLine_1 and newLine_2:

                            # Avoid the creation of lines with the same combination of m/z and sequence,
                            # which will be problematic for the matching
                            # Redundancy is solved by adding additional values to the end of the m/z values
                            if newLine_1 in unique_lines:

                                for var in range(0, 9):
                                    if newLine_1 + str(var) not in unique_lines:
                                        digest_list_out.append(newLine_1.split()[0] + str(var) + " " + " ".join(
                                            newLine_1.split()[1:]) + " " + "|".join(seq_molecules) + newLine_2 + "\n")
                                        unique_lines.append(newLine_1 + str(var))
                                        break

                            else:
                                digest_list_out.append(newLine_1 + " " + "|".join(seq_molecules) + newLine_2 + "\n")

    return digest_list_out


def find_sequence_with_x(consensus, single_sequences):
    """
    Determine the correct consensus sequence with X replacing consolidated nucleotides, based on the original
    sequences to be consolidated
    """
    output_sequence = ''

    for position, nt in enumerate(consensus):
        if nt == 'X':
            reference_base, flag = single_sequences[0][position], None
            for seq in single_sequences:
                if seq[position] != reference_base:
                    output_sequence += 'X'
                    flag = 1
                    break

            if not flag:
                output_sequence += reference_base

        else:
            output_sequence += nt

    return output_sequence
