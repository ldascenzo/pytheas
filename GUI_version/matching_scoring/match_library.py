#!/usr/bin/python3

"""
Last update: April 2021
Author: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzoluigi@gmail.com
GitHub project repository: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Core library for the Pytheas matching and scoring algorithm. The matching between the theoretical and the experimental
library is performed with the aid of an empirical scoring function, to estimate the quality of the fit. Additional
output include FDR estimation and csv files containing all the targets and decoys information.
More details about the algorithm, the scoring function and the output can be found within Pytheas paper and the manual.

***OUTPUT***
1) match_output -> txt file containing all the info on the precursor ions alongside all the MS2 ions matched and scored
2) FDR -> csv file with the FDR estimations on the given dataset
3) targets/decoys list -> csv files with all the information on the matched target and decoy sequences
"""

import re
import numpy as np
import itertools
import stats_tools as stats
from datetime import datetime
import ntpath

# Hydrogen and neutron masses
H_mass, neutron_mass = 1.007825032, 1.008665


class Match:
    def __init__(self, digest, mgf, isotope, MS1_mzmin, MS1_mzmax, MS2_mzmin, MS2_mzmax, MS1_ppm, MS2_ppm,
                 MS1_offset_ppm, MS2_offset_ppm, MS2_intensity_min, MS2_peaks_max, precursor_removal,
                 loss_removal, beta, alpha, MS2_cutoff_normalized_int, isotopologues, charges_mgf, FDR_isotopes,
                 targets_with_decoys):
        self.digest_file, self.mgf_file, self.light_heavy = digest, mgf, isotope
        self.MS1_mzlow, self.MS1_mzhigh, self.MS2_mzlow, self.MS2_mzhigh = MS1_mzmin, MS1_mzmax, MS2_mzmin, MS2_mzmax
        self.MS1_ppm, self.MS2_ppm, self.MS1_ppm_offset = round(MS1_ppm, 1), round(MS2_ppm, 1), -round(MS1_offset_ppm, 1)
        self.MS2_ppm_offset = -round(MS2_offset_ppm, 1)
        self.MS2_peak_int_min, self.MS2_peak_num_max = MS2_intensity_min, MS2_peaks_max
        self.precursor_window_removal, self.losses_window_removal = precursor_removal, loss_removal
        self.beta_increment, self.alpha, self.MS2_normint_cutoff = beta, alpha, MS2_cutoff_normalized_int
        self.precursor_isotopologues, self.ignore_charges_mgf = isotopologues, charges_mgf
        self.FDR_light_heavy, self.targets_without_decoys = FDR_isotopes, targets_with_decoys

    def threshold_MS2_int(self):
        """
        Define minimum intensity threshold for MS2 ions matching
        """
        if self.MS2_peak_int_min == 'None':
            threshold = 0
        else:
            threshold = int(self.MS2_peak_int_min)

        return threshold

    def normalize_int_MS2(self, mgf_peaks, match_peaks, digest):
        """
        Normalize the intensity of the matching ions where 100 is the most intense ion excluding M-xx series +
        their first two isotopic peaks
        Only the ion with normalized intensities >= MS2_normint_cutoff are considered for matching
        """
        output_dic = match_peaks

        for key in match_peaks:
            if match_peaks[key]:

                for t in match_peaks[key]:
                    max_int = 0

                    # Loop among all ions in the mgf ordered by intensity
                    for f in mgf_peaks[key]:

                        flag = 0
                        # Loops among the M-xx ions to check the first non M-xx to be 100 in normalization step
                        for peak in digest[t][-1]:

                            # Add the offset mass to the mz_digest
                            mz_digest, mz_mgf, charge = np.float64(peak.split(":")[-1]) + ppm_range(
                                np.float64(peak.split(":")[-1]), self.MS2_ppm_offset), np.float64(f), int(
                                re.findall(r'\d+', peak.split("(")[1].split(")")[0])[0])

                            # Check if the mz of the mgf from the highest are not peaks for M-xx ions
                            # (within ppm values)
                            # or their isotopic peaks, in case they are not set 100 for the highest and normalize the I
                            # of the matching ions based on that mz value
                            if (mz_mgf < self.MS2_mzlow or mz_mgf > self.MS2_mzhigh or
                                    (mz_digest - ppm_range(mz_digest, self.MS2_ppm) <= mz_mgf <= mz_digest +
                                     ppm_range(mz_digest, self.MS2_ppm)) or
                                    (mz_digest - ppm_range(mz_digest, self.MS2_ppm) + neutron_mass / charge <=
                                     mz_mgf <= mz_digest + ppm_range(mz_digest,
                                                                     self.MS2_ppm) + neutron_mass / charge) or
                                    (mz_digest - ppm_range(mz_digest, self.MS2_ppm) - neutron_mass / charge <=
                                     mz_mgf <= mz_digest + ppm_range(mz_digest,
                                                                     self.MS2_ppm) - neutron_mass / charge) or
                                    (mz_digest - ppm_range(mz_digest, self.MS2_ppm) + 2 * neutron_mass / charge <=
                                     mz_mgf <= mz_digest + ppm_range(mz_digest,
                                                                     self.MS2_ppm) + 2 * neutron_mass / charge)):
                                flag = 1

                        # The highest intensity mgf peak non matching to M-xx ions is set to be 100 in the normalization
                        if flag == 0:
                            max_int = np.float64(mgf_peaks[key][f])
                            break

                    normalized_list = []

                    # Normalizing the intensity of the matching ions
                    for ele in match_peaks[key][t][7:]:
                        norm_int = round(np.float64(ele[-1]) / max_int * 100, 2)
                        new_list = ele[:3] + [str(norm_int)]
                        if norm_int >= self.MS2_normint_cutoff:
                            normalized_list.append(new_list)

                    output_dic[key].update({t: match_peaks[key][t][:7] + normalized_list})

        return output_dic

    def MS2_max_peaks(self):
        """
        Return the number of most intense MS2 ions to limit the search to
        """
        if self.MS2_peak_num_max == 'all':
            max_value = 999999
        else:
            max_value = self.MS2_peak_num_max

        return int(max_value)

    def exclusion_windows_matching(self, match_peaks):
        """
        Discard the occurrences of matching and non-matching ions when they are found in the window
        (+-losses_window_removal) around M-xx or free bases ions
        """
        output_dic = match_peaks

        for key in match_peaks:
            if match_peaks[key]:

                for t in match_peaks[key]:
                    mass_losses_list, new_list = find_losses_freebases(match_peaks[key][t][7:]), []
                    for ion in match_peaks[key][t][7:]:

                        # Keep ion losses and free bases matched in the MS2_matches list
                        if 'M-' not in ion[1] and len(ion[1].split('(')[0]) != 1:

                            flag, mz_ion = 1, np.float64(ion[2])

                            for mass_loss in mass_losses_list:

                                # Add the MS2 offset
                                mass_loss_offseted = mass_loss + ppm_range(mass_loss, self.MS2_ppm_offset)

                                # Check and discard any sequencing ion is found in the M-xx exclusion window
                                if mass_loss_offseted - self.losses_window_removal <= \
                                        mz_ion <= mass_loss_offseted + self.losses_window_removal:
                                    flag = 0
                                    break

                            if flag == 1:
                                new_list.append(ion)

                        else:
                            new_list.append(ion)

                    output_dic[key].update({t: match_peaks[key][t][:7] + new_list})

        return output_dic

    def consecutive_series_5(self, list_MS2, ion_series):
        """
        Return beta for each of the 5 ion series groups.
        beta is determined as:
                if no consecutive matches: beta = 0
                if at least one consecutive match: beta = sum(beta_increment + (alpha * consecutive_matches ) *
                beta_increment) over all consecutive matches in the same ion series group (e.g. a/b, w/x...)

        Note that the series a-B, y-P and z-P require a special correction to take into account their different
        numbering scheme
        """
        beta, consec_flag = 0, 0

        for i in range(1, 20):
            if len(ion_series) == 4:
                if (((ion_series[0] + str(i)) in list_MS2 or (ion_series[1] + str(i)) in list_MS2 or (
                        ion_series[2][0] + str(i) + ion_series[2][1:]) in list_MS2 or (
                             ion_series[3][0] + str(i) + ion_series[3][1:]) in list_MS2) and
                        ((ion_series[0] + str(i + 1)) in list_MS2 or (ion_series[1] + str(i + 1)) in list_MS2 or (
                                ion_series[2][0] + str(i + 1) + ion_series[2][1:]) in list_MS2 or (
                                 ion_series[3][0] + str(i + 1) + ion_series[3][1:]) in list_MS2)):

                    beta += self.beta_increment + (self.alpha * consec_flag) * self.beta_increment
                    consec_flag += 1

                else:
                    consec_flag = 0

            else:
                if ((ion_series[0] + str(i)) in list_MS2 or (ion_series[1] + str(i)) in list_MS2) and (
                        (ion_series[0] + str(i + 1)) in list_MS2 or (ion_series[1] + str(i + 1)) in list_MS2):

                    beta += self.beta_increment + (self.alpha * consec_flag) * self.beta_increment
                    consec_flag += 1

                else:
                    consec_flag = 0

        return beta

    def consecutive_series_11(self, list_MS2, ion_series):
        """
        Returns beta for each of the 11 ion series groups.
        beta is determined as:
                if no consecutive matches: beta = 0
                if at least one consecutive match: beta = sum(beta_increment + (alpha * consecutive_matches ) *
                beta_increment) over all consecutive matches in the same ion series

        Note that the series a-B, y-P and z-P require a special correction to take into account their different
        numbering scheme
        """
        beta, consec_flag = 0, 0

        for i in range(1, 20):

            # Add the support for a-b, y-P and z-P that have the number after the first character and not at the end
            if len(ion_series) == 1:
                current, consecutive = ion_series + str(i), ion_series + str(i + 1)

            else:
                current, consecutive = ion_series[0] + str(i) + ion_series[1:], ion_series[0] + str(i + 1) + ion_series[
                                                                                                             1:]

            if current in list_MS2 and consecutive in list_MS2:
                beta += self.beta_increment + (self.alpha * consec_flag) * self.beta_increment
                consec_flag += 1

            else:
                consec_flag = 0

        return beta

    def L(self, mz_match, d, mod_bases):
        """
        Return L from the scoring function. Only sequence-defining fragment ions are considered (11 series)
        Disregard multiple charged states for the same series (e.g. a2(-1) and a2(-2) will be counted as 1)
        """
        count, unique_series = 0, []
        for ion in d[mz_match][-1]:

            # Consider for the purpose of L only the MS2 ions with m/z within the option-specified window
            if self.MS2_mzlow <= np.float64(ion.split(':')[1]) <= self.MS2_mzhigh:
                series = ion.split('(')[0]

                # Exclude all ions from neutral/charged losses and free bases
                if 'M-' not in series and series != 'G' and series != 'A' and series != 'C' and series != 'U':

                    # Exclude free bases from modified nucleotides
                    if mod_bases:
                        if series not in mod_bases and series not in unique_series:
                            count += 1
                            unique_series.append(series)

                    else:
                        if series not in unique_series:
                            count += 1
                            unique_series.append(series)

        return count

    def sumI(self, list_MS2, mod_bases):
        """
        Return the sum of intensities for a given precursor ion (first term of the scoring function)
        """
        sumI, mgf_peaks = 0, []

        # Adds recursively the intensities values to the sum
        for l in list_MS2:
            if self.MS2_mzlow < np.float64(l[2]) < self.MS2_mzhigh:
                series, mgf_match = l[1].split('(')[0], str("{0:.6f}".format(np.float64(l[2])))
                # When two series correspond to the same MS2 ion peak from the mgf,
                # the intensity is summed only once since they both participate to the same intensity
                # value being indistinguishable
                if mgf_match not in mgf_peaks:

                    mgf_peaks.append(mgf_match)

                    # Exclude all ions from neutral/charged losses and free bases
                    if 'M-' not in series and series != 'G' and series != 'A' and series != 'C' and series != 'U':

                        # Excludes free bases from modified nucleotides
                        if mod_bases:
                            if series not in mod_bases:
                                sumI += np.float64(l[-1])

                        else:
                            sumI += np.float64(l[-1])

        return sumI

    def n_calc(self, list_MS2, mod_bases):
        """
        Return the value n (number of MS2 ions) for the scoring function disregarding multiple charged species for the
        same ion series
        (e.g. a2(-1) and a2(-2) will count as 1)
        """
        n, unique_series = 0, []

        # MS2 ions following specific criteria are recursively added to the n output
        for l in list_MS2:

            if self.MS2_mzlow < np.float64(l[0]) < self.MS2_mzhigh:

                series = l[1].split('(')[0]

                # Exclude all ions from neutral/charged losses and free bases
                if 'M-' not in series and series != 'G' and series != 'A' and series != 'C' and series != 'U':

                    # Exclude free bases from modified nucleotides
                    if mod_bases:
                        if series not in mod_bases and series not in unique_series:
                            n += 1
                            unique_series.append(series)

                    else:
                        if series not in unique_series:
                            n += 1
                            unique_series.append(series)

        return n

    def digest_dic(self, digest_file):
        """
        Make a dictionary with the info from the digest file given as input
        Format -> m/z_precursor : all other info, ..., [CID fragments]
        """
        d = {}
        global enzyme
        for line in digest_file:

            if "ENZYME" in line:
                enzyme = line.split()[1]

            if line[0].isdigit():

                # Check if the lines to consider are all or only light/heavy
                if self.light_heavy == 'all':

                    x = line.split()[1:13]
                    x.append(line.split()[13:])
                    d[line.split()[0] + "_" + line.split()[7] + "_" + re.findall(r'\d+', line.split()[5])[0]] = x

                else:

                    # Only light/heavy lines (as specified by option) are considered
                    if self.light_heavy in line:
                        x = line.split()[1:13]
                        x.append(line.split()[13:])
                        d[line.split()[0] + "_" + line.split()[7] + "_" + re.findall(r'\d+', line.split()[5])[0]] = x

        return d

    def mz_freebases(self):
        """
        Extract the theoretical m/z values for all free bases, in order to include them in the matching even if they
        are outside the MS2_mz window. These ions will not be considered for scoring but will nonetheless appear in
        the output together with all the other matching ions
        """
        mz_list, digest = [], self.digest_dic(open(self.digest_file, 'r'))

        # Add the m/z values from free bases extracted from the digest
        for b in digest:
            for ion in digest[b][-1]:
                if len(ion.split("(")[0]) == 1:
                    if ion.split(':')[-1] not in mz_list:
                        mz_list.append(ion.split(':')[-1])

        return mz_list

    def mgf_dic(self):
        """
        Make a dictionary with MS1/MS2 ions from input .mgf file
        Format -> {m/z_precursor_RT : { m/z ion1 : norm_intensity1 , m/z ion2 : norm_intensity2 ...}, ... }
        """

        mgf_file, mz_freeb = open(self.mgf_file, 'r'), self.mz_freebases()
        d = {}
        # Create a dictionary with all precursor m/z_RT as keys, to be filled with info for the final output
        # (M intensity and absolute sumI value)
        global mgf_peaks, info_MS2_scans
        mgf_peaks, info_MS2_scans = 0, {}

        # Read the lines of mgf input file
        for line in mgf_file:

            # Set values for prec_mass, rt and charge in order to avoid errors when those values are missing from
            # the mgf. This means precursors with missing mass values, charge or rt will not be considered for matching

            # NOTE: charge is assigned as default to 1 in case no charge information is present in the mgf
            if "BEGIN IONS" in line:
                prec_mass, rt, charge = 0, 0, str(1)

            # m/z of the precursor ions are used as keys for the dictionary
            if "PEPMASS" in line:
                prec_mass = '{:.6f}'.format(np.float64(line.split()[0].split('=')[1]))
                mgf_peaks += 1

            # Charge of the precursor ion
            if "CHARGE" in line:
                charge = re.findall(r'\d+', line.split('=')[1])[0]

            # RT for the precursor ion
            if "RTINSECONDS" in line:
                rt, list_mz_MS2 = str(np.float64(line.split('=')[1][:-1])), []

            if line[0].isdigit():
                # Add m/z value to the list of MS2 m/z
                list_mz_MS2.append(line.split()[0] + "_" + line.split()[1])

            # Determine the intensity for the precursor ion peak (closest value within MS2 scans to the PEPMASS)
            if "END IONS" in line:
                # Control on charges and m/z boundaries for MS1 precursor ions
                if self.MS1_mzlow <= np.float64(prec_mass) <= self.MS1_mzhigh:
                    d[prec_mass + "_" + rt + "_" + charge], info_MS2_scans[prec_mass + "_" + rt] = {}, []

                    if list_mz_MS2:

                        info_MS2_scans[prec_mass + "_" + rt].append(closest_value(list_mz_MS2, prec_mass))

                        for line in list_mz_MS2:
                            # Add all the ion m/z : intensity info (excludes values around precursor ion,
                            # also if below MS2_peak_int_min)
                            l = line.split('_')

                            if (prec_mass + "_" + rt + "_" + charge in d.keys() and np.float64(
                                    l[1]) > self.threshold_MS2_int() and (
                                    np.float64(l[0]) < np.float64(prec_mass) - self.precursor_window_removal or
                                    np.float64(l[0]) > np.float64(prec_mass) + self.precursor_window_removal)):
                                d[prec_mass + "_" + rt + "_" + charge].update({l[0]: l[1]})

                            # If no MS2 ions are listed in the mgf, a single fake 1 m/z : 1 int peak is created to avoid
                            # errors in the matching process (will result in score 0)
                    else:
                        info_MS2_scans[prec_mass + "_" + rt].append('1'), d[prec_mass + "_" + rt + "_" + charge].update(
                            {'1': '1'})

        # Normalize and sort the list of MS2 m/z:intensity
        for key in d:
            lista = []

            # Create a list with all MS2 ions entries, if m/z of MS2 ions are within the specified MS2 m/z window
            for k in d[key]:
                if self.MS2_mzlow <= np.float64(k) <= self.MS2_mzhigh:
                    lista.append([k, d[key][k]])

                # Adds free bases ions even if they are outside the specified MS2 window
                else:
                    for ion in mz_freeb:
                        if (np.float64(ion) + ppm_range(np.float64(ion), self.MS2_ppm_offset) -
                                ppm_range(np.float64(ion) + ppm_range(np.float64(ion), self.MS2_ppm_offset),
                                          self.MS2_ppm) <= np.float64(k) <= np.float64(ion) + ppm_range(np.float64(ion),
                                                                                                        self.MS2_ppm_offset) + ppm_range(
                                    np.float64(ion) + ppm_range(np.float64(ion),
                                                                self.MS2_ppm_offset), self.MS2_ppm)):
                            lista.append([k, d[key][k]])
                            break

            # Check if in the mgf file there is at least a MS2 ion
            if lista:
                lista.sort(key=lambda x: np.float64(x[1]), reverse=True)

                # Add the absolute value of the highest intensity to the info dictionary (in preparation for the output)
                info_MS2_scans[key.split('_')[0] + '_' + key.split('_')[1]].append(lista[0][1])

                # Eliminate eventual eccess of MS2 ions found based on a given upper value
                if len(lista) > self.MS2_max_peaks():
                    lista = lista[:self.MS2_max_peaks()]

                # Writes the normalized and sorted MS2 scans into the output dictionary
                d[key] = {x[0]: str(x[1]) for x in lista}

        return d

    def matching(self):
        """
        Match the MS1/MS2 ions between theoretical digest and experimental .mgf file
        The 'all' addition are meant to create an additional dictionary with a list of all ions in the mgf irrespective
        of matching
        """
        digest_dic, mgf_dic = self.digest_dic(open(self.digest_file, 'r')), self.mgf_dic()
        MS2_output, MS2_output_all = {}, {}

        # Cycle among all entries in the mgf dictionary
        for key in mgf_dic:
            MS2_match, prec_mass, prec_charge, MS2_match_all = {}, np.float64(key.split('_')[0]), int(
                key.split('_')[2]), {}

            # MS1 precursor ions matching. Check and apply the search for isotopologues if selected by the user
            if not self.precursor_isotopologues:

                # The charges in the mgf file are ignored when matching if specified by the user
                if not self.ignore_charges_mgf:
                    prec_match = (dict((k, v) for k, v in digest_dic.items() if
                                       prec_mass + ppm_range(prec_mass, self.MS1_ppm_offset) - ppm_range(prec_mass +
                                                                                                         ppm_range(
                                                                                                             prec_mass,
                                                                                                                       self.MS1_ppm_offset),
                                                                                                         self.MS1_ppm) <= np.float64(
                                           k.split('_')[0]) <= prec_mass +
                                       ppm_range(prec_mass, self.MS1_ppm_offset) + ppm_range(
                                           prec_mass + ppm_range(prec_mass, self.MS1_ppm_offset),
                                           self.MS1_ppm) and prec_charge == int(k.split('_')[2])))
                else:
                    prec_match = (dict((k, v) for k, v in digest_dic.items() if
                                       prec_mass + ppm_range(prec_mass, self.MS1_ppm_offset) - ppm_range(prec_mass +
                                                                                                         ppm_range(
                                                                                                             prec_mass,
                                                                                                             self.MS1_ppm_offset),
                                                                                                         self.MS1_ppm) <=
                                       np.float64(k.split('_')[0]) <= prec_mass +
                                       ppm_range(prec_mass, self.MS1_ppm_offset) +
                                       ppm_range(prec_mass + ppm_range(prec_mass, self.MS1_ppm_offset), self.MS1_ppm)))

            else:

                # The charges in the mgf file are used for matching if specified by the user
                if not self.ignore_charges_mgf:
                    prec_match = (dict((k, v) for k, v in digest_dic.items() if prec_charge == int(k.split('_')[2])
                                       and ((prec_mass + ppm_range(prec_mass, self.MS1_ppm_offset) - ppm_range(
                        prec_mass + ppm_range(prec_mass, self.MS1_ppm_offset), self.MS1_ppm) <=
                                             np.float64(k.split('_')[0]) <= prec_mass +
                                             ppm_range(prec_mass, self.MS1_ppm_offset) + ppm_range(prec_mass +
                                                                                                   ppm_range(prec_mass,
                                                                                                             self.MS1_ppm_offset),
                                                                                                   self.MS1_ppm)) or (
                                                    prec_mass +
                                                    neutron_mass / abs(prec_charge) + ppm_range(prec_mass,
                                                                                                self.MS1_ppm_offset) -
                                                    ppm_range(prec_mass + ppm_range(prec_mass, self.MS1_ppm_offset),
                                                              self.MS1_ppm) <= np.float64(
                                                k.split('_')[0]) <= prec_mass + neutron_mass / abs(
                                                prec_charge) + ppm_range(prec_mass,
                                                                         self.MS1_ppm_offset) + ppm_range(
                                                prec_mass + ppm_range(prec_mass, self.MS1_ppm_offset),
                                                self.MS1_ppm)) or
                                            (prec_mass - neutron_mass / abs(prec_charge) +
                                             ppm_range(prec_mass, self.MS1_ppm_offset) - ppm_range(
                                                        prec_mass + ppm_range(
                                                            prec_mass, self.MS1_ppm_offset),
                                                        self.MS1_ppm) <= np.float64(k.split('_')[0])
                                             <= prec_mass - neutron_mass / abs(prec_charge) +
                                             ppm_range(prec_mass, self.MS1_ppm_offset) + ppm_range(prec_mass +
                                                                                                   ppm_range(prec_mass,
                                                                                                             self.MS1_ppm_offset),
                                                                                                   self.MS1_ppm)))))
                else:
                    prec_match = (dict((k, v) for k, v in digest_dic.items() if
                                       ((prec_mass + ppm_range(prec_mass, self.MS1_ppm_offset) -
                                         ppm_range(prec_mass + ppm_range(prec_mass, self.MS1_ppm_offset),
                                                   self.MS1_ppm) <=
                                         np.float64(k.split('_')[0]) <= prec_mass +
                                         ppm_range(prec_mass, self.MS1_ppm_offset) + ppm_range(prec_mass +
                                                                                               ppm_range(prec_mass,
                                                                                                         self.MS1_ppm_offset),
                                                                                               self.MS1_ppm))
                                        or (prec_mass + neutron_mass / abs(prec_charge) +
                                            ppm_range(prec_mass, self.MS1_ppm_offset) - ppm_range(prec_mass +
                                                                                                  ppm_range(prec_mass,
                                                                                                            self.MS1_ppm_offset),
                                                                                                  self.MS1_ppm) <=
                                            np.float64(k.split('_')[0]) <= prec_mass + neutron_mass / abs(prec_charge) +
                                            ppm_range(prec_mass, self.MS1_ppm_offset) + ppm_range(prec_mass +
                                                                                                  ppm_range(prec_mass,
                                                                                                            self.MS1_ppm_offset),
                                                                                                  self.MS1_ppm))
                                        or (prec_mass - neutron_mass / abs(prec_charge) +
                                            ppm_range(prec_mass, self.MS1_ppm_offset) -
                                            ppm_range(prec_mass + ppm_range(prec_mass, self.MS1_ppm_offset),
                                                      self.MS1_ppm)
                                            <= np.float64(k.split('_')[0]) <= prec_mass - neutron_mass / abs(
                                                           prec_charge) +
                                            ppm_range(prec_mass, self.MS1_ppm_offset) + ppm_range(prec_mass +
                                                                                                  ppm_range(prec_mass,
                                                                                                            self.MS1_ppm_offset),
                                                                                                  self.MS1_ppm)))))

            # Cycle among all precursor ions matching between digest and mgf file
            for key_match in prec_match:
                MS2_header, MS2_list = [prec_match[key_match][0], prec_match[key_match][4], prec_match[key_match][6],
                                        prec_match[key_match][7], prec_match[key_match][8], prec_match[key_match][9],
                                        prec_match[key_match][11]], []
                MS2_header_all = [prec_match[key_match][0], prec_match[key_match][4], prec_match[key_match][6],
                                  prec_match[key_match][7], prec_match[key_match][8], prec_match[key_match][9],
                                  prec_match[key_match][11]]

                # Find matches of MS2 fragments between mgf file and theoretical digest
                for x in prec_match[key_match][-1]:
                    ion_mass, ion_charge = np.float64(x.split(':')[-1]), int(x.split(')')[0][-1])
                    # Matched only ions within specified MS2 charge, and within given ppm offset
                    add = (list([k, v] for k, v in mgf_dic[key].items() if ion_mass +
                                ppm_range(ion_mass, self.MS2_ppm_offset) -
                                ppm_range(ion_mass + ppm_range(ion_mass, self.MS2_ppm_offset), self.MS2_ppm) <=
                                np.float64(k) <= ion_mass + ppm_range(ion_mass, self.MS2_ppm_offset) + ppm_range(
                        ion_mass + ppm_range(ion_mass, self.MS2_ppm_offset), self.MS2_ppm)))

                    # Select only MS2 ions that have a relative intensity above a certain threshold
                    if add:
                        MS2_max_intensity_ion = pick_max_intensity(add)

                        MS2_list.append(
                            [x.split(':')[-1], x.split(':')[0], MS2_max_intensity_ion[0], MS2_max_intensity_ion[1]])

                MS2_list.sort(key=lambda x: np.float64(x[3]))

                matched_ions = {}
                for ion in MS2_list:
                    matched_ions.update({ion[2]: ion[1]})

                all_mgf_list, flag = [], 0.01

                for a in map(list, mgf_dic[key].items()):
                    if a[0] not in matched_ions.keys():
                        all_mgf_list.append([self.MS2_mzlow + flag, ''] + a)

                    else:
                        all_mgf_list.append([self.MS2_mzlow + flag, matched_ions[a[0]]] + a)

                    flag += 0.01

                # Add the ions to the dictionary of MS2 ions
                if MS2_list:
                    MS2_header.extend(MS2_list[::-1])

                # Add the ions to the dictionary of MS2 ions
                if all_mgf_list:
                    MS2_header_all.extend(all_mgf_list[::-1])

                MS2_match.update({key_match: MS2_header})
                MS2_match_all.update({key_match: MS2_header_all})

            MS2_output.update({key: MS2_match})
            MS2_output_all.update({key: MS2_match_all})

        # Obtain a reduced version of the digest with only M-xx ions to be excluded in the normalization
        red_digest = reduced_digest(digest_dic)

        # Normalize the intensities of matching MS2 ions
        self.normalize_int_MS2(mgf_dic, MS2_output, red_digest)
        self.normalize_int_MS2(mgf_dic, MS2_output_all, red_digest)

        # Discard the occurrences of sequencing ions found in the vicinity (+-MS2_ppm) of M-xx +
        # their isotopic peak (-1, +1 and +2 peaks) and free bases
        self.exclusion_windows_matching(MS2_output)
        self.exclusion_windows_matching(MS2_output_all)

        return MS2_output, MS2_output_all

    def scoring(self, d, d_all):
        """
        Add the Sp score to the precursor ion matches
        """
        d2, dig = {}, self.digest_dic(open(self.digest_file, 'r'))
        for key_rt in d:

            d2[key_rt] = {}
            for key in d[key_rt]:

                # Sp Score value calculation with the two score_calc functions:
                mod_nts = mod_detection(d[key_rt][key][2])

                sumi, MS2_scans, gamma = self.sumI(d[key_rt][key][7:], mod_nts), d[key_rt][key][7:], 0
                sumi_all = self.sumI(d_all[key_rt][key][7:], mod_nts)

                n, l = self.n_calc(d[key_rt][key][7:], mod_nts), self.L(key, dig, mod_nts)

                betaa, betaaB, betab, betac, betad, betaw, betax, betay, betaz, betayP, betazP = (
                    self.consecutive_series_11(" ".join(str(r) for v in MS2_scans for r in v), 'a'),
                    self.consecutive_series_11(" ".join(str(r) for v in MS2_scans for r in v), 'a-B'),
                    self.consecutive_series_11(" ".join(str(r) for v in MS2_scans for r in v), 'b'),
                    self.consecutive_series_11(" ".join(str(r) for v in MS2_scans for r in v), 'c'),
                    self.consecutive_series_11(" ".join(str(r) for v in MS2_scans for r in v), 'd'),
                    self.consecutive_series_11(" ".join(str(r) for v in MS2_scans for r in v), 'w'),
                    self.consecutive_series_11(" ".join(str(r) for v in MS2_scans for r in v), 'x'),
                    self.consecutive_series_11(" ".join(str(r) for v in MS2_scans for r in v), 'y'),
                    self.consecutive_series_11(" ".join(str(r) for v in MS2_scans for r in v), 'z'),
                    self.consecutive_series_11(" ".join(str(r) for v in MS2_scans for r in v), 'y-P'),
                    self.consecutive_series_11(" ".join(str(r) for v in MS2_scans for r in v), 'z-P'))

                # Score calculation (+ addition of all single factors for debugging)
                score = (str(
                    score_calc_11(sumi, n, betaa, betaaB, betab, betac, betad, betaw, betax, betay, betaz, betayP,
                                  betazP, gamma, l, sumi_all)) + "(sumI=" + str(int(round(sumi, 0))) + ";n=" + str(
                    n) + ";ba=" + str(round(np.float64(betaa), 3)) + ";ba-B=" + str(round(np.float64(betaaB), 3)) +
                         ";bb=" + str(round(np.float64(betab), 3)) + ";bc=" + str(
                            round(np.float64(betac), 3)) + ";bd=" + str(round(np.float64(betad), 3)) + ";bw=" + str(
                            round(np.float64(betaw), 3)) + ";bx=" + str(round(np.float64(betax), 3)) + ";by=" + str(
                            round(np.float64(betay), 3))
                         + ";bz=" + str(round(np.float64(betaz), 3)) + ";by-P=" + str(
                            round(np.float64(betayP), 3)) + ";bz-P=" + str(
                            round(np.float64(betazP), 3)) + ";gamma=" + str(gamma) + ";L=" + str(l) + ";SumIall=" + str(
                            int(round(sumi_all, 0))) + ")")

                d2[key_rt][key + "_" + str(score)] = d[key_rt][key]

        return d2

    def consolidate_match(self):
        """
        All the occurrences of the same precursor with different RT are clustered together
        """
        dict, dict_all = self.matching()
        dic = self.scoring(dict, dict_all)
        lista, d = [], {}

        for key in dic:

            # Determine the "unique" precursor ions m/z values
            if key.split('_')[0] not in lista:
                lista.append(key.split('_')[0])

        # Cluster entries of precursor ions with different RT but same m/z together
        for key in lista:
            d[key] = list([k, v] for k, v in dic.items() if np.float64(k.split('_')[0]) == np.float64(key))

        return d

    def output(self):
        """
        Write the lines for the output file
        """
        match_dic = self.consolidate_match()
        global MS1_hits, MS2_hits
        MS1_hits, MS2_hits, list_ion_rt = 0, 0, []

        if self.MS1_ppm_offset == 0:
            MS1_ppm_off = 0.0
        else:
            MS1_ppm_off = - self.MS1_ppm_offset

        if self.MS2_ppm_offset == 0:
            MS2_ppm_off = 0.0
        else:
            MS2_ppm_off = - self.MS2_ppm_offset

        output_lines = ["#theoretical_digest " + self.digest_file + "\n#enzyme " + str(
            enzyme) + "\n#MS_data " + self.mgf_file + "\n#isotopic_species " + self.light_heavy +
                        "\n#MS1_ppm " + str(self.MS1_ppm) + "\n#MS2_ppm " + str(
            self.MS2_ppm) + "\n#MS1_offset_ppm " + str(MS1_ppm_off) + "\n#MS2_offset_ppm " + str(
            MS2_ppm_off) + "\n#MS1_mz_minimum " + str(self.MS1_mzlow) + "\n#MS1_mz_maximum " + str(
            self.MS1_mzhigh) + "\n#MS2_mz_minimum " + str(self.MS2_mzlow) +
                        "\n#MS2_mz_maximum " + str(self.MS2_mzhigh) + "\n#MS2_abs_peak_intensity " + str(
            self.MS2_peak_int_min) + "\n#MS2_peak_num_maximum " + str(
            self.MS2_peak_num_max) + "\n#MS2_normint_cutoff " + str(self.MS2_normint_cutoff) +
                        "\n#precursor_exclusion_window " + str(self.precursor_window_removal) +
                        "\n#precursor_losses_exclusion_window " + str(self.losses_window_removal) +
                        "\n#alpha " + str(self.alpha) + "\n#beta " + str(self.beta_increment) +
                        "\n#precursor_isotopologues " + str(self.precursor_isotopologues) + "\n"]
        output_lines.append(
            "#MATCHES_HEADER:m/z(meas) RT m/z(theo) offset Sp dSp rank"
            " #MS2_matches isotope length charge sequence sequence_mod 5'-end 3'-end molecule_location"
            " MS2_matches-->m/z_measured(offset(ppm))[norm_intensity]:m/z_theoretical[CID_ion] detailed_score\n\n")

        # Write into the output lines the info on each precursor ion and MS2 matches
        sort = sorted(list(match_dic), key=lambda x: np.float64(x))
        for key in sort:
            list_precursor = []

            # Make sure precursor that are within the score thresholds set by parameters are written in the output
            ctrl = False
            for l in match_dic[key]:
                if l[1]:
                    ctrl = True
                    break

            if ctrl:
                output_lines.append("\nPRECURSOR_ION=" + key + "\n")  # ION line for each measured precursor ion peak

                for ion_rt in match_dic[key]:
                    MS1_hits += 1

                    for match in ion_rt[1]:

                        # Calculate the number of MS2 ions excluding neutral/charged losses and free bases
                        n_ms2_ions = self.n_calc(ion_rt[1][match][7:], mod_detection(ion_rt[1][match][2]))

                        # Correct the offset for matches of isotopologues and mark them with a *
                        calculated_offset = ppm_offset(np.float64(ion_rt[0].split('_')[0]) + ppm_range(
                                                           np.float64(ion_rt[0].split('_')[0]), self.MS1_ppm_offset),
                                                       match.split('_')[0])
                        if abs(calculated_offset) > self.MS1_ppm:
                            corrected_offset_minus = ppm_offset(np.float64(ion_rt[0].split('_')[0])
                                                                 - neutron_mass / abs(int(match.split('_')[2])) +
                                                                 ppm_range(np.float64(ion_rt[0].split('_')[0]),
                                                                           self.MS1_ppm_offset),
                                                                np.float64(match.split('_')[0]))
                            corrected_offset_plus = ppm_offset(
                                                                 np.float64(ion_rt[0].split('_')[0])
                                                                 + neutron_mass / abs(int(match.split('_')[2])) +
                                                                 ppm_range(np.float64(ion_rt[0].split('_')[0]),
                                                                           self.MS1_ppm_offset),
                                                                 np.float64(match.split('_')[0]), )

                            if abs(corrected_offset_plus) > abs(corrected_offset_minus):
                                corrected_offset = corrected_offset_minus
                            else:
                                corrected_offset = corrected_offset_plus

                            corrected_offset = str(corrected_offset) + '*'

                        else:
                            corrected_offset = str(calculated_offset)

                        line_list = [ion_rt[0].split('_')[0], str(np.float64(ion_rt[0].split('_')[1])),
                                     match.split('_')[0],
                                     match.split('_')[-1],
                                     corrected_offset, str(n_ms2_ions), " ".join(ion_rt[1][match][:7])]

                        if n_ms2_ions != 0 and ion_rt[0].split('_')[0] + " RT=" + str(
                                np.float64(ion_rt[0].split('_')[1])) not in list_ion_rt:
                            MS2_hits += 1
                            list_ion_rt.append(
                                ion_rt[0].split('_')[0] + " RT=" + str(np.float64(ion_rt[0].split('_')[1])))

                        # Add all info on MS2 matches
                        for i, MS2 in enumerate(ion_rt[1][match][6:-1]):
                            measured_mass = np.float64(ion_rt[1][match][i + 7][2])
                            theoretical_mass = np.float64(ion_rt[1][match][i + 7][0]) + \
                                               ppm_range(np.float64(ion_rt[1][match][i + 7][0]), self.MS2_ppm_offset)
                            line_list.append(str("{0:.6f}".format(np.float64(ion_rt[1][match][i + 7][2]))) + "(" + str(
                                -ppm_offset(theoretical_mass, measured_mass)) +
                                             "ppm)[" + str(
                                int(round(np.float64(ion_rt[1][match][i + 7][3]), 0))) + "]:" + str(
                                "{0:.6f}".format(np.float64(ion_rt[1][match][i + 7][0]))) + "[" +
                                             ion_rt[1][match][i + 7][
                                                 1] + "]")

                        list_precursor.append(line_list)

            list_prec_sorted = sorted(list_precursor, key=lambda x: np.float64(x[1]))

            # Group the precursor ions based on the same retention time (keeping the original order
            # based on descending score and adding an ascending order by abs(MS1_ppm offset))
            groups = []
            for key, group in itertools.groupby(list_prec_sorted, lambda x: x[1]):
                groups.append(sorted(list(group), key=lambda x: (np.float64(x[3].split('_')[0].split('(')[0]),
                                                                 -abs(int(re.findall(r'\d+', x[4])[0]))), reverse=True))

            list_prec_sorted = sorted(groups, key=lambda x: np.float64(x[0][3].split('_')[0].split('(')[0]),
                                      reverse=True)

            for m in list_prec_sorted:
                h_score, flag, counter = 0, 1, 1
                for l in m:
                    score = np.float64(l[3].split('(')[0])
                    if score > h_score:
                        h_score = score

                    if score == 0 and h_score == 0:
                        delta_sp = 0

                    else:
                        delta_sp = (h_score - score) / h_score

                    # Add a star to the mass of isotopologues matches and correct the MS1_ppm error to the value
                    # centered on the actual isotopologue
                    if '*' in l[4]:
                        matched_mass, MS1_ppm = l[0] + '*', l[4][:-1]
                    else:
                        matched_mass, MS1_ppm = l[0], l[4]

                    isotope, charge, sequence, sequence_mod, chem5, chem3, sequence_loc = l[6].split()
                    line = ("{} RT={} TH_MATCH={} {}ppm Sp={} dSp={} rank={} #MS2={} {} {} {} {} {} {} {} {} {} "
                            "SCORE={}\n".format(matched_mass, round(np.float64(l[1]) / 60, 3), l[2], MS1_ppm, score,
                                                round(delta_sp, 2), counter, l[5], isotope, len(sequence), charge,
                                                sequence,
                                                sequence_mod, chem5, chem3, sequence_loc, " ".join(l[7:]), l[3]))

                    # Control on the number of top scoring matches to output
                    if h_score == 0 or score / h_score >= 0:
                        output_lines.append(line)

                    flag += 1
                    counter += 1

        return output_lines

    def final_output(self):
        """
        Write the output files
        """
        time = datetime.now()

        out = self.output()

        # Determine the name for all the output files, from the input mgf
        mgf_name = filename_from_path(self.mgf_file)[:-4]

        # Write the lines in the output file "match_output.txt"
        open("./match_output_" + mgf_name + ".txt", 'w').writelines(out)

        # Write a log file with statistics on matching hits
        open("./log.txt", 'w').writelines(log_file())

        if not self.targets_without_decoys:
            targets_with_decoys = 'y'
        else:
            targets_with_decoys = 'n'

        # Output specific for the statistical analysis
        stats.csv_output(
            stats.filter_data(stats.input_data("match_output_" + mgf_name + ".txt", 'n', 'y'), 'y', 'Sp'),
            'all', 'n', 0, mgf_name, targets_with_decoys,
            self.FDR_light_heavy)

        print("\nDone! Output file(s) -> {} {} {} {} {}".format("match_output_" + mgf_name + ".txt",
                                                                "log.txt",
                                                                "targets_{}.csv".format(mgf_name),
                                                                "decoys_{}.csv".format(mgf_name),
                                                                "FDR_{}.csv".format(mgf_name)))

        print("start: {} end: {}".format(time, datetime.now()))


def score_calc_5(sumI, n, beta_a_b, beta_aB, beta_c_d, beta_w_x, beta_y_z_P, gamma, L, sumI_all):
    """
    Scoring function with 5 ion series: sumI(matched) / sumI(all) * n / L * ( 1 + beta[a + b] + beta[a-B] + beta[c + d]
     + beta[w + x] + beta[y + z + y-P + z-P])

    n = total number of matching MS2 ions for the given precursor, excluding non sequence-defining ions
    (only 11 main series) within the specified m/z window, charge independent
    beta is determined as:
            if no consecutive matches in the same ion series: beta = 0
            if at least one consecutive match in the same ion series: beta = sum(beta_increment +
            (x_beta * consecutive_matches ) * beta_increment) over all consecutive matches in the same
            ion series (a,a-b,b..) or
            ion series group (e.g. a/b, w/x...). beta_increment and x_beta are tunable via input options

     L = normalization factor base on total number of theoretical predicted ions for a precursor ion excluding
     non sequence-defining ions (only 11 main series), wihin the specified m/z window
    """
    return round(sumI * n * (1 + beta_a_b + beta_aB + beta_c_d + beta_w_x + beta_y_z_P) / (L * sumI_all),
                 3)


def score_calc_11(sumI, n, beta_a, beta_aB, beta_b, beta_c, beta_d, beta_w, beta_x, beta_y, beta_z, beta_yP, beta_zP,
                  gamma, L, sumI_all):
    """
    Scoring function with 11 ion series: sumI(matched) / sumI(all) * n / L * ( 1 + beta[a] + beta[a-B] + beta[b] +
    beta[c] + beta[d] + beta[w] + beta[x] + beta[y] + beta[z] + beta[y-P] + beta[z-P])

     n = total number of matching MS2 ions for the given precursor, excluding non sequence-defining ions
     (only 11 main series) within the specified m/z window, charge independent
     beta is determined as:
            if no consecutive matches in the same ion series: beta = 0
            if at least one consecutive match in the same ion series: beta = sum(beta_increment +
            (x_beta * consecutive_matches ) * beta_increment) over all consecutive matches in the same
            ion series (a,a-b,b..) or
            ion series group (e.g. a/b, w/x...). beta_increment and x_beta are tunable via input options
     
     L = normalization factor base on total number of theoretical predicted ions for a precursor ion excluding
     non sequence-defining ions (only 11 main series), wihin the specified m/z window
     """
    return round(sumI * n * (
            1 + beta_a + beta_aB + beta_b + beta_c + beta_d + beta_w + beta_x + beta_y + beta_z + beta_yP + beta_zP)
                 / (L * sumI_all), 3)


def mod_detection(sequence):
    """
    Find the list of modified bases present in the matching MS2 ions, to be used for sumI, n, L calculations
    """
    mod_bases = []
    for nt in sequence:
        if nt != "A" and nt != "C" and nt != "G" and nt != "U" and nt != "X":
            mod_bases.append(nt)

    return mod_bases


def ppm_range(value, ppm):
    """
    Calculate incertitude on MS1/MS2 masses equivalent to given ppm
    """
    return round(ppm / 1000000 * value, 5)


def ppm_offset(measured_mass, theoretical_mass):
    """
    Calculate the ppm offset between matching m/z values for MS1/MS2
    """
    difference = np.float64(measured_mass) - np.float64(theoretical_mass)
    ppm_offset = difference / np.float64(measured_mass) * 1000000

    return round(ppm_offset, 1)


def reduced_digest(digest_peaks):
    """
    Create a reduced version of the digest with only the M-xx ions, to be used in the step of normalization of
    intensities
    """
    output_dic = {}
    for key in digest_peaks:
        a = [x for x in digest_peaks[key][-1] if 'M-' in x.split("(")[0] or len(x.split("(")[0]) == 1]
        output_dic[key] = digest_peaks[key][:-1] + [a]

    return output_dic


def find_losses_freebases(MS2_ions):
    """
    Find the M-xx ions and free bases among the MS2 matches ions per each precursor ion
    """
    output_loss_masses = []
    for x in MS2_ions:
        if 'M-' in x[1].split("(")[0] or len(x[1].split("(")[0]) == 1:
            output_loss_masses.append(np.float64(x[2]))

    return output_loss_masses


def closest_value(mz_list, mz_precursor):
    """
    Return the closest mz value to the precursor m/z within the MS2 list in the mgf file
    This value will be used to output the absolute intensity of the precursor ion
    """

    return min(mz_list, key=lambda x: abs(np.float64(x.split('_')[0]) - np.float64(mz_precursor))).split('_')[1]


def pick_max_intensity(ion_list):
    """
    Return the MS2 ion to pick in case of multiple matching by the same theoretical ion (the most intense one)
    """
    m = sorted(ion_list, key=lambda x: np.float64(x[1]), reverse=True)
    return m[0]


def log_file():
    """
    Write a log file with info about total MS1 matches and 'productive' matches (with at least one MS2 match)
    """
    log_lines = ["# Peaks in mgf file " + str(mgf_peaks) + "\n# Matching precursor ions " + str(
        MS1_hits) + "\n#'Productive' precursor ions (with at least one MS2 hit) " + str(MS2_hits)]

    return log_lines


def csv_output(match_file):
    """
    Output an additional csv file that makes consultation and browsing of the output easier
    """
    csv_lines = []
    with open(match_file, 'r') as infile:
        for line in infile:
            if "PRECURSOR_ION" in line:
                csv_lines.extend(
                    ["\n", "Precursor m/z", "RT(min)", "MS1_offset", "", "Rank", "Score(Sp)", "dSp",
                     "MS2 peaks matches #", "sumI", "sumIall", "L", "sumBeta", "", "Length (nts)", "Species",
                     "Charge", "5'-chemistry", "3'-chemistry", "", "Sequence", "Sequence_ext", "Molecule\n"])

            if line[0].isdigit():
                s = line.split()
                sumb = 0
                for x in s[7].split(';'):
                    if x[0] == 'b':
                        sumb += np.float64(x.split('=')[1])
                csv_lines.extend([s[0], s[1].split('=')[1], s[3], "", s[-2], str(np.float64(s[8])),
                                  str(np.float64(s[6])), s[10].split('=')[1], s[-1].split(';')[0].split('=')[2],
                                  s[-1].split(';')[-1].split('=')[1][:-1], s[-1].split(';')[-2].split('=')[1],
                                  str(round(sumb, 3)),
                                  "", s[9], s[11], s[12], s[15], s[16], "", s[13], s[14], s[17],
                                  " ".join(s[17:-2]) + "\n"])

    return ";".join(csv_lines)


def filename_from_path(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)
