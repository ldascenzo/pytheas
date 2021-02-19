#!/usr/bin/python3

"""
Last update: January 2021
Credits to: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzo@scripps.edu

DESCRIPTION
This script visualizes as annotated spectra the results of theoretical digest files and of the matching/scoring
algorithm.
For the digest all theoretical peaks are represented with an arbitrary intensity.
For the matching/scoring output the spectrum from .mgf file is plotted with annotations on the matching MS2 peaks.
The final output is a html file that contains an output table with links to the spectra in the same page
"""

import os, sys, re
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pyteomics.mgf as pymgf
import multiprocessing
import shutil
import ntpath

params = {
    'figure.dpi': 300,
    'axes.labelsize': 9,
    'font.size': 10,
    'legend.fontsize': 8,
    'legend.frameon': True,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'font.family': 'serif',
    'axes.linewidth': 0.5,
    'xtick.major.size': 4,  # major tick size in points
    'xtick.minor.size': 2,  # minor tick size in points
    'xtick.direction': 'out',
    'ytick.major.size': 4,  # major tick size in points
    'ytick.minor.size': 2,  # minor tick size in points
    'ytick.direction': 'out',
}
plt.rcParams.update(params)

# Color code for MS2 ion series in the output spectra
color_MS2 = {'a': 'green', 'a-B': 'green', 'w': 'green', 'b': 'blue', 'x': 'blue', 'c': 'magenta', 'y': 'magenta',
             'd': 'orange', 'z': 'orange'}

MS2_ion_series = ['a', 'b', 'c', 'd', 'w', 'x', 'y', 'z', 'a-B', 'y-P', 'z-P']

# Parameter used to exclude Ion + 1/2 Na from being a reference peak for scaling intensities
NA_mass, H_mass = 22.989769282, 1.007825032


class Visualize:
    def __init__(self, digest, mgf, match, high_peaks, x_min, x_max, spectra_digest, Sp_cutoff, dSp_cutoff,
                 visualize_decoys, only_modified, rank, unique_positions, remove_redundant_x):
        self.digest_file, self.mgf_file, self.score_file, self.highest_peaks = digest, mgf, match, high_peaks
        self.x_min, self.x_max, self.digest_spectra = x_min, x_max, spectra_digest
        self.Sp_cutoff, self.dSp_cutoff = Sp_cutoff, dSp_cutoff
        self.visualize_decoys, self.only_modified, self.rank = visualize_decoys, only_modified, rank
        self.only_unique_positions, self.redundant_x = unique_positions, remove_redundant_x

    def digest_peaks(self, input_file):
        """
        Creates a dictionary with the info on m/z:intensities from the digest file input
        Args:
            input_file (text file): the digest input file
        Returns:
            peaks (dictionary with peaks stored in pandas dataframes)
        """
        ######NOTE!!!!
        # The first part of this function is necessary to generate the output txt file for the digest peaks,
        ###########

        dic, outlist = {}, []

        # Iterates within the lines of the digest file
        for line in open(input_file, 'r'):

            # Only lines starting with a number are considered (ones with precursor ions)
            if line[0].isdigit():
                # A dictionary key is created with { M_sequence : all the info }
                dic[line.split()[0] + "_" + line.split()[7]] = line.split()[1:]

        # Order the precursor ions by crescent m/z value
        precursors = sorted(dic.keys(), key=lambda x: np.float64(x.split('_')[0]))

        # Creates a dictionary with dataframes containing  all the ions from the digest
        digest_df = digest_table(precursors, dic)

        # First line of the output file carries the info on the precursor ion
        outlist.append(self.digest_file + "\n\n")

        # Preparing the lines for the final output, with header and m/z intensity ion
        for ion in precursors:
            outlist.append(
                "BEGIN IONS\nprecursor= {} seq= {} mod= {} charge= {}\n".format(ion.split('_')[0], ion.split('_')[1],
                                                                                dic[ion][7], dic[ion][4]))

            MS2_ions = sorted(dic[ion][12:], key=lambda x: np.float64(x.split(':')[-1]))

            for i in MS2_ions:
                outlist.append(str(round(np.float64(i.split(':')[-1]), 6)) + "\t100\t" + i.split(':')[0] + "\n")

            outlist.append("END IONS\n")

        flag, out_dic = 0, {}

        for line in outlist:

            if "precursor" in line:
                flag, prec_mass, seq, charge = 1, line.split('=')[1].split()[0], line.split('=')[2].split()[0], \
                                               line.split('=')[4].split()[0]
                out_dic[prec_mass + "_" + seq + "_" + charge] = {'m/z': [], 'intensity': [], 'ion': []}

            if "END IONS" in line:
                flag = 0

            if flag == 1 and line[0].isdigit():
                (out_dic[prec_mass + "_" + seq + "_" + charge]['m/z'].append(np.float64(line.split("\t")[0])),
                 out_dic[prec_mass + "_" + seq + "_" + charge]['intensity'].append(np.float64(line.split("\t")[1])),
                 out_dic[prec_mass + "_" + seq + "_" + charge]['ion'].append(line.split("\t")[2][:-1]))

        # Transform the input data in pandas dataframes
        for key in out_dic.keys():
            df = pd.DataFrame(out_dic[key])
            out_dic[key] = df

        return out_dic, outlist, digest_df

    def scored_info(self, input_file):
        """
        Creates a dictionary from the scored file given as input

        Args:
            input_file: .txt file output of the scoring/matching script MS2_match_seqscore.py

        Returns:
            dic: dictionary with all the matching precursor ions with keys = m/z precursor_rt_seq
        """
        dic, flag, = {}, 0

        sequences_with_x = []

        # Prepare the variables for the test on the options about decoys and modified matches
        if self.visualize_decoys:
            decoys_string = '@#$%'
        else:
            decoys_string = 'decoy'

        if self.only_modified:
            mod_string = '['
        else:
            mod_string = ''

        if self.only_unique_positions:
            unique_length = 1
        else:
            unique_length = 999

        for line in open(input_file, 'r'):

            if "#precursor_window_removal" in line:
                global prec_window_removal
                prec_window_removal = np.float64(line.split()[1])

            if "PRECURSOR_ION" in line:
                prec_ion, flag = line.split()[0].split('=')[1], 1
                dic[prec_ion] = {}

            if flag == 1 and line[0].isdigit():

                split = line.split()
                # Check Sp, dSp and Rank cutoffs if applied
                score, dsp = split[8], split[6]
                molecule, mod = split[17], split[14]
                unique_list = split[17].split(';')
                rank = int(split[-2])

                # Apply all the cutoffs and conditions from user-selected options
                if float(score) >= self.Sp_cutoff and float(dsp) <= self.dSp_cutoff and decoys_string not in molecule\
                        and mod_string in mod and len(unique_list) <= unique_length and rank <= self.rank:

                    prec, rt, th_match = split[0], split[1].split('=')[1], split[2].split('=')[1]
                    MS2_scans, charge, sequence, positions = split[10].split('=')[1], split[12], split[13], \
                                                             split[17]

                    chem3, chem5 = split[15], split[16]

                    # Make sure that only the top scoring/lowest MS1 offset sequence containing 'X' is visualized
                    flag_x = False

                    if self.redundant_x:
                        if 'X' in sequence:
                            if '{}_{}_{}'.format(sequence, prec_ion, rt) not in sequences_with_x:
                                sequences_with_x.append('{}_{}_{}'.format(sequence, prec_ion, rt))
                            else:
                                flag_x = True

                    if not flag_x:

                        sequence_mod = split[14]
                        unique_name = "{}_{}_{}".format(prec, rt, sequence)

                        dic[prec_ion][unique_name] = {}

                        dic[prec_ion][unique_name]['th_match'] = th_match
                        dic[prec_ion][unique_name]['score'] = score
                        dic[prec_ion][unique_name]['dSp'] = dsp
                        dic[prec_ion][unique_name]['MS2_scans'] = MS2_scans
                        dic[prec_ion][unique_name]['charge'] = charge
                        dic[prec_ion][unique_name]['sequence'] = sequence
                        dic[prec_ion][unique_name]['sequence_mod'] = sequence_mod
                        dic[prec_ion][unique_name]['positions'] = positions
                        dic[prec_ion][unique_name]['MS2_match'] = []
                        dic[prec_ion][unique_name]['3chemistry'] = chem3
                        dic[prec_ion][unique_name]['5chemistry'] = chem5

                        for scan in split[18:-2]:
                            mgf_peak, relative_intensity, th_peak, ion = scan.split(':')[0].split('(')[0], \
                                                                         scan.split(']')[0].split('[')[1], \
                                                                         scan.split(':')[1].split('[')[0], \
                                                                         scan.split('[')[2].split(']')[0]
                            dic[prec_ion][unique_name]['MS2_match'].append((mgf_peak, relative_intensity, th_peak, ion))

        return dic

    def df_output_html(self):
        """
         Exports the data from the matching output file as a html table for the final output
         """
        infile = self.score_file
        # Define a dictionary with all the column to be included in the dataframe
        d = {'m/z': [], 'Rank': [], 'RT': [], 'MS1_offset(ppm)': [], 'length': [], "5'-end": [], 'sequence': [],
             "3'-end": [], 'charge': [], 'isotope': [], 'Score (Sp)': [], 'dSp': [], 'n/L': [], 'sequence_mods': [],
             'molecule': []}

        sequences_with_x = []

        # Prepare the variables for the test on the options about decoys and modified matches
        if self.visualize_decoys:
            decoys_string = '@#$%'
        else:
            decoys_string = 'decoy'

        if self.only_modified:
            mod_string = '['
        else:
            mod_string = ''

        if self.only_unique_positions:
            unique_length = 1
        else:
            unique_length = 999

        with open(infile, 'r') as input_file:
            for line in input_file:

                if 'PRECURSOR_ION' in line:
                    unique_precursor = line.split('=')[1]

                if line[0].isdigit():
                    sp = line.split()

                    score, dSp, molecule, mod = sp[8], sp[6], sp[17], sp[14]
                    unique_list = sp[17].split(';')
                    rank = int(sp[-2])

                    # Only include data within requested cutoffs for Sp and dSp
                    if float(score) >= self.Sp_cutoff and float(dSp) <= self.dSp_cutoff \
                            and decoys_string not in molecule and mod_string in mod and \
                            len(unique_list) <= unique_length and rank <= self.rank:
                        mz, rt, seq = sp[0], sp[1].split('=')[1], sp[13]

                        # Make sure that only the top scoring/lowest MS1 offset sequence containing 'X' is
                        # used for output
                        flag_x = False

                        if self.redundant_x:
                            if 'X' in seq:
                                if '{}_{}_{}'.format(seq, unique_precursor, rt) not in sequences_with_x:
                                    sequences_with_x.append('{}_{}_{}'.format(seq, unique_precursor, rt))
                                else:
                                    flag_x = True

                        if not flag_x:
                            # Insert the m/z value with an anchor tag for the html table output directing
                            # towards the respective spectrum
                            if '*' in mz:
                                mz_isotopologues_correction = mz[:-1]
                            else:
                                mz_isotopologues_correction = mz

                            d['m/z'].append('<a href="#{}_{}_{}_{}">{}</a>'.format(mz_isotopologues_correction, rt,
                                                                                   score, seq, mz))

                            # Insert all the other info for all ions in preparation for the table output
                            d['Rank'].append(sp[-2])
                            d['RT'].append(rt), d['Score (Sp)'].append(score), d['dSp'].append(dSp), d['charge'].append(
                                sp[12]),
                            d['isotope'].append(sp[11]), d['molecule'].append(sp[17]), d['sequence'].append(seq), d[
                                'sequence_mods'].append(sp[14]),
                            d["5'-end"].append(sp[16]), d["3'-end"].append(sp[15]), d['MS1_offset(ppm)'].append(
                                sp[3].split('p')[0]), d['length'].append(len(seq)), \
                            d['n/L'].append(round(int(sp[-1].split(';')[1].split('=')[1]) /
                                                  int(sp[-1].split(';')[-2].split('=')[1]), 3))

        df = pd.DataFrame(data=d)
        df.index += 1
        pd.set_option('display.max_colwidth', 999)

        return df.to_html(classes='sortable', justify='center', escape=False, table_id='output_table')

    def digest_spectra_multicpu(self, key):
        """
        Generate the spectra plots for the digest
        """
        input_values = self.digest_peaks(self.digest_file)[0]

        # Axis limits
        x_min = plot_limits(input_values, key)[0]
        x_max = plot_limits(input_values, key)[1]
        y_max = plot_limits(input_values, key)[2]

        fig = plt.figure(1, figsize=(conv_inch(self.figwidth_mm), conv_inch(self.figheight_mm)))
        ax = fig.add_subplot(111, autoscale_on=False, xlim=(x_min, x_max), ylim=(0, y_max))

        bars = ax.bar(input_values[key]['m/z'], input_values[key]['intensity'], color='grey', width=self.bars_width,
                      alpha=.5)
        plt.xticks(np.arange(x_min, x_max + 1, 25), rotation=90)
        ax.set_xlabel('m/z')
        ax.set_ylabel('Intensity')
        ax.set_title('precursor=' + key.split('_')[0] + '  seq=' + key.split('_')[1] + '  charge=' + key.split('_')[2])

        ###### Annotation
        # Iterates among all the rows of the dataframe with m/z, intensities and ion
        texts = []
        for index, row in input_values[key].iterrows():

            if row['ion'][0] in color_MS2.keys():
                texts.append(plt.text(row['m/z'], 100, row['ion'] + " " + str(round(row['m/z'], 3)), rotation=90,
                                      color=color_MS2[row['ion'][0]], size=5))
                plt.axvline(np.float64(row['m/z']), ymax=100 / 150, color=color_MS2[row['ion'][0]],
                            linewidth=self.bars_width,
                            linestyle='dashed')

            else:
                texts.append(
                    plt.text(row['m/z'], 100, row['ion'] + " " + str(round(row['m/z'], 3)), rotation=90, color='black',
                             size=5))
                plt.axvline(np.float64(row['m/z']), ymax=100 / 150, color='black', linewidth=self.bars_width,
                            linestyle='dashed')

        # adjust_text(texts, add_objects = bars, expand_text = (0.75, 0.75),
        # arrowprops=dict(arrowstyle="simple, head_width=0.15, tail_width=0.05", color = 'grey', lw=0.3, alpha=0.5))

        fig.canvas.draw()
        fig.tight_layout()
        plt.savefig("./digest_spectra/" + '_'.join(key.split('_')[:2]) + '.png')
        # plt.savefig("./digest_spectra/" + '_'.join(key.split('_')[:2]) + '.pdf')
        plt.close()

    def scored_spectra(self, match):
        """
        Generate the spectra plots for the scored/matched results
        """
        # Only matches with at least one MS2 ion and Sp score > 0 can be plotted
        if match[1]['MS2_match'] and np.float64(match[1]['score']) > 0:

            input_mgf = mgf_peaks(self.mgf_file)
            key_mgf = re.findall("[\d]+[.]?[\d]+", match[0].split('_')[0])[0] + '_' + match[0].split('_')[1]

            MS2_matches, charge, prec, mz_highest_ion = match[1]['MS2_match'], int(
                re.findall(r'\d+', match[1]['charge'])[0]), np.float64(
                re.findall("[\d]+[.]?[\d]+", match[0].split('_')[0])[0]), 0

            # Determines the sequence ion with highest intensity, used to scale the y axis of the spectra
            for single_ion in MS2_matches:
                if single_ion[3][0] != 'M' and len(single_ion[3].split('(')[0]) > 1:
                    mz_highest_ion = np.float64(single_ion[0])
                    break

            # AXIS LIMITS #
            # Limit on Y defined as the max value of peaks intensity without precursor + a set amount
            peak_max = input_mgf[key_mgf].loc[input_mgf[key_mgf]['m/z'] == mz_highest_ion]['intensity'] * 1.175
            if peak_max.empty:
                y_max = round_number(max(input_mgf[key_mgf]['intensity']) * 1.175)
            else:
                y_max = round_number(input_mgf[key_mgf].loc[input_mgf[key_mgf]['m/z'] ==
                                                            mz_highest_ion]['intensity'] * 1.175)

            # Limits on X are specified from the minimum and maximum m/z values for matching ions
            sorted_MS2 = sorted(MS2_matches, key=lambda x: np.float64(x[0]))

            if self.x_min == -1:
                x_min = round_number(sorted_MS2[0][0]) - 50

                if x_min >= prec:
                    x_min = prec - 50
            else:

                if self.x_min >= prec:
                    x_min = prec - 50

                else:
                    x_min = self.x_min

            if self.x_max == -1:

                # Set the x axis maximum as the m/z value of the highest ion in the mgf scan
                x_max = round_number(input_mgf[key_mgf]['m/z'].max()) + 50

                if x_max <= prec:
                    x_max = prec + 50
            else:
                if self.x_max <= prec:
                    x_max = prec + 50

                else:
                    x_max = self.x_max

            # Generates the spectrum axis with given plot dimensions
            fig = plt.figure(1, figsize=(conv_inch(self.figwidth_mm), conv_inch(self.figheight_mm)))
            ax = fig.add_subplot(111, autoscale_on=False, xlim=(x_min, x_max), ylim=(0, y_max))

            # Select only the highest n (option highest_peaks) peaks to be shown
            if self.highest_peaks != 'all':
                red_dic = input_mgf[key_mgf].nlargest(self.highest_peaks, 'intensity')
                limit_peaks = self.highest_peaks
            else:
                red_dic = input_mgf[key_mgf].nlargest(len(input_mgf[key_mgf]), 'intensity')
                limit_peaks = 99999

            # Reset the index of the dataframe (avoid an error with ax.bar missing the index
            red_dic = red_dic.reset_index(drop=True)
            bars = ax.bar(red_dic['m/z'], red_dic['intensity'], color='grey', width=self.bars_width, alpha=.7)

            plt.xticks(np.arange(x_min, x_max + 1, 25), rotation=90)
            ax.set_xlabel('m/z')
            ax.set_ylabel('Intensity')

            # Add title with info on the matching ion in case the match is not a decoy
            if 'decoy' in match[1]['positions']:
                molecule, resn = 'decoy', '-'

            else:
                molecule, resn = match[1]['positions'].split(",")[0], "{}-{}".format(
                    match[1]['positions'].split(",")[1],
                    match[1]['positions'].split(",")[2])

            if match[1]['sequence_mod'] != '-':
                ax.set_title(
                    "M={}  RT={}  seq=$\it{}$-{}-$\it{}$  seq_mod=$\it{}$-{}-$\it{}$  Sp={}  dSp={}  charge={}  "
                    "molecule={}  resn={}".format(
                        match[0].split('_')[0], match[0].split('_')[1], match[1]['5chemistry'], match[1]['sequence'],
                        match[1]['3chemistry'], match[1]['5chemistry'], match[1]['sequence_mod'],
                        match[1]['3chemistry'],
                        match[1]['score'], match[1]['dSp'], match[1]['charge'], molecule, resn))

            else:
                ax.set_title(
                    "M={}  RT={}  seq=$\it{}$-{}-$\it{}$  Sp={}  dSp={}  charge={}  molecule={}  "
                    "resn={}".format(match[0].split('_')[0],
                                     match[0].split('_')[1],
                                     match[1]['5chemistry'],
                                     match[1]['sequence'],
                                     match[1]['3chemistry'],
                                     match[1]['score'],
                                     match[1]['dSp'],
                                     match[1]['charge'],
                                     molecule, resn))

            ##### Annotation
            # Iterates among all the rows of the dataframe with m/z, intensities and ion
            texts = []

            ### Steps to update the MS2 ions table with matches
            for index, i in enumerate(MS2_matches):

                if index <= limit_peaks:

                    # Add the labels to matching cells in the output table with corresponding series, changing
                    # the value with the theoretical  m/z
                    m = re.search(r'\d+$', i[3].split("(")[0].split('-')[0])
                    if m:
                        ion_series = i[3].split("(")[0].replace(m.group(), '') + "(" + i[3].split(':')[0].split("(")[1]
                        pd.set_option('display.max_columns', None)

                        if i[3][0] == 'a' or i[3][0] == 'b' or i[3][0] == 'c' or i[3][0] == 'd':

                            dig_tab[match[1]['th_match'] + '_' + match[1]['sequence']].loc[
                                int(m.group()) - 1, ion_series] = str(round(np.float64(i[2]), 4)) + '_' + i[3][
                                0] + 'clr'

                        else:
                            dig_tab[match[1]['th_match'] + '_' + match[1]['sequence']].loc[
                                len(match[1]['sequence']) - int(m.group()), ion_series] = str(
                                round(np.float64(i[2]), 4)) + '_' + i[3][0] + 'clr'

                    # Add the labels to matching cells in the output table for free bases and losses, changing
                    # the value with the theoretical  m/z
                    if re.search(r'^(?![abcdzyxw])', i[3].split("(")[0].split('-')[0]):
                        dig_tab[match[1]['th_match'] + '_' + match[1]['sequence']].loc[0, i[3]] = str(
                            round(np.float64(i[2]), 4)) + '_Mclr'

                    # Annotate each matching bar with info on ion + m/z,
                    # using a color based on the ion series it belongs
                    prec_fl = np.float64(i[0])

                    # Make sure the precursor ion is in the MS2 matches m/z values to avoid errors
                    if not input_mgf[key_mgf].loc[input_mgf[key_mgf]['m/z'] == prec_fl].empty and np.float64(
                            i[0]) < x_max and np.float64(i[0]) > x_min:

                        y_lim_annotation = int(
                            input_mgf[key_mgf].loc[input_mgf[key_mgf]['m/z'] == prec_fl, 'intensity'].iloc[0])

                        # Set the annotation for ions to be shown even when they are offscale
                        if y_lim_annotation > y_max:
                            y_coord_ions, x_coord_ions = y_max * 0.95, np.float64(i[0]) * 1.02 + self.bars_width * 4
                        else:
                            y_coord_ions, x_coord_ions = y_lim_annotation * 1.01, np.float64(i[0])

                        if i[3][0] in color_MS2.keys() and re.search(r'\d+$', i[3].split("(")[0].split('-')[0]):

                            texts.append(
                                plt.text(x_coord_ions, y_coord_ions, i[3], color=color_MS2[i[3][0]], size=5, rotation=0,
                                         ha='center'))

                            # Add dashed colored bar to matching ions
                            plt.axvline(np.float64(i[0]), ymax=y_lim_annotation / y_max, color=color_MS2[i[3][0]],
                                        linewidth=self.bars_width, linestyle='dashed')

                        else:
                            texts.append(
                                plt.text(x_coord_ions, y_coord_ions, i[3], color='black', size=5, rotation=0,
                                         ha='center'))

                            # Add dashed colored bar to matching ions
                            plt.axvline(np.float64(i[0]), ymax=y_lim_annotation / y_max, color='black',
                                        linewidth=self.bars_width, linestyle='dashed')

            # Add precursor ion dashed line
            prec_x = input_mgf[key_mgf].iloc[(input_mgf[key_mgf]['m/z'] - prec).abs().argsort()[:1]]['m/z']
            prec_y = input_mgf[key_mgf].loc[input_mgf[key_mgf]['m/z'] == np.float64(prec_x), 'intensity']

            plt.axvline(np.float64(prec_x), ymax=np.float64(prec_y) / y_max, color='black', linewidth=self.bars_width,
                        linestyle='-')

            # Add a label for the precursor ion
            if np.float64(prec_y) < y_max:
                ax.annotate("M({})".format(match[1]['charge']), xy=(np.float64(prec_x), 1.01 * np.float64(prec_y)),
                            xycoords='data',
                            xytext=(np.float64(prec_x) + self.bars_width * 4, np.float64(prec_y)), textcoords='data',
                            size=7,
                            ha='center')
            else:
                ax.annotate("M({})".format(match[1]['charge']), xy=(np.float64(prec_x) * 1.02, y_max * 0.95),
                            xycoords='data',
                            xytext=(np.float64(prec_x) * 1.01 + self.bars_width * 4, y_max * 0.95), textcoords='data',
                            size=7,
                            ha='center')

            plt.savefig(
                "./scored_spectra_" + filename_from_path(self.mgf_file)[:-4] + "/" + key_mgf + '_' + match[1][
                    'score'] + '_' + match[1][
                    'sequence'] + '.png')
            plt.close()

            # Creates the tabled for the final output from the pandas dataframes, separating ion series from losses and
            # free bases
            df_ions = dig_tab[match[1]['th_match'] + '_' + match[1]['sequence']].filter(regex='^[abcdzyxw]').set_index(
                keys=np.array(list(match[1]['sequence'])))
            df_losses = dig_tab[match[1]['th_match'] + '_' + match[1]['sequence']].filter(regex='^(?![abcdzyxw])').iloc[
                [0]]

            # Order the columns alphabetically based on the header row
            df_ions = df_ions.reindex(sorted(df_ions.columns), axis=1)
            df_losses = df_losses.reindex(sorted(df_losses.columns, key=len), axis=1)

            # Create the html lines for the spectra + tables
            self.html_lines.append("""<div class="wrapper">
                                      <img class="center-fit" id="{}_{}_{}" />
                                  </div>
                                  """.format(key_mgf, match[1]['score'], match[1]['sequence']))

            self.html_lines.extend(
                (df_ions.to_html(justify='center', classes='scrollable') + "\n\n" + df_losses.to_html(
                    index=False, justify='center', classes='scrollable'),
                 '<br>\n', '<br>\n', '<hr>\n', '<br>\n', '<br>\n'))

    def html_css_header(self):
        """
        Writes the header for the output html  file with all the visualized spectra
        """
        return """<style>
                * {
                margin: 0;
                padding: 0;
            }

            .center-fit {
                max-width: 100%;
                max-height: 100vh;
                margin: auto;
            }

            #output_table {
            border:2px solid black;
            width:95%;
            border-collapse:collapse;
            margin:5;
        }

            table, th, td {

            padding:4px;
            }

        th {
            background: #dddddd;
            position: sticky;
            top: 0;
            }
            
        table tbody tr:hover {
            background-color: #dddddd;
            }

        hr {
                border: 10px solid green;
                border-radius: 5px;
            }
        p { 
        padding: 15px; }

        table.sortable thead {
        background-color:#eee;
        color:#666666;
        font-weight: bold;
        cursor: default;
        }	

        table.scrollable {
        display: block;
        overflow-x: auto;
        white-space: nowrap;
        }
        
        #ResultTable {
        max-height: 1000 px;
        overflow: auto;
        }
        </style>

       <head>
       <title>""" + filename_from_path(self.mgf_file)[:-4] + """</title>
        <p><strong>Dataset: """ + filename_from_path(self.mgf_file)[:-4] + """</strong></p>
       </head>
       
       <div id=ResultTable">
        """

    def final_output(self):
        if self.x_min > self.x_max:
            print("x_min value must be smaller than x_max!!!!")
            sys.exit()

        # Dimensions of the plots
        self.figwidth_mm = 400
        self.figheight_mm = self.figwidth_mm * 0.3
        self.bars_width = 0.5

        output_name = filename_from_path(self.mgf_file)[:-4]

        create_directory('scored_spectra_{}'.format(output_name))

        # Table with all the MS2 ions from the digest, for final output
        global dig_tab
        dig_tab = self.digest_peaks(self.digest_file)[2]

        input_scored = self.scored_info(self.score_file)

        self.html_lines = [self.html_css_header(), self.df_output_html(), "</div>"]

        for key in input_scored:

            # Sort the input entries by score
            d = input_scored[key]
            score_sorted = sorted(d.items(), key=lambda x: np.float64(x[1]['score']), reverse=True)

            for match in score_sorted:
                self.scored_spectra(match)

        # Adds the final lines to call the Javascript script to color the matching tables
        self.html_lines.append("""

                <script src='js/visualize_spectra.js'></script>       
                <script src='js/color_tables.js'> </script>
                <script src='js/sorttable.js'></script>
                """)

        with open("visualization_{}.html".format(output_name), 'w') as _file:
            _file.writelines(self.html_lines)

        if self.digest_spectra:
            # Plot both digest and scored spectra in case all three optional files are given
            open("./parsed_digest.txt", 'w').writelines(self.digest_peaks(self.digest_file)[1])

            create_directory('digest_spectra')

            key_list = []

            for key in self.digest_peaks(self.digest_file)[0].keys():
                key_list.append(key)

                p = multiprocessing.Pool()
                p.map(self.digest_spectra_multicpu, key_list)
                p.close()


def round_number(x, base=25):
    """
    Round a number to its closest integer value multiple of 25. Used to determine values for the axis
    on the spectra plots
    """
    return int(base * round(np.float64(x) / base))


def conv_inch(length_mm):
    """
    Converts a length from millimeters to inch
    """
    return length_mm / 25.4


def plot_limits(xy_dataframe, dic_key):
    """
    Determines the limits for x and y axes, defined based on the min and max values of m/z (x) and intensity (y)
    from the input file
    """
    return round_number(min(xy_dataframe[dic_key]['m/z'])) - 25, round_number(
        max(xy_dataframe[dic_key]['m/z'])) + 25, round_number(max(xy_dataframe[dic_key]['intensity'])) + 50


def annotation(row, ax, bars_width):
    """
    Annotate the spectrum with ion series and a dashed line over the bar
    """
    ax.annotate(row['ion'] + " " + str(round(row['m/z'], 3)), xy=(row['m/z'], 105), xycoords='data',
                xytext=(row['m/z'] - bars_width * 4, 130), textcoords='data', rotation=90,
                color=color_MS2[row['ion'][0]], size=5)

    plt.axvline(np.float64(row['m/z']), ymax=100 / 150, color=color_MS2[row['ion'][0]], linewidth=bars_width,
                linestyle='dashed')


def create_directory(name):
    """
    Create the directory in which to output the spectra
    """
    if os.path.exists("./{}/".format(name)):
        shutil.rmtree("./{}".format(name))

    os.mkdir("./{}".format(name))


def digest_table(p, dic):
    """
    Generates a dataframe with all the info to be output as a final table with all matching ions
    """
    out_dic = {}
    for precursor in p:

        seq_length, d = len(dic[precursor][6]), {}

        # Creates an empty dictionary for all the MS2 ion series 
        for ion in dic[precursor][12:]:

            m = re.search(r'\d+$', ion.split("(")[0].split('-')[0])
            if m:
                ion_series = ion.split("(")[0].replace(m.group(), '') + "(" + ion.split(':')[0].split("(")[1]

                if ion_series not in d.keys():
                    d[ion_series] = [''] * seq_length

            else:
                ion_series = ion.split(':')[0]
                if ion_series not in d.keys():
                    d[ion_series] = [''] * seq_length

        # Fills the dictionary with m/z values on  
        for ion in dic[precursor][12:]:

            m = re.search(r'\d+$', ion.split("(")[0].split('-')[0])
            if m:
                ion_series = ion.split("(")[0].replace(m.group(), '') + "(" + ion.split(':')[0].split("(")[1]
                if ion[0] == 'a' or ion[0] == 'b' or ion[0] == 'c' or ion[0] == 'd':
                    d[ion_series][int(m.group()) - 1] = str(round(np.float64(ion.split(':')[1]), 4))

                else:
                    d[ion_series][seq_length - int(m.group())] = str(round(np.float64(ion.split(':')[1]), 4))

            else:
                ion_series = ion.split(':')[0]

                # Add a control to account if some of the ion series names are used as modified bases
                if ion_series.split("(")[0] in MS2_ion_series:
                    ion_series = "B" + ion_series

                d[ion_series] = str(round(np.float64(ion.split(':')[1]), 4))

        # Convert the dictionary with all theoretical MS2 ions into a dataframe
        df = pd.DataFrame(data=d)

        out_dic[precursor] = df

    return out_dic


def mgf_peaks(input_file):
    """
    Creates a dictionary with the info on m/z:intensities from the mgf file input

    Args:
        input_file (text file): the mgf input file

    Returns:
        peaks (dictionary with peaks stored in pandas dataframes)
    """
    dic, out_dic = {}, {}

    # Correct a problem with the charges being non int which crashes pyteomics parser
    mgf_lines, outlines = open(input_file, 'r').readlines(), []
    for line in mgf_lines:
        if "CHARGE" in line:
            charge = re.findall(r'\d+', line.split('=')[1])[0]
            outlines.append("CHARGE={}\n".format(charge))
        else:
            outlines.append(line)

    name_output = filename_from_path(input_file)
    open("temp_{}".format(name_output), 'w').writelines(outlines)

    # Create a dictionary with all the precursor ions_rt as key
    for scan in pymgf.read("temp_{}".format(name_output)):
        prec_mass, rt = '{:.6f}'.format(scan['params']['pepmass'][0]), round(
            np.float64(scan['params']['rtinseconds']) / 60, 3)
        dic[str(prec_mass) + "_" + str(rt)] = scan

    # Transform the m/z and intensity data in pandas dataframe format
    for key in dic:
        out_dic[key] = pd.DataFrame({'m/z': np.round_(dic[key]['m/z array'], decimals=6),
                                     'intensity': np.round_(dic[key]['intensity array'], decimals=6)},
                                    dtype=np.float64)

    try:
        os.remove("temp_{}".format(input_file.split('/')[-1]))
    except OSError:
        pass

    return out_dic


def filename_from_path(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)
