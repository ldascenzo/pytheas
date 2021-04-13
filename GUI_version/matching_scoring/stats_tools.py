#!/usr/bin/python3

"""
Last update: April 2021
Author: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzoluigi@gmail.com
GitHub project repository: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
A collection of tools used for FDR estimation, generation of the targets and decoys output files and
the generation of the statistical plots in the Pytheas workflow
"""

import sys, os, re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import statistics

np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


def read_excel_input(nts_file):
    """
     Produce a dataframe with all the info on the nucleobases from the input file nts_alphabet_light
     """
    # Checking that the nts_alphabet_light file given in argument exists
    if not os.path.exists(nts_file):
        print("ERROR! File " + nts_file + " does not exist. Execution terminated without generating any output."
                                          "This file is generated during the Pytheas in silico digestion")
        sys.exit(1)

    # Create a dataframe with info from Excel spreadsheet
    df = pd.read_excel(nts_file, sheet_name=[0, 1], header=0)

    # Drop rows with additional legend indications
    df[0] = df[0][df[0].length != 'l']
    df[1] = df[1][df[1].length != 'l']

    # Transform all ID values in string
    df[0] = df[0].astype({"m/z obs": str, "RT": str, "z": str})
    df[1] = df[1].astype({"m/z obs": str, "RT": str, "z": str})

    return df


def input_data(infile, excel, sumi_all):
    """
    Extract the data from the input match file creating a pandas dataframe
    """
    # Define a dictionary with all the column to be included in the dataframe
    if sumi_all == 'n':
        d = {'m/z': [], 'RT': [], 'm/z_RT': [], 'rank': [], 'Theoretical_m/z': [], 'isotopologue': [], 'length': [],
             'Sp': [], 'dSp': [], 'n': [], 'sumI': [], 'L': [], 'n/L': [], 'Beta': [], 'charge': [], 'isotope': [],
             'molecule_location': [], 'sequence': [], 'sequence_mods': [],
             'MS1_ppm': [], 'MS2_ppm': [], 'MS2_matches': []}

    else:
        d = {'m/z': [], 'RT': [], 'm/z_RT': [], 'rank': [], 'Theoretical_m/z': [], 'isotopologue': [], 'length': [],
             'Sp': [], 'dSp': [], 'n': [], 'sumI': [], 'sumI_all': [], 'L': [], 'n/L': [], 'Beta': [], 'charge': [],
             'isotope': [], 'molecule_location': [], 'sequence': [],
             'sequence_mods': [], 'MS1_ppm': [], 'MS2_ppm': [], 'MS2_matches': []}

    with open("./" + infile, 'rU') as input_file:
        for line in input_file:

            if line[0].isdigit():
                sp = line.split()

                if '*' in sp[0]:
                    d['isotopologue'].append('y')
                else:
                    d['isotopologue'].append('')

                prec_M = re.findall("[\d]+[.]?[\d]+", sp[0])[0]
                d['m/z'].append(prec_M), d['length'].append(int(sp[9])), d['Sp'].append(
                    np.float64(sp[4].split('=')[1])),
                d['dSp'].append(sp[5].split('=')[1]), d['charge'].append(sp[10]),
                d['isotope'].append(sp[8]), d['molecule_location'].append(sp[15]), d['sequence'].append(sp[11]),
                d['sequence_mods'].append(sp[12]), d['RT'].append(sp[1].split('=')[1]),
                d['m/z_RT'].append(prec_M + '_' + sp[1].split('=')[1]), d['MS1_ppm'].append(sp[3][:-3]),
                d['MS2_matches'].append(",".join(sp[16:-1])),
                d['Theoretical_m/z'].append(round(np.float64(sp[2].split('=')[1]), 5))
                d['rank'].append(int(sp[6].split('=')[1]))

                # Obtain a  list of MS2 matches ppm offsets
                MS2_ppm_list = []
                for match in sp[16:-1]:
                    MS2_ppm_list.append(match.split(":")[0].split(")")[0].split("(")[1][:-3])

                d['MS2_ppm'].append(",".join(MS2_ppm_list))

                terms = sp[-1].split("(")[1].split(";")
                b = 0

                for t in terms:

                    if t[0] == 'n':
                        n = int(t.split("=")[1])
                        d['n'].append(n)

                    if t[0] == 's':
                        d['sumI'].append(int(t.split("=")[1]))

                    if sumi_all == 'y':
                        if t[0] == 'S':
                            d['sumI_all'].append(int(t.split("=")[1][:-1]))

                    if t[0] == 'L':
                        if t[-1] == ')':
                            L = int(t.split("=")[1][:-1])

                        else:
                            L = int(t.split("=")[1])

                        d['L'].append(L)

                    if t[0] == 'b':
                        b += np.float64(t.split("=")[1])

                d['Beta'].append(np.float64(b))
                d['n/L'].append(round(n / L, 3))

    df = pd.DataFrame(data=d)

    # Add info about the maximum intensity of precursor ions.
    # Also clusters the ions based on length + charge values
    if excel != 'n':
        df_excel = read_excel_input(excel)
        for index, row in df.iterrows():

            max_int = \
                df_excel[0][
                    (df_excel[0]['m/z obs'].str.contains(row['m/z'])) & (df_excel[0]['RT'].str.contains(row['RT']))][
                    'Prec_max_int']
            for m in max_int:
                df.loc[df.index[index], 'Prec_max_int'] = m

            max_int = \
                df_excel[1][
                    (df_excel[1]['m/z obs'].str.contains(row['m/z'])) & (df_excel[1]['RT'].str.contains(row['RT']))][
                    'Prec_max_int']
            for m in max_int:
                df.loc[df.index[index], 'Prec_max_int'] = m

    z = [str(row.length) + "_" + str(row.charge) for index, row in df.iterrows()]
    df['length_z'] = z

    df = df.astype({"m/z": np.float64, 'dSp': np.float64})

    return df


def filter_data(df, top, yl):
    """
    Filter the data in preparation for the final graphic output
    """
    if top == 'y':
        filt_df = df.sort_values(['m/z', yl], ascending=[True, False]).reset_index(drop=True)

        # Terminate the script if no matches are found
        if filt_df.empty:
            print("ERROR!!!!!! No target matches found. Execution terminated without output")
            sys.exit(1)

        # Extract the lines with decoys
        decoy_df = filt_df[filt_df['molecule_location'].str.contains('decoy')]

        # Keep only the highest ranking decoy per unique combination m/z_RT
        decoy_df = decoy_df.sort_values(['Sp', 'rank'], ascending=[False, True]).drop_duplicates(subset=['m/z', 'RT'])
        decoy_df = decoy_df.sort_values(['m/z', yl], ascending=[True, False]).reset_index(drop=True)

        # Remove the decoys from the targets filtered dataframe
        filt_df = filt_df[filt_df['molecule_location'].str.contains('decoy') == False]

        # Create a dataframe with all the "nontop" scoring targets
        dup_df = filt_df.sort_values(['m/z', 'Sp'], ascending=[True, False])
        dup_df = dup_df.drop_duplicates(subset=['m/z_RT', 'sequence'])
        dup_df = dup_df[dup_df.duplicated(subset=['m/z', 'RT'])]
        dup_df = dup_df.sort_values('Sp', ascending=False).drop_duplicates(subset=['m/z', 'RT'])
        dup_df['dSp2'] = dup_df['dSp']
        dup_df = dup_df[['m/z_RT', 'dSp2']]

        # Keep only the highest ranking target per unique combination m/z_RT
        filt_df = filt_df.sort_values(['Sp', 'rank'], ascending=[False, True]).drop_duplicates(subset=['m/z', 'RT'])
        filt_df = filt_df.sort_values(['m/z', yl], ascending=[True, False]).reset_index(drop=True)

        # Add a column to the targets dataframe with the values of the second best dSp
        filt_df = pd.merge(filt_df, dup_df, on='m/z_RT', how='outer')
        filt_df = filt_df[
            ['m/z', 'RT', 'm/z_RT', 'Theoretical_m/z', 'isotopologue', 'length', 'sequence', 'sequence_mods',
             'rank', 'Sp', 'dSp', 'dSp2', 'MS1_ppm', 'n', 'sumI', 'sumI_all', 'L', 'n/L', 'Beta', 'charge',
             'isotope', 'molecule_location', 'length_z', 'MS2_ppm', 'MS2_matches']]

    else:
        filt_df = df

        # Terminate the script if no matches are found
        if filt_df.empty:
            print("ERROR!!!!! No target matches found. Execution terminated without output")
            sys.exit(1)

        decoy_df = filt_df[filt_df['molecule_location'].str.contains('decoy')]

    # Add a column with info on targets having at least one competing decoy
    filt_df.insert(8, 'has_decoy', '')
    filt_df.loc[filt_df['m/z_RT'].isin(decoy_df['m/z_RT']), 'has_decoy'] = 'y'

    return filt_df, decoy_df


def csv_output(df, lengths, per, Sp_cutoff, name_output, targets_with_decoys, isotopes):
    """
     Generate the output .csv file with the values used to create the graphs
     &
     FDR file
     """
    in_plot_df, in_decoy_df = df

    # Select the matches of the isotopic channel selected by the user for the FDR and csv output of targets and decoys
    if isotopes == 'light':
        in_plot_df, in_decoy_df = in_plot_df[in_plot_df['isotope'] == 'light'], \
                                  in_decoy_df[in_decoy_df['isotope'] == 'light']
    elif isotopes == 'heavy':
        in_plot_df, in_decoy_df = in_plot_df[in_plot_df['isotope'] == 'heavy'], \
                                  in_decoy_df[in_decoy_df['isotope'] == 'heavy']
    else:
        pass

    plot_df = in_plot_df.loc[in_plot_df['Sp'] > Sp_cutoff]

    # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
    # are always kept
    decoy_df = in_decoy_df[in_decoy_df['m/z_RT'].isin(plot_df['m/z_RT'])]

    decoy_df = decoy_df[
        ['m/z', 'RT', 'm/z_RT', 'Theoretical_m/z', 'isotopologue', 'length', 'sequence', 'sequence_mods',
         'rank', 'Sp', 'dSp', 'MS1_ppm', 'n', 'sumI', 'sumI_all', 'L', 'n/L', 'Beta', 'charge',
         'isotope', 'molecule_location', 'length_z', 'MS2_ppm', 'MS2_matches']]

    plot_df.to_csv("./targets_{}.csv".format(name_output), index=False)
    decoy_df.to_csv("./decoys_{}.csv".format(name_output), index=False)

    # Output a csv with considerations on how many target matches are above a certain level that is xxth percentile
    # of decoys Sp
    if per == 'y':
        p99, p95 = decoy_df.quantile(0.99)['Sp'], decoy_df.quantile(0.95)['Sp']
        p90 = decoy_df.quantile(0.9)['Sp']
        p80 = decoy_df.quantile(0.8)['Sp']
        p75 = decoy_df.quantile(0.75)['Sp']
        p50 = decoy_df.quantile(0.5)['Sp']
        p0 = decoy_df.quantile(0)['Sp']

        nts, d = ['all'] + sorted(plot_df['length'].unique()), {'length': [], '99th': [], '95th': [], '90th': [],
                                                                '80th': [],
                                                                '75th': [], '50th': [], '0th': []}
        for i in nts:
            d['length'].append(i)
            if i != 'all':
                d['99th'].append(len(plot_df.loc[(plot_df['Sp'] > p99) & (plot_df['length'] == i)]['m/z'])), d[
                    '95th'].append(len(plot_df.loc[(plot_df['Sp'] > p95) & (plot_df['length'] == i)]['m/z']))
                d['90th'].append(len(plot_df.loc[(plot_df['Sp'] > p90) & (plot_df['length'] == i)]['m/z'])), d[
                    '80th'].append(len(plot_df.loc[(plot_df['Sp'] > p80) & (plot_df['length'] == i)]['m/z']))
                d['75th'].append(len(plot_df.loc[(plot_df['Sp'] > p75) & (plot_df['length'] == i)]['m/z'])), d[
                    '50th'].append(len(plot_df.loc[(plot_df['Sp'] > p50) & (plot_df['length'] == i)]['m/z']))
                d['0th'].append(len(plot_df.loc[(plot_df['Sp'] > p0) & (plot_df['length'] == i)]['m/z']))

            else:
                d['99th'].append(len(plot_df.loc[plot_df['Sp'] > p99]['M'])), d['95th'].append(
                    len(plot_df.loc[plot_df['Sp'] > p95]['m/z'])), d['90th'].append(
                    len(plot_df.loc[plot_df['Sp'] > p90]['m/z']))
                d['80th'].append(len(plot_df.loc[plot_df['Sp'] > p80]['m/z'])), d['75th'].append(
                    len(plot_df.loc[plot_df['Sp'] > p75]['m/z'])), d['50th'].append(
                    len(plot_df.loc[plot_df['Sp'] > p50]['m/z']))
                d['0th'].append(len(plot_df.loc[plot_df['Sp'] > p0]['m/z']))

        d['95th(' + str(round(p95, 1)) + ')'] = d.pop('95th')
        d['90th(' + str(round(p90, 1)) + ')'] = d.pop('90th')
        d['80th(' + str(round(p80, 1)) + ')'] = d.pop('80th')
        d['75th(' + str(round(p75, 1)) + ')'] = d.pop('75th')
        d['50th(' + str(round(p50, 1)) + ')'] = d.pop('50th')
        d['0th(' + str(round(p0, 1)) + ')'] = d.pop('0th')
        d['99th(' + str(round(p99, 1)) + ')'] = d.pop('99th')
        percentile_df = pd.DataFrame(data=d)
        percentile_df.to_csv("./Sp_vs_percentileSpdecoys_{}.csv".format(name_output), index=False)

    label_lengths = 'all'
    # Apply the filter on sequence lengths for FDR estimation
    if lengths != 'all':
        plot_df = in_plot_df.loc[(in_plot_df['length'].isin(lengths)) & (in_plot_df['Sp'] > Sp_cutoff)]

        # Decoys that share the same mz_RT with targets above cutoff are always kept
        decoy_df = in_decoy_df[in_decoy_df['m/z_RT'].isin(plot_df['m/z_RT'])]
        label_lengths = "/".join([str(x) for x in lengths])

    # Create a dataframe of targets for FDR estimations, with only targets that have at least one decoy
    if targets_with_decoys == 'y':
        targets_FDR_df = plot_df[plot_df['m/z_RT'].isin(decoy_df['m/z_RT'])]
        targets_w_decoys_df = targets_FDR_df
    else:
        targets_FDR_df = plot_df
        targets_w_decoys_df = plot_df[plot_df['m/z_RT'].isin(decoy_df['m/z_RT'])]

    # Create a dataframe concatenating all the targets and decoys with rank = 1
    merged_df = pd.concat([targets_FDR_df.loc[targets_FDR_df['rank'] == 1], decoy_df.loc[decoy_df['rank'] == 1]],
                          ignore_index=True, sort=True)
    merged_df = merged_df.drop_duplicates(subset='m/z_RT', keep=False)

    # Creates a dataframe with only the top match targets
    toptargets_df = targets_FDR_df.loc[targets_FDR_df['rank'] == 1]

    decoy_max = np.float64(merged_df[merged_df['molecule_location'].str.contains('decoy') == True]['Sp'].max())

    if not np.isnan(decoy_max):
        step = round(decoy_max / 100, 3)

        d = {'Sp_cutoff': [], 'top_match_targets(t)': [], 'top_unique_targets': [], 'top_match_decoys(D)': [],
             'FDR(D/t)': [], 'FDR(%)': []}

        # Loop through various Sp_cutoffs to calculate FDR
        for i in np.arange(0, round(decoy_max, 3) + step, step):
            tot = len(toptargets_df[toptargets_df['Sp'] >= i]['m/z'])
            if tot > 0:
                d['Sp_cutoff'].append(i), d['top_match_targets(t)'].append(tot)
                top_decoy = len(
                    merged_df[(merged_df['molecule_location'].str.contains('decoy') == True) & (merged_df['Sp'] >= i)
                              & (merged_df['rank'] == 1)]['m/z'])
                d['top_match_decoys(D)'].append(top_decoy), d['FDR(D/t)'].append(round(top_decoy / tot, 4)), d[
                    'FDR(%)'].append(round(top_decoy * 100 / tot, 2))

                # Define unique the targets with unique sequences, disregarding everything else. Ignore
                # targets with one or more X in the sequence
                d['top_unique_targets'].append(toptargets_df[(toptargets_df['Sp'] > i)
                                                             & (toptargets_df['sequence'].str.contains('X') == False)][
                                                   'sequence'].nunique())

        out_df = pd.DataFrame(data=d)
        outcsv = out_df.to_csv(index=False)

    else:
        d = {'Sp_cutoff': [0], 'top_match_targets(t)': [len(merged_df[merged_df['Sp'] > 0]['m/z'])],
             'top_unique_targets': [toptargets_df['sequence'].nunique()], 'top_match_decoys(D)': [0],
             'FDR(D/t)': [0], 'FDR(%)': [0]}

        out_df = pd.DataFrame(data=d)
        outcsv = out_df.to_csv(index=False)

    # Outputs the csv file with the FDR estimations
    targets_w_decoys, targets_without_decoys = len(targets_w_decoys_df.index), len(targets_FDR_df.index) - len(
        targets_w_decoys_df.index)
    open("FDR_{}.csv".format(name_output), 'w').writelines(
        ["# targets with competing decoys used for FDR,{}\n # targets without competing decoys used for FDR,{}\n"
         "Sequence lengths targets/decoys,{}\n".format(targets_w_decoys, targets_without_decoys, label_lengths),
         outcsv])


def subplots_number(nts):
    """
    Determine the number and distribution of the subplots when representing the scatterplot Sp_targetvs Sp_decoys
    based on precursor nucleotide lengths
    """
    lengths = len(nts)
    if lengths <= 4:
        rows, columns = 2, 2
    elif 5 <= lengths <= 6:
        rows, columns = 2, 3
    elif 7 <= lengths <= 8:
        rows, columns = 2, 4
    elif 9 <= lengths <= 12:
        rows, columns = 3, 4
    elif 13 <= lengths <= 16:
        rows, columns = 4, 4
    else:
        rows, columns = 5, 5

    outvalues = []

    for r in range(0, rows):
        for c in range(0, columns):
            outvalues.append([r, c])

    return outvalues


def scatter_Sp_vs_Spdecoy(df, lengths, Sp_cutoff, info_box):
    """
    Generate a scatter plot considering only the best scoring candidate for each combination of precursor/RT
    """
    inplot_df, indecoy_df = df
    # Apply the Sp cutoff
    plot_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff]

    # The Sp cutoff is applied only to targets,
    # decoys that share the same mz_RT with targets above cutoff are always kept
    decoy_df = indecoy_df[indecoy_df['m/z_RT'].isin(plot_df['m/z_RT'])]

    plot_df = plot_df[plot_df['m/z_RT'].isin(decoy_df['m/z_RT'])]

    # Merge the targets and decoy dataframe based on common m/z_RT values to print Sp vs Sp
    df_merged = pd.merge(plot_df, decoy_df, how='inner', on='m/z_RT')
    df_merged.dropna(inplace=True, subset=['Sp_x', 'Sp_y'])

    if 'all' in lengths:
        x1 = df_merged['Sp_x']
        y1 = df_merged['Sp_y']
        label_lengths = 'all'

    elif 'analysis' in lengths:
        # Determine the max values for x and y axis in the graphs
        x1 = df_merged['Sp_x']
        y1 = df_merged['Sp_y']

        max_value = [y1.max(), x1.max()]

        nts = sorted(df_merged['length_x'].unique().tolist())

        # Create the scatter plot with 16 subplots divided by length
        fig, axs = plt.subplots(4, 4, sharex=True, sharey=True)
        fig.suptitle('Best targets vs corresponding decoys Sp scores at different precursor lengths', y=0.95,
                     fontsize=7, fontweight='bold')

        # Determine the number of subplots for the scatterplots of Sp_targets vs Sp_decoys on different precursor
        # Nucleotide lengths
        graphs = subplots_number(nts)

        for i, l in enumerate(nts):

            x = df_merged['Sp_x'].loc[df_merged['length_x'] == l]
            y = df_merged['Sp_y'].loc[df_merged['length_x'] == l]

            plt.ylim(0 - max(max_value) / 50, max(max_value) + max(max_value) / 25)
            plt.xlim(0 - max(max_value) / 50, max(max_value) + max(max_value) / 25)

            axs[graphs[i][0], graphs[i][1]].scatter(x, y, s=10, facecolors='dimgray', edgecolors='black', marker=".",
                                                    alpha=0.7)
            axs[graphs[i][0], graphs[i][1]].set_title('{}-mer[{}]'.format(l, x.count()), pad=-5, fontsize=6)
            for ax in axs.flat:
                ax.set(xlabel='Targets [Sp]', ylabel='Decoys [Sp]')

            # Hide x labels and tick labels for top plots and y ticks for right plots.
            for ax in axs.flat:
                ax.label_outer()

            # Plot the line x = y
            axs[graphs[i][0], graphs[i][1]].plot(list(range(100)), list(range(100)), color='black', linewidth=0.1,
                                                 linestyle="-.")



    else:
        x1 = df_merged.loc[df_merged['length_x'].isin(lengths)]['Sp_x']
        y1 = df_merged.loc[df_merged['length_y'].isin(lengths)]['Sp_y']
        label_lengths = ",".join([str(x) for x in lengths])

    if 'analysis' not in lengths:
        max_value = [y1.max(), x1.max()]
        plt.ylim(0 - max(max_value) / 50, max(max_value) + max(max_value) / 25)
        plt.xlim(0 - max(max_value) / 50, max(max_value) + max(max_value) / 25)

        plt.scatter(x1, y1, facecolors='dimgray', edgecolors='black', marker=".", alpha=0.7)

        plt.xlabel('Targets [Sp]'), plt.ylabel('Decoys [Sp]')
        plt.title('Best targets vs corresponding decoys Sp score')

        # Place a box with info on the graph about the total number of points and parameters
        if info_box == 'y':
            textstr = '\n'.join(('N = {}'.format(len(plot_df.index)),
                                 'Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                 'sequence lengths = {}'.format(label_lengths)))

            props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
            plt.text(0.03, 0.96, textstr, transform=plt.gca().transAxes, fontsize=4.5, verticalalignment='top',
                     bbox=props,
                     linespacing=1.75)

        # Plot the line x = y
        plt.plot(list(range(100)), list(range(100)), color='black', linewidth=0.2, linestyle="-.")

    else:
        pass


def scatter_dSp_vs_Sp(df, lengths, Sp_cutoff, info_box, targets_with_decoys):
    """
     Generate a scatter plot considering only thebest scoring targets & competing decoys
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff]

    # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoffs
    # are always kept
    decoy_df = indecoy_df[indecoy_df['m/z_RT'].isin(plot_df['m/z_RT'])]

    if targets_with_decoys == 'y':
        plot_df = plot_df[plot_df['m/z_RT'].isin(decoy_df['m/z_RT'])]

    if 'all' in lengths:
        x1 = plot_df['dSp']
        y1 = plot_df['Sp']
        x2 = decoy_df['dSp']
        y2 = decoy_df['Sp']
        label_lengths = 'all'

    else:
        x1 = plot_df.loc[plot_df['length'].isin(lengths)]['dSp']
        y1 = plot_df.loc[plot_df['length'].isin(lengths)]['Sp']
        x2 = decoy_df.loc[decoy_df['length'].isin(lengths)]['dSp']
        y2 = decoy_df.loc[decoy_df['length'].isin(lengths)]['Sp']
        label_lengths = ",".join([str(x) for x in lengths])

    # Determine the y axis min and max values
    max_val = [y1.max(), y2.max()]
    plt.ylim(0 - max(max_val) / 50, max(max_val) + max(max_val) / 25)

    plt.scatter(x1, y1, facecolors='royalblue', edgecolors='blue', marker=".", alpha=0.75)
    plt.scatter(x2, y2, facecolors='lightcoral', edgecolors='red', marker=".", alpha=0.75)

    # Place a box with info on the graph about the total number of points and parameters
    if info_box == 'y':
        textstr = ('\n'.join(('Sp cutoff [targets] = {}'.format(Sp_cutoff),
                              'Targets w/ decoys = {}'.format(targets_with_decoys),
                              'sequence lengths = {}'.format(label_lengths))))

        props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
        plt.text(0.85, 0.9, textstr, transform=plt.gca().transAxes, fontsize=4, verticalalignment='top', bbox=props,
                 linespacing=1.75)

    # Label axis and title of the plot
    plt.xlabel('dSp'), plt.ylabel('Sp')
    plt.title('Sp vs dSp scores for the best targets & competing decoys')

    # Create a legend
    dec_legend = plt.scatter([], [], facecolors='lightcoral', edgecolors='red', marker=".", alpha=0.75)
    tar_legend = plt.scatter([], [], facecolors='royalblue', edgecolors='blue', marker=".", alpha=0.75)
    plt.legend((tar_legend, dec_legend), ('targets ({})'.format(len(x1.index)), 'decoys ({})'.format(len(x2.index))))


def scatter_nts_vs_score(df, yl, y_bot, y_top, lengths, Sp_cutoff, info_box, targets_with_decoys):
    """
     Scatter plot of target Sp scores vs sequence length
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff]

    # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
    # are always kept
    decoy_df = indecoy_df[indecoy_df['m/z_RT'].isin(plot_df['m/z_RT'])]

    if targets_with_decoys == 'y':
        plot_df = plot_df[plot_df['m/z_RT'].isin(decoy_df['m/z_RT'])]

    if not plot_df.empty:
        if 'all' in lengths:
            x1 = plot_df['length']
            y1 = plot_df[yl]
            x2 = decoy_df['length']
            y2 = decoy_df[yl]
            label_lengths = 'all'

        else:
            x1 = plot_df.loc[plot_df['length'].isin(lengths)]['length']
            y1 = plot_df.loc[plot_df['length'].isin(lengths)][yl]
            x2 = decoy_df.loc[decoy_df['length'].isin(lengths)]['length']
            y2 = decoy_df.loc[decoy_df['length'].isin(lengths)][yl]
            label_lengths = ",".join([str(x) for x in lengths])

        plot_df = plot_df.astype({"length": int})
        # Set the minimum value of x to be 1
        plt.xticks(list(range(plot_df['length'].min(), plot_df['length'].max() + 1)))

        if y_bot and y_top != 0:
            plt.ylim(y_bot, y_top)

        else:
            if y_top != 0:
                plt.ylim(top=y_top)
            else:
                plt.ylim(y_bot - plot_df[yl].max() / 50, plot_df[yl].max() + plot_df[yl].max() / 25)

        x1 -= 0.1
        x2 += 0.1

        plt.scatter(x1, y1, facecolors='royalblue', edgecolors='blue', marker=".", alpha=0.75)
        plt.scatter(x2, y2, facecolors='lightcoral', edgecolors='red', marker=".", alpha=0.75)

        # Place a box with info on the graph about the total number of points and parameters
        if info_box == 'y':
            textstr = ('\n'.join(('Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  'sequence lengths = {}'.format(label_lengths),
                                  'Targets w/ decoys = {}'.format(targets_with_decoys))))

            props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
            plt.text(0.855, 0.79, textstr, transform=plt.gca().transAxes, fontsize=4, verticalalignment='top',
                     bbox=props, linespacing=1.75)

        plt.xlabel('sequence_length'), plt.ylabel(yl)
        plt.title('Sp score vs sequence length for best targets & competing decoys')

        # Create a legend
        dec_legend = plt.scatter([], [], facecolors='lightcoral', edgecolors='red', marker=".", alpha=0.75)
        tar_legend = plt.scatter([], [], facecolors='royalblue', edgecolors='blue', marker=".", alpha=0.75)
        plt.legend((tar_legend, dec_legend),
                   ('targets ({})'.format(len(x1.index)), 'decoys ({})'.format(len(x2.index))))


def box_nts_vs_score(df, y, y_bot, y_top, lengths, Sp_cutoff, targets_with_decoys):
    """
     Box plot length vs scores of targets and decoys
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff]

    # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
    # are always kept
    decoy_df = indecoy_df[indecoy_df['m/z_RT'].isin(plot_df['m/z_RT'])]

    if targets_with_decoys == 'y':
        plot_df = plot_df[plot_df['m/z_RT'].isin(decoy_df['m/z_RT'])]

    if not plot_df.empty:
        plot_df = plot_df.astype({"length": int})
        max_nts = plot_df['length'].max()

        data1, data2 = [], []

        if 'all' in lengths:
            for i in range(3, max_nts + 1):
                data1.append((list(plot_df.loc[plot_df['length'] == i, y])))
                data2.append((list(decoy_df.loc[decoy_df['length'] == i, y])))

            # Delete the entries of shorter length from the list, empty values cause problems with boxplots
            value = -1
            for i, x in enumerate(data1):
                if x:
                    break
                else:
                    value = i

            data1 = data1[(value + 1):]
            data2 = data2[(value + 1):]

        else:
            lengths.sort()
            for i in lengths:
                data1.append((list(plot_df.loc[plot_df['length'] == int(i), y])))
                data2.append((list(decoy_df.loc[decoy_df['length'] == int(i), y])))

        plt.xlabel('sequence_length'), plt.ylabel(y)

        if y_bot and y_top != 0:
            plt.ylim((y_bot, y_top))

        else:

            if y_top != 0:
                plt.ylim(top=y_top)

            else:
                plt.ylim((y_bot - plot_df['Sp'].max() / 50, plot_df['Sp'].max() + plot_df['Sp'].max() / 25))

        if 'all' in lengths:

            lengths_labels = list(np.arange(plot_df['length'].min(), plot_df['length'].max() + 1))

            target = plt.boxplot(data1, positions=[np.float64(x) - 0.15 for x in lengths_labels], manage_ticks=False,
                                 widths=0.25, sym='b.')
            decoys = plt.boxplot(data2, positions=[np.float64(x) + 0.15 for x in lengths_labels], manage_ticks=False,
                                 widths=0.25, sym='r.')

            for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
                plt.setp(target[element], color='blue')
                plt.setp(decoys[element], color='red')

            plt.xticks(lengths_labels, lengths_labels)

        else:
            lengths.sort()
            target = plt.boxplot(data1, positions=[np.float64(x) - 0.15 for x in lengths],
                                 manage_ticks=False, widths=0.25, sym='b.')
            decoys = plt.boxplot(data2, positions=[np.float64(x) + 0.15 for x in lengths],
                                 manage_ticks=False, widths=0.25, sym='r.')

            for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
                plt.setp(target[element], color='blue')
                plt.setp(decoys[element], color='red')

            plt.xticks(lengths, lengths)

        # Create a legend
        plt.plot([], c='blue', label='targets')
        plt.plot([], c='red', label='decoys')
        plt.legend()


def hist_Sp(df, lengths, Sp_cutoff, targets_with_decoys, info_box):
    """
     Histogram with all Sp values for targets and decoys
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff]

    # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
    # are always kept
    decoy_df = indecoy_df[indecoy_df['m/z_RT'].isin(plot_df['m/z_RT'])]

    if targets_with_decoys == 'y':
        plot_df = plot_df[plot_df['m/z_RT'].isin(decoy_df['m/z_RT'])]

    # Apply the option about lengths if specified
    if 'all' not in lengths:
        label_lengths = ",".join([str(x) for x in lengths])
        plot_df = plot_df.loc[plot_df['length'].isin(lengths)]
        decoy_df = decoy_df.loc[decoy_df['length'].isin(lengths)]
    else:
        label_lengths = 'all'

    # Determine statistical descriptors for targets and decoys
    median_t, mean_t, std_t = plot_df['Sp'].median(), plot_df['Sp'].mean(), plot_df['Sp'].std()
    median_d, mean_d, std_d = decoy_df['Sp'].median(), decoy_df['Sp'].mean(), decoy_df['Sp'].std()

    plt.ylabel('Frequency')
    plt.title('Best targets & competing decoys Sp scores')

    if not plot_df.empty:
        maxvalue = plot_df['Sp'].max()
        if decoy_df['Sp'].max() > plot_df['Sp'].max():
            maxvalue = decoy_df['Sp'].max()

        plt.xlim(left=0, right=maxvalue + 0.1)

        plt.hist(plot_df['Sp'], bins=50, range=(0, maxvalue + 0.1), alpha=0.5,
                 label='targets ({})'.format(len(plot_df.index)), edgecolor='black', color='blue')

        plt.hist(decoy_df['Sp'], bins=50, range=(0, maxvalue + 0.1), alpha=0.5,
                 label='decoys ({})'.format(len(decoy_df.index)), edgecolor='black', color='red')

        # Place a box with info on the graph about the total number of points and parameters
        if info_box == 'y':
            textstr = ('\n'.join(('Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  'Targets w/ decoys = {}'.format(targets_with_decoys),
                                  'sequence lengths = {}'.format(label_lengths),
                                  r'$\mu(t)=%.2f$ M(t)=%.2f $\sigma(t)=%.2f$' % (mean_t, median_t, std_t,),
                                  r'$\mu(d)=%.2f$ M(d)=%.2f $\sigma(d)=%.2f$' % (mean_d, median_d, std_d,))))
            props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
            plt.text(0.79, 0.78, textstr, transform=plt.gca().transAxes, fontsize=4, verticalalignment='top',
                     bbox=props, linespacing=1.75)

        plt.legend(loc='upper right')


def hist_top_Sp(df, lengths, Sp_cutoff, targets_with_decoys, info_box):
    """
     Histogram with Sp values (rank = 1) for top (rank = 1) targets and decoys
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff]

    # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
    # are always kept
    decoy_df = indecoy_df[indecoy_df['m/z_RT'].isin(plot_df['m/z_RT'])]

    if targets_with_decoys == 'y':
        plot_df = plot_df[plot_df['m/z_RT'].isin(decoy_df['m/z_RT'])]

    # Select only the the targets and decoys with rank = 1, top targets and top decoys
    top_plot_df, top_decoy_df = plot_df[plot_df['rank'] == 1], decoy_df[decoy_df['rank'] == 1]

    # Apply the option about lengths if specified
    if 'all' not in lengths:
        label_lengths = ",".join([str(x) for x in lengths])
        top_plot_df = top_plot_df.loc[top_plot_df['length'].isin(lengths)]
        top_decoy_df = top_decoy_df.loc[top_decoy_df['length'].isin(lengths)]
    else:
        label_lengths = 'all'

    # Determine statistical descriptors for targets and decoys
    median_t, mean_t, std_t = top_plot_df['Sp'].median(), top_plot_df['Sp'].mean(), top_plot_df['Sp'].std()
    median_d, mean_d, std_d = top_decoy_df['Sp'].median(), top_decoy_df['Sp'].mean(), top_decoy_df['Sp'].std()

    plt.xlabel('Sp'), plt.ylabel('Frequency')
    plt.title('Sp scores of top (rank = 1) targets & decoys')

    if not plot_df.empty:
        maxvalue = plot_df['Sp'].max()
        if decoy_df['Sp'].max() > plot_df['Sp'].max():
            maxvalue = decoy_df['Sp'].max()

        plt.xlim(left=0, right=maxvalue + 0.1)

        plt.hist(top_plot_df['Sp'], bins=50, range=(0, maxvalue + 0.1), alpha=0.5,
                 label='targets ({})'.format(len(top_plot_df.index)), edgecolor='black', color='blue')
        plt.hist(top_decoy_df['Sp'], bins=50, range=(0, maxvalue + 0.1), alpha=0.5,
                 label='decoys ({})'.format(len(top_decoy_df.index)), edgecolor='black', color='red')

        # Place a box with info on the graph about the total number of points and parameters
        if info_box == 'y':
            textstr = ('\n'.join(('Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  'Targets w/ decoys = {}'.format(targets_with_decoys),
                                  'sequence lengths = {}'.format(label_lengths),
                                  r'$\mu(t)=%.2f$ M(t)=%.2f $\sigma(t)=%.2f$' % (mean_t, median_t, std_t,),
                                  r'$\mu(d)=%.2f$ M(d)=%.2f $\sigma(d)=%.2f$' % (mean_d, median_d, std_d,))))
            props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
            plt.text(0.79, 0.78, textstr, transform=plt.gca().transAxes, fontsize=4, verticalalignment='top',
                     bbox=props,
                     linespacing=1.75)

        plt.legend(loc='upper right')


def hist_second_dSp(df, lengths, Sp_cutoff, targets_with_decoys, info_box):
    """
     Histogram with dSp2 values for top targets (rank = 1) with at least another competing target
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff]

    # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
    # are always kept
    decoy_df = indecoy_df[indecoy_df['m/z_RT'].isin(plot_df['m/z_RT'])]

    if targets_with_decoys == 'y':
        plot_df = plot_df[plot_df['m/z_RT'].isin(decoy_df['m/z_RT'])]

    # Select only the the targets with rank = 1
    sec_plot_df = plot_df[(plot_df['dSp2'] >= 0) & (plot_df['rank'] == 1)]

    # Apply the option about lengths if specified
    if 'all' not in lengths:
        label_lengths = ",".join([str(x) for x in lengths])
        sec_plot_df = sec_plot_df.loc[sec_plot_df['length'].isin(lengths)]
    else:
        label_lengths = 'all'

    # Determine statistical descriptors for targets and decoys
    median_t, mean_t, std_t = sec_plot_df['dSp2'].median(), sec_plot_df['dSp2'].mean(), sec_plot_df[
        'dSp2'].std()

    plt.xlabel('dSp'), plt.ylabel('Frequency')
    plt.title('dSp2 scores distribution of top targets (rank = 1) with at least another competing target')

    if not plot_df.empty:
        plt.hist(sec_plot_df['dSp2'], bins=20, alpha=0.5, label='targets ({})'.format(len(sec_plot_df.index)),
                 edgecolor='black', color='blue')

        # Place a box with info on the graph about the total number of points and parameters
        if info_box == 'y':
            textstr = ('\n'.join(('Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  'sequence lengths = {}'.format(label_lengths),
                                  r'$\mu(t)=%.2f$ M(t)=%.2f $\sigma(t)=%.2f$' % (mean_t, median_t, std_t,))))
            props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
            plt.text(0.79, 0.93, textstr, transform=plt.gca().transAxes, fontsize=4, verticalalignment='top',
                     bbox=props,
                     linespacing=1.75)

        plt.legend(loc='upper right')


def ppm_errors_histogram(df, Sp_cutoff, label, match_file, info_box):
    """
     Histogram with the distribution of ppm incertitude on MS1 and MS2 matches for top Sp values (rank = 1)
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff and select only the top targets (rank = 1)
    plot_df = inplot_df.loc[(inplot_df['Sp'] > Sp_cutoff) & (inplot_df['rank'] == 1)]

    # Obtain the info on MS1/MS2 ppm offsets used in matching&scoring
    with open(match_file) as infile:
        head = [next(infile) for x in range(20)]

    for element in head:
        if "#MS1_ppm" in element:
            MS1_ppm = float(element.split()[1])
        if "#MS2_ppm" in element:
            MS2_ppm = float(element.split()[1])
        if "#MS1_offset" in element:
            MS1_offset = float(element.split()[1])
        if "#MS2_offset" in element:
            MS2_offset = float(element.split()[1])

    # Create a unique list with all the MS2 offset ppm values
    MS2_ppm_list = []

    MS2_ppm_series = plot_df['MS2_ppm'].dropna()
    for x in MS2_ppm_series:
        for i in x.split(","):
            MS2_ppm_list.append(float(i))

    # Determine meand and std dev for MS1 and MS2 offset values
    median_MS1, mean_MS1, std_MS1 = plot_df['MS1_ppm'].median(), plot_df['MS1_ppm'].mean(), plot_df['MS1_ppm'].std()
    median_MS2, mean_MS2, std_MS2 = statistics.median(MS2_ppm_list), statistics.mean(MS2_ppm_list), statistics.stdev(
        MS2_ppm_list)

    # Add the mean and std dev values in the plot
    if label == 'MS1' and info_box == 'y':
        textstr = ('\n'.join(('Sp cutoff [targets] = {}'.format(Sp_cutoff),
                              'MS1 input offset [ppm] = {}'.format(MS1_offset),
                              r'$\mu=%.1f$ M=%.1f $\sigma=%.1f$' % (mean_MS1, median_MS1, std_MS1,),)))

        props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
        plt.text(0.03, 0.93, textstr, transform=plt.gca().transAxes, fontsize=4, verticalalignment='top', bbox=props,
                 linespacing=1.75)
        plt.ylabel('Frequency')

    elif label == 'MS2' and info_box == 'y':
        textstr = ('\n'.join(('Sp cutoff [targets] = {}'.format(Sp_cutoff),
                              'MS2 input offset [ppm] = {}'.format(MS2_offset),
                              r'$\mu=%.1f$ M=%.1f $\sigma=%.1f$' % (mean_MS2, median_MS2, std_MS2,),)))
        props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
        plt.text(0.03, 0.93, textstr, transform=plt.gca().transAxes, fontsize=4, verticalalignment='top', bbox=props,
                 linespacing=1.75)
        plt.xlabel('ppm_offset'), plt.ylabel('Frequency')

    plt.title('m/z offsets (ppm) for the top {} targets (centered around the input offset)'.format(label))

    if label == 'MS1':
        plt.hist(plot_df['MS1_ppm'], bins=35, label='targets ({})'.format(len(plot_df.index)),
                 range=[-MS1_ppm, MS1_ppm], edgecolor='black', color='blue')
    else:
        plt.hist(MS2_ppm_list, bins=50, label='MS2 ions ({})'.format(len(MS2_ppm_list)),
                 range=[-MS2_ppm, MS2_ppm], edgecolor='black', color='red')

    plt.legend(loc='upper right')


def scatter_nts_z_vs_score(df, yl, y_bot, y_top, lengths, Sp_cutoff, targets_with_decoys, info_box):
    """
     Generate a scatter plot considering only the top scoring candidate for each combination of precursor/RT
     """
    # Keep only the top scoring entry for each combination of precursor + RT

    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df, decoy_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff], indecoy_df.loc[indecoy_df['Sp'] > Sp_cutoff]

    # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
    # are always kept
    decoy_df = indecoy_df[indecoy_df['m/z_RT'].isin(plot_df['m/z_RT'])]

    if targets_with_decoys == 'y':
        plot_df = plot_df[plot_df['m/z_RT'].isin(decoy_df['m/z_RT'])]

    plot_df = plot_df.sort_values(['length', 'charge'], ascending=[True, True])

    if not plot_df.empty:
        if 'all' in lengths:
            x1 = plot_df['length_z']
            y1 = plot_df[yl]
            x2 = decoy_df['length_z']
            y2 = decoy_df[yl]
            label_lengths = 'all'

        else:
            x1 = plot_df.loc[plot_df['length'].isin(lengths)]['length_z']
            y1 = plot_df.loc[plot_df['length'].isin(lengths)][yl]
            x2 = decoy_df.loc[decoy_df['length'].isin(lengths)]['length_z']
            y2 = decoy_df.loc[decoy_df['length'].isin(lengths)][yl]
            label_lengths = ",".join([str(x) for x in lengths])

        if y_bot and y_top != 0:
            plt.ylim(y_bot, y_top)

        else:

            if y_top != 0:
                plt.ylim(top=y_top)

            else:
                plt.ylim(y_bot - plot_df[yl].max() / 50, plot_df[yl].max() + plot_df[yl].max() / 25)

        plt.scatter(x1, y1, facecolors='royalblue', edgecolors='blue', marker=".", alpha=0.75)
        plt.scatter(x2, y2, facecolors='lightcoral', edgecolors='red', marker=".", alpha=0.75)

        plt.margins(x=0.025)

        # Place a box with info on the graph about the total number of points and parameters
        if info_box == 'y':
            textstr = ('\n'.join(('Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  'sequence lengths = {}'.format(label_lengths),
                                  'Targets w/ decoys = {}'.format(targets_with_decoys))))

            props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
            plt.text(0.855, 0.79, textstr, transform=plt.gca().transAxes, fontsize=4, verticalalignment='top',
                     bbox=props, linespacing=1.75)

        plt.xlabel('sequence_length_charge'), plt.ylabel(yl)
        plt.title('Sp score vs sequence length & charge for best targets & competing decoys')


        # Create a legend
        dec_legend = plt.scatter([], [], facecolors='lightcoral', edgecolors='red', marker=".", alpha=0.75)
        tar_legend = plt.scatter([], [], facecolors='royalblue', edgecolors='blue', marker=".", alpha=0.75)
        plt.legend((tar_legend, dec_legend),
                   ('targets ({})'.format(len(x1.index)), 'decoys ({})'.format(len(x2.index))))


def box_nts_z_vs_score(df, y, y_bot, y_top, lengths, Sp_cutoff, targets_with_decoys):
    """
     Box plot with Sp scores vs target/decoy sequences length grouped by charge state
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df, decoy_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff], indecoy_df.loc[indecoy_df['Sp'] > Sp_cutoff]

    # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
    # are always kept
    decoy_df = indecoy_df[indecoy_df['m/z_RT'].isin(plot_df['m/z_RT'])]

    if targets_with_decoys == 'y':
        plot_df = plot_df[plot_df['m/z_RT'].isin(decoy_df['m/z_RT'])]

    plot_df = plot_df.sort_values(['length', 'charge'], ascending=[True, True])

    if not plot_df.empty:
        if 'all' not in lengths:
            plot_df = plot_df[plot_df['length'].isin(lengths)]

        nts_z = plot_df['length_z'].unique()

        if y_bot and y_top != 0:
            plt.ylim((y_bot, y_top))

        else:

            if y_top != 0:
                plt.ylim(top=y_top)

            else:
                plt.ylim((y_bot - plot_df['Sp'].max() / 50, plot_df['Sp'].max() + plot_df['Sp'].max() / 25))

        targets, decoys = [], []
        for i in nts_z:
            targets.append((list(plot_df.loc[plot_df['length_z'] == i, y])))
            decoys.append((list(decoy_df.loc[decoy_df['length_z'] == i, y])))

        plt.xlabel('sequence_length_charge'), plt.ylabel(y)
        target = plt.boxplot(targets, labels=nts_z.tolist(), sym='b.')
        decoy = plt.boxplot(decoys, labels=nts_z.tolist(), sym='b.')

        for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
            plt.setp(target[element], color='blue')
            plt.setp(decoy[element], color='red')

        # Create a legend
        plt.plot([], c='blue', label='targets')
        plt.plot([], c='red', label='decoys')
        plt.legend()


def FDR_update(df, lengths, Sp_cutoff, name_output):
    """
     Update the FDR table if target without decoys are selected for statistics
     """
    in_plot_df, in_decoy_df = df

    if 'all' in lengths:
        plot_df = in_plot_df.loc[in_plot_df['Sp'] > Sp_cutoff]

        # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
        # are always kept
        decoy_df = in_decoy_df[in_decoy_df['m/z_RT'].isin(plot_df['m/z_RT'])]
        label_lengths = 'all'

    else:
        plot_df = in_plot_df.loc[(in_plot_df['length'].isin(lengths)) & (in_plot_df['Sp'] > Sp_cutoff)]

        # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
        # are always kept
        decoy_df = in_decoy_df[in_decoy_df['m/z_RT'].isin(plot_df['m/z_RT'])]
        label_lengths = "/".join([str(x) for x in lengths])

    # Creates a dataframe of targets for FDR estimations, with targets with and without decoys
    targets_FDR_df = plot_df
    targets_w_decoys_df = plot_df[plot_df['m/z_RT'].isin(decoy_df['m/z_RT'])]

    # Create a dataframe with all values of rank = 1
    merged_df = pd.concat([targets_FDR_df.loc[targets_FDR_df['rank'] == 1], decoy_df.loc[decoy_df['rank'] == 1]],
                          ignore_index=True, sort=True)
    merged_df = merged_df.drop_duplicates(subset='m/z_RT', keep=False)

    # Create a dataframe with only the top match targets
    toptargets_df = targets_FDR_df.loc[targets_FDR_df['rank'] == 1]

    decoy_max = np.float64(merged_df[merged_df['molecule_location'].str.contains('decoy') == True]['Sp'].max())

    if not np.isnan(decoy_max):
        step = round(decoy_max / 100, 3)

        d = {'Sp_cutoff': [], 'top_match_targets(t)': [], 'top_unique_targets': [], 'top_match_decoys(D)': [],
             'FDR(D/t)': [], 'FDR(%)': []}

        # Cycle through various Sp_cutoffs to calculate FDR
        for i in np.arange(0, round(decoy_max, 3) + step, step):
            tot = len(toptargets_df.loc[toptargets_df['Sp'] >= i]['m/z'])
            if tot > 0:
                d['Sp_cutoff'].append(i), d['top_match_targets(t)'].append(tot)
                top_decoy = len(
                    merged_df.loc[(merged_df['molecule_location'].str.contains('decoy') == True) &
                                  (merged_df['Sp'] >= i) & (merged_df['rank'] == 1)]['m/z'])
                d['top_match_decoys(D)'].append(top_decoy), d['FDR(D/t)'].append(round(top_decoy / tot, 4)), d[
                    'FDR(%)'].append(round(top_decoy * 100 / tot, 2))

                # Define unique the targets with unique sequences, disregarding everything else. Ignore
                # targets with one or more X in the sequence
                d['top_unique_targets'].append(toptargets_df[(toptargets_df['Sp'] > i)
                                                             & (toptargets_df['sequence'].str.contains('X') == False)][
                                                   'sequence'].nunique())

        out_df = pd.DataFrame(data=d)
        outcsv = out_df.to_csv(index=False)

    else:
        d = {'Sp_cutoff': [0], 'top_match_targets(t)': [len(merged_df[merged_df['Sp'] > 0]['m/z'])],
             'top_unique_targets': [toptargets_df['sequence'].nunique()], 'top_match_decoys(D)': [0],
             'FDR(D/t)': [0], 'FDR(%)': [0]}

        out_df = pd.DataFrame(data=d)
        outcsv = out_df.to_csv(index=False)

    # Output the csv file with the FDR estimations
    targets_with_decoys, targets_without_decoys = len(targets_w_decoys_df.index), len(targets_FDR_df.index) - len(
        targets_w_decoys_df.index)
    open("FDR_{}_edited.csv".format(name_output), 'w').writelines([
        "Targets with competing decoys used for FDR,{}\nTargets without competing decoys used for FDR,{}\n"
         "Sequence lengths targets/decoys,{}\n".format(
            targets_with_decoys, targets_without_decoys,
            label_lengths), outcsv])
