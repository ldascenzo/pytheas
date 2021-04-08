#!/usr/bin/python3

"""
Last update: December 2020
Credits to: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzo@scripps.edu
"""

import sys, os, re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import statistics

np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

def read_excel_input(nts_file):
    """
     Produces a dataframe with all the info on the nucleobases from the input file nts_alphabet_light
     """
    # Checking that the nts_alphabet_light file given in argument exists
    if not os.path.exists(nts_file):
        print("ERROR! File " + nts_file + " does not exist. Execution terminated without generating any output")
        sys.exit(1)

    # Creates a dataframe with info from Excel spreadsheet
    df = pd.read_excel(nts_file, sheet_name=[0, 1], header=0)

    # Drops rows with NaN values
    # df = df[pd.notnull(df['length'])]

    # Drops rows with additional legend indications
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
        d = {'M': [], 'RT': [], 'M_RT': [], 'rank': [], 'Theo_mono_M': [], 'isotopologue': [], 'nts': [], 'Sp': [],
             'dSp': [], 'n': [], 'sumI': [], 'L': [], 'n/L': [], 'Beta': [], 'charge': [], 'isotope': [],
             'molecule': [], 'sequence': [], 'sequence_ext': [],
             'MS1_ppm': [], 'MS2_ppm': [], 'MS2_matches': []}

    else:
        d = {'M': [], 'RT': [], 'M_RT': [], 'rank': [], 'Theo_mono_M': [], 'isotopologue': [], 'nts': [], 'Sp': [],
             'dSp': [], 'n': [], 'sumI': [], 'sumI_all': [], 'L': [], 'n/L': [], 'Beta': [], 'charge': [],
             'isotope': [], 'molecule': [], 'sequence': [],
             'sequence_ext': [], 'MS1_ppm': [], 'MS2_ppm': [], 'MS2_matches': []}

    with open("./" + infile, 'rU') as input_file:
        for line in input_file:

            if line[0].isdigit():
                sp = line.split()
                prec_M = re.findall("[\d]+[.]?[\d]+", sp[0])[0]
                d['M'].append(prec_M), d['nts'].append(int(sp[9])), d['Sp'].append(np.float64(sp[4].split('=')[1])),
                d['dSp'].append(sp[5].split('=')[1]), d['charge'].append(sp[10]),
                d['isotope'].append(sp[8]), d['molecule'].append(sp[15]), d['sequence'].append(sp[11]),
                d['sequence_ext'].append(sp[12]), d['RT'].append(sp[1].split('=')[1]),
                d['M_RT'].append(prec_M + '_' + sp[1].split('=')[1]), d['MS1_ppm'].append(sp[3][:-3]),
                d['MS2_matches'].append(",".join(sp[16:-1])),
                d['Theo_mono_M'].append(round(np.float64(sp[2].split('=')[1]), 5))
                d['rank'].append(int(sp[6].split('=')[1]))

                if '*' in sp[0]:
                    d['isotopologue'].append('y')
                else:
                    d['isotopologue'].append('')

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
    # Also clusters the ions based on nts + charge values
    if excel != 'n':
        df_excel = read_excel_input(excel)
        for index, row in df.iterrows():

            max_int = \
                df_excel[0][
                    (df_excel[0]['m/z obs'].str.contains(row['M'])) & (df_excel[0]['RT'].str.contains(row['RT']))][
                    'Prec_max_int']
            for m in max_int:
                df.loc[df.index[index], 'Prec_max_int'] = m

            max_int = \
                df_excel[1][
                    (df_excel[1]['m/z obs'].str.contains(row['M'])) & (df_excel[1]['RT'].str.contains(row['RT']))][
                    'Prec_max_int']
            for m in max_int:
                df.loc[df.index[index], 'Prec_max_int'] = m

    z = [str(row.nts) + "_" + str(row.charge) for index, row in df.iterrows()]
    df['nts_z'] = z

    df = df.astype({"M": np.float64, 'dSp': np.float64})

    return df


def filter_data(df, top, yl):
    """
    Filter the data in preparation for the final graphic output
    """
    if top == 'y':
        filt_df = df.sort_values(['M', yl], ascending=[True, False]).reset_index(drop=True)

        # Terminates the script if no matches are found
        if filt_df.empty:
            print("ERROR!!!!!! No target matches found. Execution terminated without output")
            sys.exit(1)

        # Extracts the lines with decoys
        decoy_df = filt_df[filt_df['molecule'].str.contains('decoy')]

        # Keep only the highest scoring decoy per unique combination M_RT       
        decoy_df = decoy_df.sort_values('Sp', ascending=False).drop_duplicates(subset=['M', 'RT'])
        decoy_df = decoy_df.sort_values(['M', yl], ascending=[True, False]).reset_index(drop=True)

        # Removes the decoys from the targets filtered dataframe
        filt_df = filt_df[filt_df['molecule'].str.contains('decoy') == False]

        # Creates a dataframe with all the "nontop" scoring targets
        dup_df = filt_df.sort_values(['M', 'Sp'], ascending=[True, False])
        dup_df = dup_df.drop_duplicates(subset=['M_RT', 'sequence'])
        dup_df = dup_df[dup_df.duplicated(subset=['M', 'RT'])]
        dup_df = dup_df.sort_values('Sp', ascending=False).drop_duplicates(subset=['M', 'RT'])
        dup_df['dSp_2nd_target'] = dup_df['dSp']
        dup_df = dup_df[['M_RT', 'dSp_2nd_target']]

        # Keep only the highest scoring target per unique combination M_RT
        filt_df = filt_df.sort_values('Sp', ascending=False).drop_duplicates(subset=['M', 'RT'])
        filt_df = filt_df.sort_values(['M', yl], ascending=[True, False]).reset_index(drop=True)

        # Adds a column to the targets dataframe with the values of the second best dSp
        filt_df = pd.merge(filt_df, dup_df, on='M_RT', how='outer')
        filt_df = filt_df[
            ['M', 'RT', 'M_RT' , 'rank' , 'Theo_mono_M', 'isotopologue', 'nts', 'Sp', 'dSp', 'dSp_2nd_target',
             'MS1_ppm', 'n', 'sumI', 'sumI_all', 'L', 'n/L', 'Beta', 'charge', 'isotope', 'molecule', 'sequence',
             'sequence_ext', 'nts_z', 'MS2_ppm', 'MS2_matches']]

    else:
        filt_df = df

        # Terminates the script if no matches are found
        if filt_df.empty:
            print("ERROR!!!!! No target matches found. Execution terminated without output")
            sys.exit(1)

        decoy_df = filt_df[filt_df['molecule'].str.contains('decoy')]

    # Add a column with info on targets having at least one competing decoy
    filt_df.insert(9, 'has_decoy', '')
    filt_df.loc[filt_df.M_RT.isin(decoy_df.M_RT), 'has_decoy'] = 'y'

    return filt_df, decoy_df


def csv_output(df, l, per, Sp_cutoff, name_output, targets_with_decoys, isotopes):
    """
     Generates the output .csv file with the values used to create the graphs
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

    if 'all' in l:
        plot_df = in_plot_df.loc[in_plot_df['Sp'] > Sp_cutoff]

        # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
        # are always kept
        decoy_df = in_decoy_df[in_decoy_df.M_RT.isin(plot_df.M_RT)]

    else:
        plot_df = in_plot_df.loc[(in_plot_df['nts'].isin(l)) & (in_plot_df['Sp'] > Sp_cutoff)]

        # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
        # are always kept
        decoy_df = in_decoy_df[in_decoy_df.M_RT.isin(plot_df.M_RT)]

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

        nts, d = ['all'] + sorted(plot_df['nts'].unique()), {'nts': [], '99th': [], '95th': [], '90th': [], '80th': [],
                                                             '75th': [], '50th': [], '0th': []}
        for i in nts:
            d['nts'].append(i)
            if i != 'all':
                d['99th'].append(len(plot_df.loc[(plot_df['Sp'] > p99) & (plot_df['nts'] == i)]['M'])), d[
                    '95th'].append(len(plot_df.loc[(plot_df['Sp'] > p95) & (plot_df['nts'] == i)]['M']))
                d['90th'].append(len(plot_df.loc[(plot_df['Sp'] > p90) & (plot_df['nts'] == i)]['M'])), d[
                    '80th'].append(len(plot_df.loc[(plot_df['Sp'] > p80) & (plot_df['nts'] == i)]['M']))
                d['75th'].append(len(plot_df.loc[(plot_df['Sp'] > p75) & (plot_df['nts'] == i)]['M'])), d[
                    '50th'].append(len(plot_df.loc[(plot_df['Sp'] > p50) & (plot_df['nts'] == i)]['M']))
                d['0th'].append(len(plot_df.loc[(plot_df['Sp'] > p0) & (plot_df['nts'] == i)]['M']))

            else:
                d['99th'].append(len(plot_df.loc[plot_df['Sp'] > p99]['M'])), d['95th'].append(
                    len(plot_df.loc[plot_df['Sp'] > p95]['M'])), d['90th'].append(
                    len(plot_df.loc[plot_df['Sp'] > p90]['M']))
                d['80th'].append(len(plot_df.loc[plot_df['Sp'] > p80]['M'])), d['75th'].append(
                    len(plot_df.loc[plot_df['Sp'] > p75]['M'])), d['50th'].append(
                    len(plot_df.loc[plot_df['Sp'] > p50]['M']))
                d['0th'].append(len(plot_df.loc[plot_df['Sp'] > p0]['M']))

        d['95th(' + str(round(p95, 1)) + ')'] = d.pop('95th')
        d['90th(' + str(round(p90, 1)) + ')'] = d.pop('90th')
        d['80th(' + str(round(p80, 1)) + ')'] = d.pop('80th')
        d['75th(' + str(round(p75, 1)) + ')'] = d.pop('75th')
        d['50th(' + str(round(p50, 1)) + ')'] = d.pop('50th')
        d['0th(' + str(round(p0, 1)) + ')'] = d.pop('0th')
        d['99th(' + str(round(p99, 1)) + ')'] = d.pop('99th')
        percentile_df = pd.DataFrame(data=d)
        percentile_df.to_csv("./Sp_vs_percentileSpdecoys_{}.csv".format(name_output), index=False)

    # Creates a dataframe of targets for FDR estimations, with only targets that have at least one decoy
    if targets_with_decoys == 'y':
        targets_FDR_df = plot_df[plot_df.M_RT.isin(decoy_df.M_RT)]
        targets_w_decoys_df = targets_FDR_df
    else:
        targets_FDR_df = plot_df
        targets_w_decoys_df = plot_df[plot_df.M_RT.isin(decoy_df.M_RT)]

    # Creates a dataframe with all values of dSp = 0
    merged_df = pd.concat([targets_FDR_df.loc[targets_FDR_df['dSp'] == 0], decoy_df.loc[decoy_df['dSp'] == 0]],
                          ignore_index=True, sort=True)
    merged_df = merged_df.drop_duplicates(subset='M_RT', keep=False)
    toptargets_df = merged_df[merged_df['molecule'].str.contains('decoy') == False]

    # Creates a dataframe with only the top match targets
    toptargets_df = targets_FDR_df.loc[targets_FDR_df['rank'] == 1]

    decoy_max = np.float64(merged_df[merged_df['molecule'].str.contains('decoy') == True]['Sp'].max())

    if not np.isnan(decoy_max):
        step = round(decoy_max / 100, 3)

        d = {'Sp_cutoff': [], 'top_match_targets(t)': [], 'unique_sequences_targets': [], 'top_match_decoys(D)': [],
             'FDR(D/t)': [], 'FDR(%)': []}

        # Cycle through various Sp_cutoffs to calculate FDR
        for i in np.arange(0, round(decoy_max, 3) + step, step):
            tot = len(toptargets_df[toptargets_df['Sp'] > i]['M'])
            if tot > 0:
                d['Sp_cutoff'].append(i), d['top_match_targets(t)'].append(tot)
                top_decoy = len(
                    merged_df[(merged_df['molecule'].str.contains('decoy') == True) & (merged_df['Sp'] > i)
                    & (merged_df['rank'] == 1)]['M'])
                d['top_match_decoys(D)'].append(top_decoy), d['FDR(D/t)'].append(round(top_decoy / tot, 4)), d[
                    'FDR(%)'].append(round(top_decoy * 100 / tot, 2))
                d['unique_sequences_targets'].append(toptargets_df[toptargets_df['Sp'] > i]['sequence'].nunique())

        out_df = pd.DataFrame(data=d)
        outcsv = out_df.to_csv(index=False)

    else:
        d = {'Sp_cutoff': [0], 'top_match_targets(t)': [len(merged_df[merged_df['Sp'] > 0]['M'])],
             'unique_sequences_targets': [toptargets_df['sequence'].nunique()], 'top_match_decoys(D)': [0],
             'FDR(D/t)': [0], 'FDR(%)': [0]}

        out_df = pd.DataFrame(data=d)
        outcsv = out_df.to_csv(index=False)

    # Outputs the csv file with the FDR estimations
    targets_w_decoys, targets_without_decoys = len(targets_w_decoys_df.index), len(targets_FDR_df.index) - len(
        targets_w_decoys_df.index)
    open("FDR_{}.csv".format(name_output), 'w').writelines(
        ["Targets_with_decoys,{}\nTargets_without_decoys,{}\n\n".format(targets_w_decoys, targets_without_decoys),
         outcsv])


def subplots_number(nts):
    """
    Determines the number and distribution of the subplots when representing the scatterplot Sp_targetvsSp_decoys
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


def scatter_Sp_vs_Spdecoy(df, lengths, Sp_cutoff, MS1_offset_min, MS1_offset_max, info_box):
    """
     Generate a scatter plot considering only the top scoring candidate for each combination of precursor/RT
     """
    inplot_df, indecoy_df = df
    # Apply the Sp cutoff
    plot_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff]

    # Apply the offset cutoff if requested
    if isinstance(MS1_offset_min, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] >= MS1_offset_min]
    if isinstance(MS1_offset_max, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] <= MS1_offset_max]

        # The Sp cutoff is applied only to targets,
        # decoys that share the same mz_RT with targets above cutoff are always kept
    decoy_df = indecoy_df[indecoy_df.M_RT.isin(plot_df.M_RT)]

    plot_df = plot_df[plot_df.M_RT.isin(decoy_df.M_RT)]

    # Merge the targets and decoy dataframe based on common M_RT values to print Sp vs Sp
    df_merged = pd.merge(plot_df, decoy_df, how='inner', on='M_RT')
    df_merged.dropna(inplace=True, subset=['Sp_x', 'Sp_y'])

    if 'all' in lengths:
        x1 = df_merged['Sp_x']
        y1 = df_merged['Sp_y']

    elif 'analysis' in lengths:
        # Determine the max values for x and y axis in the graphs
        x1 = df_merged['Sp_x']
        y1 = df_merged['Sp_y']

        max_value = [y1.max(), x1.max()]

        nts = sorted(df_merged['nts_x'].unique().tolist())

        # Create the scatter plot with 16 subplots divided by nts length
        fig, axs = plt.subplots(4, 4, sharex=True, sharey=True)
        fig.suptitle('Targets score (Sp) vs corresponding decoys score (Sp) at different precursor lengths', y=0.95,
                     fontsize=7, fontweight='bold')

        # Determines the number of subplots for the scatterplots of Sp_targets vs Sp_decoys on different precursor
        # Nucleotide lengths
        graphs = subplots_number(nts)

        for i, l in enumerate(nts):

            x = df_merged['Sp_x'].loc[df_merged['nts_x'] == l]
            y = df_merged['Sp_y'].loc[df_merged['nts_x'] == l]

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
        x1 = df_merged.loc[df_merged['nts_x'].isin(lengths)]['Sp_x']
        y1 = df_merged.loc[df_merged['nts_y'].isin(lengths)]['Sp_y']

    if 'analysis' not in lengths:
        max_value = [y1.max(), x1.max()]
        plt.ylim(0 - max(max_value) / 50, max(max_value) + max(max_value) / 25)
        plt.xlim(0 - max(max_value) / 50, max(max_value) + max(max_value) / 25)

        plt.scatter(x1, y1, facecolors='dimgray', edgecolors='black', marker=".", alpha=0.7)

        plt.xlabel('Targets [Sp]'), plt.ylabel('Decoys [Sp]')
        plt.title('Targets score (Sp) vs corresponding decoys score (Sp)')

        # Place a box with info on the graph about the total number of points and parameters
        if info_box == 'y':
            if isinstance(MS1_offset_min, int) or isinstance(MS1_offset_max, int):
                textstr = '\n'.join(('N = {}'.format(len(plot_df.index)),
                                     'Nts. lengths = {}'.format(','.join([str(x) for x in lengths])),
                                     'MS1_offset = [{},{}]'.format(MS1_offset_min, MS1_offset_max),
                                     'Sp cutoff [targets] = {}'.format(Sp_cutoff)))
            else:
                textstr = '\n'.join(('N = {}'.format(len(plot_df.index)),
                                     'Nts. lengths = {}'.format(','.join([str(x) for x in lengths])),
                                     'Sp cutoff [targets] = {}'.format(Sp_cutoff)))

            props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
            plt.text(0.03, 0.96, textstr, transform=plt.gca().transAxes, fontsize=4.5, verticalalignment='top',
                     bbox=props,
                     linespacing=1.75)

        # Plot the line x = y
        plt.plot(list(range(100)), list(range(100)), color='black', linewidth=0.2, linestyle="-.")

    else:
        pass


def barchart_global_normalized_scores(df, lengths, Sp_cutoff, MS1_offset_min, MS1_offset_max, info_box):
    """
     Generate a barchart with all the normalized Sp scores for the highest scoring target/decoy per each unique
     sequence
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff]

    # Apply the offset cutoff if requested
    if isinstance(MS1_offset_min, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] >= MS1_offset_min]
    if isinstance(MS1_offset_max, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] <= MS1_offset_max]

    # The Sp cutoff is applied only to targets,
    # decoys that share the same mz_RT with targets above cutoff are always kept
    decoy_df = indecoy_df[indecoy_df.M_RT.isin(plot_df.M_RT)]

    plot_df = plot_df[plot_df.M_RT.isin(decoy_df.M_RT)]

    # Merge the targets and decoy dataframe based on common M_RT values to print Sp vs Sp
    df_merged = pd.merge(plot_df, decoy_df, how='inner', on='M_RT')
    df_merged.dropna(inplace=True, subset=['Sp_x', 'Sp_y'])

    # Create a smaller dataframe with only the top target/decoy per precursor ion
    # Multiple occurrences of matches with the same M but different RT are consolidated and only the highest
    # scoring target/decoy duo is used
    df_noduplicates_MRT = df_merged.drop_duplicates(subset='M_x')

    df_onlyvalues = df_noduplicates_MRT[['Sp_x', 'Sp_y']]

    df_onlyvalues = df_onlyvalues.astype(np.float64)
    df_onlyvalues = df_onlyvalues[(df_onlyvalues['Sp_x'] > 0.15) & (df_onlyvalues['Sp_y'] > 0.15)]

    df_onlyvalues['Sp_max'] = df_onlyvalues.max(axis=1)

    # Normalize the Sp scores for targets and decoys based on the highest scoring of the couple which will have a value
    # of 1
    df_normalized = df_onlyvalues[['Sp_x', 'Sp_y']].div(df_onlyvalues['Sp_max'], axis=0)
    df_normalized.reset_index(drop=True, inplace=True)

    labels = df_normalized.index.tolist()
    targets_sp = df_normalized['Sp_x']
    decoys_sp = df_normalized['Sp_y']

    x = np.arange(len(labels))  # the label locations
    width = 0.3  # the width of the bars

    # Determines the dimensions of the graph based on how many pairs have to be shown
    width_graph, height_graph = len(labels) / 5, len(labels) / 30

    fig, ax = plt.subplots(figsize=(width_graph, height_graph))
    ax.bar(x - width / 2, targets_sp, width, alpha=0.5,
           label='targets ({})'.format(len(targets_sp)), edgecolor='black', color='blue')
    ax.bar(x + width / 2, decoys_sp, width, alpha=0.5,
           label='decoys ({})'.format(len(targets_sp)), edgecolor='black', color='red')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Normalized Sp', fontsize='large')
    ax.set_title('Normalized Sp scores for highest scoring targets/decoys', fontsize='xx-large')
    ax.set_xticks(x)
    ax.set_xticklabels([x + 1 for x in labels])

    plt.legend(loc='upper right', fontsize='large')


def scatter_dSp_vs_Sp(df, lengths, per, fdr_lines, Sp_cutoff, name_output, fdr_df_input, info_box, targets_with_decoys,
                      MS1_offset_min, MS1_offset_max):
    """
     Generate a scatter plot considering only the top scoring candidate for each combination of precursor/RT
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff]

    # Apply the offset cutoff if requested
    if isinstance(MS1_offset_min, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] >= MS1_offset_min]
    if isinstance(MS1_offset_max, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] <= MS1_offset_max]

    # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoffs
    # are always kept
    decoy_df = indecoy_df[indecoy_df.M_RT.isin(plot_df.M_RT)]

    if targets_with_decoys == 'y':
        plot_df = plot_df[plot_df.M_RT.isin(decoy_df.M_RT)]

    if 'all' in lengths:
        x1 = plot_df['dSp']
        y1 = plot_df['Sp']
        x2 = decoy_df['dSp']
        y2 = decoy_df['Sp']

        if per == 'y':
            # Add horizontal lines with info on percentiles populations

            frac_value_95 = round(
                plot_df[plot_df > (decoy_df.quantile(0.95)['Sp'])].count()['Sp'] / plot_df['Sp'].count(), 2)
            frac_value_90 = round(
                plot_df[plot_df > (decoy_df.quantile(0.9)['Sp'])].count()['Sp'] / plot_df['Sp'].count(), 2)

            plt.axhline(y=decoy_df.quantile(0.95)['Sp'], linestyle='--', linewidth=1, color='blue')
            plt.text(x2.max() - x2.max() / 7, decoy_df.quantile(0.95)['Sp'] + decoy_df.quantile(0.95)['Sp'] / 100,
                     "95th p = " + str(frac_value_95), fontsize=5, color='blue')

            plt.axhline(y=decoy_df.quantile(0.9)['Sp'], linestyle='--', linewidth=1, color='cyan')
            plt.text(x2.max() - x2.max() / 7, decoy_df.quantile(0.9)['Sp'] + decoy_df.quantile(0.95)['Sp'] / 100,
                     "90th p = " + str(frac_value_90), fontsize=5, color='cyan')

    else:
        x1 = plot_df.loc[plot_df['nts'].isin(lengths)]['dSp']
        y1 = plot_df.loc[plot_df['nts'].isin(lengths)]['Sp']
        x2 = decoy_df.loc[decoy_df['nts'].isin(lengths)]['dSp']
        y2 = decoy_df.loc[decoy_df['nts'].isin(lengths)]['Sp']

        if per == 'y':
            # Add horizontal lines with info on percentiles populations

            red_plot_df, red_decoy_df = plot_df.loc[plot_df['nts'].isin(lengths)], decoy_df.loc[
                decoy_df['nts'].isin(lengths)]
            frac_value_95 = round(
                red_plot_df[red_plot_df > (red_decoy_df.quantile(0.95)['Sp'])].count()['Sp'] / red_plot_df[
                    'Sp'].count(), 2)
            frac_value_90 = round(
                red_plot_df[red_plot_df > (red_decoy_df.quantile(0.9)['Sp'])].count()['Sp'] / red_plot_df['Sp'].count(),
                2)

            plt.axhline(y=red_decoy_df.quantile(0.95)['Sp'], linestyle='--', linewidth=1, color='blue')
            plt.text(x2.max() - x2.max() / 7,
                     red_decoy_df.quantile(0.95)['Sp'] + red_decoy_df.quantile(0.95)['Sp'] / 100,
                     "95th p = " + str(frac_value_95), fontsize=5, color='blue')

            plt.axhline(y=red_decoy_df.quantile(0.9)['Sp'], linestyle='--', linewidth=1, color='cyan')
            plt.text(x2.max() - x2.max() / 7,
                     red_decoy_df.quantile(0.9)['Sp'] + red_decoy_df.quantile(0.95)['Sp'] / 100,
                     "90th p = " + str(frac_value_90), fontsize=5, color='cyan')

    # Adds the FDR line, if given as input
    if fdr_lines != 0:

        fdr_df = pd.read_csv(fdr_df_input, skiprows=3)

        # Create a dataframe with only entries where top_matches > 100
        fdr_df_100 = fdr_df.loc[fdr_df['top_match_targets(t)'] > 100]

        for fdr in fdr_lines:

            if not fdr_df_100.empty:
                # Draw FDR lines at given Sp cutoffs
                Sp_line = fdr_df_100.loc[(fdr_df_100['FDR(%)'] - fdr).abs().argsort()[:]]['Sp_cutoff'].iloc[0]
                plt.axhline(y=Sp_line, linestyle='--', linewidth=0.75, color='blue')
                plt.text(x2.max() - x2.max() / 20, Sp_line, "{}% FDR".format(int(round(fdr, 0))), fontsize=5,
                         color='blue')

            else:
                print(
                    "Warning!! The FDR assessment has been made on less than 100 top targets, therefore it is not significant and will not be shown in the graphs")

    # Determine the y axis min and max values
    max_val = [y1.max(), y2.max()]
    plt.ylim(0 - max(max_val) / 50, max(max_val) + max(max_val) / 25)

    plt.scatter(x1, y1, facecolors='royalblue', edgecolors='blue', marker=".", alpha=0.75)
    plt.scatter(x2, y2, facecolors='lightcoral', edgecolors='red', marker=".", alpha=0.75)

    # Place a box with info on the graph about the total number of points and parameters
    if info_box == 'y':
        if isinstance(MS1_offset_min, int) or isinstance(MS1_offset_max, int):
            textstr = ('\n'.join(('Nts. lengths = {}'.format(','.join([str(x) for x in lengths])),
                                  'Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  'MS1_offset = [{},{}]'.format(MS1_offset_min, MS1_offset_max),
                                  'Targets w/ decoys = {}'.format(targets_with_decoys))))
        else:
            textstr = ('\n'.join(('Nts. lengths = {}'.format(','.join([str(x) for x in lengths])),
                                  'Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  'Targets w/ decoys = {}'.format(targets_with_decoys))))

        props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
        plt.text(0.85, 0.9, textstr, transform=plt.gca().transAxes, fontsize=4, verticalalignment='top', bbox=props,
                 linespacing=1.75)

    # Label axis and title of the plot
    plt.xlabel('dSp'), plt.ylabel('Sp')
    plt.title('Sp vs dSp for the best match targets & decoys')

    # Create a legend
    dec_legend = plt.scatter([], [], facecolors='lightcoral', edgecolors='red', marker=".", alpha=0.75)
    tar_legend = plt.scatter([], [], facecolors='royalblue', edgecolors='blue', marker=".", alpha=0.75)
    plt.legend((tar_legend, dec_legend), ('targets ({})'.format(len(x1.index)), 'decoys ({})'.format(len(x2.index))))


def scatter_nts_vs_score(df, yl, y_bot, y_top, lengths, per, fdr_lines, Sp_cutoff, name_output, fdr_df_input, info_box,
                         targets_with_decoys, MS1_offset_min, MS1_offset_max):
    """
     Generate a scatter plot considering only the top scoring candidate for each combination of precursor/RT
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff]

    # Apply the offset cutoff if requested
    if isinstance(MS1_offset_min, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] >= MS1_offset_min]
    if isinstance(MS1_offset_max, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] <= MS1_offset_max]

        # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
        # are always kept
    decoy_df = indecoy_df[indecoy_df.M_RT.isin(plot_df.M_RT)]

    if targets_with_decoys == 'y':
        plot_df = plot_df[plot_df.M_RT.isin(decoy_df.M_RT)]

    if 'all' in lengths:
        x1 = plot_df['nts']
        y1 = plot_df[yl]
        x2 = decoy_df['nts']
        y2 = decoy_df[yl]

        if per == 'y':
            # Add horizontal lines with info on percentiles populations

            frac_value_95 = round(
                plot_df[plot_df > (decoy_df.quantile(0.95)['Sp'])].count()['Sp'] / plot_df['Sp'].count(), 2)
            frac_value_90 = round(
                plot_df[plot_df > (decoy_df.quantile(0.9)['Sp'])].count()['Sp'] / plot_df['Sp'].count(), 2)

            plt.axhline(y=decoy_df.quantile(0.95)['Sp'], linestyle='--', linewidth=1, color='blue')
            plt.text(x2.max() - x2.max() / 20, decoy_df.quantile(0.95)['Sp'] + decoy_df.quantile(0.95)['Sp'] / 100,
                     "95th p = " + str(frac_value_95), fontsize=5, color='blue')

            plt.axhline(y=decoy_df.quantile(0.9)['Sp'], linestyle='--', linewidth=1, color='cyan')
            plt.text(x2.max() - x2.max() / 20, decoy_df.quantile(0.9)['Sp'] + decoy_df.quantile(0.95)['Sp'] / 100,
                     "90th p = " + str(frac_value_90), fontsize=5, color='cyan')

    else:
        x1 = plot_df.loc[plot_df['nts'].isin(lengths)]['nts']
        y1 = plot_df.loc[plot_df['nts'].isin(lengths)][yl]
        x2 = decoy_df.loc[decoy_df['nts'].isin(lengths)]['nts']
        y2 = decoy_df.loc[decoy_df['nts'].isin(lengths)][yl]

        if per == 'y':
            # Add horizontal lines with info on percentiles populations

            red_plot_df, red_decoy_df = plot_df.loc[plot_df['nts'].isin(lengths)], decoy_df.loc[
                decoy_df['nts'].isin(lengths)]
            frac_value_95 = round(
                red_plot_df[red_plot_df > (red_decoy_df.quantile(0.95)['Sp'])].count()['Sp'] / red_plot_df[
                    'Sp'].count(), 2)
            frac_value_90 = round(
                red_plot_df[red_plot_df > (red_decoy_df.quantile(0.9)['Sp'])].count()['Sp'] / red_plot_df['Sp'].count(),
                2)

            plt.axhline(y=red_decoy_df.quantile(0.95)['Sp'], linestyle='--', linewidth=1, color='blue')
            plt.text(x2.max() - x2.max() / 20,
                     red_decoy_df.quantile(0.95)['Sp'] + red_decoy_df.quantile(0.95)['Sp'] / 100,
                     "95th p = " + str(frac_value_95), fontsize=5, color='blue')

            plt.axhline(y=red_decoy_df.quantile(0.9)['Sp'], linestyle='--', linewidth=1, color='cyan')
            plt.text(x2.max() - x2.max() / 20,
                     red_decoy_df.quantile(0.9)['Sp'] + red_decoy_df.quantile(0.95)['Sp'] / 100,
                     "90th p = " + str(frac_value_90), fontsize=5, color='cyan')

    plot_df = plot_df.astype({"nts": int})
    # Set the minimum value of x to be 1
    plt.xticks(list(range(plot_df['nts'].min(), plot_df['nts'].max() + 1)))

    if y_bot and y_top != 0:
        plt.ylim(y_bot, y_top)

    else:

        if y_top != 0:
            plt.ylim(top=y_top)

        else:
            plt.ylim(y_bot - plot_df[yl].max() / 50, plot_df[yl].max() + plot_df[yl].max() / 25)

    x1 -= 0.1
    x2 += 0.1

    # Adds the FDR line, if given as input
    if fdr_lines != 0:
        fdr_df = pd.read_csv(fdr_df_input, skiprows=3)

        # Create a dataframe with only entries where top_matches > 100
        fdr_df_100 = fdr_df.loc[fdr_df['top_match_targets(t)'] > 100]

        for fdr in fdr_lines:

            if not fdr_df_100.empty:
                # Draw FDR lines at given Sp cutoffs
                Sp_line = fdr_df_100.loc[(fdr_df_100['FDR(%)'] - fdr).abs().argsort()[:]]['Sp_cutoff'].iloc[0]
                plt.axhline(y=Sp_line, linestyle='--', linewidth=0.75, color='blue')
                plt.text(x2.max() - x2.max() / 20, Sp_line, "{}% FDR".format(round(fdr, 0)), fontsize=5, color='blue')

            else:
                print(
                    "Warning!! The FDR assessment has been made on less than 100 top targets, therefore it is not significant and will not be shown in the graphs")

    plt.scatter(x1, y1, facecolors='royalblue', edgecolors='blue', marker=".", alpha=0.75)
    plt.scatter(x2, y2, facecolors='lightcoral', edgecolors='red', marker=".", alpha=0.75)

    # Place a box with info on the graph about the total number of points and parameters
    if info_box == 'y':
        if isinstance(MS1_offset_min, int) or isinstance(MS1_offset_max, int):
            textstr = ('\n'.join(('Nts. lengths = {}'.format(','.join([str(x) for x in lengths])),
                                  'Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  'MS1_offset = [{},{}]'.format(MS1_offset_min, MS1_offset_max),
                                  'Targets w/ decoys = {}'.format(targets_with_decoys))))
        else:
            textstr = ('\n'.join(('Nts. lengths = {}'.format(','.join([str(x) for x in lengths])),
                                  'Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  'Targets w/ decoys = {}'.format(targets_with_decoys))))

        props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
        plt.text(0.855, 0.79, textstr, transform=plt.gca().transAxes, fontsize=4, verticalalignment='top', bbox=props,
                 linespacing=1.75)

    plt.xlabel('sequence_length'), plt.ylabel(yl)
    plt.title('Score (Sp) vs length for the best match targets & decoys')

    # Create a legend
    dec_legend = plt.scatter([], [], facecolors='lightcoral', edgecolors='red', marker=".", alpha=0.75)
    tar_legend = plt.scatter([], [], facecolors='royalblue', edgecolors='blue', marker=".", alpha=0.75)
    plt.legend((tar_legend, dec_legend), ('targets ({})'.format(len(x1.index)), 'decoys ({})'.format(len(x2.index))))


def box_nts_vs_score(df, y, y_bot, y_top, lengths, Sp_cutoff, targets_with_decoys, MS1_offset_min, MS1_offset_max):
    """
     Box plot nts vs scores of targets and decoys
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff]

    # Apply the offset cutoff if requested
    if isinstance(MS1_offset_min, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] >= MS1_offset_min]
    if isinstance(MS1_offset_max, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] <= MS1_offset_max]

    # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
    # are always kept
    decoy_df = indecoy_df[indecoy_df.M_RT.isin(plot_df.M_RT)]

    if targets_with_decoys == 'y':
        plot_df = plot_df[plot_df.M_RT.isin(decoy_df.M_RT)]

    plot_df = plot_df.astype({"nts": int})
    max_nts = plot_df['nts'].max()

    data1, data2 = [], []

    if 'all' in lengths:
        for i in range(3, max_nts + 1):
            data1.append((list(plot_df.loc[plot_df['nts'] == i, y])))
            data2.append((list(decoy_df.loc[decoy_df['nts'] == i, y])))

            # Delete the entries of shorter nts from the list, empty values cause problems with boxplots
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
            data1.append((list(plot_df.loc[plot_df['nts'] == int(i), y])))
            data2.append((list(decoy_df.loc[decoy_df['nts'] == int(i), y])))

    plt.xlabel('sequence_length'), plt.ylabel(y)

    if y_bot and y_top != 0:
        plt.ylim((y_bot, y_top))

    else:

        if y_top != 0:
            plt.ylim(top=y_top)

        else:
            plt.ylim((y_bot - plot_df['Sp'].max() / 50, plot_df['Sp'].max() + plot_df['Sp'].max() / 25))

    if 'all' in lengths:

        lenghts_labels = list(np.arange(plot_df['nts'].min(), plot_df['nts'].max() + 1))

        target = plt.boxplot(data1, positions=[np.float64(x) - 0.15 for x in lenghts_labels], manage_ticks=False,
                             labels=lenghts_labels, widths=0.25, sym='b.')
        decoys = plt.boxplot(data2, positions=[np.float64(x) + 0.15 for x in lenghts_labels], manage_ticks=False,
                             labels=lenghts_labels, widths=0.25, sym='r.')

        for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
            plt.setp(target[element], color='blue')
            plt.setp(decoys[element], color='red')


    else:
        lengths.sort()
        target = plt.boxplot(data1, labels=lengths, positions=[np.float64(x) - 0.15 for x in lengths],
                             manage_ticks=False, widths=0.25, sym='b.')
        decoys = plt.boxplot(data2, labels=lengths, positions=[np.float64(x) + 0.15 for x in lengths],
                             manage_ticks=False, widths=0.25, sym='r.')

        for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
            plt.setp(target[element], color='blue')
            plt.setp(decoys[element], color='red')

    # Create a legend
    plt.plot([], c='blue', label='targets')
    plt.plot([], c='red', label='decoys')
    plt.legend()


def hist_Sp(df, lengths, Sp_cutoff, targets_with_decoys, MS1_offset_min, MS1_offset_max, info_box):
    """
     Histogram with all Sp values for targets and decoys
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff]

    # Apply the offset cutoff if requested
    if isinstance(MS1_offset_min, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] >= MS1_offset_min]
    if isinstance(MS1_offset_max, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] <= MS1_offset_max]

    # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
    # are always kept
    decoy_df = indecoy_df[indecoy_df.M_RT.isin(plot_df.M_RT)]

    if targets_with_decoys == 'y':
        plot_df = plot_df[plot_df.M_RT.isin(decoy_df.M_RT)]

    # Apply the option about lengths if specified
    if 'all' not in lengths:
        plot_df = plot_df.loc[plot_df['nts'].isin(lengths)]
        decoy_df = decoy_df.loc[decoy_df['nts'].isin(lengths)]

    # Determine statistical descriptors for targets and decoys
    median_t, mean_t, std_t = plot_df['Sp'].median(), plot_df['Sp'].mean(), plot_df['Sp'].std()
    median_d, mean_d, std_d = decoy_df['Sp'].median(), decoy_df['Sp'].mean(), decoy_df['Sp'].std()

    plt.ylabel('Frequency')
    plt.title('Top targets & decoys scores [Sp]')

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
        if isinstance(MS1_offset_min, int) or isinstance(MS1_offset_max, int):
            textstr = ('\n'.join(('Nts. lengths = {}'.format(','.join([str(x) for x in lengths])),
                                  'Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  'Targets w/ decoys = {}'.format(targets_with_decoys),
                                  'MS1_offset = [{},{}]'.format(MS1_offset_min, MS1_offset_max),
                                  r'$\mu(t)=%.2f$ M(t)=%.2f $\sigma(t)=%.2f$' % (mean_t, median_t, std_t,),
                                  r'$\mu(d)=%.2f$ M(d)=%.2f $\sigma(d)=%.2f$' % (mean_d, median_d, std_d,))))
        else:
            textstr = ('\n'.join(('Nts. lengths = {}'.format(','.join([str(x) for x in lengths])),
                                  'Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  'Targets w/ decoys = {}'.format(targets_with_decoys),
                                  r'$\mu(t)=%.2f$ M(t)=%.2f $\sigma(t)=%.2f$' % (mean_t, median_t, std_t,),
                                  r'$\mu(d)=%.2f$ M(d)=%.2f $\sigma(d)=%.2f$' % (mean_d, median_d, std_d,))))
        props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
        plt.text(0.79, 0.78, textstr, transform=plt.gca().transAxes, fontsize=4, verticalalignment='top', bbox=props,
                 linespacing=1.75)

    plt.legend(loc='upper right')


def hist_top_Sp(df, lengths, Sp_cutoff, targets_with_decoys, MS1_offset_min, MS1_offset_max, info_box):
    """
     Histogram with top Sp values (dSp = 0) for targets and decoys
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff]

    # Apply the offset cutoff if requested
    if isinstance(MS1_offset_min, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] >= MS1_offset_min]
    if isinstance(MS1_offset_max, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] <= MS1_offset_max]

        # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
        # are always kept
    decoy_df = indecoy_df[indecoy_df.M_RT.isin(plot_df.M_RT)]

    if targets_with_decoys == 'y':
        plot_df = plot_df[plot_df.M_RT.isin(decoy_df.M_RT)]

    # Select only the the targets and decoys with dSp = 0, top targets and top decoys
    top_plot_df, top_decoy_df = plot_df[plot_df['dSp'] == 0], decoy_df[decoy_df['dSp'] == 0]

    # Apply the option about lengths if specified
    if 'all' not in lengths:
        top_plot_df = top_plot_df.loc[top_plot_df['nts'].isin(lengths)]
        top_decoy_df = top_decoy_df.loc[top_decoy_df['nts'].isin(lengths)]

    # Determine statistical descriptors for targets and decoys
    median_t, mean_t, std_t = top_plot_df['Sp'].median(), top_plot_df['Sp'].mean(), top_plot_df['Sp'].std()
    median_d, mean_d, std_d = top_decoy_df['Sp'].median(), top_decoy_df['Sp'].mean(), top_decoy_df['Sp'].std()

    plt.xlabel('Sp'), plt.ylabel('Frequency')
    plt.title('Scores (Sp) of target & decoys with dSp = 0')

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
        if isinstance(MS1_offset_min, int) or isinstance(MS1_offset_max, int):
            textstr = ('\n'.join(('Nts. lengths = {}'.format(','.join([str(x) for x in lengths])),
                                  'Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  'Targets w/ decoys = {}'.format(targets_with_decoys),
                                  'MS1_offset = [{},{}]'.format(MS1_offset_min, MS1_offset_max),
                                  r'$\mu(t)=%.2f$ M(t)=%.2f $\sigma(t)=%.2f$' % (mean_t, median_t, std_t,),
                                  r'$\mu(d)=%.2f$ M(d)=%.2f $\sigma(d)=%.2f$' % (mean_d, median_d, std_d,))))
        else:
            textstr = ('\n'.join(('Nts. lengths = {}'.format(','.join([str(x) for x in lengths])),
                                  'Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  'Targets w/ decoys = {}'.format(targets_with_decoys),
                                  r'$\mu(t)=%.2f$ M(t)=%.2f $\sigma(t)=%.2f$' % (mean_t, median_t, std_t,),
                                  r'$\mu(d)=%.2f$ M(d)=%.2f $\sigma(d)=%.2f$' % (mean_d, median_d, std_d,))))
        props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
        plt.text(0.79, 0.78, textstr, transform=plt.gca().transAxes, fontsize=4, verticalalignment='top', bbox=props,
                 linespacing=1.75)

    plt.legend(loc='upper right')


def hist_second_dSp(df, lengths, Sp_cutoff, targets_with_decoys, MS1_offset_min, MS1_offset_max, info_box):
    """
     Histogram with top Sp values (dSp = 0) for targets and decoys
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff]

    # Apply the offset cutoff if requested
    if isinstance(MS1_offset_min, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] >= MS1_offset_min]
    if isinstance(MS1_offset_max, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] <= MS1_offset_max]

    # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
    # are always kept
    decoy_df = indecoy_df[indecoy_df.M_RT.isin(plot_df.M_RT)]

    if targets_with_decoys == 'y':
        plot_df = plot_df[plot_df.M_RT.isin(decoy_df.M_RT)]

    # Select only the the targets and decoys with dSp = 0, top targets and top decoys
    sec_plot_df = plot_df[(plot_df['dSp_2nd_target'] >= 0) & (plot_df['dSp'] == 0)]

    # Apply the option about lengths if specified
    if 'all' not in lengths:
        sec_plot_df = sec_plot_df.loc[sec_plot_df['nts'].isin(lengths)]

    # Determine statistical descriptors for targets and decoys
    median_t, mean_t, std_t = sec_plot_df['dSp_2nd_target'].median(), sec_plot_df['dSp_2nd_target'].mean(), sec_plot_df[
        'dSp_2nd_target'].std()

    plt.xlabel('dSp'), plt.ylabel('Frequency')
    plt.title('dSp scores of the 2nd best match targets')

    plt.hist(sec_plot_df['dSp_2nd_target'], bins=20, alpha=0.5, label='targets ({})'.format(len(sec_plot_df.index)),
             edgecolor='black', color='blue')

    # Place a box with info on the graph about the total number of points and parameters
    if info_box == 'y':
        if isinstance(MS1_offset_min, int) or isinstance(MS1_offset_max, int):
            textstr = ('\n'.join(('Nts. lengths = {}'.format(','.join([str(x) for x in lengths])),
                                  'Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  'MS1_offset = [{},{}]'.format(MS1_offset_min, MS1_offset_max),
                                  r'$\mu(t)=%.2f$ M(t)=%.2f $\sigma(t)=%.2f$' % (mean_t, median_t, std_t,))))
        else:
            textstr = ('\n'.join(('Nts. lengths = {}'.format(','.join([str(x) for x in lengths])),
                                  'Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  r'$\mu(t)=%.2f$ M(t)=%.2f $\sigma(t)=%.2f$' % (mean_t, median_t, std_t,))))
        props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
        plt.text(0.79, 0.93, textstr, transform=plt.gca().transAxes, fontsize=4, verticalalignment='top', bbox=props,
                 linespacing=1.75)

    plt.legend(loc='upper right')


def ppm_errors_histogram(df, lengths, Sp_cutoff, label, match_file, MS1_offset_min, MS1_offset_max, info_box):
    """
     Histogram with the distribution of ppm incertitude on MS1 and MS2 matches for top Sp values (dSp = 0)
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff and select only the winning targets (dSp = 0)
    plot_df = inplot_df.loc[(inplot_df['Sp'] > Sp_cutoff) & (inplot_df['dSp'] == 0)]

    # Apply the offset cutoff if requested
    if isinstance(MS1_offset_min, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] >= MS1_offset_min]
    if isinstance(MS1_offset_max, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] <= MS1_offset_max]

     # Obtain the info on MS1/MS2 ppm offsets used in matching&scoring
    with open(match_file) as infile:
        head = [next(infile) for x in range(20)]

    for element in head:
        if "#MS1_ppm" in element:
            MS1_ppm = int(element.split()[1])
        if "#MS2_ppm" in element:
            MS2_ppm = int(element.split()[1])
        if "#MS1_offset" in element:
            MS1_offset = int(element.split()[1])
        if "#MS2_offset" in element:
            MS2_offset = int(element.split()[1])

    # Create a unique list with all the MS2 offset ppm values
    MS2_ppm_list = []

    MS2_ppm_series = plot_df['MS2_ppm'].dropna()
    for x in MS2_ppm_series:
        for i in x.split(","):
            MS2_ppm_list.append(np.int64(i))

    # Determine meand and std dev for MS1 and MS2 offset values
    median_MS1, mean_MS1, std_MS1 = plot_df['MS1_ppm'].median(), plot_df['MS1_ppm'].mean(), plot_df['MS1_ppm'].std()
    median_MS2, mean_MS2, std_MS2 = statistics.median(MS2_ppm_list), statistics.mean(MS2_ppm_list), statistics.stdev(
        MS2_ppm_list)

    # Add the mean and std dev values in the plot
    if label == 'MS1' and info_box == 'y':
        if isinstance(MS1_offset_min, int) or isinstance(MS1_offset_max, int):
            textstr = ('\n'.join(('Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  'MS1 set offset [ppm] = {}'.format(MS1_offset),
                                  'MS1_offset = [{},{}]'.format(MS1_offset_min, MS1_offset_max),
                                  r'$\mu=%.1f$ M=%.1f $\sigma=%.1f$' % (mean_MS1, median_MS1, std_MS1,),)))
        else:
            textstr = ('\n'.join(('Sp cutoff [targets] = {}'.format(Sp_cutoff),
                                  'MS1 set offset [ppm] = {}'.format(MS1_offset),
                                  r'$\mu=%.1f$ M=%.1f $\sigma=%.1f$' % (mean_MS1, median_MS1, std_MS1,),)))

        props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
        plt.text(0.03, 0.93, textstr, transform=plt.gca().transAxes, fontsize=4, verticalalignment='top', bbox=props,
                 linespacing=1.75)
        plt.ylabel('Frequency')

    elif label == 'MS2' and info_box == 'y':
        textstr = ('\n'.join(('Sp cutoff [targets] = {}'.format(Sp_cutoff),
                              'MS2 set offset [ppm] = {}'.format(MS2_offset),
                              r'$\mu=%.1f$ M=%.1f $\sigma=%.1f$' % (mean_MS2, median_MS2, std_MS2,),)))
        props = dict(boxstyle='round, pad = 1', facecolor='palegreen', edgecolor='green', alpha=0.5)
        plt.text(0.03, 0.93, textstr, transform=plt.gca().transAxes, fontsize=4, verticalalignment='top', bbox=props,
                 linespacing=1.75)
        plt.xlabel('ppm_offset'), plt.ylabel('Frequency')

    plt.title('ppm m/z offsets for the top {} matches (centered around the input offset)'.format(label))

    if label == 'MS1':
        plt.hist(plot_df['MS1_ppm'], bins=30, label='targets', range=[-MS1_ppm, MS1_ppm], edgecolor='black',
                 color='blue')
    else:
        plt.hist(MS2_ppm_list, bins=50, label='targets', range=[-MS2_ppm, MS2_ppm], edgecolor='black', color='red')


def scatter_nts_z_vs_score(df, yl, y_bot, y_top, lengths, Sp_cutoff, MS1_offset_min, MS1_offset_max,
                           targets_with_decoys):
    """
     Generate a scatter plot considering only the top scoring candidate for each combination of precursor/RT
     """
    # Keep only the top scoring entry for each combination of precursor + RT

    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df, decoy_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff], indecoy_df.loc[indecoy_df['Sp'] > Sp_cutoff]

    # Apply the offset cutoff if requested
    if isinstance(MS1_offset_min, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] >= MS1_offset_min]
        decoy_df = decoy_df.loc[decoy_df['MS1_ppm'] >= MS1_offset_min]
    if isinstance(MS1_offset_max, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] <= MS1_offset_max]
        decoy_df = decoy_df.loc[decoy_df['MS1_ppm'] <= MS1_offset_max]

    if targets_with_decoys == 'y':
        plot_df = plot_df[plot_df.M_RT.isin(decoy_df.M_RT)]

    plot_df = plot_df.sort_values(['nts', 'charge'], ascending=[True, True])

    if 'all' in lengths:
        x1 = plot_df['nts_z']
        y1 = plot_df[yl]
        x2 = decoy_df['nts_z']
        y2 = decoy_df[yl]

    else:
        x1 = plot_df.loc[plot_df['nts'].isin(lengths)]['nts_z']
        y1 = plot_df.loc[plot_df['nts'].isin(lengths)][yl]
        x2 = decoy_df.loc[decoy_df['nts'].isin(lengths)]['nts_z']
        y2 = decoy_df.loc[decoy_df['nts'].isin(lengths)][yl]

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

    plt.xlabel('nts_z'), plt.ylabel(yl)

    # Create a legend
    dec_legend = plt.scatter([], [], facecolors='lightcoral', edgecolors='red', marker=".", alpha=0.75)
    tar_legend = plt.scatter([], [], facecolors='royalblue', edgecolors='blue', marker=".", alpha=0.75)
    plt.legend((tar_legend, dec_legend), ('targets ({})'.format(len(x1.index)), 'decoys ({})'.format(len(x2.index))))


def box_nts_z_vs_score(df, y, y_bot, y_top, lengths, Sp_cutoff, MS1_offset_min, MS1_offset_max, targets_with_decoys):
    """
     Box plot
     """
    inplot_df, indecoy_df = df

    # Apply the Sp cutoff
    plot_df, decoy_df = inplot_df.loc[inplot_df['Sp'] > Sp_cutoff], indecoy_df.loc[indecoy_df['Sp'] > Sp_cutoff]

    # Apply the offset cutoff if requested
    if isinstance(MS1_offset_min, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] >= MS1_offset_min]
        decoy_df = decoy_df.loc[decoy_df['MS1_ppm'] >= MS1_offset_min]
    if isinstance(MS1_offset_max, int):
        plot_df = plot_df.loc[plot_df['MS1_ppm'] <= MS1_offset_max]
        decoy_df = decoy_df.loc[decoy_df['MS1_ppm'] <= MS1_offset_max]

    if targets_with_decoys == 'y':
        plot_df = plot_df[plot_df.M_RT.isin(decoy_df.M_RT)]

    plot_df = plot_df.sort_values(['nts', 'charge'], ascending=[True, True])

    if 'all' not in lengths:
        plot_df = plot_df[plot_df['nts'].isin(lengths)]

    # max_nts = plot_df['nts'].max()
    nts_z = plot_df['nts_z'].unique()

    if y_bot and y_top != 0:
        plt.ylim((y_bot, y_top))

    else:

        if y_top != 0:
            plt.ylim(top=y_top)

        else:
            plt.ylim((y_bot - plot_df['Sp'].max() / 50, plot_df['Sp'].max() + plot_df['Sp'].max() / 25))

    data = []
    for i in nts_z:
        data.append((list(plot_df.loc[plot_df['nts_z'] == i, y])))

    plt.xlabel('nts_z'), plt.ylabel(y)
    target = plt.boxplot(data, labels=nts_z.tolist(), sym='b.')

    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(target[element], color='blue')


def FDR_update(df, lengths, Sp_cutoff, name_output, MS1_offset_min, MS1_offset_max, targets_with_decoys):
    """
     Update the FDR table with the options selected while elaborating statistics and plotting the graphs
     """
    in_plot_df, in_decoy_df = df

    if 'all' in lengths:
        plot_df = in_plot_df.loc[in_plot_df['Sp'] > Sp_cutoff]

        # Apply the offset cutoff if requested
        if isinstance(MS1_offset_min, int):
            plot_df = plot_df.loc[plot_df['MS1_ppm'] >= MS1_offset_min]
        if isinstance(MS1_offset_max, int):
            plot_df = plot_df.loc[plot_df['MS1_ppm'] <= MS1_offset_max]

            # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
            # are always kept
        decoy_df = in_decoy_df[in_decoy_df.M_RT.isin(plot_df.M_RT)]

    else:
        plot_df = in_plot_df.loc[(in_plot_df['nts'].isin(lengths)) & (in_plot_df['Sp'] > Sp_cutoff)]

        # Apply the offset cutoff if requested
        if isinstance(MS1_offset_min, int):
            plot_df = plot_df.loc[plot_df['MS1_ppm'] >= MS1_offset_min]
        if isinstance(MS1_offset_max, int):
            plot_df = plot_df.loc[plot_df['MS1_ppm'] <= MS1_offset_max]

            # The Sp cutoff is applied only to targets, decoys that share the same mz_RT with targets above cutoff
            # are always kept
        decoy_df = in_decoy_df[in_decoy_df.M_RT.isin(plot_df.M_RT)]

        # Creates a dataframe of targets for FDR estimations, with only targets that have at least one decoy
    if targets_with_decoys == 'y':
        targets_FDR_df = plot_df[plot_df.M_RT.isin(decoy_df.M_RT)]
        targets_w_decoys_df = targets_FDR_df

    else:
        targets_FDR_df = plot_df
        targets_w_decoys_df = plot_df[plot_df.M_RT.isin(decoy_df.M_RT)]

    # Creates a dataframe with all values of dSp = 0
    merged_df = pd.concat([targets_FDR_df.loc[targets_FDR_df['dSp'] == 0], decoy_df.loc[decoy_df['dSp'] == 0]],
                          ignore_index=True, sort=True)
    merged_df = merged_df.drop_duplicates(subset='M_RT', keep=False)
    toptargets_df = merged_df[merged_df['molecule'].str.contains('decoy') == False]

    # Creates a dataframe with only the top match targets
    toptargets_df = targets_FDR_df.loc[targets_FDR_df['dSp'] == 0]

    decoy_max = np.float64(merged_df[merged_df['molecule'].str.contains('decoy') == True]['Sp'].max())

    if not np.isnan(decoy_max):
        step = round(decoy_max / 100, 3)

        d = {'Sp_cutoff': [], 'top_match_targets(t)': [], 'unique_sequences_targets': [], 'top_match_decoys(D)': [],
             'FDR(D/t)': [], 'FDR(%)': []}

        # Cycle through various Sp_cutoffs to calculate FDR
        for i in np.arange(0, round(decoy_max, 3) + step, step):
            tot = len(toptargets_df.loc[toptargets_df['Sp'] > i]['M'])
            if tot > 0:
                d['Sp_cutoff'].append(i), d['top_match_targets(t)'].append(tot)
                top_decoy = len(
                    merged_df.loc[(merged_df['molecule'].str.contains('decoy') == True) & (merged_df['Sp'] > i)]['M'])
                d['top_match_decoys(D)'].append(top_decoy), d['FDR(D/t)'].append(round(top_decoy / tot, 4)), d[
                    'FDR(%)'].append(round(top_decoy * 100 / tot, 2))
                d['unique_sequences_targets'].append(toptargets_df.loc[toptargets_df['Sp'] > i]['sequence'].nunique())

        out_df = pd.DataFrame(data=d)
        outcsv = out_df.to_csv(index=False)

    else:
        d = {'Sp_cutoff': [0], 'top_match_targets(t)': [len(merged_df[merged_df['Sp'] > 0]['M'])],
             'unique_sequences_targets': [toptargets_df['sequence'].nunique()], 'top_match_decoys(D)': [0],
             'FDR(D/t)': [0], 'FDR(%)': [0]}

        out_df = pd.DataFrame(data=d)
        outcsv = out_df.to_csv(index=False)

    # Outputs the csv file with the FDR estimations
    targets_with_decoys, targets_without_decoys = len(targets_w_decoys_df.index), len(targets_FDR_df.index) - len(
        targets_w_decoys_df.index)
    open("FDR_{}_edited.csv".format(name_output), 'w').writelines([
        "Targets_with_decoys,{}\nTargets_without_decoys,{}\nLengths,{}\n".format(
            targets_with_decoys, targets_without_decoys,
            ",".join(x for x in lengths)), outcsv])