import pandas as pd
import os
import sys
import pickle
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
# from matplotlib.pyplot import *
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import numpy as np
from datetime import datetime
import random
import matplotlib.gridspec as gridspec
from matplotlib import ticker
import statistics
import subprocess
from collections import defaultdict
from multiprocessing import Queue, Process
import itertools
from scipy.spatial.distance import braycurtis
from skbio.stats.ordination import pcoa
from mpl_toolkits.mplot3d import Axes3D

# plot that will be the three species pcoas with the its2 correlations in the plot below
def plot_pcoa_spp_18s_its2(is_three_d = False):
    info_df = generate_info_df_for_samples()
    fig_info_df = generate_fig_indo_df(info_df)


    # get the pcoa df
    pcoa_df_dict = generate_bray_curtis_distance_and_pcoa_spp()
    # setup figure
    spp_list = ['Porites', 'Pocillopora', 'Millepora']
    axarr=[]
    fig = plt.figure(figsize=(18, 10))

    # here we can read in the its2 data that we are going to need
    # read in the SymPortal relative abundance output

    path_to_tab_delim_count_DIV = '/home/humebc/projects/tara/initial_its2_processing/2018-10-21_08-59-37.620726.DIVs.relative.txt'

    path_to_tab_delim_count_type = '/home/humebc/projects/tara/initial_its2_processing/34_init_tara_standalone_all_samps_151018_2018-10-21_08-45-56.507454.profiles.relative.txt'

    smp_id_to_smp_name_dict, smp_name_to_smp_id_dict, sp_output_df_div = process_div_df(path_to_tab_delim_count_DIV)

    # now read in a process the type df too.

    colour_dict_type, sp_output_df_type, sorted_type_prof_names_by_local_abund, max_n_cols_type, max_n_rows_type, num_leg_cells_type, = process_type_df(
        path_to_tab_delim_count_type)

    colour_dict_div, max_n_cols_div, max_n_rows_div, num_leg_cells_div, ordered_list_of_seqs = get_div_colour_dict_and_ordered_list_of_seqs(
        sp_output_df_div)

    # the ordered_list_of_seqs can also be used for the plotting order

    # we should consider doing a plot per clade but for the time being lets start by doing a single plot that will
    # contain all of the clades

    # if we are plotting this in companion with an ITS2 type profile output then we will be passed a
    # sample_order_list. It is very useful to have the ITS2 type profile output figure and the seq figure
    # in the same sample order for direct comparison

    ordered_sample_list = sp_output_df_type.index.values.tolist()
    # pickle.dump(ordered_sample_list, open('ordered_sample_list_for_18s_work.pickle', 'wb'))
    # let's reorder the columns and rows of the sp_output_df according to the sequence sample and sequence
    # order so that plotting the data is easier

    sp_output_df_div = sp_output_df_div[ordered_list_of_seqs]



    # this is going to be quite a complicated setup but should be worth it
    #
    gs = plt.GridSpec(8, 21, figure=fig, height_ratios=[3, 0.2, 0.2, 1, 0.2, 1, 0.2, 1], width_ratios=[0.2, 0.2, 1, 0.2, 1, 0.2, 1, 0.9, 1, 0.2, 1, 0.2, 1, 0.9,1, 0.2, 1, 0.2, 1, 0.2, 0.2])
    # we can have several axes lists so that we can populate them one by one
    # first lets make the pcoa plot axes list
    pcoa_axes_list = []
    #porites ax
    pcoa_axes_list.append(plt.subplot(gs[0, 2:7]))
    # pocillopora ax
    # test = plt.subplot(gs[0, 7])
    pcoa_axes_list.append(plt.subplot(gs[0, 8:13]))
    # millepora ax
    pcoa_axes_list.append(plt.subplot(gs[0, 14:19]))

    # now make the its2 matrix of plots
    its2_axes_list = []
    for start_y, start_x in [(3,2),(3,8),(3,14)]:
        for i in range(0, 5, 2):
            for j in range(0, 5, 2):
                ind_x = j + start_x
                ind_y = i + start_y
                its2_axes_list.append(plt.subplot(gs[ind_y, ind_x]))

    #todo we can add the rest of the axes e.g. for the headers later on


    for spp in spp_list:

        ax = pcoa_axes_list[spp_list.index(spp)]


        df_of_spp = pcoa_df_dict[spp]

        # x values
        x_values = df_of_spp['PC1'].values.tolist()[:-1]

        # y values
        y_values = df_of_spp['PC2'].values.tolist()[:-1]



        samples_of_spp = df_of_spp.index.values.tolist()[:-1]
        # get list of colours and list of markers
        # colours can designate islands
        # TODO we will need to look up the colour according to the its2 type profile designation of the sample


        island_colour_dict = {'ISLAND06':'#C0C0C0', 'ISLAND10':'#808080', 'ISLAND15':'#000000'}
        island_colour_list = [island_colour_dict[fig_info_df.loc[smp, 'island']] for smp in samples_of_spp]

        # shapes can designate sites
        site_marker_dict = {'SITE01': '^', 'SITE02': 'o', 'SITE03': 's'}
        site_marker_list = [site_marker_dict[fig_info_df.loc[smp, 'site']] for smp in samples_of_spp]


        # plot the points
        for x_val, y_val, col, mark in zip(x_values, y_values, island_colour_list, site_marker_list):
            ax.scatter(x_val, y_val, c=col, marker=mark)

        # add axes labels
        ax.set_xlabel('PC1; explained = {}'.format('%.3f' % df_of_spp['PC1'][-1]))
        ax.set_ylabel('PC2; explained = {}'.format('%.3f' % df_of_spp['PC2'][-1]))
        # set axis title
        ax.set_title('{}'.format(spp))




    fig.show()
    if not is_three_d:
        plt.savefig('spp_pcoa_with_pc3_pc4.png'.format())
        plt.savefig('spp_pcoa_with_pc3_pc4.svg'.format())

    return

def get_div_colour_dict_and_ordered_list_of_seqs(sp_output_df_div):
    colour_palette_div = get_colour_list()
    grey_palette_div = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']
    # get a list of the sequences in order of their abundance and use this list to create the colour dict
    # the abundances can be got by simply summing up the columns making sure to ommit the last columns
    abundance_dict = {}
    for col in list(sp_output_df_div):
        abundance_dict[col] = sum(sp_output_df_div[col])
    # get the names of the sequences sorted according to their totalled abundance
    ordered_list_of_seqs = [x[0] for x in sorted(abundance_dict.items(), key=lambda x: x[1], reverse=True)]
    # create the colour dictionary that will be used for plotting by assigning a colour from the colour_palette
    # to the most abundant seqs first and after that cycle through the grey_pallette assigning colours
    # If we aer only going to have a legend that is cols x rows as shown below, then we should only use
    # that many colours in the plotting.
    max_n_cols = 8
    max_n_rows = 7
    num_leg_cells = max_n_cols * max_n_rows
    colour_dict_div = {}
    for i in range(len(ordered_list_of_seqs)):
        if i < num_leg_cells:
            colour_dict_div[ordered_list_of_seqs[i]] = colour_palette_div[i]
        else:
            grey_index = i % len(grey_palette_div)
            colour_dict_div[ordered_list_of_seqs[i]] = grey_palette_div[grey_index]
    return colour_dict_div, max_n_cols, max_n_rows, num_leg_cells, ordered_list_of_seqs

def process_type_df(path_to_tab_delim_count_type):
    sp_output_df_type = pd.read_csv(path_to_tab_delim_count_type, sep='\t', lineterminator='\n',
                                    skiprows=[0, 1, 2, 3, 5],
                                    header=None)
    # get a list of tups that are the seq names and the abundances zipped together
    type_profile_to_abund_tup_list = [(name, int(abund)) for name, abund in
                                      zip(sp_output_df_type.iloc[1][2:].values.tolist(),
                                          sp_output_df_type.iloc[0][2:].values.tolist())]
    # convert the names that are numbers into int strings rather than float strings.
    int_temp_list = []
    for name_abund_tup in type_profile_to_abund_tup_list:
        try:
            int_temp_list.append((str(int(name_abund_tup[0])), int(name_abund_tup[1])))
        except:
            int_temp_list.append((name_abund_tup[0], int(name_abund_tup[1])))
    type_profile_to_abund_tup_list = int_temp_list
    # need to drop the rows that contain the sequence accession and species descriptions
    for i, row_name in enumerate(sp_output_df_type.iloc[:, 0]):
        if 'Sequence accession' in row_name:
            # then we want to drop all rows from here until the end
            index_to_drop_from = i
            break
    sp_output_df_type = sp_output_df_type.iloc[:index_to_drop_from]
    # now drop the sample name columns
    sp_output_df_type.drop(columns=1, inplace=True)
    # make headers
    sp_output_df_type.columns = ['sample_id'] + [a[0] for a in type_profile_to_abund_tup_list]
    # now drop the local abund row and promote the its2_type_prof names to columns headers.
    sp_output_df_type.drop(index=[0, 1], inplace=True)
    sp_output_df_type = sp_output_df_type.set_index(keys='sample_id', drop=True).astype('float')

    max_n_cols = 4
    max_n_rows = 7
    num_leg_cells = max_n_cols * max_n_rows

    # we will use the col headers as the its2 type profile order for plotting but we
    # we should colour according to the abundance of the its2 type profiles
    # as we don't want to run out of colours by the time we get to profiles that are very abundant.
    # The sorted_type_prof_names_by_local_abund object has the names of the its2 type profile in order of abundance
    # we will use the index order as the order of samples to plot
    # create the colour dictionary that will be used for plotting by assigning a colour from the colour_palette
    # to the most abundant seqs first and after that cycle through the grey_pallette assigning colours
    sorted_type_prof_names_by_local_abund = [a[0] for a in
                                             sorted(type_profile_to_abund_tup_list, key=lambda x: x[1], reverse=True)]

    # we should plot sample by sample and its2 type by its2 type in the order of the output
    # the problem with doing he convert_to_pastel is that the colours become very similar
    # colour_palette = convert_to_pastel(get_colour_list())
    # Rather, I will attempt to generate a quick set of colours that are pastel and have a minimum distance
    # rule for any colours that are generated from each other.
    # let's do this for 50 colours to start with and see how long it takes.
    # turns out it is very quick. Easily quick enough to do dynamically.
    # When working with pastel colours (i.e. mixing with 255,255,255 it is probably best to work with a smaller dist cutoff

    colour_dict_type = pickle.load(open('/home/humebc/projects/tara/initial_its2_processing/colour_dict_type.pickle'.format(os.getcwd()), 'rb'))

    return colour_dict_type, sp_output_df_type, sorted_type_prof_names_by_local_abund, max_n_cols, max_n_rows, num_leg_cells

def process_div_df(path_to_tab_delim_count_DIV):
    sp_output_df = pd.read_csv(path_to_tab_delim_count_DIV, sep='\t', lineterminator='\n', header=0, index_col=0)

    # In order to be able to drop the DIV row at the end and the meta information rows, we should
    # drop all rows that are after the DIV column. We will pass in an index value to the .drop
    # that is called here. To do this we need to work out which index we are working with
    index_values_as_list = sp_output_df.index.values.tolist()
    for i in range(-1, -(len(index_values_as_list)), -1):
        if index_values_as_list[i].startswith('DIV'):
            # then this is the index (in negative notation) that we need to cut from
            meta_index_to_cut_from = i
            break
    sp_output_df = sp_output_df.iloc[:meta_index_to_cut_from]

    # create sample id to sample name dict
    smp_id_to_smp_name_dict = {ID: '_'.join(nm.split('_')[:3]) for ID, nm in
                               zip(sp_output_df.index.values.tolist(), sp_output_df['sample_name'].values.tolist())}
    smp_name_to_smp_id_dict = {'_'.join(nm.split('_')[:3]): ID for ID, nm in
                               zip(sp_output_df.index.values.tolist(), sp_output_df['sample_name'].values.tolist())}

    # now lets drop the QC columns from the SP output df and also drop the clade summation columns
    # we will be left with just clumns for each one of the sequences found in the samples
    sp_output_df.drop(columns=['sample_name', 'noName Clade A', 'noName Clade B', 'noName Clade C', 'noName Clade D',
                               'noName Clade E', 'noName Clade F', 'noName Clade G', 'noName Clade H',
                               'noName Clade I', 'raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                               'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                               'post_taxa_id_absolute_non_symbiodinium_seqs',
                               'post_taxa_id_unique_non_symbiodinium_seqs',
                               'size_screening_violation_absolute', 'size_screening_violation_unique',
                               'post_med_absolute', 'post_med_unique'
                               ], inplace=True)
    sp_output_df = sp_output_df.astype('float')
    return smp_id_to_smp_name_dict, smp_name_to_smp_id_dict, sp_output_df

# plotting of the pcoa coordinates for each species
# I'm thinking of using shapes and colours for doing sites and islands
def plot_pcoa_spp(is_three_d = False):
    info_df = generate_info_df_for_samples()
    fig_info_df = generate_fig_indo_df(info_df)

    # for each species
    # get the pcoa df
    pcoa_df_dict = generate_bray_curtis_distance_and_pcoa_spp()
    # setup figure
    spp_list = ['Porites', 'Pocillopora', 'Millepora']
    axarr=[]
    fig = plt.figure(figsize=(18, 10))

    gs = plt.GridSpec(3, 6, figure=fig, width_ratios=[1, 0.2, 1, 0.2, 1, 0.5], height_ratios=[1,0.2,1])
    for j in range(0,3,2):
        for i in range(0, 5, 2):
            if is_three_d:
                if j == 1:
                    axarr.append(fig.add_subplot(gs[j,i], projection='3d'))
                else:
                    axarr.append(plt.subplot(gs[j,i]))
            else:
                axarr.append(plt.subplot(gs[j, i]))

    legend_ax = plt.subplot(gs[:, 5])
    remove_axes_but_allow_labels(legend_ax)
    for spp in spp_list:

        ax = axarr[spp_list.index(spp)]
        ax_second = axarr[spp_list.index(spp) + 3]

        df_of_spp = pcoa_df_dict[spp]

        # x values
        x_values = df_of_spp['PC1'].values.tolist()[:-1]

        # y values
        y_values = df_of_spp['PC2'].values.tolist()[:-1]

        # z values
        z_values = df_of_spp['PC3'].values.tolist()[:-1]

        # pc4 values
        pc4_values = df_of_spp['PC4'].values.tolist()[:-1]

        samples_of_spp = df_of_spp.index.values.tolist()[:-1]
        # get list of colours and list of markers
        # colours can designate islands
        island_colour_dict = {'ISLAND06':'#C0C0C0', 'ISLAND10':'#808080', 'ISLAND15':'#000000'}
        island_colour_list = [island_colour_dict[fig_info_df.loc[smp, 'island']] for smp in samples_of_spp]

        # shapes can designate sites
        site_marker_dict = {'SITE01': '^', 'SITE02': 'o', 'SITE03': 's'}
        site_marker_list = [site_marker_dict[fig_info_df.loc[smp, 'site']] for smp in samples_of_spp]


        # plot the points
        for x_val, y_val, col, mark in zip(x_values, y_values, island_colour_list, site_marker_list):
            ax.scatter(x_val, y_val, c=col, marker=mark)

        # add axes labels
        ax.set_xlabel('PC1; explained = {}'.format('%.3f' % df_of_spp['PC1'][-1]))
        ax.set_ylabel('PC2; explained = {}'.format('%.3f' % df_of_spp['PC2'][-1]))
        # set axis title
        ax.set_title('{}'.format(spp))
        if is_three_d:
            # then lets plot the 3d equivalent below the 3d figs
            # plot the points
            for x_val, y_val, z_val, col, mark in zip(x_values, y_values, z_values, island_colour_list, site_marker_list):
                ax_second.scatter(x_val, y_val, z_val, c=col, marker=mark)

            # add axes labels
            ax_second.set_xlabel('PC1; explained = {}'.format('%.3f' % df_of_spp['PC1'][-1]))
            ax_second.set_ylabel('PC2; explained = {}'.format('%.3f' % df_of_spp['PC2'][-1]))
            ax_second.set_zlabel('PC3; explained = {}'.format('%.3f' % df_of_spp['PC3'][-1]))


        else:
            # else lets just plot out the 3rd and 4th PCs below
            # plot the points
            for z_val, pc4_val, col, mark in zip(z_values, pc4_values, island_colour_list, site_marker_list):
                ax_second.scatter(z_val, pc4_val, c=col, marker=mark)

            # add axes labels
            ax_second.set_xlabel('PC3; explained = {}'.format('%.3f' % df_of_spp['PC3'][-1]))
            ax_second.set_ylabel('PC4; explained = {}'.format('%.3f' % df_of_spp['PC4'][-1]))

    # here we should put together the legend axis
    vert_leg_axis(island_colour_dict, legend_ax, site_marker_dict)

    fig.show()
    if not is_three_d:
        plt.savefig('spp_pcoa_with_pc3_pc4.png'.format())
        plt.savefig('spp_pcoa_with_pc3_pc4.svg'.format())

    return

def plot_pcoa_spp_island():
    info_df = generate_info_df_for_samples()
    fig_info_df = generate_fig_indo_df(info_df)

    # for each species
    # get the pcoa df
    pcoa_df_dict = generate_bray_curtis_distance_and_pcoa_spp_and_island()
    # setup figure
    spp_list = ['PORITES', 'POCILLOPORA', 'MILLEPORA']

    island_list = ['ISLAND06', 'ISLAND10', 'ISLAND15']

    axarr = []

    fig = plt.figure(figsize=(18, 10))

    gs = plt.GridSpec(5, 7, figure=fig, width_ratios=[1, 0.2, 1, 0.2, 1,0.2, 0.5], height_ratios=[1, 0.2, 1, 0.2, 1])
    for j in range(0, 5, 2):
        for i in range(0, 5, 2):
            axarr.append(plt.subplot(gs[j, i]))

    legend_ax = plt.subplot(gs[:, 6])
    remove_axes_but_allow_labels(legend_ax)
    ax_count = 0
    for spp in spp_list:
        for island in island_list:

            ax = axarr[ax_count]

            df_of_spp_island = pcoa_df_dict['{}_{}'.format(spp, island)]

            # x values
            x_values = df_of_spp_island['PC1'].values.tolist()[:-1]

            # y values
            y_values = df_of_spp_island['PC2'].values.tolist()[:-1]

            samples_of_spp = df_of_spp_island.index.values.tolist()[:-1]
            # get list of colours and list of markers
            # colours can designate islands
            island_colour_dict = {'ISLAND06': '#C0C0C0', 'ISLAND10': '#808080', 'ISLAND15': '#000000'}
            island_colour_list = [island_colour_dict[fig_info_df.loc[smp, 'island']] for smp in samples_of_spp]

            # shapes can designate sites
            site_marker_dict = {'SITE01': '^', 'SITE02': 'o', 'SITE03': 's'}
            site_marker_list = [site_marker_dict[fig_info_df.loc[smp, 'site']] for smp in samples_of_spp]

            # plot the points
            for x_val, y_val, col, mark in zip(x_values, y_values, island_colour_list, site_marker_list):
                ax.scatter(x_val, y_val, c=col, marker=mark)

            # add axes labels
            ax.set_xlabel('PC1; explained = {}'.format('%.3f' % df_of_spp_island['PC1'][-1]))
            ax.set_ylabel('PC2; explained = {}'.format('%.3f' % df_of_spp_island['PC2'][-1]))
            # set axis title
            if ax_count in [0,1,2]:
                ax.set_title('{}'.format(island) , fontsize='large', fontweight='bold')
            if ax_count in [2,5,8]:
                ax2 = ax.twinx()
                ax2.yaxis.set_ticklabels([])
                ax2.yaxis.set_ticks_position('none')
                ax2.set_ylabel(spp, fontsize='large', fontweight='bold')

            ax_count += 1

    # here we should put together the legend axis
    vert_leg_axis(island_colour_dict, legend_ax, site_marker_dict)

    fig.show()
    plt.savefig('spp_island_pcoa.png'.format())
    plt.savefig('spp_island_pcoa.svg'.format())
    return


def vert_leg_axis(colour_dict, legend_ax, site_dict):
    legend_ax.set_ylim(0, 1)
    legend_ax.set_xlim(0, 1)
    legend_ax.invert_yaxis()
    number_of_icons = 6
    island_list = ['ISLAND06', 'ISLAND10', 'ISLAND15']
    site_list = ['SITE01', 'SITE02', 'SITE03']
    icon_list = []
    # first populate the island icons
    for island in island_list:
        icon_list.append((colour_dict[island], site_dict['SITE01']))
    # then populate the site icons
    for site in site_list:
        icon_list.append((colour_dict['ISLAND06'], site_dict[site]))
    # lets assume that the axis is divided into 20 spaces for icons
    max_number_icons = 20
    # the  icon position should be mid way so max_number_icons
    # first icon position
    first_icon_position = int((max_number_icons - number_of_icons) / 2)
    pos_counter = first_icon_position
    for i in range(len(icon_list)):
        y_val_for_icon_and_text = (1 / max_number_icons) * pos_counter
        x_val_for_icon = 0.1
        x_val_for_text = x_val_for_icon + 0.2
        legend_ax.scatter(x=x_val_for_icon, y=y_val_for_icon_and_text, c=icon_list[i][0], marker=icon_list[i][1], s=100)

        if int(i / 3) == 0:
            legend_ax.text(s=island_list[i], x=x_val_for_text, y=y_val_for_icon_and_text)
        elif int(i / 3) == 1:
            legend_ax.text(s=site_list[i % 3], x=x_val_for_text, y=y_val_for_icon_and_text)
        pos_counter += 1


# I want to produce distance matrices for each of the coral spp.
# I then will plot a PCOA for each of the coral species.
# I will try to reuse as much of the SP code as possible for doing this.
# To start with I will do a bray curtis distance as this doesn't require the sequences to be aligned
def generate_bray_curtis_distance_and_pcoa_spp():
    # Read in the minor div dataframe which should have normalised abundances in them
    # For each sample we have a fasta that we can read in which has the normalised (to 1000) sequences
    # For feeding into med. We can use a default dict to collect the sequences and abundances from this fairly
    # simply.
    # This is likely best done on for each sample outside of the pairwise comparison to save on redoing the same
    # collection of the sequences.
    # get info_df
    info_df = generate_info_df_for_samples()

    # get fig_info_df
    fig_info_df = generate_fig_indo_df(info_df)

    if os.path.isfile('{}/minor_div_abundance_dict.pickle'.format(os.getcwd())):
        minor_div_abundance_dict = pickle.load(open('{}/minor_div_abundance_dict.pickle'.format(os.getcwd()), 'rb'))

    else:

        # this dict will have sample name as key and then a dict as value with seq to abundance values
        # we can then work with this for the pairwise comparison
        minor_div_abundance_dict = {}


        for ind in fig_info_df.index.values.tolist():
            sample_dir = fig_info_df.loc[ind, 'sample_dir']
            with open('{}/fasta_for_med.fasta'.format(sample_dir), 'r') as f:
                sample_fasta = [line.rstrip() for line in f]

            sample_minor_abundance_dict = defaultdict(int)
            for i in range(0, len(sample_fasta), 2):
                sample_minor_abundance_dict[sample_fasta[i+1]] += 1

            # here we have the dict popoulated for the sample
            # we can now add this to the minor_div_abundace_dict
            minor_div_abundance_dict[ind] = sample_minor_abundance_dict

        # we should now pickle out this sample_minor_abundance_dict
        pickle.dump(minor_div_abundance_dict, open('{}/minor_div_abundance_dict.pickle'.format(os.getcwd()), 'wb'))


    # For each of the spp.
    spp_pcoa_df_dict = {}
    for spp in ['Porites', 'Pocillopora', 'Millepora']:
        if os.path.isfile('{}/spp_pcoa_df_{}.pickle'.format(os.getcwd(), spp)):
            spp_pcoa_df = pickle.load(open('{}/spp_pcoa_df_{}.pickle'.format(os.getcwd(), spp), 'rb'))
        else:
            # Get a list of the samples that we should be working with
            sample_names_of_spp = fig_info_df.loc[fig_info_df['genus'] == spp.upper()].index.values.tolist()

            #TODO remov the two porites species form this that seem to be total outliers
            if spp == 'Porites':
                sample_names_of_spp.remove('CO0001674')
                sample_names_of_spp.remove('CO0001669')


            if os.path.isfile('{}/spp_distance_dict_{}.pickle'.format(os.getcwd(), spp)):
                spp_distance_dict = pickle.load(open('{}/spp_distance_dict_{}.pickle'.format(os.getcwd(), spp), 'rb'))

            else:
                # Create a dictionary that will hold the distance between the two samples
                spp_distance_dict = {}

                # For pairwise comparison of each of these sequences
                for smp_one, smp_two in itertools.combinations(sample_names_of_spp, 2):
                    print('Calculating distance for {}_{}'.format(smp_one, smp_two))
                    # Get a set of the sequences found in either one of the samples
                    smp_one_abund_dict = minor_div_abundance_dict[smp_one]
                    smp_two_abund_dict = minor_div_abundance_dict[smp_two]
                    list_of_seqs_of_pair = []
                    list_of_seqs_of_pair.extend(list(smp_one_abund_dict.keys()))
                    list_of_seqs_of_pair.extend(list(smp_two_abund_dict.keys()))
                    list_of_seqs_of_pair = list(set(list_of_seqs_of_pair))


                    # then create a list of abundances for sample one by going through the above list and checking
                    sample_one_abundance_list = []
                    for seq_name in list_of_seqs_of_pair:
                        if seq_name in smp_one_abund_dict.keys():
                            sample_one_abundance_list.append(smp_one_abund_dict[seq_name])
                        else:
                            sample_one_abundance_list.append(0)

                    # then create a list of abundances for sample two by going through the above list and checking
                    sample_two_abundance_list = []
                    for seq_name in list_of_seqs_of_pair:
                        if seq_name in smp_two_abund_dict.keys():
                            sample_two_abundance_list.append(smp_two_abund_dict[seq_name])
                        else:
                            sample_two_abundance_list.append(0)


                    # Do the Bray Curtis.
                    distance = braycurtis(sample_one_abundance_list, sample_two_abundance_list)

                    # Add the distance to the dictionary using both combinations of the sample names
                    spp_distance_dict['{}_{}'.format(smp_one, smp_two)] = distance
                    spp_distance_dict['{}_{}'.format(smp_two, smp_one)] = distance

                # doing this takes a bit of time so let's pickle it out
                pickle.dump(spp_distance_dict, open('{}/spp_distance_dict_{}.pickle'.format(os.getcwd(), spp), 'wb'))

            # Generate the distance out file from this dictionary
            # from this dict we can produce the distance file that can be passed into the generate_PCoA_coords method
            distance_out_file = [len(sample_names_of_spp)]
            for sample_outer in sample_names_of_spp:
                # The list that will hold the line of distance. This line starts with the name of the sample
                temp_sample_dist_string = [sample_outer]

                for sample_inner in sample_names_of_spp:
                    if sample_outer == sample_inner:
                        temp_sample_dist_string.append(0)
                    else:
                        temp_sample_dist_string.append(spp_distance_dict[
                                                          '{}_{}'.format(sample_outer, sample_inner)])
                distance_out_file.append('\t'.join([str(distance_item) for distance_item in temp_sample_dist_string]))

            # from here we can hopefully rely on the rest of the methods as they already are. The .dist file should be
            # written out

            dist_out_path = '{}/bray_curtis_within_spp_sample_distances_{}.dist'.format(os.getcwd(), spp)

            with open(dist_out_path, 'w') as f:
                for line in distance_out_file:
                    f.write('{}\n'.format(line))

            # Feed this into the generate_PCoA_coords method
            spp_pcoa_df = generate_PCoA_coords(distance_out_file, spp)
            pickle.dump(spp_pcoa_df, open('{}/spp_pcoa_df_{}.pickle'.format(os.getcwd(), spp), 'wb'))
        spp_pcoa_df_dict[spp] = spp_pcoa_df
    return spp_pcoa_df_dict


# I also want to generate a pcoa df set for each island within each species instead of with all of the islands
# considered per species.
# that is what we will generate here.
def generate_bray_curtis_distance_and_pcoa_spp_and_island():
    # Read in the minor div dataframe which should have normalised abundances in them
    # For each sample we have a fasta that we can read in which has the normalised (to 1000) sequences
    # For feeding into med. We can use a default dict to collect the sequences and abundances from this fairly
    # simply.
    # This is likely best done on for each sample outside of the pairwise comparison to save on redoing the same
    # collection of the sequences.
    # get info_df
    info_df = generate_info_df_for_samples()

    # get fig_info_df
    fig_info_df = generate_fig_indo_df(info_df)

    if os.path.isfile('{}/minor_div_abundance_dict.pickle'.format(os.getcwd())):
        minor_div_abundance_dict = pickle.load(open('{}/minor_div_abundance_dict.pickle'.format(os.getcwd()), 'rb'))

    else:

        # this dict will have sample name as key and then a dict as value with seq to abundance values
        # we can then work with this for the pairwise comparison
        minor_div_abundance_dict = {}


        for ind in fig_info_df.index.values.tolist():
            sample_dir = fig_info_df.loc[ind, 'sample_dir']
            with open('{}/fasta_for_med.fasta'.format(sample_dir), 'r') as f:
                sample_fasta = [line.rstrip() for line in f]

            sample_minor_abundance_dict = defaultdict(int)
            for i in range(0, len(sample_fasta), 2):
                sample_minor_abundance_dict[sample_fasta[i+1]] += 1

            # here we have the dict popoulated for the sample
            # we can now add this to the minor_div_abundace_dict
            minor_div_abundance_dict[ind] = sample_minor_abundance_dict

        # we should now pickle out this sample_minor_abundance_dict
        pickle.dump(minor_div_abundance_dict, open('{}/minor_div_abundance_dict.pickle'.format(os.getcwd()), 'wb'))


    # For each of the spp.
    spp_island_pcoa_df_dict = {}
    for spp in ['PORITES', 'POCILLOPORA', 'MILLEPORA']:
        for island in ['ISLAND06', 'ISLAND10', 'ISLAND15']:
            if os.path.isfile('{}/spp_island_pcoa_df_{}_{}.pickle'.format(os.getcwd(), spp, island)):
                spp_island_pcoa_df = pickle.load(open('{}/spp_island_pcoa_df_{}_{}.pickle'.format(os.getcwd(), spp, island), 'rb'))
            else:
                # Get a list of the samples that we should be working with
                sample_names_of_spp = fig_info_df.loc[(fig_info_df['genus'] == spp.upper()) & (fig_info_df['island'] == island.upper())].index.values.tolist()

                #TODO remov the two porites species form this that seem to be total outliers
                if spp == 'PORITES':
                    if 'CO0001674' in sample_names_of_spp:
                        sample_names_of_spp.remove('CO0001674')
                        sample_names_of_spp.remove('CO0001669')


                if os.path.isfile('{}/spp_island_distance_dict_{}_{}.pickle'.format(os.getcwd(), spp, island)):
                    spp_island_distance_dict = pickle.load(open('{}/spp_island_distance_dict_{}_{}.pickle'.format(os.getcwd(), spp, island), 'rb'))

                else:
                    # Create a dictionary that will hold the distance between the two samples
                    spp_island_distance_dict = {}

                    # For pairwise comparison of each of these sequences
                    for smp_one, smp_two in itertools.combinations(sample_names_of_spp, 2):
                        print('Calculating distance for {}_{}'.format(smp_one, smp_two))
                        # Get a set of the sequences found in either one of the samples
                        smp_one_abund_dict = minor_div_abundance_dict[smp_one]
                        smp_two_abund_dict = minor_div_abundance_dict[smp_two]
                        list_of_seqs_of_pair = []
                        list_of_seqs_of_pair.extend(list(smp_one_abund_dict.keys()))
                        list_of_seqs_of_pair.extend(list(smp_two_abund_dict.keys()))
                        list_of_seqs_of_pair = list(set(list_of_seqs_of_pair))


                        # then create a list of abundances for sample one by going through the above list and checking
                        sample_one_abundance_list = []
                        for seq_name in list_of_seqs_of_pair:
                            if seq_name in smp_one_abund_dict.keys():
                                sample_one_abundance_list.append(smp_one_abund_dict[seq_name])
                            else:
                                sample_one_abundance_list.append(0)

                        # then create a list of abundances for sample two by going through the above list and checking
                        sample_two_abundance_list = []
                        for seq_name in list_of_seqs_of_pair:
                            if seq_name in smp_two_abund_dict.keys():
                                sample_two_abundance_list.append(smp_two_abund_dict[seq_name])
                            else:
                                sample_two_abundance_list.append(0)


                        # Do the Bray Curtis.
                        distance = braycurtis(sample_one_abundance_list, sample_two_abundance_list)

                        # Add the distance to the dictionary using both combinations of the sample names
                        spp_island_distance_dict['{}_{}'.format(smp_one, smp_two)] = distance
                        spp_island_distance_dict['{}_{}'.format(smp_two, smp_one)] = distance

                    # doing this takes a bit of time so let's pickle it out
                    pickle.dump(spp_island_distance_dict, open('{}/spp_island_distance_dict_{}_{}.pickle'.format(os.getcwd(), spp, island), 'wb'))

                # Generate the distance out file from this dictionary
                # from this dict we can produce the distance file that can be passed into the generate_PCoA_coords method
                distance_out_file = [len(sample_names_of_spp)]
                for sample_outer in sample_names_of_spp:
                    # The list that will hold the line of distance. This line starts with the name of the sample
                    temp_sample_dist_string = [sample_outer]

                    for sample_inner in sample_names_of_spp:
                        if sample_outer == sample_inner:
                            temp_sample_dist_string.append(0)
                        else:
                            temp_sample_dist_string.append(spp_island_distance_dict[
                                                              '{}_{}'.format(sample_outer, sample_inner)])
                    distance_out_file.append('\t'.join([str(distance_item) for distance_item in temp_sample_dist_string]))

                # from here we can hopefully rely on the rest of the methods as they already are. The .dist file should be
                # written out

                dist_out_path = '{}/bray_curtis_within_spp_island_sample_distances_{}_{}.dist'.format(os.getcwd(), spp, island)

                with open(dist_out_path, 'w') as f:
                    for line in distance_out_file:
                        f.write('{}\n'.format(line))

                # Feed this into the generate_PCoA_coords method
                spp_island_pcoa_df = generate_PCoA_coords(distance_out_file, spp)
                pickle.dump(spp_island_pcoa_df, open('{}/spp_island_pcoa_df_{}_{}.pickle'.format(os.getcwd(), spp, island), 'wb'))
            spp_island_pcoa_df_dict['{}_{}'.format(spp, island)] = spp_island_pcoa_df
    return spp_island_pcoa_df_dict

def generate_PCoA_coords(raw_dist_file, spp):
    # simultaneously grab the sample names in the order of the distance matrix and put the matrix into
    # a twoD list and then convert to a numpy array
    temp_two_D_list = []
    sample_names_from_dist_matrix = []
    for line in raw_dist_file[1:]:
        temp_elements = line.split('\t')
        sample_names_from_dist_matrix.append(temp_elements[0])
        temp_two_D_list.append([float(a) for a in temp_elements[1:]])
    uni_frac_dist_array = np.array(temp_two_D_list)
    sys.stdout.write('\rcalculating PCoA coordinates')
    pcoA_full_path = '{}/pcoa_coords_{}.csv'.format(os.getcwd(), spp)
    pcoa_df = pcoa(uni_frac_dist_array)

    # rename the dataframe index as the sample names
    pcoa_df.samples['sample'] = sample_names_from_dist_matrix
    renamed_dataframe = pcoa_df.samples.set_index('sample')

    # now add the variance explained as a final row to the renamed_dataframe

    renamed_dataframe = renamed_dataframe.append(pcoa_df.proportion_explained.rename('proportion_explained'))


    return renamed_dataframe




# This code will create MED node profiles for each of the samples, disregarding the most abundant sequence
# NB after creating the MED sequences it was not much different to the raw sequences. This is likely becauase
# all of the sequences were found at such low and even abundances. I think, moving forwards we should just stick with
# working with the raw sequencs rather than the MED nodes.
def create_MED_node_sample_abundance_df_for_minor_intras():
    # we should go sample by sample using the abundance_df and write out a fasta for the MED
    # this will be a subsample of 1000 sequences of the minor DIVs, which will be run with a dynamic
    # m value on MED.
    # we will then read in the MED and create a dict for each sample that is MED node sequence to relabund
    # we will then come out of MP and go serial to collect all of the MED nodes and get a cumulative abundance
    # across all samples. We will order this and use the MED nodes as columns for the new minor_intras_abund_df
    # we will then go back through all of the samples in serial and populate the MED df
    # at this point we can pickle this df out and we can put it directly into plotting code of below
    # to create the same minor intra figure that we had before but with having done the MED.

    if os.path.isfile('{}/minor_intra_med_abund_df.pickle'.format(os.getcwd())):
        minor_intra_med_abund_df = pickle.load(open('{}/minor_intra_med_abund_df.pickle'.format(os.getcwd()), 'rb'))
    else:
        # get info_df
        info_df = generate_info_df_for_samples()

        # get fig_info_df
        fig_info_df = generate_fig_indo_df(info_df)

        # get the sample_abundance_df
        sample_abundance_df = generate_seq_abundance_df(fig_info_df)

        #MP sample by sample to do the MED
        input_q = Queue()

        # put in a tup of series. The first series being the info, the second being the current raw seq abundance df series
        for ind in fig_info_df.index.values.tolist():
            input_q.put((fig_info_df.loc[ind], sample_abundance_df.loc[ind]))

        numProc = 20
        for n in range(numProc):
            input_q.put('STOP')

        all_procs = []
        for n in range(numProc):
            p = Process(target=MED_worker, args=(input_q,))
            all_procs.append(p)
            p.start()

        for p in all_procs:
            p.join()

        apples = 'asdf'
        # here we have pickled out dictionaries of the MED nodes and their abundnaaces in each of the sample direcotires
        # these abundances are relative abundances
        # now have a master default dict and collect the abundances of the med nodes
        master_abund_med_dict = defaultdict(float)
        for ind in fig_info_df.index.values.tolist():
            print('adding med nodes to master dict for {}'.format(ind))
            sample_dir = fig_info_df.loc[ind, 'sample_dir']

            # read in the pickled abundance dict
            sample_med_abundance_dict = pickle.load(open('{}/sample_node_dict_rel_abund.pickle'.format(sample_dir), 'rb'))
            for seq_key, seq_rel_abund in sample_med_abundance_dict.items():
                master_abund_med_dict[seq_key] += seq_rel_abund

        # here we have the master dict updated with all of the med node sequences
        # sort and get list of med_node_seqs
        sorted_med_node_seqs_tups = sorted(master_abund_med_dict.items(), key=lambda x:x[1], reverse=True)
        med_node_seqs_sorted = [a[0] for a in sorted_med_node_seqs_tups]

        # now create a df that will hold all of the sample abundances
        minor_intra_med_abund_df = pd.DataFrame(columns=med_node_seqs_sorted)

        # now go back through the samples again and populate the df by making a temp series
        for ind in fig_info_df.index.values.tolist():
            print('populating minor_intra_med_abund_df with {}'.format(ind))
            sample_dir = fig_info_df.loc[ind, 'sample_dir']
            sample_med_abundance_dict = pickle.load(open('{}/sample_node_dict_rel_abund.pickle'.format(sample_dir), 'rb'))
            sample_data_in_order = []
            for seq in med_node_seqs_sorted:
                if seq in sample_med_abundance_dict.keys():
                    sample_data_in_order.append(sample_med_abundance_dict[seq])
                else:
                    sample_data_in_order.append(float(0))

            minor_intra_med_abund_df = minor_intra_med_abund_df.append(pd.Series(name=ind, data=sample_data_in_order, index=med_node_seqs_sorted))

        # here we have the minor_intra_med_abund_df populated and we can pickle it out
        pickle.dump(minor_intra_med_abund_df, open('{}/minor_intra_med_abund_df.pickle'.format(os.getcwd()), 'wb'))
    return minor_intra_med_abund_df

def MED_worker(input_q):
    for info_series, abund_series in iter(input_q.get, 'STOP'):
        # check to see if the MED has already been completed for the sample. If so then simply pass onto next sample
        print('MED for {}'.format(abund_series.name))
        sample_dir = info_series['sample_dir']
        if os.path.isfile('{}/sample_node_dict_rel_abund.pickle'.format(sample_dir)):
            continue

        sample_name = abund_series.name
        # now create a fasta for writing out
        sample_fasta_for_MED = []
        # https://pandas.pydata.org/pandas-docs/stable/generated/pandas.Series.nonzero.html
        non_zero_series = abund_series[abund_series.nonzero()[0]]

        seq_list_in_order = []
        relabund_list_in_order = []

        # we will skip the first sequence which is the most abundant
        for seq in non_zero_series.index.values.tolist()[1:]:
            rel_abund_div = abund_series.loc[seq]
            seq_list_in_order.append(seq)
            relabund_list_in_order.append(rel_abund_div)

        # now re normalise the abundances
        seq_counter = 0
        tot = sum(relabund_list_in_order)
        normalised_seq_abund_list = []
        for i in range(len(seq_list_in_order)):
            normalised_seq_abund = relabund_list_in_order[i] / tot
            normalised_seq_abund_list.append(normalised_seq_abund)
            absolute_abund_div = normalised_seq_abund * 1000
            if int(absolute_abund_div) > 0:
                for j in range(int(absolute_abund_div)):
                    sample_fasta_for_MED.append('>raw_seq_{}'.format(seq_counter))
                    sample_fasta_for_MED.append(seq_list_in_order[i])
                    seq_counter += 1
        check_tot = sum(normalised_seq_abund_list)

        # here we have a fasta populated that we will use for running the MED
        # now write it out, pad it, MED it
        path_to_fasta_for_med = '{}/fasta_for_med.fasta'.format(sample_dir)
        with open(path_to_fasta_for_med, 'w') as f:
            for line in sample_fasta_for_MED:
                f.write('{}\n'.format(line))

        subprocess.run([r'o-pad-with-gaps', r'{}'.format(path_to_fasta_for_med)])

        # Now run MED
        listOfFiles = []
        for (dirpath, dirnames, filenames) in os.walk(sample_dir):
            listOfFiles.extend(filenames)
            break
        for file in listOfFiles:
            if 'PADDED' in file:
                pathToFile = '{}/{}'.format(sample_dir, file)
                break
        MEDOutDir = '{}/{}/'.format(sample_dir, 'MEDOUT')
        os.makedirs(MEDOutDir, exist_ok=True)
        sys.stdout.write('{}: running MED\n'.format(sample_name))
        # Here we need to make sure that the M value is defined dynamically
        # the M value is a cutoff that looks at the abundance of the most abundant unique sequence in a node
        # if the abundance is lower than M then the node is discarded
        # we have been working recently with an M that equivaltes to 0.4% of 0.004. This was
        # calculated when working with a modelling project where I was subsampling to 1000 sequences. In this
        # scenario the M was set to 4.
        # We should also take care that M doesn't go below 4, so we should use a max choice for the M
        M_value = max(4, int(0.004 * (len(sample_fasta_for_MED) / 2)))
        subprocess.run(
            [r'decompose', '-M', str(M_value), '--skip-gexf-files', '--skip-gen-figures', '--skip-gen-html',
             '--skip-check-input', '-o', MEDOutDir, pathToFile])
        sys.stdout.write('{}: MED complete\n'.format(sample_name))

        # here the MED has been conducted and we can now read in the nodes and create a dict for the sample
        # read in the node file
        with open('{}/NODE-REPRESENTATIVES.fasta'.format(MEDOutDir), 'r') as f:
            node_file = [line.rstrip() for line in f]

        # get
        sample_node_dict = {}
        for i in range(0, len(node_file), 2):
            sample_node_dict[node_file[i+1]] = int(node_file[i].split(':')[1])

        # now do the rel abunds again and then pickle out
        tot = sum(sample_node_dict.values())
        sample_node_dict_rel_abund = {k: v/tot for k,v in sample_node_dict.items()}

        pickle.dump(sample_node_dict_rel_abund, open('{}/sample_node_dict_rel_abund.pickle'.format(sample_dir), 'wb'))



    return



# Code for generating the figure that will let us examine what our profile look like.
# this first go at this will not take into account any generated type profiles but rather just plot the sequences
# in each sample.



def generate_fig(plot_type):
    # this is the sample order from the its2 work
    sample_order = pickle.load(open('ordered_sample_names_from_its2_work.pickle', 'rb'))
    # to do this we'll have to read this in from the its2 work.
    # plot_type should be either 'full', 'low', or 'med'
    # full means all of the sequences including the maj
    # low means without the maj
    # med means without the maj and having been through med
    # qc_taxa_rel_abund means plot the relative abundance of the post qc taxa categories
    # qc_absolute means plot the absolute abundance of the post-qc sequences (all tax categories)
    info_df = generate_info_df_for_samples()

    fig_info_df = generate_fig_indo_df(info_df)

    # here we have the fig_info_df generated. We can use this for making the figure
    # SETUP AXES
    ax_list, fig = setup_axes()
    # here we have the axes set up.
    # before we move onto the actual plotting of the samples we should create an sp_output_df_div equivalent
    # for the samples. This will mean going through each of the coral sample directories and collecting relative
    # abundances. We will do this using the genus specific fastas.
    if plot_type == 'full' or plot_type == 'low':
        sample_abundance_df = generate_seq_abundance_df(fig_info_df)
        colour_dict = generate_colour_dict(sample_abundance_df=sample_abundance_df, is_med=False, is_qc=False)
    elif plot_type == 'med':
        sample_abundance_df = create_MED_node_sample_abundance_df_for_minor_intras()
        colour_dict = generate_colour_dict(sample_abundance_df=sample_abundance_df, is_med=True, is_qc=False)
    elif 'qc' in plot_type:
        sample_abundance_df = generate_post_qc_taxa_df()
        colour_dict = generate_colour_dict(sample_abundance_df=sample_abundance_df, is_med=True, is_qc=True)

    if plot_type == 'low':
        plot_data_axes_18s(ax_list, colour_dict, fig_info_df, sample_abundance_df, sample_order, minor_DIV=True)
    elif plot_type == 'full' or plot_type == 'med':
        plot_data_axes_18s(ax_list, colour_dict, fig_info_df, sample_abundance_df, sample_order, minor_DIV=False)
    elif 'qc' in plot_type:
        plot_data_axes_18s(ax_list, colour_dict, fig_info_df, sample_abundance_df, sample_order, minor_DIV=False,
                       qc=True, plot_type=plot_type)

    add_labels(ax_list)

    if plot_type == 'full':
        plt.savefig('{}/raw_seqs_abund_stacked.png'.format(os.getcwd()))
        plt.savefig('{}/raw_seqs_abund_stacked.svg'.format(os.getcwd()))
    elif plot_type == 'low':
        plt.savefig('{}/raw_seqs_abund_minor_div_only_stacked.png'.format(os.getcwd()))
        plt.savefig('{}/raw_seqs_abund_minor_div_only_stacked.svg'.format(os.getcwd()))
    elif plot_type == 'med':
        plt.savefig('{}/raw_seqs_abund_minor_div_only_stacked_with_med.png'.format(os.getcwd()))
        plt.savefig('{}/raw_seqs_abund_minor_div_only_stacked_with_med.svg'.format(os.getcwd()))
    elif plot_type == 'qc_taxa_rel_abund':
        plt.savefig('{}/post_qc_taxa_rel_abund.png'.format(os.getcwd()))
        plt.savefig('{}/post_qc_taxa_rel_abund.svg'.format(os.getcwd()))
    elif plot_type == 'qc_absolute':
        plt.savefig('{}/post_qc_absolute.png'.format(os.getcwd()))
        plt.savefig('{}/post_qc_absolute.svg'.format(os.getcwd()))



    return

# the first step for this is to generate an info df.
# we can generate this from the current info_df
# we should call it fig_info_df
def generate_fig_indo_df(info_df):
    if os.path.isfile('{}/fig_info_df.pickle'.format(os.getcwd())):
        return pickle.load(open('{}/fig_info_df.pickle'.format(os.getcwd()), 'rb'))
    else:
        fig_info_df = pd.DataFrame(columns=['island', 'site', 'genus', 'individual', 'sample_dir'])
        for ind in info_df.index.values.tolist():
            if 'CORAL' in info_df.loc[ind, 'fastq_fwd_file_path']:
                fastq_string = info_df.loc[ind, 'fastq_fwd_file_path']
                components = fastq_string.split('/')
                smp_site = components[-6]
                smp_island = components[-7]
                smp_individual = components[-3]
                smp_genus = components[-4]
                smp_dir = '/'.join(components[:-1])
                fig_info_df = fig_info_df.append(
                    pd.Series(name=ind, data=[smp_island, smp_site, smp_genus, smp_individual, smp_dir],
                              index=['island', 'site', 'genus', 'individual', 'sample_dir']))
        pickle.dump(fig_info_df, open('{}/fig_info_df.pickle'.format(os.getcwd()), 'wb'))
        return fig_info_df

def setup_axes():
    # https://matplotlib.org/users/gridspec.html
    fig = plt.figure(figsize=(14, 8))
    # the bottom row will be for the legend
    # the second to last will just be invisible to give a space between the legend and the other plots
    # we also want to include a gridspec plot after each of the main three. These will hold the csw and surface
    # samples
    gs = plt.GridSpec(5, 3, figure=fig, height_ratios=[1, 0.3, 1, 0.3, 1])
    # within each of the GrdiSpec subplots we will make a subplotspec which is three plots on one row
    ax_list = []
    grid_spec_subplot_list = []
    increaser = 0
    # make the main 3x3 axis
    for row_ind in range(3):
        for col_ind in range(3):
            # put in the main data 3 plots
            temp_grid_spec_subplot = gridspec.GridSpecFromSubplotSpec(1, 3,
                                                                      subplot_spec=gs[row_ind + increaser, col_ind])
            grid_spec_subplot_list.append(temp_grid_spec_subplot)
            for i in range(3):
                # NB this might be a 2d array, lets see.
                ax = plt.Subplot(fig, temp_grid_spec_subplot[i])
                ax_list.append(ax)
                fig.add_subplot(ax)
        # now put in the spacer row
        if increaser < 2:
            ax_space = plt.subplot(gs[row_ind + 1 + increaser, :])
            remove_axes_but_allow_labels(ax_space)
            increaser += 1
    return ax_list, fig


def plot_data_axes_18s(ax_list, colour_dict, fig_info_df, sample_abundance_df, sample_order, minor_DIV=False, qc=False, plot_type=None):

    ax_count = 0
    for site in ['SITE01', 'SITE02', 'SITE03']:
        for location in ['ISLAND06', 'ISLAND10', 'ISLAND15']:
            for spp in ['PORITES', 'POCILLOPORA', 'MILLEPORA']:
                ax = ax_list[ax_count]
                patches_list = []
                ind = 0
                colour_list = []

                # for each set of location, site and spp, we basically want to get a list of the samples
                # that meet the set criteria, we then want to plot samples according to the ordered_sample_list
                # order which will be in IDs. As such we will have to convert the sample_name in the info_df
                # to a sample ID using the smp_name_to_smp_id_dict.

                # get sample_names that fit the requirements
                sample_names_of_set = fig_info_df.loc[
                    (fig_info_df['island'] == location) &
                    (fig_info_df['site'] == site) &
                    (fig_info_df['genus'] == spp)
                    ].index.values.tolist()

                # get the above sample_names_of_set in order of the sample_order list
                ordered_sample_names_of_set = []
                for samp_name in sample_order:
                    if samp_name in sample_names_of_set:
                        ordered_sample_names_of_set.append(samp_name)

                # there may be some samples in the sample_names_of_set that aren't in the sample_order
                # add these here
                for sample_name in sample_names_of_set:
                    if sample_name not in sample_order:
                        ordered_sample_names_of_set.append(sample_name)


                num_smp_in_this_subplot = len(sample_names_of_set)
                x_tick_label_list = []

                for smple_id_to_plot in ordered_sample_names_of_set:
                    # General plotting
                    sys.stdout.write('\rPlotting sample: {}'.format(smple_id_to_plot))
                    x_tick_label_list.append(smple_id_to_plot)
                    # for each sample we will start at 0 for the y and then add the height of each bar to this

                    # PLOT DIVs
                    if minor_DIV:
                        plot_div_over_type_minor_div_only(colour_dict, colour_list, ind, patches_list, smple_id_to_plot, sample_abundance_df)
                    elif qc:
                        plot_div_over_type_qc_abund(colour_dict=colour_dict, colour_list=colour_list, ind=ind,
                                                    patches_list=patches_list, smple_id_to_plot=smple_id_to_plot,
                                                    sample_abundance_df=sample_abundance_df, plot_type=plot_type)
                    else:
                        plot_div_over_type(colour_dict, colour_list, ind, patches_list, smple_id_to_plot, sample_abundance_df)


                    ind += 1
                if qc:
                    paint_rect_to_axes_div_and_type(ax=ax, colour_list=colour_list,
                                                    num_smp_in_this_subplot=num_smp_in_this_subplot,
                                                    patches_list=patches_list,
                                                    x_tick_label_list=x_tick_label_list,
                                                    max_num_smpls_in_subplot=10, qc=True, plot_type=plot_type,
                                                    ax_count=ax_count)
                else:
                    paint_rect_to_axes_div_and_type(ax=ax, colour_list=colour_list,
                                                    num_smp_in_this_subplot=num_smp_in_this_subplot,
                                                    patches_list=patches_list,
                                                    x_tick_label_list=x_tick_label_list,
                                                    max_num_smpls_in_subplot=10)

                ax_count += 1


def add_labels(ax_list):
    ax_list[1].set_title('ISLAND06')
    ax_list[4].set_title('ISLAND10')
    ax_list[7].set_title('ISLAND15')

    ax_list[0].set_ylabel('SITE 1', fontsize='x-large')
    ax_list[9].set_ylabel('SITE 2', fontsize='x-large')
    ax_list[18].set_ylabel('SITE 3', fontsize='x-large')

    ax_list[18].set_xlabel('porites', fontsize='medium')
    ax_list[19].set_xlabel('pocillopora', fontsize='medium')
    ax_list[20].set_xlabel('millepora', fontsize='medium')

    ax_list[21].set_xlabel('porites', fontsize='medium')
    ax_list[22].set_xlabel('pocillopora', fontsize='medium')
    ax_list[23].set_xlabel('millepora', fontsize='medium')

    ax_list[24].set_xlabel('porites', fontsize='medium')
    ax_list[25].set_xlabel('pocillopora', fontsize='medium')
    ax_list[26].set_xlabel('millepora', fontsize='medium')


def paint_rect_to_axes_div_and_type(ax, colour_list, num_smp_in_this_subplot,  patches_list, x_tick_label_list=None,
                                    max_num_smpls_in_subplot=10, qc=False, plot_type=None, ax_count=None):
    # We can try making a custom colour map
    # https://matplotlib.org/api/_as_gen/matplotlib.colors.ListedColormap.html
    this_cmap = ListedColormap(colour_list)
    # here we should have a list of Rectangle patches
    # now create the PatchCollection object from the patches_list
    patches_collection = PatchCollection(patches_list, cmap=this_cmap)
    patches_collection.set_array(np.arange(len(patches_list)))
    # if n_subplots is only 1 then we can refer directly to the axarr object
    # else we will need ot reference the correct set of axes with i
    # Add the pathces to the axes
    ax.add_collection(patches_collection)
    ax.autoscale_view()
    ax.figure.canvas.draw()
    # also format the axes.
    # make it so that the x axes is constant length
    ax.set_xlim(0 - 0.5, max_num_smpls_in_subplot - 0.5)
    if plot_type != 'qc_absolute':
        ax.set_ylim(0, 1)
    ax.set_xticks(range(num_smp_in_this_subplot))
    ax.set_xticklabels(x_tick_label_list, rotation='vertical', fontsize=6)

    if plot_type == 'qc_absolute':
        ax.set_ylim(0, 2000000)
        ax.set_yscale('symlog')
        if ax_count % 9 != 0:
            remove_axes_but_allow_labels(ax, x_tick_label_list)
        else:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)



    else:
        remove_axes_but_allow_labels(ax, x_tick_label_list)

    # as well as getting rid of the top and right axis splines
    # I'd also like to restrict the bottom spine to where there are samples plotted but also
    # maintain the width of the samples
    # I think the easiest way to do this is to hack a bit by setting the x axis spines to invisible
    # and then drawing on a line at y = 0 between the smallest and largest ind (+- 0.5)
    # ax.spines['bottom'].set_visible(False)
    if plot_type != 'qc_absolute':
        ax.add_line(Line2D((0 - 0.5, num_smp_in_this_subplot - 0.5), (0, 0), linewidth=2, color='black'))
    elif plot_type == 'qc_absolute':
        ax.add_line(Line2D((0 - 0.5, num_smp_in_this_subplot - 0.5), (0.1, 0.1), linewidth=2, color='black'))


def plot_div_over_type(colour_dict, colour_list, ind, patches_list, smple_id_to_plot, sample_abundance_df):
    bottom_div = 0
    # for each sequence, create a rect patch
    # the rect will be 1 in width and centered about the ind value.

    sample_series = sample_abundance_df.loc[smple_id_to_plot]
    # https://pandas.pydata.org/pandas-docs/stable/generated/pandas.Series.nonzero.html
    non_zero_series = sample_series[sample_series.nonzero()[0]]

    for seq in non_zero_series.index.values.tolist():
        # class matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
        rel_abund_div = sample_abundance_df.loc[smple_id_to_plot, seq]
        if rel_abund_div > 0:
            patches_list.append(Rectangle((ind - 0.5, bottom_div), 1, rel_abund_div, color=colour_dict[seq]))
            # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
            colour_list.append(colour_dict[seq])
            bottom_div += rel_abund_div

def plot_div_over_type_minor_div_only(colour_dict, colour_list, ind, patches_list, smple_id_to_plot, sample_abundance_df):
    # so it appears that there is a predominant sequence that occupies about 93% of each sample
    # so the intragenomic diverstiy is very compressed in the top of the plots. To expand this a bit so that we can
    # look at the intragenomic, I will leave out the predomiant sequence
    bottom_div = 0
    # for each sequence, create a rect patch
    # the rect will be 1 in width and centered about the ind value.
    # have a list that holds the seq and another that holds the abundance in order.
    # that way we can re normalise the list
    seq_list_in_order = []
    relabund_list_in_order = []
    # we will skip the first three sequences which were the most abundant
    # this is relatively slow and I think we can use the nonzero function of a series to speed this up
    sample_series = sample_abundance_df.loc[smple_id_to_plot]
    # https://pandas.pydata.org/pandas-docs/stable/generated/pandas.Series.nonzero.html
    non_zero_series = sample_series[sample_series.nonzero()[0]]

    for seq in non_zero_series.index.values.tolist()[1:]:
        rel_abund_div = sample_abundance_df.loc[smple_id_to_plot, seq]
        if rel_abund_div > 0 and rel_abund_div < 0.5:
            seq_list_in_order.append(seq)
            relabund_list_in_order.append(rel_abund_div)
    # now re normalise the abundances
    tot = sum(relabund_list_in_order)
    for i in range(len(seq_list_in_order)):
        normalised_seq_abund = relabund_list_in_order[i]/tot
        patches_list.append(Rectangle((ind - 0.5, bottom_div), 1, normalised_seq_abund, color=colour_dict[seq_list_in_order[i]]))
        # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
        colour_list.append(colour_dict[seq_list_in_order[i]])
        bottom_div += normalised_seq_abund


def plot_div_over_type_qc_abund(colour_dict, colour_list, ind, patches_list,
                                smple_id_to_plot, sample_abundance_df, plot_type=None):

    bottom_div = 0
    # for each sequence, create a rect patch
    # the rect will be 1 in width and centered about the ind value.

    sample_series = sample_abundance_df.loc[smple_id_to_plot]
    # https://pandas.pydata.org/pandas-docs/stable/generated/pandas.Series.nonzero.html
    non_zero_series = sample_series[sample_series.nonzero()[0]]

    if plot_type == 'qc_absolute':
        for tax_category in non_zero_series.index.values.tolist():
            # class matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
            abs_abund = sample_abundance_df.loc[smple_id_to_plot, tax_category]
            patches_list.append(Rectangle((ind - 0.5, bottom_div), 1, abs_abund, color='#808080'))
            colour_list.append('#808080')
            bottom_div += abs_abund

    elif plot_type == 'qc_taxa_rel_abund':
        cat_abund_list = []
        cat_name_list = []
        for tax_category in non_zero_series.index.values.tolist():
            cat_abs_abund = sample_abundance_df.loc[smple_id_to_plot, tax_category]
            cat_abund_list.append(cat_abs_abund)
            cat_name_list.append(tax_category)
        tot = sum(cat_abund_list)
        rel_abunds_list = [abund/tot for abund in cat_abund_list]
        check = sum(rel_abunds_list)
        for i, cat_name in enumerate(cat_name_list):
            patches_list.append(Rectangle((ind - 0.5, bottom_div), 1, rel_abunds_list[i], color=colour_dict[cat_name]))
            colour_list.append(colour_dict[cat_name])
            bottom_div += rel_abunds_list[i]



def generate_colour_dict(sample_abundance_df, is_med=False, is_qc=False):
    # the purpose of this is to return a seuence to colour dictionary that we will pickle out to maintain
    # continuity
    # make the colour list and the grey list
    # then get the sequence order from the sample_abundance_df
    # then associate to the colours until we are out of colours
    colour_dict = {}
    if is_qc:
        return {'Porites':'#FFFF00', 'Pocillopora':'#87CEFA', 'Millepora':'#FF6347',
         'other_coral':'#C0C0C0', 'Symbiodiniaceae':'#00FF00', 'other_taxa':'#696969'}

    if is_med:
        colour_list = get_colour_list()[3:]
    else:
        colour_list = get_colour_list()
    grey_palette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']
    for i, seq_name in enumerate(list(sample_abundance_df)):
        if i < len(colour_list):
            colour_dict[seq_name] = colour_list[i]
        else:
            colour_dict[seq_name] = grey_palette[i%6]
    return colour_dict


def generate_seq_abundance_df(fig_info_df, numProc=20):
    # the purpose of this will be to get a dataframe that contains the sequence abundances for each of the samples
    # this will only contain the sequence abundances for the specific genus coral sequences
    # so it will not contain any zooxs seqs or non-scleractinian zooxs.
    # One thing that will be tough is getting a global collection of the sequences that have been found.
    # I think the best way to do this is to MP this and return a dictionary for each sample
    # where the key will be the sample name
    # and the vlaue will be a dictionary which will be key of actual sequences and relabund of the seq
    # once we have all of these we can get a gloab diversity of sequence and we can also get a glbal abundance
    # of the seuqencs so that we can plot the most abundant sequencs first.

    # this df only holds info for the coral samples so there is no need to filter

    if os.path.isfile('{}/sample_abundance_df.pkl'.format(os.getcwd())):
        sample_abundance_df = pd.read_pickle('{}/sample_abundance_df.pkl'.format(os.getcwd()))
    else:
        input_q = Queue()

        for ind in fig_info_df.index.values.tolist():
            input_q.put((ind, fig_info_df.loc[ind, 'sample_dir']))

        for n in range(numProc):
            input_q.put('STOP')


        all_procs = []
        for n in range(numProc):
            p = Process(target=abundance_worker, args=(input_q, ))
            all_procs.append(p)
            p.start()

        for p in all_procs:
            p.join()

        # at this stage we have pickled out abundance dictionaries for each of the samples in their
        # respective directories
        # now go through each of samples and check to see if each of the seqs is found in the master seq dict
        # or if it is a subset or superset of any of the strings. This has the potential to be very slow
        # so we could maybe chunk this and do it MP. But for the time being lets try to do it serial.
        # first we will go through each of the samples and generate a set of sequence for the master_abund_dict
        # that will be the largest superset sequences of all the sequences found in the samples
        # once we've done that, then we will go back through the sequences and see which sequences the sample's sequecnes
        # should be associated to. We should pickle out along the way to save us time if we have to re-run

        # just out of interest I want to see how many times we have a match list greater than 1
        if os.path.isfile('{}/master_abund_dict.pickle'.format(os.getcwd())):
            master_abund_dict = pickle.load(open('{}/master_abund_dict.pickle'.format(os.getcwd()), 'rb'))
        else:
            big_match_counter = 0
            total = 0
            master_abund_dict = defaultdict(float)
            for ind in fig_info_df.index.values.tolist():
                print('Counting {}'.format(ind))
                sample_dir = fig_info_df.loc[ind, 'sample_dir']
                sample_abund_dict = pickle.load(open('{}/seq_abund_dict.pickle'.format(sample_dir), 'rb'))
                for sample_key_sequence, v_sample in sample_abund_dict.items():
                    total += 1
                    # for each of the sequences, see if it is found within one of the master_abund_dicts
                    # check to see that it is not found in > than one of the sequences too. If it is then this is a problem
                    # but there may be a simple solution, i.e. consolidating the two sequences.

                    # list that holds the sequences that are either a subset or superset of the sample sequence
                    match_list = []
                    for master_key_sequence, v_master in master_abund_dict.items():

                        if sample_key_sequence in master_key_sequence or master_key_sequence in sample_key_sequence:
                            match_list.append(master_key_sequence)

                    # here we have checked each of the master sequences
                    if not match_list:
                        # if the match list is empty then this is a new sequence and we can add it to the master_abund_dict
                        master_abund_dict[sample_key_sequence] += v_sample
                    elif len(match_list) == 1:
                        # then we have a single match.
                        if sample_key_sequence in match_list[0]:
                            # if the sample is a subset of one of the master_abund_dicts then its
                            # abundance should be attributed to this sequence
                            master_abund_dict[match_list[0]] += v_sample
                        elif match_list[0] in sample_key_sequence:
                            # else if the sample sequences contains the master sequence, then the master sequence
                            # should be updated to the seq in question, the current abundance should be associated
                            # to this new sequence and the abundance of the seq in q should also be attributed (added)
                            temp_abundance = master_abund_dict[match_list[0]]
                            del master_abund_dict[match_list[0]]
                            master_abund_dict[sample_key_sequence] = temp_abundance
                            master_abund_dict[sample_key_sequence] += v_sample
                    elif len(match_list) > 1:
                        # then we have a bit of a problem here and we need to work out what to do from looking at the
                        # situation
                        # this is a particular situation and has happened when there is a sequence that is shorter than the
                        # majority and one of the sequences has a differntiating snp in the difference in sequence length.
                        # for the time being I will attribute the abundance of this sequence to the most abundant match
                        most_abund_seq = [a[0] for a in sorted([(seq, master_abund_dict[seq]) for seq in match_list], key=lambda x: x[1], reverse=True)][0]
                        master_abund_dict[most_abund_seq] += v_sample
                        big_match_counter += 1
                        continue
            # at this point we have the master_abund_dict populated and we should pickle it out
            pickle.dump(master_abund_dict, open('{}/master_abund_dict.pickle'.format(os.getcwd()), 'wb'))
        apples = 'asdf'

        # now we once again go sample by sample and then seq by seq and check all of the seqs
        # thar are in the master dict. Then when we find an association we use that to populate to the dataframe
        sorted_tups = sorted(list(master_abund_dict.items()), key=lambda x:x[1], reverse=True)
        header_seq_s = [a[0] for a in sorted_tups]
        sample_abundance_df = pd.DataFrame(columns=header_seq_s)

        total_matches_evaluated = 0
        multi_match_list = 0
        no_match = 0
        for ind in fig_info_df.index.values.tolist():

            # create an empty sequence that will be full of float 0s to start with
            sample_abundance_series = pd.Series(name=ind, data=[float(0) for i in header_seq_s], index=header_seq_s)
            print('Populating df with {}'.format(ind))
            sample_dir = fig_info_df.loc[ind, 'sample_dir']
            if os.path.isfile('{}/sample_abundance_series.pkl'.format(sample_dir)):
                sample_abundance_series = pd.read_pickle('{}/sample_abundance_series.pkl'.format(sample_dir))
            else:
                sample_abund_dict = pickle.load(open('{}/seq_abund_dict.pickle'.format(sample_dir), 'rb'))
                # for every sequence in the sample
                for sample_key_sequence, v_sample in sample_abund_dict.items():
                    total_matches_evaluated += 1
                    # list that holds the sequences that are either a subset or superset of the sample sequence
                    match_list = []
                    exact_match = False
                    for master_key_sequence, v_master in master_abund_dict.items():
                        if sample_key_sequence == master_key_sequence:
                        # then we have an exact match and this is the sequence that the abundance should be associated
                        # to.
                        # else continue to work with the match_list etc.
                            sample_abundance_series[master_key_sequence] += v_sample
                            exact_match = True
                            break
                        if sample_key_sequence in master_key_sequence or master_key_sequence in sample_key_sequence:
                            match_list.append(master_key_sequence)
                    if exact_match:
                        continue
                    # here we have not acheived an exact match and we should work with the list
                    # here we have checked each of the master sequences
                    if not match_list:
                        # this should never happen
                        no_match += 1
                        continue
                    elif len(match_list) == 1:
                        # then we have a single match and this is the column of the df that we should associate the
                        # abundance to
                        sample_abundance_series[match_list[0]] += v_sample

                    elif len(match_list) > 1:
                        # then we have a bit of a problem here and we need to work out what to do from looking at the
                        # situation
                        # this is a particular situation and has happened when there is a sequence that is shorter than the
                        # majority and one of the sequences has a differntiating snp in the difference in sequence length.
                        # for the time being I will attribute the abundance of this sequence to the most abundant match
                        multi_match_list += 1
                        sorted_tups_temp = sorted([(seq, master_abund_dict[seq]) for seq in match_list], key=lambda x: x[1],
                                                 reverse=True)
                        most_abund_seq = [a[0] for a in sorted_tups_temp][0]
                        sample_abundance_series[most_abund_seq] += v_sample
                # here we have the series for the sample populated.
                # now add this to the dataframe
                # pickle out the series to save time
                sample_abundance_series.to_pickle('{}/sample_abundance_series.pkl'.format(sample_dir))
            sample_abundance_df = sample_abundance_df.append(sample_abundance_series)
        # at this point we have fully populated the sample_abundance_df and we should pickle it out
        sample_abundance_df.to_pickle('{}/sample_abundance_df.pkl'.format(os.getcwd()))

    return sample_abundance_df

def abundance_worker(input_q):
    # within this method we will aim to create a dict and pickle it out so that it can be read in outside
    # of this MP.
    for sample_name, sample_dir in iter(input_q.get, 'STOP'):
        # check to see if the sample's abundance dictionary has already been created
        if os.path.isfile('{}/seq_abund_dict.pickle'.format(sample_dir)):
            continue

        # I want this dictionary to have key as the raw sequence as value
        sample_dict = defaultdict(int)

        sys.stdout.write('\nSample {}\n'.format(sample_dir))


        # for each sample we will only need two files
        # we will need the name file inorder to create an abundance dict
        # then we will need the coral_genus_only.fasta file so that we know which sequences we should be counting

        #read in the name file
        with open('{}/stability.trim.contigs.good.abund.pcr.names'.format(sample_dir), 'r') as f:
            name_file = [line.rstrip() for line in f]

        abund_dict = {line.split('\t')[0]:len(line.split('\t')[1].split(',')) for line in name_file}

        # now read in the coral_genus_only.fasta
        with open('{}/coral_genus_only.fasta'.format(sample_dir), 'r') as f:
            coral_genus_fasta = [line.rstrip() for line in f]

        # here we have the abund dict and the coral_genus fasta
        # now we can go sesquence by sequence through the fasta and look up the abund and add to the sample_dict
        # Becuase we didn't do a find seqs.unique it is quite possible that there are sequences that are the same
        # and so we should use a default dict to make the abundance dictionary for the sample
        for i in range(0, len(coral_genus_fasta), 2):
            sample_dict[coral_genus_fasta[i + 1]] += abund_dict[coral_genus_fasta[i][1:]]

        # now normalise the seq abundances
        tot_seqs = sum(sample_dict.values())
        normalised_sample_dict = {k: v/tot_seqs for k, v in sample_dict.items()}

        # here we have the samples abundance dictionary populated. Now pickle out
        pickle.dump(normalised_sample_dict, open('{}/seq_abund_dict.pickle'.format(sample_dir), 'wb'))


    return

def remove_axes_but_allow_labels(ax, x_tick_label_list=None):
    ax.set_frame_on(False)
    if not x_tick_label_list:
        ax.set_xticks([])
    ax.set_yticks([])
    ax.minorticks_off()


def generate_post_qc_taxa_df():

    if os.path.isfile('{}/tax_absolute_abund_df.pickle'.format(os.getcwd())):
        tax_absolute_abund_df = pickle.load(open('{}/tax_absolute_abund_df.pickle'.format(os.getcwd()), 'rb'))
    else:
        info_df = generate_info_df_for_samples()

        fig_info_df = generate_fig_indo_df(info_df)

        input_q = Queue()

        for ind in fig_info_df.index.values.tolist():
            input_q.put((ind, fig_info_df.loc[ind, 'sample_dir']))

        numProc = 20
        for n in range(numProc):
            input_q.put('STOP')

        all_procs = []
        for n in range(numProc):
            p = Process(target=qc_abundance_worker, args=(input_q,))
            all_procs.append(p)
            p.start()

        for p in all_procs:
            p.join()

        # at this point we have the sample_tax_category_dict items for each sample pickled out
        # we can now go through each of the sample_dir again and use them to populate the master df
        # which we will then return
        tax_columns = ['Porites', 'Millepora', 'Pocillopora', 'other_coral', 'Symbiodiniaceae', 'other_taxa']
        tax_absolute_abund_df = pd.DataFrame(columns=tax_columns)
        for ind in fig_info_df.index.values.tolist():
            sample_dir = fig_info_df.loc[ind, 'sample_dir']

            # load the sample_tax_category_dict
            sample_tax_category_dict = pickle.load(open('{}/sample_tax_category_dict.pickle'.format(sample_dir), 'rb'))

            sample_data = []
            for cat in tax_columns:
                if cat in sample_tax_category_dict.keys():
                    sample_data.append(sample_tax_category_dict[cat])
                else:
                    sample_data.append(0)
            tax_absolute_abund_df = tax_absolute_abund_df.append(pd.Series(data=sample_data, index=tax_columns, name=ind))

        # here we have the tax_absolute_abund_df populated
        # now pickle it out
        pickle.dump(tax_absolute_abund_df, open('{}/tax_absolute_abund_df.pickle'.format(os.getcwd()), 'wb'))
    return tax_absolute_abund_df

def qc_abundance_worker(input_q):
    # the aim of this worker will be to pickle out a abund_tax_dict that will have the keys
    # pocillopora, millepora, porites, other_host, symbiodiniacea, other_taxa
    for sample_name, sample_dir in iter(input_q.get, 'STOP'):
        if os.path.isfile('{}/sample_tax_category_dict.pickle'.format(sample_dir)):
            continue
        sys.stdout.write('\nSample {}\n'.format(sample_name))
        # read in names file
        with open('{}/stability.trim.contigs.good.abund.pcr.names'.format(sample_dir), 'r') as f:
            name_file = [line.rstrip() for line in f]

        # from this we can get the other taxa as any taxa that's value isn't 'Scleractinia' or 'Symbiodiniacea'
        sample_tax_dict = pickle.load(open('{}/sample_tax_dict.pickle'.format(sample_dir), 'rb'))
        # from this the value is the genus. So we can search for 'Porites', 'Millepora', 'Pocillopora'
        coral_dict = pickle.load(open('{}/coral_dict.pickle'.format(sample_dir), 'rb'))
        # each of the sequence names in this is symbiodiniaceae so we can simply count these
        symbiodiniaceae_dict = pickle.load(open('{}/symbiodiniaceae_dict.pickle'.format(sample_dir), 'rb'))

        sample_tax_category_dict = defaultdict(int)
        # go through the name file line by line identifying what it belongs to
        for line in name_file:
            abundance = len(line.split('\t')[1].split(','))
            seq_name = line.split('\t')[0]
            if seq_name in coral_dict.keys():
                genus_name = coral_dict[seq_name]
                if genus_name in ['Porites', 'Millepora', 'Pocillopora']:
                    sample_tax_category_dict[genus_name] += abundance
                else:
                    # then this another coral
                    sample_tax_category_dict['other_coral'] += abundance
            elif seq_name in symbiodiniaceae_dict.keys():
                sample_tax_category_dict['Symbiodiniaceae'] += abundance
            elif seq_name in sample_tax_dict.keys():
                sample_tax_category_dict['other_taxa'] += abundance

        # here we have the sample_tax_category_dict populated and we can now pickle it out
        pickle.dump(sample_tax_category_dict, open('{}/sample_tax_category_dict.pickle'.format(sample_dir), 'wb'))





# For each sample I will generate a scleractinian sequences fasta and a fasta that contains only the sequencs
# for the coral that the sample is supposed to be, i.e. porites, pocillopora or millepora
# From here we should then be able to make the figures as we have before for the ITS2
def generate_sclerac_fastas(numProc=20):
    # we are only interested in the coral samples.
    # as such we should filter these out when doing the MPing.
    info_df = generate_info_df_for_samples()

    input_q_for_fasta_gen = Queue()

    for ind in info_df.index.values.tolist():
        if 'CORAL' in info_df.loc[ind, 'fastq_fwd_file_path']:
            input_q_for_fasta_gen.put(info_df.loc[ind, 'fastq_fwd_file_path'])

    for n in range(numProc):
        input_q_for_fasta_gen.put('STOP')


    all_procs = []
    for n in range(numProc):
        p = Process(target=scleractinia_fasta_worker, args=(input_q_for_fasta_gen, ))
        all_procs.append(p)
        p.start()

    for p in all_procs:
        p.join()

    return


def scleractinia_fasta_worker(input_q):
    # here we want to read in the pickled out dicts and lists
    # we will also need to read in the fasta that holds the actual sequences
    for path in iter(input_q.get, 'STOP'):
        sys.stdout.write('\nSample {}\n'.format(path.split('/')[-1]))
        genus_of_sample = path.split('/')[-4]

        sample_dir = '/'.join(path.split('/')[:-1])


        # for the time being, while we still have the blast worker going we should check to see if
        # the sample has the dictionaries already completed.
        if not os.path.isfile('{}/sample_tax_dict.pickle'.format(sample_dir)) or not os.path.isfile('{}/coral_dict.pickle'.format(sample_dir)) or not os.path.isfile('{}/symbiodiniaceae_dict.pickle'.format(sample_dir)):
            continue

        # next check to see if we have already generated the fastas for this sample.
        if os.path.isfile('{}/coral.fasta'.format(sample_dir)) and os.path.isfile('{}/coral_genus_only.fasta'.format(sample_dir)):
            continue

        coral_dict = pickle.load(open('{}/coral_dict.pickle'.format(sample_dir), 'rb'))

        # need to read in the fasta that holds the sequences and create a dict
        with open('{}/stability.trim.contigs.good.unique.abund.pcr.fasta'.format(sample_dir), 'r') as f:
            fasta_file = [line.rstrip() for line in f]
        fasta_dict = {fasta_file[i][1:].split('\t')[0] :fasta_file[i+1] for i in range(0, len(fasta_file),2)}

        coral_fasta = []
        genus_fasta = []
        for seq_name, seq_genus in coral_dict.items():
            coral_fasta.append('>{}'.format(seq_name))
            coral_fasta.append(fasta_dict[seq_name])
            if seq_genus.lower() == genus_of_sample.lower():
                genus_fasta.append('>{}'.format(seq_name))
                genus_fasta.append(fasta_dict[seq_name])

        # here we have created both of the fastas and we can now write them out
        with open('{}/coral.fasta'.format(sample_dir), 'w') as f:
            for line in coral_fasta:
                f.write('{}\n'.format(line))

        # here we have created both of the fastas and we can now write them out
        with open('{}/coral_genus_only.fasta'.format(sample_dir), 'w') as f:
            for line in genus_fasta:
                f.write('{}\n'.format(line))

    return


def generate_summary_figure():
    # I think it would be good to generate a summary figure the same as we did for the
    # tara data
    # I will try to recycle the code as best I can
    # would be good to have: raw contigs, post-qc, scleractinian and non-scleractinian
    # However, to work with the same code as we had with the tara data we will need to already have the
    # quantitative stats collected. We will have to use the directories listed
    # in the info_df and look for the summary files for each of the steps and extract what we need from them.

    info_df = generate_info_df_for_samples()

    # a dict that will hold the sample_name and the values for each of the qc levels
    for ind in info_df.index.values.tolist():
        # for each of the directories, we should

        # run a summary.seqs on the stability.trim.contigs.good.unique.abund.pcr.fasta

        # once we have run the summary, read in this summary which is the post-qc
        # also read in the summary of the contigs for the number of contigs pre-qc
        # then we can read in the sclearactinian fasta to get the number of

        # read in the dictionaries to get the taxa counts for both coral and non-coral reads
        apples = 'asdf'

    # here we can add the new data as columns to the info_df and pickle this out so that it is automatically
    # loaded in future

    # then we should be able make the figure using the tara code.
    return


# the sequence_QC method has taken care of the basic mothur qc for us. The next step will be to run a blast
# analysis for each of the samples that we have
# In each directory we should have a stability.trim.contigs.good.unique.abund.pcr.fasta paired with a
# stability.trim.contigs.good.abund.pcr.names. We should work with these pairs and produce a dictionary
# that will give sequence name to taxonomic id.

def generate_tax_dictionaries(numProc=20):
    # we are only interested in the coral samples.
    # as such we should filter these out when doing the MPing.
    info_df = generate_info_df_for_samples()

    ncbircFile = []
    db_path = '/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/ntdbdownload'
    ncbircFile.extend(["[BLAST]", "BLASTDB={}".format(db_path)])

    # write out the ncbircFile
    with open('{}/.ncbirc'.format(os.getcwd()), 'w') as f:
        for line in ncbircFile:
            f.write('{}\n'.format(line))


    input_q_for_blast = Queue()

    for ind in info_df.index.values.tolist():
        if 'CORAL' in info_df.loc[ind, 'fastq_fwd_file_path']:
            input_q_for_blast.put(info_df.loc[ind, 'fastq_fwd_file_path'])

    for n in range(numProc):
        input_q_for_blast.put('STOP')

    # generate the dictionaries used for taxa_designation
    # send a copy to each of the processes
    node_dict, name_dict = generate_taxa_designation_from_staxid_dicts()

    all_procs = []
    for n in range(numProc):
        p = Process(target=blast_worker, args=(input_q_for_blast, node_dict, name_dict))
        all_procs.append(p)
        p.start()

    for p in all_procs:
        p.join()

    return

def generate_taxa_designation_from_staxid_dicts(taxdump_dir = '/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/ntdbdownload/taxdump'):
    # read in the .nodes file. This file tells us which tax level the node is and which node is the parent level
    with open('{}/nodes.dmp'.format(taxdump_dir), 'r') as f:
        node_file = [line.rstrip() for line in f]
    # now make a dict from this where key is the tax id and the value is a tup where 0 = parent 1 = tax level
    node_dict = {line.split('\t|\t')[0]: (line.split('\t|\t')[1], line.split('\t|\t')[2]) for line in node_file}

    # next read in the names file. This file hold the name of the node.
    with open('{}/names.dmp'.format(taxdump_dir), 'r') as f:
        name_file = [line.rstrip() for line in f]

    # now make a dict from the names file where the key is the staxid and the value is the name
    name_dict = {line.split('\t|\t')[0]: line.split('\t|\t')[1] for line in name_file if line.split('\t|\t')[3].replace('\t|', '') == 'scientific name'}

    return node_dict, name_dict

def get_taxa_designation_from_staxid(staxid, tax_level_list, node_dict=None, name_dict=None, taxdump_dir = '/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/ntdbdownload/taxdump'):
    #I have set this up so that you can either feed in the already made node_dict and name_dict or
    # you they can be generated by this method. If you are running this method many times it will make
    # sense to run the above generate_taxa_designation_from_staxid_dicts methods
    # to generate the dicts to feed into this
    # This will take in a list as its argument and find the levels that are listed in the list
    list_to_return = [False for i in tax_level_list]
    if node_dict==None or name_dict==None:
        node_dict, name_dict = generate_taxa_designation_from_staxid_dicts()

    # now we can do the searching
    # go through staxid nodes until we get to the tax_level required
    # then grab the name
    current_staxid = staxid
    while True:
        if current_staxid == '1':
            return list_to_return

        current_tax_level = node_dict[current_staxid][1]

        if current_tax_level in tax_level_list:
            # then this is the taxonomic level we want, grab the name and return
            list_to_return[tax_level_list.index(current_tax_level)] = name_dict[current_staxid]
            if False not in list_to_return:
                # then we have found all of the taxa_levels
                return list_to_return
            else:
                # else we have not and we should continue the procession through the tax ids

                current_staxid = node_dict[current_staxid][0]

        else:
            current_staxid = node_dict[current_staxid][0]

def blast_worker(input_q, node_dict, name_dict):
    # this is the mp worker. We will aim to go to each of the directories and perform a blast and
    # and produce a dictionary that will be sequence id to taxonomic designation

    for path in iter(input_q.get, 'STOP'):
        sys.stdout.write('\nSample {}\n'.format(path.split('/')[-1]))

        sample_dir = '/'.join(path.split('/')[:-1])
        # we only need to be doing this sample if it hasn't already been done.
        if os.path.isfile('{}/sample_tax_dict.pickle'.format(sample_dir)) and os.path.isfile('{}/coral_dict.pickle'.format(sample_dir)) and os.path.isfile('{}/symbiodiniaceae_dict.pickle'.format(sample_dir)):
            continue

        # write out the screened fasta so that it can be read in to the blast
        # make sure to reference the sequence support and the iteration
        input_fasta_path = '{}/stability.trim.contigs.good.unique.abund.pcr.fasta'.format(sample_dir)


        # Set up environment for running local blast
        blastOutputPath = '{}/blast.out'.format(sample_dir)
        outputFmt = "6 qseqid sseqid staxids evalue pident qcovs staxid stitle ssciname"

        # no need to run blast if it has already been run and results have been pickled out
        if os.path.isfile('{}.pickle'.format(blastOutputPath)):
            blast_output_file = pickle.load(open('{}.pickle'.format(blastOutputPath), 'rb'))
        else:
            # Run local blast
            # completedProcess = subprocess.run([blastnPath, '-out', blastOutputPath, '-outfmt', outputFmt, '-query', inputPath, '-db', 'symbiodinium.fa', '-max_target_seqs', '1', '-num_threads', '1'])
            subprocess.run(
                ['blastn', '-out', blastOutputPath, '-outfmt', outputFmt, '-query', input_fasta_path, '-db', 'nt',
                 '-max_target_seqs', '10', '-num_threads', '20'])

            # Read in blast output
            with open(blastOutputPath, 'r') as f:
                blast_output_file = [line.rstrip() for line in f]

            # pickle out the blast results here to possibly save us time in the future
            pickle.dump(blast_output_file, open('{}.pickle'.format(blastOutputPath), 'wb'))

        # now create a dict that is the list of 10 results for linked to a key of the sample name
        blast_output_dict = defaultdict(list)
        for output_line in blast_output_file:
            components = output_line.split('\t')
            blast_output_dict[components[0]].append(components)

        # at this point we should have the blast output result read in
        # we can now make the taxa dictionary
        # this dict will hold the taxa for all samples
        sample_tax_dict = {}
        # this dict will hold the matches for a subset which are the coral samples
        coral_dict = {}
        # this dict will hold the matches for a subset which are the symbiodiniaceae samples
        symbiodiniaceae_dict = {}
        for blast_key, comp_list_list in blast_output_dict.items():
            sys.stdout.write('\rsequence {}'.format(blast_key))
            symbiodiniaceae_count = 0
            symbiodiniaceae_genus_dict = defaultdict(int)
            coral_count = 0
            # a dictionary that will keep track of which genus was the most abundant within the scleractinians
            coral_genus_dict = defaultdict(int)
            other_genus_count_dict = defaultdict(int)
            for comp_list in comp_list_list:

                # first check to see if the family is Symbiodiniacea
                try:
                    genus_level, family_level, order_level = get_taxa_designation_from_staxid(staxid=comp_list[6], tax_level_list=['genus', 'family', 'order'], node_dict=node_dict, name_dict=name_dict)
                except:
                    # there are some tax_ids that we don't seem to be able to find
                    # for the time being we will just continue over this seqeunce
                    continue
                if False in [genus_level, family_level, order_level]:
                    continue

                if family_level == 'Symbiodiniaceae':
                    symbiodiniaceae_count += 1
                    symbiodiniaceae_genus_dict[genus_level] += 1

                elif  order_level == 'Scleractinia' or order_level== 'Anthoathecata':
                    # check to see if we have a coral here and if we do, sore the genus
                    coral_count += 1
                    coral_genus_dict[genus_level] += 1
                else:
                    # if this is not coral and not symbiodinium then we populate the orther_genus_count_dict
                    other_genus_count_dict[genus_level] += 1

            # here we have been through each of the 10 hits for the sample and we should populate the
            # sample_tax_dict and the other dicts if the sample is either coral or Symbiodinium
            # see which is the biggest the Sym count or the scler or the other
            # first check to see if there were some annotations found

            if len(other_genus_count_dict.items()) > 0:
                max_other_taxa_count = max(other_genus_count_dict.values())
            else:
                max_other_taxa_count = 0
                if coral_count == 0 and symbiodiniaceae_count == 0:
                    # then we have had no valid annotations for this seqeunce
                    continue
            if max_other_taxa_count > max(coral_count, symbiodiniaceae_count):
                # then this is a not one of the subsets
                sample_tax_dict[blast_key] = [a[0] for a in sorted(other_genus_count_dict.items(), key=lambda x:x[1], reverse=True)][0]
            elif coral_count > max(max_other_taxa_count, symbiodiniaceae_count):
                # then this is a scleratinian and we should store the order
                sample_tax_dict[blast_key] = 'Scleractinia'
                # for the scleractinian dictionary we should associate to the most abundant genus
                most_abundant_genus_sclerac = [a[0] for a in sorted(coral_genus_dict.items(), key=lambda x:x[1], reverse=True)][0]
                coral_dict[blast_key] = most_abundant_genus_sclerac
            elif symbiodiniaceae_count > max(max_other_taxa_count, coral_count):
                # then this is a Symbiodiniaceae and we should store the genus in the specific dictionary
                sample_tax_dict[blast_key] = 'Symbiodiniaceae'
                # for the symbiodiniaceae dictinary we should associate to the most abundant genus
                most_abundant_genus_symbiodin = \
                [a[0] for a in sorted(symbiodiniaceae_genus_dict.items(), key=lambda x: x[1], reverse=True)][0]
                symbiodiniaceae_dict[blast_key] = most_abundant_genus_symbiodin
        # here we have populated the taxonomy dictionaries for the sample in question
        # we can now pickle out the dictionaries
        pickle.dump(sample_tax_dict, open('{}/sample_tax_dict.pickle'.format(sample_dir), 'wb'))
        pickle.dump(coral_dict, open('{}/coral_dict.pickle'.format(sample_dir), 'wb'))
        pickle.dump(symbiodiniaceae_dict, open('{}/symbiodiniaceae_dict.pickle'.format(sample_dir), 'wb'))






# this function will be responsible for creating the processed name and fasta paris from the fastq files
# that can be found in the info df
def sequence_QC(numProc=20):
    # for the time being I will run the samples through essentially the same processing as the basic
    # SymPortal processing.
    # this processing should be MPed
    info_df = generate_info_df_for_samples()

    # lets make an iput queue that is going to be each of the samples
    mothur_qc_input_queue = Queue()

    # for each sample in the info_df I will add a tuple that is a pair of the fwd and rev directories to the fastq.gz
    for ind in info_df.index.values.tolist():
        mothur_qc_input_queue.put((info_df.loc[ind, 'fastq_fwd_file_path'], info_df.loc[ind, 'fastq_rev_file_path']))

    for n in range(numProc):
        mothur_qc_input_queue.put('STOP')

    all_procs = []
    for n in range(numProc):
        p = Process(target=mothur_worker, args=(mothur_qc_input_queue,))
        all_procs.append(p)
        p.start()

    for p in all_procs:
        p.join()

    return

def mothur_worker(input_q):

    for path_tup in iter(input_q.get, 'STOP'):

        sample_dir = '/'.join(path_tup[0].split('/')[:-1])

        # the 18S V9 region primers
        primerFwdSeq = 'TTGTACACACCGCCC'  # Written 5'-->3'
        primerRevSeq = 'CCTTCYGCAGGTTCACCTAC'  # Written 5'-->3'

        oligoFile = [
            r'#1389F',
            'forward\t{0}'.format(primerFwdSeq),
            r'#1510R',
            'reverse\t{0}'.format(primerRevSeq)
        ]

        oligo_file_path = '{}/primers.oligos'.format(sample_dir)



        stability_file = ['{}\t{}'.format(path_tup[0], path_tup[1])]
        stability_file_path = '{}/stability.files'.format(sample_dir)
        root_name = 'stability'

        # write out stability file
        with open(stability_file_path, 'w') as f:
            for line in stability_file:
                f.write('{}\n'.format(line))

        # write out oligo file
        with open(oligo_file_path, 'w') as f:
            for line in oligoFile:
                f.write('{}\n'.format(line))

        # here we have the stability file written out and the oligo file

        # The mothur batch file that will be run by mothur.
        # todo if you find yourself reusing this you should add in a final unique.seqs at the end
        mBatchFile = [
            r'set.dir(input={0})'.format(sample_dir),
            r'set.dir(output={0})'.format(sample_dir),
            r'make.contigs(file={}, processors=20)'.format(stability_file_path),
            r'summary.seqs(fasta={}/{}.trim.contigs.fasta)'.format(sample_dir, root_name),
            r'screen.seqs(fasta={0}/{1}.trim.contigs.fasta, maxambig=0, maxhomop=5)'.format(
                sample_dir, root_name),
            r'summary.seqs(fasta={0}/{1}.trim.contigs.good.fasta)'.format(sample_dir, root_name),
            r'unique.seqs(fasta={0}/{1}.trim.contigs.good.fasta)'.format(sample_dir, root_name),
            r'summary.seqs(fasta={0}/{1}.trim.contigs.good.unique.fasta, name={0}/{1}.trim.contigs.good.names)'.format(
                sample_dir, root_name),
            r'split.abund(cutoff=2, fasta={0}/{1}.trim.contigs.good.unique.fasta, name={0}/{1}.trim.contigs.good.names)'.format(
                sample_dir, root_name),
            r'summary.seqs(fasta={0}/{1}.trim.contigs.good.unique.abund.fasta, name={0}/{1}.trim.contigs.good.abund.names)'.format(
                sample_dir, root_name),
            r'summary.seqs(fasta={0}/{1}.trim.contigs.good.unique.rare.fasta, name={0}/{1}.trim.contigs.good.rare.names)'.format(
                sample_dir, root_name),
            r'pcr.seqs(fasta={0}/{1}.trim.contigs.good.unique.abund.fasta, name={0}/{1}.trim.contigs.good.abund.names, oligos={0}/primers.oligos, pdiffs=2, rdiffs=2)'.format(
                sample_dir, root_name)
        ]

        mBatchFile_path = '{}/mBatchFile'.format(sample_dir)

        # write out batch file
        with open(mBatchFile_path, 'w') as f:
            for line in mBatchFile:
                f.write('{}\n'.format(line))

        # run the mothur processing
        # subprocess.run(['mothur', r'{0}'.format(mBatchFile_path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        subprocess.run(['mothur', r'{0}'.format(mBatchFile_path)])

        sys.stdout.write('mothur complete for {}'.format(sample_dir))

    return

# This is a function to parse over the directory stucture that was on the TARA/Genescope ftp and create
# an information dataframe
# This method used the next three methods as well generate_info_collection_dict, create_sample_dict and
# create_info_df_from_info_collection_dict
def generate_info_df_for_samples():

    if os.path.isfile('{}/info_df.pickle'.format(os.getcwd())):
        info_df = pickle.load(open('{}/info_df.pickle'.format(os.getcwd()), 'rb'))
    else:
        tara_data_dir = '/home/humebc/projects/tara/18s_data/18S_V9_1389F_1510R'

        # lets create a dict where the key will be the sample_name and the value will be dict with each of the above values
        info_collection_dict = {}

        # now lets parse through the directories using them to get some of the information facts above
        generate_info_collection_dict(info_collection_dict, tara_data_dir)

        # here we should have the info_collection_dict populated. We can now turn each of these into series and then
        # add them to the info_df
        columns_for_df = ['sample_name', 'fastq_fwd_file_path', 'fastq_rev_file_path', 'coral_plankton', 'spp_water', 'location', 'site', 'size_fraction']

        info_df = create_info_df_from_info_collection_dict(columns_for_df, info_collection_dict)

        pickle.dump(info_df, open('{}/info_df.pickle'.format(os.getcwd()), 'wb'))

        # here we have a dataframe with all of the samples in it
        # we can now move on to do the analysis of them.
        # this hsould be done in a separate method
    return info_df

def generate_info_collection_dict(info_collection_dict, tara_data_dir):
    for location in os.listdir(tara_data_dir):
        if 'ISLAND' in location:
            parsing_dir_loc = '{}/{}'.format(tara_data_dir, location)
            for site in os.listdir(parsing_dir_loc):
                parsing_dir_site = '{}/{}'.format(parsing_dir_loc, site)
                for sample_type in os.listdir(parsing_dir_site):
                    parsing_dir_sample_type = '{}/{}'.format(parsing_dir_site, sample_type)
                    if sample_type == 'CORAL':
                        for species in os.listdir(parsing_dir_sample_type):
                            parsing_dir_species = '{}/{}'.format(parsing_dir_sample_type, species)
                            for individual in os.listdir(parsing_dir_species):
                                parsing_dir_indi = '{}/{}/CS4L'.format(parsing_dir_species, individual)
                                # now we are in the directory that contains the actual paired fastq.gz files for a
                                # given coral individual
                                # collect the information we need
                                create_sample_dict(location, parsing_dir_indi, sample_type, site,
                                                                              species, info_collection_dict)

                    elif sample_type == 'PLANKTON':
                        for water_type in os.listdir(parsing_dir_sample_type):
                            if water_type == 'CSW':
                                parsing_dir_water_type = '{}/{}'.format(parsing_dir_sample_type, water_type)
                                for individual in os.listdir(parsing_dir_water_type):
                                    parsing_dir_indi = '{}/{}'.format(parsing_dir_water_type, individual)
                                    for size_fraction_indi in os.listdir(parsing_dir_indi):
                                        parsing_dir_indi_size = '{}/{}'.format(parsing_dir_indi, size_fraction_indi)
                                        # now we are in the directory that contains the actual paired fastq.gz files for a
                                        # given water sample

                                        # collect the information we need
                                        create_sample_dict(location, parsing_dir_indi_size, sample_type, site,
                                                                                      water_type, info_collection_dict, size_fraction_indi)


                            elif water_type == 'SURFACE':
                                # then this is a SURFACE sample and there are no individuals
                                parsing_dir_water_type = '{}/{}/S320'.format(parsing_dir_sample_type, water_type)

                                # collect the information we need
                                create_sample_dict(location, parsing_dir_water_type, sample_type, site, water_type, info_collection_dict, size_fraction_indi)


        elif 'OA' in location:
            parsing_dir_loc = '{}/{}/PLANKTON/SURFACE'.format(tara_data_dir, location)
            for size_fraction_indi in os.listdir(parsing_dir_loc):
                parsing_dir_indi = '{}/{}'.format(parsing_dir_loc, size_fraction_indi)

                sample_name = os.listdir(parsing_dir_indi)[0].split('_')[0]
                # FWD and REV PATHS
                files = os.listdir(parsing_dir_indi)
                if len(files) != 2:
                    print('more than 2 files in individual\'s directory {}'.format(parsing_dir_indi))
                    sys.exit(1)
                fwd_found = False
                rev_found = False
                for file_name in files:
                    if 'R1' in file_name:
                        fwd_path = '{}/{}'.format(parsing_dir_indi, file_name)
                        fwd_found = True
                    elif 'R2' in file_name:
                        rev_path = '{}/{}'.format(parsing_dir_indi, file_name)
                        rev_found = True
                # make sure that both the fwd and rev paths have been identified
                if not fwd_found or not rev_found:
                    print('fwd or rev read not found')
                    sys.exit(1)
                sample_dict = {'sample_name': sample_name,
                               'fastq_fwd_file_path': fwd_path,
                               'fastq_rev_file_path': rev_path,
                               'coral_plankton': 'OA',
                               'spp_water': 'PLANKTON',
                               'location': location,
                               'site': 'OA',
                               'size_fraction': size_fraction_indi
                               }
                info_collection_dict[sample_name] = sample_dict
    return

def create_sample_dict(location, parsing_dir, sample_type, site, water_type, info_collection_dict, size_fraction='coral'):
    # SAMPLE NAME
    sample_name = os.listdir(parsing_dir)[0].split('_')[0]
    # FWD and REV PATHS
    files = os.listdir(parsing_dir)
    if len(files) != 2:
        print('mergin fastqs in {}'.format(parsing_dir))
        # this method will merge the fastqs together so that 4 turn into 2.
        merge_fastqs(files, parsing_dir)

        # finally re-read in the files
        files = os.listdir(parsing_dir)
        if len(files) != 2:
            print('we still gotta problem {}'.format(parsing_dir))
            sys.exit(0)


    fwd_found = False
    rev_found = False
    for file_name in files:
        if 'R1' in file_name:
            fwd_path = '{}/{}'.format(parsing_dir, file_name)
            fwd_found = True
        elif 'R2' in file_name:
            rev_path = '{}/{}'.format(parsing_dir, file_name)
            rev_found = True
    # make sure that both the fwd and rev paths have been identified
    if not fwd_found or not rev_found:
        print('fwd or rev read not found')
        sys.exit(1)
    sample_dict = {'sample_name': sample_name,
                   'fastq_fwd_file_path': fwd_path,
                   'fastq_rev_file_path': rev_path,
                   'coral_plankton': sample_type,
                   'spp_water': water_type,
                   'location': location,
                   'site': site,
                   'size_fraction': size_fraction}
    info_collection_dict[sample_name] = sample_dict
    return


def merge_fastqs(files, parsing_dir):
    # If we get here then there have been multiple sequencing runs for the same sample
    # we will aim to simply merge the fastqs together insitu
    # first get the name that we want to use (this is the one containing BG
    one_found = False
    two_found = False
    for file_name in files:
        if 'BG' in file_name and 'R1' in file_name:
            R1_name = file_name
            one_found = True
        if 'BG' in file_name and 'R2' in file_name:
            R2_name = file_name
            two_found = True
    if not one_found or not two_found:
        # then this may be due to there not being a BG file
        # so we just take the first name
        for file_name in files:
            if 'R1' in file_name:
                R1_name = file_name
                R2_name = file_name.replace('R1', 'R2')
        print('couldnt find the right files {}'.format(parsing_dir))
    # second unzip all of the files
    for file_name in files:
        subprocess.run(['gunzip', '{}/{}'.format(parsing_dir, file_name)])
    # now we have all of the files unzipped
    # now look at each of the files and add them to a master R1 and R2.
    # it doesn't matter whether we add the R1 or R2 first. So long as we add the pair at the same time
    un_files = os.listdir(parsing_dir)
    master_fastq_R1 = []
    master_fastq_R2 = []
    for un_file in un_files:
        # we should go into this if twice once for each pair
        if 'R1' in un_file:
            # then we can read this and its pair
            rone_path = '{}/{}'.format(parsing_dir, un_file)
            with open(rone_path, 'r') as f:
                rone_file = [line.rstrip() for line in f]
            master_fastq_R1.extend(rone_file)

            # now do the same for the corresponing R2 file
            rtwo_path = '{}/{}'.format(parsing_dir, un_file.replace('R1', 'R2'))
            with open(rtwo_path, 'r') as f:
                rtwo_file = [line.rstrip() for line in f]
            master_fastq_R2.extend(rtwo_file)

            # now delte the files
            os.remove(rone_path)
            os.remove(rtwo_path)
    # here we have a master file for the R1 and R2. We can now write it out to the name that we have
    rone_path_to_write = '{}/{}'.format(parsing_dir, R1_name.replace('.gz', ''))
    with open(rone_path_to_write, 'w') as f:
        for line in master_fastq_R1:
            f.write('{}\n'.format(line))
    rtwo_path_to_write = '{}/{}'.format(parsing_dir, R2_name.replace('.gz', ''))
    with open(rtwo_path_to_write, 'w') as f:
        for line in master_fastq_R2:
            f.write('{}\n'.format(line))
    # now we simply need to recompress the files
    un_files = os.listdir(parsing_dir)
    for un_file in un_files:
        subprocess.run(['gzip', '{}/{}'.format(parsing_dir, un_file)])


def create_info_df_from_info_collection_dict(columns_for_df, info_collection_dict):
    series_list = []
    for sample_key, sample_info_dict in info_collection_dict.items():
        data_for_series = [sample_info_dict[ind] for ind in columns_for_df]
        temp_series = pd.Series(data_for_series, index=columns_for_df, name=sample_key)
        series_list.append(temp_series)
    # now we can populate the info df using the series list
    info_df = pd.DataFrame.from_items([(s.name, s) for s in series_list]).T
    return info_df

def get_colour_list():
    colour_list = ["#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059", "#FFDBE5",
                  "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693",
                  "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900",
                  "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
                  "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744",
                  "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68",
                  "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED",
                  "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                  "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00", "#7900D7",
                  "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F",
                  "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55",
                  "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C", "#83AB58", "#001C1E", "#D1F7CE", "#004B28",
                  "#C8D0F6", "#A3A489", "#806C66", "#222800", "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59",
                  "#8ADBB4", "#1E0200", "#5B4E51", "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94",
                  "#7ED379", "#012C58", "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393",
                  "#943A4D", "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
                  "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5", "#E773CE",
                  "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4", "#00005F", "#A97399",
                  "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01", "#6B94AA", "#51A058", "#A45B02",
                  "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966", "#64547B", "#97979E", "#006A66", "#391406",
                  "#F4D749", "#0045D2", "#006C31", "#DDB6D0", "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9",
                  "#FFFFFE", "#C6DC99", "#203B3C", "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527",
                  "#8BB400", "#797868", "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C",
                  "#B88183", "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
                  "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F", "#003109",
                  "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E", "#1A3A2A", "#494B5A",
                  "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F", "#BDC9D2", "#9FA064", "#BE4700",
                  "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00", "#061203", "#DFFB71", "#868E7E", "#98D058",
                  "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66", "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F",
                  "#545C46", "#866097", "#365D25", "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"]
    return colour_list
plot_pcoa_spp_18s_its2()