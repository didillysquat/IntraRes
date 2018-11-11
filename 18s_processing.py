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



# Code for generating the figure that will let us examine what our profile look like.
# this first go at this will not take into account any generated type profiles but rather just plot the sequences
# in each sample.

# the first step for this is to generate an info df.
# we can generate this from the current info_df
# we should call it fig_info_df
def generate_fig_info_df():
    info_df = generate_info_df_for_samples()

    fig_info_df = pd.DataFrame(columns = ['island', 'site', 'genus', 'individual', 'sample_dir'])
    for ind in info_df.index.values.tolist():
        if 'CORAL' in info_df.loc[ind, 'fastq_fwd_file_path']:
            fastq_string = info_df.loc[ind, 'fastq_fwd_file_path']
            components = fastq_string.split('/')
            smp_site = components[-6]
            smp_island = components[-7]
            smp_individual = components[-3]
            smp_genus = components[-4]
            smp_dir = '/'.join(components[:-1])
            fig_info_df = fig_info_df.append(pd.Series(name = ind, data = [smp_island, smp_site, smp_genus, smp_individual, smp_dir], index=['island', 'site', 'genus', 'individual', 'sample_dir']))

    # here we have the fig_info_df generated. We can use this for making the figure
    # SETUP AXES
    # https://matplotlib.org/users/gridspec.html
    fig = plt.figure(figsize=(14, 10))

    # the bottom row will be for the legend
    # the second to last will just be invisible to give a space between the legend and the other plots
    # we also want to include a gridspec plot after each of the main three. These will hold the csw and surface
    # samples
    gs = plt.GridSpec(5, 3, figure=fig, height_ratios=[1, 1, 1, 0.2, 1])
    # within each of the GrdiSpec subplots we will make a subplotspec which is three plots on one row

    ax_list = []

    grid_spec_subplot_list = []
    # make the main 3x3 axis
    for row_ind in range(3):
        for col_ind in range(3):
            # put in the main data 3 plots
            temp_grid_spec_subplot = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[row_ind, col_ind])
            grid_spec_subplot_list.append(temp_grid_spec_subplot)
            for i in range(3):
                # NB this might be a 2d array, lets see.
                ax = plt.Subplot(fig, temp_grid_spec_subplot[i])
                ax_list.append(ax)
                fig.add_subplot(ax)


    # now do the invisible row that will give us the space we want
    ax_space = plt.subplot(gs[3, :])
    remove_axes_but_allow_labels(ax_space)

    leg_axes = plt.subplot(gs[4, :])

    # here we have the axes set up.
    # before we move onto the actual plotting of the samples we should create an sp_output_df_div equivalent
    # for the samples. This will mean going through each of the coral sample directories and collecting relative
    # abundances. We will do this using the genus specific fastas.
    sample_seqs_df = generate_seq_abundance_df()
    # so first lets finish these fastas
    # plot_data_axes(ax_list, extra_ax_list, colour_dict_div, colour_dict_type, info_df, ordered_sample_list,
    #                smp_id_to_smp_name_dict,
    #                smp_name_to_smp_id_dict, sp_output_df_div, sp_output_df_type, individual_info)

    return

def generate_seq_abundance_df(numProc=20):
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
    fig_info_df = generate_info_df_for_samples()

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
    return

def abundance_worker(input_q):
    for sample_name, sample_dir in iter(input_q.get, 'STOP'):
        sys.stdout.write('\nSample {}\n'.format(sample_dir))
        genus_of_sample = sample_dir.split('/')[-3]


    return

# def plot_data_axes(ax_list, extra_ax_list, colour_dict_div, colour_dict_type, info_df, ordered_sample_list,
#                    smp_id_to_smp_name_dict,
#                    smp_name_to_smp_id_dict, sp_output_df_div, sp_output_df_type, indi_indo):
#     ax_count = 0
#     extra_ax_count = 0
#     for site in ['SITE01', 'SITE02', 'SITE03']:
#         for location in ['ISLAND06', 'ISLAND10', 'ISLAND15']:
#             for spp in ['PORITES', 'POCILLOPORA', 'MILLEPORA']:
#                 ax = ax_list[ax_count]
#                 patches_list = []
#                 ind = 0
#                 colour_list = []
#
#                 # for each set of location, site and spp, we basically want to get a list of the samples
#                 # that meet the set criteria, we then want to plot samples according to the ordered_sample_list
#                 # order which will be in IDs. As such we will have to convert the sample_name in the info_df
#                 # to a sample ID using the smp_name_to_smp_id_dict.
#
#                 # get sample_names that fit the requirements
#                 sample_names_of_set = info_df.loc[
#                     (info_df['location'] == location) &
#                     (info_df['site'] == site) &
#                     (info_df['spp_water'] == spp)
#                     ].index.values.tolist()
#
#                 # temporarily remove CO0002044, CO0002041
#                 # from the above list
#                 # if 'CO0002044' or 'CO0002041' in sample_names_of_set:
#                 #     sample_names_of_set = [name for name in sample_names_of_set if
#                 #                            name not in ['CO0002044', 'CO0002041']]
#
#                 # convert these to sample IDs
#                 # The sample names in symportal are actually the full file names version rather than
#                 # the shorter versions in the info_df. As such we should we will have to do a conversion here
#                 full_sample_names = [
#                     '_'.join(info_df.loc[smp_name]['fastq_fwd_file_path'].split('/')[-1].split('_')[:3]) for
#                     smp_name in
#                     sample_names_of_set]
#                 try:
#                     smple_ids_of_set = [smp_name_to_smp_id_dict[smp_name] for smp_name in full_sample_names]
#                 except:
#                     apples = 'asdf'
#                 # now we want to plot in the order of the ordered_sample_list
#                 ordered_smple_ids_of_set = [smpl_id for smpl_id in ordered_sample_list if
#                                             smpl_id in smple_ids_of_set]
#
#                 coral_csw_x_val_tup_list = None
#                 if spp == 'POCILLOPORA':
#                     # here we need to work out if any of these are csw associated samples
#                     # if so then we need to follow through and get the id of them
#                     coral_csw_id_one = None
#                     coral_csw_id_ten = None
#                     for smp_name in sample_names_of_set:
#                         if smp_name in indi_indo.keys():
#                             # use the same logic below to get the id of this sample
#                             temp_full_name = '_'.join(
#                                 info_df.loc[smp_name]['fastq_fwd_file_path'].split('/')[-1].split('_')[:3])
#                             # then this is one of the csw associated samples
#                             if indi_indo[smp_name] == 'INDIVIDUAL1':
#                                 coral_csw_id_one = smp_name_to_smp_id_dict[temp_full_name]
#                             elif indi_indo[smp_name] == 'INDIVIDUAL10':
#                                 coral_csw_id_ten = smp_name_to_smp_id_dict[temp_full_name]
#
#                     # now we need to get the x values of where to draw the line to and from for the csw associated corals
#                     # we can work this out by seeing what the index of the sample ids are in the ordered_smple list
#                     # as the width we use is 1, the x values will then be the index -+ 0.5
#                     coral_csw_x_val_tup_list = []
#                     if coral_csw_id_one:
#                         if coral_csw_id_one in ordered_smple_ids_of_set:
#                             temp_index_of_sample_in_list = ordered_smple_ids_of_set.index(coral_csw_id_one)
#                             coral_csw_x_value_tup = (
#                             temp_index_of_sample_in_list - 0.5, temp_index_of_sample_in_list + 0.5, 'brown')
#                             coral_csw_x_val_tup_list.append(coral_csw_x_value_tup)
#                     if coral_csw_id_ten:
#                         if coral_csw_id_ten in ordered_smple_ids_of_set:
#                             temp_index_of_sample_in_list = ordered_smple_ids_of_set.index(coral_csw_id_ten)
#                             coral_csw_x_value_tup = (
#                             temp_index_of_sample_in_list - 0.5, temp_index_of_sample_in_list + 0.5, 'gray')
#                             coral_csw_x_val_tup_list.append(coral_csw_x_value_tup)
#
#                 num_smp_in_this_subplot = len(ordered_smple_ids_of_set)
#                 x_tick_label_list = []
#
#                 for smple_id_to_plot in ordered_smple_ids_of_set:
#                     # General plotting
#                     sys.stdout.write('\rPlotting sample: {}'.format(smple_id_to_plot))
#                     x_tick_label_list.append(smp_id_to_smp_name_dict[smple_id_to_plot].split('_')[0])
#                     # for each sample we will start at 0 for the y and then add the height of each bar to this
#
#                     # PLOT DIVs
#                     plot_div_over_type(colour_dict_div, colour_list, ind, patches_list, smple_id_to_plot,
#                                        sp_output_df_div)
#
#                     # PLOT type
#                     plot_type_under_div(colour_dict_type, colour_list, ind, patches_list, smple_id_to_plot,
#                                         sp_output_df_type)
#                     ind += 1
#
#                 paint_rect_to_axes_div_and_type(ax=ax, colour_list=colour_list,
#                                                 num_smp_in_this_subplot=num_smp_in_this_subplot,
#                                                 patches_list=patches_list,
#                                                 coral_csw_x_val_tup_list=coral_csw_x_val_tup_list,
#                                                 x_tick_label_list=x_tick_label_list,
#                                                 max_num_smpls_in_subplot=10)
#
#                 ax_count += 1
#
#             # PLOT csw and surface
#
#             # first identify which the csw sample is
#             csw_samples = info_df.loc[
#                 (info_df['location'] == location) &
#                 (info_df['site'] == site) &
#                 (info_df['spp_water'] == 'CSW')
#                 ].index.values.tolist()
#
#             # convert these to sample IDs
#             # The sample names in symportal are actually the full file names version rather than
#             # the shorter versions in the info_df. As such we should we will have to do a conversion here
#             full_sample_names_csw = [
#                 '_'.join(info_df.loc[smp_name]['fastq_fwd_file_path'].split('/')[-1].split('_')[:3]) for
#                 smp_name in csw_samples]
#             smple_ids_of_set_csw = [smp_name_to_smp_id_dict[smp_name] for smp_name in full_sample_names_csw]
#
#             # here we need to work out which sample is which individual for the csw association to corals
#
#             coral_csw_id_one = None
#             coral_csw_id_ten = None
#             for smp_name in csw_samples:
#                 if smp_name in indi_indo.keys():
#                     # use the same logic below to get the id of this sample
#                     temp_full_name = '_'.join(
#                         info_df.loc[smp_name]['fastq_fwd_file_path'].split('/')[-1].split('_')[:3])
#                     # then this is one of the csw associated samples
#                     if indi_indo[smp_name] == 'INDIVIDUAL1':
#                         coral_csw_id_one = smp_name_to_smp_id_dict[temp_full_name]
#                     elif indi_indo[smp_name] == 'INDIVIDUAL10':
#                         coral_csw_id_ten = smp_name_to_smp_id_dict[temp_full_name]
#
#             # now we need to get the x values of where to draw the line to and from for the csw associated corals
#             # we can work this out by seeing what the index of the sample ids are in the ordered_smple list
#             # as the width we use is 1, the x values will then be the index -+ 0.5
#             coral_csw_x_val_tup_list = []
#
#             if coral_csw_id_one and coral_csw_id_ten:
#                 # then both exist
#                 # check to see that they are in the order 'one' then 'ten'. if not, reverse
#                 if smple_ids_of_set_csw.index(coral_csw_id_one) == 1:
#                     # then we need to reverse
#                     smple_ids_of_set_csw = list(reversed(smple_ids_of_set_csw))
#                 # now add the tuples to the coral_csw_x_val_tup_list
#                 for i, colour in enumerate(['brown', 'gray']):
#                     coral_csw_x_value_tup = (i - 0.5, i + 0.5, colour)
#                     coral_csw_x_val_tup_list.append(coral_csw_x_value_tup)
#
#             num_smp_in_this_subplot = len(ordered_smple_ids_of_set)
#             x_tick_label_list = []
#
#             colour_list = []
#             ind = 0
#             patches_list = []
#             for smple_id_to_plot in smple_ids_of_set_csw:
#                 bottom_div = 0
#                 # for each sequence, create a rect patch
#                 # the rect will be 1 in width and centered about the ind value.
#                 for seq in list(sp_output_df_div):
#                     # class matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
#                     rel_abund_div = sp_output_df_div.loc[smple_id_to_plot, seq]
#                     if rel_abund_div > 0:
#                         patches_list.append(
#                             Rectangle((ind - 0.5, bottom_div), 1, rel_abund_div, color=colour_dict_div[seq]))
#                         # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
#                         colour_list.append(colour_dict_div[seq])
#                         bottom_div += rel_abund_div
#                 ind += 1
#
#             paint_rect_to_axes_div_and_type(ax=extra_ax_list[extra_ax_count], colour_list=colour_list,
#                                             num_smp_in_this_subplot=2,
#                                             coral_csw_x_val_tup_list=coral_csw_x_val_tup_list,
#                                             patches_list=patches_list, max_num_smpls_in_subplot=2)
#             extra_ax_count += 1
#
#             # now get the surface sample
#             surface_sample = info_df.loc[
#                 (info_df['location'] == location) &
#                 (info_df['site'] == site) &
#                 (info_df['spp_water'] == 'SURFACE')
#                 ].index.values.tolist()[0]
#
#             # convert these to sample IDs
#             # The sample names in symportal are actually the full file names version rather than
#             # the shorter versions in the info_df. As such we should we will have to do a conversion here
#             full_sample_name_surface = '_'.join(
#                 info_df.loc[surface_sample]['fastq_fwd_file_path'].split('/')[-1].split('_')[:3])
#
#             smple_id_of_surface_sample = smp_name_to_smp_id_dict[full_sample_name_surface]
#
#             colour_list = []
#             bottom_div = 0
#             ind = 0
#             patches_list = []
#             # for each sequence, create a rect patch
#             # the rect will be 1 in width and centered about the ind value.
#             for seq in list(sp_output_df_div):
#                 # class matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
#                 rel_abund_div = sp_output_df_div.loc[smple_id_of_surface_sample, seq]
#                 if rel_abund_div > 0:
#                     patches_list.append(
#                         Rectangle((ind - 0.5, bottom_div), 1, rel_abund_div, color=colour_dict_div[seq]))
#                     # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
#                     colour_list.append(colour_dict_div[seq])
#                     bottom_div += rel_abund_div
#
#             paint_rect_to_axes_div_and_type(ax=extra_ax_list[extra_ax_count], colour_list=colour_list,
#                                             num_smp_in_this_subplot=1,
#                                             patches_list=patches_list, max_num_smpls_in_subplot=2)
#             extra_ax_count += 1

def remove_axes_but_allow_labels(ax, x_tick_label_list=None):
    ax.set_frame_on(False)
    if not x_tick_label_list:
        ax.set_xticks([])
    ax.set_yticks([])

def setup_axes(max_n_cols, max_n_rows):
    # https://matplotlib.org/users/gridspec.html
    fig = plt.figure(figsize=(14, 10))
    # the bottom row will be for the legend
    # the second to last will just be invisible to give a space between the legend and the other plots
    # we also want to include a gridspec plot after each of the main three. These will hold the csw and surface
    # samples
    # we will put in an invisible spacer so as an extra row
    gs = plt.GridSpec(max_n_rows + 1, max_n_cols)
    ax_list = []
    # first make an ax for each of the main data plots
    for row in range(max_n_rows - 1):
        for col in range(max_n_cols):
            ax = plt.Subplot(fig, gs[row, col])
            ax_list.append(ax)
            fig.add_subplot(ax)
    # now make the ax for the space
    blank_ax = plt.subplot(gs[max_n_rows - 1, :])
    # make the axes invisible for the space ax
    remove_axes_but_allow_labels(blank_ax)
    # now split up the final row to put the legend in. One for DIVs and one for TYPEs
    temp_grid_spec_subplot_leg = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[-1, :])
    leg_axes = []
    for i in range(2):
        ax = plt.Subplot(fig, temp_grid_spec_subplot_leg[i])
        leg_axes.append(ax)
        fig.add_subplot(ax)
    return ax_list, leg_axes, fig

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
        if not os.path.isfile('{}/sample_tax_dict.pickle'.format(sample_dir)) or not os.path.isfile('{}/scleractinian_dict.pickle'.format(sample_dir)) or not os.path.isfile('{}/symbiodiniaceae_dict.pickle'.format(sample_dir)):
            continue

        # next check to see if we have already generated the fastas for this sample.
        if os.path.isfile('{}/scleractinian.fasta'.format(sample_dir)) and os.path.isfile('{}/scleractinian_genus_only.fasta'.format(sample_dir)):
            continue

        scleractinian_dict = pickle.load(open('{}/scleractinian_dict.pickle'.format(sample_dir), 'rb'))

        # need to read in the fasta that holds the sequences and create a dict
        with open('{}/stability.trim.contigs.good.unique.abund.pcr.fasta'.format(sample_dir), 'r') as f:
            fasta_file = [line.rstrip() for line in f]
        fasta_dict = {fasta_file[i][1:] :fasta_file[i+1] for i in range(0, len(fasta_file),2)}

        scleractinian_fasta = []
        genus_fasta = []
        for seq_name, seq_genus in scleractinian_dict.keys():
            scleractinian_fasta.append('>{}'.format(seq_name))
            scleractinian_fasta.append(fasta_dict[seq_name])
            if seq_genus == genus_of_sample:
                genus_fasta.append('>{}'.format(seq_name))
                genus_fasta.append(fasta_dict[seq_name])

        # here we have created both of the fastas and we can now write them out
        with open('{}/scleractinian.fasta'.format(sample_dir), 'w') as f:
            for line in scleractinian_fasta:
                f.write('{}\n'.format(line))

        # here we have created both of the fastas and we can now write them out
        with open('{}/scleractinian_genus_only.fasta'.format(sample_dir), 'w') as f:
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
        if os.path.isfile('{}/sample_tax_dict.pickle'.format(sample_dir)) and os.path.isfile('{}/scleractinian_dict.pickle'.format(sample_dir)) and os.path.isfile('{}/symbiodiniaceae_dict.pickle'.format(sample_dir)):
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
        scleractinian_dict = {}
        # this dict will hold the matches for a subset which are the symbiodiniaceae samples
        symbiodiniaceae_dict = {}
        for blast_key, comp_list_list in blast_output_dict.items():
            sys.stdout.write('\rsequence {}'.format(blast_key))
            symbiodiniaceae_count = 0
            symbiodiniaceae_genus_dict = defaultdict(int)
            scleractinian_count = 0
            # a dictionary that will keep track of which genus was the most abundant within the scleractinians
            scleractinian_genus_dict = defaultdict(int)
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

                elif  order_level == 'Scleractinia':
                    # check to see if we have a scleractinian here and if we do, sore the genus
                    scleractinian_count += 1
                    scleractinian_genus_dict[genus_level] += 1
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
                if scleractinian_count == 0 and symbiodiniaceae_count == 0:
                    # then we have had no valid annotations for this seqeunce
                    continue
            if max_other_taxa_count > max(scleractinian_count, symbiodiniaceae_count):
                # then this is a not one of the subsets
                sample_tax_dict[blast_key] = [a[0] for a in sorted(other_genus_count_dict.items(), key=lambda x:x[1], reverse=True)][0]
            elif scleractinian_count > max(max_other_taxa_count, symbiodiniaceae_count):
                # then this is a scleratinian and we should store the order
                sample_tax_dict[blast_key] = 'Scleractinia'
                # for the scleractinian dictionary we should associate to the most abundant genus
                most_abundant_genus_sclerac = [a[0] for a in sorted(scleractinian_genus_dict.items(), key=lambda x:x[1], reverse=True)][0]
                scleractinian_dict[blast_key] = most_abundant_genus_sclerac
            elif symbiodiniaceae_count > max(max_other_taxa_count, scleractinian_count):
                # then this is a Symbiodiniaceae and we should store the genus in the specific dictionary
                sample_tax_dict[blast_key] = 'Symbiodiniaceae'
                # for the symbiodiniaceae dictinary we should associate to the most abundant genus
                most_abundant_genus_symbiodin = \
                [a[0] for a in sorted(symbiodiniaceae_genus_dict.items(), key=lambda x: x[1], reverse=True)][0]
                symbiodiniaceae_dict[blast_key] = most_abundant_genus_symbiodin
        # here we have populated the taxonomy dictionaries for the sample in question
        # we can now pickle out the dictionaries
        pickle.dump(sample_tax_dict, open('{}/sample_tax_dict.pickle'.format(sample_dir), 'wb'))
        pickle.dump(scleractinian_dict, open('{}/scleractinian_dict.pickle'.format(sample_dir), 'wb'))
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


generate_tax_dictionaries()