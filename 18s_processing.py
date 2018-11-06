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


generate_info_df_for_samples()