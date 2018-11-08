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
        # go to the directory and read in each of the
        # stability.trim.contigs.fasta for the contig summary

        # run a summary.seqs on the stability.trim.contigs.good.unique.abund.pcr.fasta

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

    input_q_for_blast = Queue()

    for ind in info_df.index.values.tolist():
        if 'CORAL' in info_df[ind, 'fastq_fwd_file_path']:
            input_q_for_blast.put(info_df.loc[ind, 'fastq_fwd_file_path'])

    for n in range(numProc):
        input_q_for_blast.put('STOP')

    all_procs = []
    for n in range(numProc):
        p = Process(target=mothur_worker, args=(input_q_for_blast,))
        all_procs.append(p)
        p.start()

    for p in all_procs:
        p.join()

    return

def blast_worker(input_q):
    # this is the mp worker. We will aim to go to each of the directories and perform a blast and
    # and produce a dictionary that will be sequence id to taxonomic designation

    for path_tup in iter(input_q.get, 'STOP'):

        sample_dir = '/'.join(path_tup[0].split('/')[:-1])

        # Write out the hidden file that points to the ncbi database directory.
        ncbircFile = []
        # db_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB'))

        db_path = '/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/ntdbdownload'
        ncbircFile.extend(["[BLAST]", "BLASTDB={}".format(db_path)])

        # write out the ncbircFile
        with open(sample_dir, 'w') as f:
            for line in ncbircFile:
                f.write('{}\n'.format(line))

        # write out the screened fasta so that it can be read in to the blast
        # make sure to reference the sequence support and the iteration
        input_fasta_path = '{}/stability.trim.contigs.good.unique.abund.pcr.fasta'.format(sample_dir)


        # Set up environment for running local blast
        blastOutputPath = '{}/blast.out'.format(sample_dir)
        outputFmt = "6 qseqid sseqid staxids evalue pident qcovs staxid stitle ssciname"


        # Run local blast
        # completedProcess = subprocess.run([blastnPath, '-out', blastOutputPath, '-outfmt', outputFmt, '-query', inputPath, '-db', 'symbiodinium.fa', '-max_target_seqs', '1', '-num_threads', '1'])
        completedProcess = subprocess.run(
            ['blastn', '-out', blastOutputPath, '-outfmt', outputFmt, '-query', input_fasta_path, '-db', 'nt',
             '-max_target_seqs', '10', '-num_threads', '20'])

        # Read in blast output
        with open(blastOutputPath, 'r') as f:
            blast_output_file = [line.rstrip() for line in f]

        # at this point we should have the blast output result read in
        # we can now make the taxa dictionary
        # this dict will hold the taxa for all samples
        sample_tax_dict = {}

        # this dict will hold the matches for a subset which are the coral samples
        scleractinian_reads = {}
        for output_line in blast_output_file:
            components = output_line.split('\t')
            seq_name = 'asdfples'
            seq_taxa_match = 'asdfas'




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


sequence_QC()