import pandas as pd
import os
import sys
import pickle
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
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
from skbio.tree import TreeNode
from skbio.diversity import beta_diversity
from mpl_toolkits.mplot3d import Axes3D
from plumbum import local
import cropping, general, fasta2net
import hashlib


class EighteenSAnalysis:
    def __init__(self):
        self.root_dir = os.path.abspath(os.path.dirname(__file__))
        self.data_root_dir = os.path.join(self.root_dir, '18S_V9_1389F_1510R')
        # This is the directory where we will pickle out objects to create caches of them
        self.cache_dir = os.path.join(self.root_dir, 'cache')
        # This is the directory from which we will read in input files
        self.input_dir = os.path.join(self.root_dir, 'input')
        # This is the directory in which we will save the figure outputs
        self.figure_output_dir = os.path.join(self.root_dir, 'figures')
        # This is the directory in which we will save the .dist files that will be used in
        # calculating the PCoAs.
        self.dist_output_dir = os.path.join(self.root_dir, 'dist')
        os.makedirs(self.dist_output_dir, exist_ok=True)
        os.makedirs(self.figure_output_dir, exist_ok=True)
        os.makedirs(self.input_dir, exist_ok=True)
        os.makedirs(self.cache_dir,exist_ok=True)

        # Two dfs that contain info for samples
        # This one is all samples
        self.all_samples_info_df = None
        # This is only coral samples and is used in the plotting
        self.coral_info_df_for_figures = None
        self._generate_info_df_for_samples()

        # This is the sample order from the its2 work
        self.sample_order = pickle.load(open(os.path.join(self.input_dir, 'ordered_sample_names_from_its2_work.pickle'), 'rb'))

        self.blast_nt_db_path = '/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/ntdbdownload'

        # if running on the server we will try these first
        self.symportal_seq_output_relative_path = '/home/humebc/projects/tara/initial_its2_processing/2018-10-21_08-59-37.620726.DIVs.relative.txt'
        self.symportal_profile_output_relative_path = '/home/humebc/projects/tara/initial_its2_processing/34_init_tara_standalone_all_samps_151018_2018-10-21_08-45-56.507454.profiles.relative.txt'
        # else look to see if they are in the input directory
        self.symportal_seq_output_relative_path_input_dir = os.path.join(self.input_dir, '2018-10-21_08-59-37.620726.DIVs.relative.txt')
        self.symportal_profile_output_relative_path_input_dir = os.path.join(self.input_dir, '34_init_tara_standalone_all_samps_151018_2018-10-21_08-45-56.507454.profiles.relative.txt')
    def plot_pcoa_spp_18s_its2(self, distance_method='bray_curtis'):
        es_its2_plotter = self.E18S_ITS2_PCOA_FIGURE(parent=self, distance_method=distance_method)
        es_its2_plotter.plot()

    def plot_pcoa_spp(self, distance_method='braycurtis'):
        spp_pcoa_plotter = self.PlotPCoASpp(parent=self, distance_method=distance_method)
        spp_pcoa_plotter.plot()

    def plot_pcoa_spp_island(self, distance_method='braycurtis'):
        spp_island_pcoa_plotter = self.PlotPCoASppIsland(parent=self, distance_method=distance_method)
        spp_island_pcoa_plotter.plot()

    def do_taxa_annotations(self):
        annotater = self.TaxaAnnotation(parent=self)
        annotater.annotate()

    def plot_seq_stacked_bar_plots(self, plot_type):
        """This method produced stacked bar charts. It can produce different charts depending on the plot_type.
                'full' means all of the sequences including the maj
                'low' means without the maj
                'med' means without the maj and having been through med
                'qc_taxa_rel_abund' means plot the relative abundance of the post qc taxa categories
                'qc_absolute means' plot the absolute abundance of the post-qc sequences (all tax categories)"""
        seq_stacked_bar_plotter = self.SeqStackedBarPlotter(plot_type=plot_type, parent=self)
        seq_stacked_bar_plotter.plot()

    def do_qc(self):
        sqc = self.SequenceQC(parent = self)
        sqc.do_seq_qc()

    class Generic_PCOA_DIST_Methods:
        # TODO implement UniFrac
        """A class that will hold all of the methods that are required by the various PCoA-plotting classes"""
        def __init__(self, parent):
            self.parent = parent

            # For unifrac only
            self.abundance_df = None
            self.unaligned_fasta_path = os.path.join(self.parent.dist_output_dir, 'seqs_for_unifrac.unaligned.fasta')
            self.aligned_fasta_path = self.unaligned_fasta_path.replace('unaligned', 'aligned')
            self.tree_out_path_unrooted = self.aligned_fasta_path + '.treefile'
            self.tree_out_path_rooted = self.tree_out_path_unrooted.replace('.treefile', '.rooted.treefile')
            self.rooted_tree = None

        def _generate_unifrac_distance_and_pcoa_spp(self):
            """
            This is code for generating PCOAs for each of the coral species.
            Read in the minor div dataframe which should have normalised abundances in them
            For each sample we have a fasta that we can read in which has the normalised (to 1000) sequences
            For feeding into med. We can use a default dict to collect the sequences and abundances from this fairly
            simply.
            This is likely best done one for each sample outside of the pairwise comparison to save on redoing the same
            collection of the sequences.
            We are now implementing this to use Unifrac rather than BrayCurtis. For the unifrac calculation, we will
            need:
            An abundance dataframe that is sample as rows and seqs as columns. We should be able to use the
            minor_div_abundance_dict to make this.
            We will also need a list of sequences that we will need to build a tree for.
            The key for the dictionary is the actual nucleotide sequences so we can get the sequences from there.
            If so then we will also be able to get this from the below methods.
            We can then work with implementing the code form SymPortal to do most of the heavy
            lifting fro making the distance matrix and the pcoa.
            We will create the tree outside of the species for loop. That way we only have to compute the tree
            once and we should be able to use the same tree for each of the species that we do. We should also
            be able to use the same tree for the other unifrac work we do.
            """

            if os.path.isfile(os.path.join(self.parent.cache_dir, 'minor_div_abundance_dict.p')):
                minor_div_abundance_dict = pickle.load(
                    open(os.path.join(self.parent.cache_dir, 'minor_div_abundance_dict.p'), 'rb'))
            else:
                minor_div_abundance_dict = self._generate_minor_div_abundance_dict_from_scratch()

            self._create_df_from_minor_div_dict(minor_div_abundance_dict)

            columns, seq_fasta_list = self._make_all_seq_fasta_and_hash_names()

            self._set_df_cols_as_hashes(columns)

            self._align_seqs(seq_fasta_list)

            self._make_and_root_tree()

            spp_unifrac_pcoa_df_dict = {}
            for spp in ['Porites', 'Pocillopora', 'Millepora']:
                if os.path.isfile(os.path.join(self.parent.cache_dir, f'spp_unifrac_pcoa_df_{spp}.p')):
                    spp_unifrac_pcoa_df = pickle.load(
                        open(os.path.join(self.parent.cache_dir, f'spp_unifrac_pcoa_df_{spp}.p'), 'rb'))
                else:
                    spp_unifrac_pcoa_df = self._make_spp_unifrac_pcoa_df_from_scratch(spp)
                spp_unifrac_pcoa_df_dict[spp] = spp_unifrac_pcoa_df
            return spp_unifrac_pcoa_df_dict

        def _make_spp_unifrac_pcoa_df_from_scratch(self, spp):
            sample_names_of_spp, spp_df = self._get_subset_spp_df_for_unifrac(spp)
            # perform unifrac
            print('Performing unifrac calculations')
            wu = beta_diversity(
                metric='weighted_unifrac', counts=spp_df.to_numpy(),
                ids=[str(_) for _ in list(spp_df.index)],
                tree=self.rooted_tree, otu_ids=[str(_) for _ in list(spp_df.columns)])
            spp_unifrac_pcoa_df = self._do_spp_pcoa_unifrac(sample_names_of_spp, spp, wu)
            return spp_unifrac_pcoa_df

        def _get_subset_spp_df_for_unifrac(self, spp):
            # Get a list of the samples that we should be working with
            sample_names_of_spp = self.parent.coral_info_df_for_figures.loc[
                self.parent.coral_info_df_for_figures['genus'] == spp.upper()].index.values.tolist()
            # remove the two porites species form this that seem to be total outliers
            if spp == 'Porites':
                sample_names_of_spp.remove('CO0001674')
                sample_names_of_spp.remove('CO0001669')
            # This is a subset of the main df that contains only the samples of the species in question
            spp_df = self.abundance_df.loc[sample_names_of_spp]
            spp_df = spp_df.loc[:, (spp_df != 0).any(axis=0)]
            return sample_names_of_spp, spp_df

        def _do_spp_pcoa_unifrac(self, sample_names_of_spp, spp, wu):
            # compute the pcoa
            pcoa_output = pcoa(wu.data)
            pcoa_output.samples['sample'] = sample_names_of_spp
            spp_unifrac_pcoa_df = pcoa_output.samples.set_index('sample')
            # now add the variance explained as a final row to the renamed_dataframe
            spp_unifrac_pcoa_df = spp_unifrac_pcoa_df.append(
                pcoa_output.proportion_explained.rename('proportion_explained'))
            pickle.dump(spp_unifrac_pcoa_df,
                        open(os.path.join(self.parent.cache_dir, f'spp_unifrac_pcoa_df_{spp}.p'), 'wb'))
            return spp_unifrac_pcoa_df

        def _make_and_root_tree(self):
            # make the tree
            print('Testing models and making phylogenetic tree')
            print('This could take some time...')
            if not os.path.exists(self.tree_out_path_unrooted):
                subprocess.run(
                    ['iqtree', '-nt', 'AUTO', '-m', 'TIM3e+R3', '-s', f'{self.aligned_fasta_path}'],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            else:
                print('Tree already exists. Using existing tree.')
            # root the tree
            print('Tree creation complete')
            print('Rooting the tree at midpoint')
            tree = TreeNode.read(self.tree_out_path_unrooted)
            self.rooted_tree = tree.root_at_midpoint()
            self.rooted_tree.write(self.tree_out_path_rooted)

        def _align_seqs(self, seq_fasta_list):
            general.write_list_to_destination(destination=self.unaligned_fasta_path, list_to_write=seq_fasta_list)
            # here we have a fasta ready for alignment
            if not os.path.exists(self.aligned_fasta_path):
                general.mafft_align_fasta(input_path=self.unaligned_fasta_path, output_path=self.aligned_fasta_path,
                                          method='auto', num_proc=6)

        def _set_df_cols_as_hashes(self, columns):
            # change the columns of the df
            self.abundance_df.columns = columns

        def _make_all_seq_fasta_and_hash_names(self):
            # here we have a set of all of the sequences
            seq_fasta_list = []
            # we will change the df columns so that they match the seq names in the tree
            columns = []
            for seq in list(self.abundance_df):
                hash_of_seq = hashlib.md5(seq.encode()).hexdigest()
                if hash_of_seq in columns:
                    sys.exit('non-unique hash')
                seq_fasta_list.extend([f'>{hash_of_seq}', seq])
                columns.append(hash_of_seq)
            return columns, seq_fasta_list

        def _create_df_from_minor_div_dict(self, minor_div_abundance_dict):
            self.abundance_df = pd.DataFrame.from_dict(minor_div_abundance_dict, orient='index')
            self.abundance_df[pd.isna(self.abundance_df)] = 0

        def _generate_bray_curtis_distance_and_pcoa_spp(self):
            """
            This is code for generating PCOAs for each of the coral species.
            Read in the minor div dataframe which should have normalised abundances in them
            For each sample we have a fasta that we can read in which has the normalised (to 1000) sequences
            For feeding into med. We can use a default dict to collect the sequences and abundances from this fairly
            simply.
            This is likely best done on for each sample outside of the pairwise comparison to save on redoing the same
            collection of the sequences.
            """

            if os.path.isfile(os.path.join(self.parent.cache_dir, 'minor_div_abundance_dict.p')):
                minor_div_abundance_dict = pickle.load(
                    open(os.path.join(self.parent.cache_dir, 'minor_div_abundance_dict.p'), 'rb'))

            else:
                minor_div_abundance_dict = self._generate_minor_div_abundance_dict_from_scratch()

            # For each of the spp.
            spp_pcoa_df_dict = {}
            for spp in ['Porites', 'Pocillopora', 'Millepora']:
                if os.path.isfile(os.path.join(self.parent.cache_dir, f'spp_pcoa_df_{spp}.p')):
                    spp_pcoa_df = pickle.load(open(os.path.join(self.parent.cache_dir, f'spp_pcoa_df_{spp}.p'), 'rb'))
                else:
                    # Get a list of the samples that we should be working with
                    sample_names_of_spp = self.parent.coral_info_df_for_figures.loc[
                        self.parent.coral_info_df_for_figures['genus'] == spp.upper()].index.values.tolist()

                    # remove the two porites species form this that seem to be total outliers
                    if spp == 'Porites':
                        sample_names_of_spp.remove('CO0001674')
                        sample_names_of_spp.remove('CO0001669')

                    spp_distance_dict = self._get_spp_sample_distance_dict(minor_div_abundance_dict,
                                                                           sample_names_of_spp, spp)

                    distance_out_file = self._make_and_write_out_spp_dist_file(sample_names_of_spp, spp_distance_dict, spp)

                    # Feed this into the generate_PCoA_coords method
                    spp_pcoa_df = self._generate_PCoA_coords(distance_out_file, spp)
                    pickle.dump(spp_pcoa_df, open(os.path.join(self.parent.cache_dir, f'spp_pcoa_df_{spp}.p'), 'wb'))
                spp_pcoa_df_dict[spp] = spp_pcoa_df
            return spp_pcoa_df_dict

        def _generate_PCoA_coords(self, raw_dist_file, spp):
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

        def _generate_minor_div_abundance_dict_from_scratch(self):
            # this dict will have sample name as key and then a dict as value with seq to abundance values
            # we can then work with this for the pairwise comparison
            minor_div_abundance_dict = {}
            for ind in self.parent.coral_info_df_for_figures.index.values.tolist():
                sample_dir = self.parent.coral_info_df_for_figures.loc[ind, 'sample_dir']
                with open('{}/fasta_for_med.fasta'.format(sample_dir), 'r') as f:
                    sample_fasta = [line.rstrip() for line in f]

                sample_minor_abundance_dict = defaultdict(int)
                for i in range(0, len(sample_fasta), 2):
                    sample_minor_abundance_dict[sample_fasta[i + 1]] += 1

                # here we have the dict popoulated for the sample
                # we can now add this to the minor_div_abundace_dict
                minor_div_abundance_dict[ind] = sample_minor_abundance_dict
            # we should now pickle out this sample_minor_abundance_dict
            pickle.dump(
                minor_div_abundance_dict, open(os.path.join(self.parent.cache_dir, 'minor_div_abundance_dict.p'), 'wb'))
            return minor_div_abundance_dict

        def _make_and_write_out_spp_dist_file(self, sample_names_of_spp, spp_distance_dict, spp):
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
                distance_out_file.append(
                    '\t'.join([str(distance_item) for distance_item in temp_sample_dist_string]))
            # from here we can hopefully rely on the rest of the methods as they already are. The .dist file should be
            # written out
            dist_out_path = os.path.join(self.parent.dist_output_dir, f'bray_curtis_within_spp_sample_distances_{spp}.dist')
            with open(dist_out_path, 'w') as f:
                for line in distance_out_file:
                    f.write('{}\n'.format(line))
            return distance_out_file

        def _get_spp_sample_distance_dict(self, minor_div_abundance_dict, sample_names_of_spp, spp):
            if os.path.isfile(os.path.join(self.parent.cache_dir, f'spp_distance_dict_{spp}.p')):
                spp_distance_dict = pickle.load(
                    open(os.path.join(self.parent.cache_dir, f'spp_distance_dict_{spp}.p'), 'rb'))

            else:
                spp_distance_dict = self._make_spp_sample_distance_dict_from_scratch(
                    minor_div_abundance_dict, sample_names_of_spp, spp)
            return spp_distance_dict

        def _make_spp_sample_distance_dict_from_scratch(self, minor_div_abundance_dict, sample_names_of_spp, spp):
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
            pickle.dump(spp_distance_dict,
                        open(os.path.join(self.parent.cache_dir, f'spp_distance_dict_{spp}.p'), 'wb'))
            return spp_distance_dict

    class GenericPlottingMethods:
        """A class that will hold the methods shared by all of the plotting methods"""

        @staticmethod
        def _remove_axes_but_allow_labels(ax, x_tick_label_list=None):
            ax.set_frame_on(False)
            if not x_tick_label_list:
                ax.set_xticks([])
            ax.set_yticks([])
            ax.minorticks_off()

        @staticmethod
        def _vert_leg_axis(colour_dict, legend_ax, marker_dict):
            legend_ax.set_ylim(0, 1)
            legend_ax.set_xlim(0, 1)
            legend_ax.invert_yaxis()
            number_of_icons = 6
            island_list = ['ISLAND06', 'ISLAND10', 'ISLAND15']
            site_list = ['SITE01', 'SITE02', 'SITE03']
            icon_list = []
            # first populate the island icons
            for site in site_list:
                icon_list.append((colour_dict[site], marker_dict['ISLAND06']))
            # then populate the site icons
            for island in island_list:
                icon_list.append((colour_dict['SITE01'], marker_dict[island]))
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
                legend_ax.scatter(x=x_val_for_icon, y=y_val_for_icon_and_text, c=icon_list[i][0],
                                  marker=icon_list[i][1],
                                  s=100)

                if int(i / 3) == 0:
                    legend_ax.text(s=site_list[i], x=x_val_for_text, y=y_val_for_icon_and_text)
                elif int(i / 3) == 1:
                    legend_ax.text(s=island_list[i % 3], x=x_val_for_text, y=y_val_for_icon_and_text)
                pos_counter += 1

    class E18S_ITS2_PCOA_FIGURE(Generic_PCOA_DIST_Methods, GenericPlottingMethods):
        """A class  for holding the methods specific to plotting the figure that links the 18S ordinations with the
        zooxs its2 information.
        """
        def __init__(self, parent, distance_method):
            EighteenSAnalysis.Generic_PCOA_DIST_Methods.__init__(self, parent)
            EighteenSAnalysis.GenericPlottingMethods.__init__(self)
            self.foo = 'asdf'
            self.distance_method = distance_method
            if self.distance_method == 'braycurtis':
                self.pcoa_df_dict = self._generate_bray_curtis_distance_and_pcoa_spp()
            elif self.distance_method == 'unifrac':
                self.pcoa_df_dict = self._generate_unifrac_distance_and_pcoa_spp()
            self.spp_list = ['Porites', 'Pocillopora', 'Millepora']
            self.marker_dict = {'ISLAND06': '^', 'ISLAND10': 'o', 'ISLAND15': 's'}

            self._setup_fig_layout()

            # These three objects are all populated in _process_div_df()
            self.smp_id_to_smp_name_dict = None
            self.smp_name_to_smp_id_dict = None
            self.sp_output_df_div = None
            self._process_div_df()
            self.smp_name_to_smp_id_dict_short = {k.split('_')[0]: v for k, v in self.smp_name_to_smp_id_dict.items()}

            try:
                self.colour_dict_type = pickle.load(
                    open('/home/humebc/projects/tara/initial_its2_processing/colour_dict_type.pickle',
                         'rb'))
            except FileNotFoundError:
                self.colour_dict_type = pickle.load(
                    open(os.path.join(self.parent.input_dir, 'colour_dict_type.pickle'), 'rb'))


            # These two objects are created in _process_type_df()
            self.sp_output_df_type = None
            self.sorted_type_prof_names_by_local_abund = None
            self._process_type_df()

            # Thse two objects are created in _get_div_colour_dict_and_ordered_list_of_seqs()
            # The ordered_list_of_seqs will be used for the plotting order
            self.ordered_list_of_seqs = None
            self.colour_dict_div = None
            self._get_div_colour_dict_and_ordered_list_of_seqs()

            self.ordered_sample_list = self.sp_output_df_type.index.values.tolist()

            # Reorder the columns and rows of the sp_output_df according to the sequence sample and sequence
            # order so that plotting the data is easier
            self.sp_output_df_div = self.sp_output_df_div[self.ordered_list_of_seqs]

        def _setup_fig_layout(self):
            # Setup the figure
            self.fig = plt.figure(figsize=(18, 10))
            self.gs = plt.GridSpec(8, 21, figure=self.fig, height_ratios=[3, 0.2, 0.5, 1, 0.1, 1, 0.1, 1],
                                   width_ratios=[0.2, 0.2, 1, 0.2, 1, 0.2, 1, 1.4, 1, 0.2, 1, 0.2, 1, 1.4, 1, 0.2, 1,
                                                 0.2, 1,
                                                 0.2, 0.4])
            # we can have several axes lists so that we can populate them one by one
            # first lets make the pcoa p lot axes list
            self.pcoa_axes_list = []
            # porites ax
            self.pcoa_axes_list.append(plt.subplot(self.gs[0, 2:7]))
            # pocillopora ax
            # test = plt.subplot(gs[0, 7])
            self.pcoa_axes_list.append(plt.subplot(self.gs[0, 8:13]))
            # millepora ax
            self.pcoa_axes_list.append(plt.subplot(self.gs[0, 14:19]))
            # now make the its2 matrix of plots
            self.its2_axes_list = []
            for start_y, start_x in [(3, 2), (3, 8), (3, 14)]:
                for i in range(0, 5, 2):
                    for j in range(0, 5, 2):
                        ind_x = j + start_x
                        ind_y = i + start_y
                        self.its2_axes_list.append(plt.subplot(self.gs[ind_y, ind_x]))
            self.max_n_cols_profiles = 4
            self.max_n_rows_profiles = 7
            self.num_leg_cells_profiles = self.max_n_cols_profiles * self.max_n_rows_profiles
            self.max_n_cols_seqs = 4
            self.max_n_rows_seqs = 7
            self.num_leg_cells_seqs = self.max_n_cols_seqs * self.max_n_rows_seqs

        def plot(self):

            # plotting of the PCOAs

            for spp in self.spp_list:

                ax = self.pcoa_axes_list[self.spp_list.index(spp)]

                df_of_spp = self.pcoa_df_dict[spp]

                # x values
                x_values = df_of_spp['PC1'].values.tolist()[:-1]

                # y values
                y_values = df_of_spp['PC2'].values.tolist()[:-1]

                samples_of_spp = df_of_spp.index.values.tolist()[:-1]
                # get list of colours and list of markers
                # colours can designate islands
                # we will need to look up the colour according to the its2 type profile designation of the sample

                # here we have to see which of the sample names in samples_of_spp are available in the sp_output_df_div
                island_colour_list = []
                for smpl_name in samples_of_spp:
                    try:
                        max_type = self.sp_output_df_type.loc[self.smp_name_to_smp_id_dict_short[smpl_name]].idxmax()
                        island_colour_list.append(self.colour_dict_type[max_type])
                    except:
                        island_colour_list.append('#000000')


                # shapes can designate sites

                marker_list = [self.marker_dict[self.parent.coral_info_df_for_figures.loc[smp, 'island']] for smp in samples_of_spp]

                # plot the points
                for x_val, y_val, col, mark in zip(x_values, y_values, island_colour_list, marker_list):
                    ax.scatter(x_val, y_val, c=col, marker=mark)

                # add axes labels
                ax.set_xlabel('PC1; explained = {}'.format('%.3f' % df_of_spp['PC1'][-1]))
                ax.set_ylabel('PC2; explained = {}'.format('%.3f' % df_of_spp['PC2'][-1]))
                # set axis title
                ax.set_title('{}'.format(spp), fontsize='large', fontweight='bold')

            self._add_labels_its2()

            # here we should start to take on the plotting of the ITS2 data
            self._plot_data_axes_its2()

            legend_ax = plt.subplot(self.gs[0, 20])
            self._remove_axes_but_allow_labels(legend_ax)
            self._vert_leg_axis_with_its2(legend_ax)

            # now put in a labelling for the its2 plots to show the its2 sequences, vs its2 type profiles
            for i in range(3, 8, 2):
                ax = plt.subplot(self.gs[i, 20])
                self._remove_axes_but_allow_labels(ax)
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.text(s='its2 seqs', x=0, y=0.5)
                ax.text(s='its2 type profiles', x=0, y=0)

            self.fig.show()

            plt.savefig(os.path.join(self.parent.figure_output_dir, f'spp_{self.distance_method}_pcoa_with_its2.png'), dpi=1200)
            plt.savefig(os.path.join(self.parent.figure_output_dir, f'spp_{self.distance_method}_pcoa_with_its2.svg'))

            return

        def _plot_data_axes_its2(self):
            ax_count = 0
            for spp in ['PORITES', 'POCILLOPORA', 'MILLEPORA']:
                for site in ['SITE01', 'SITE02', 'SITE03']:
                    for location in ['ISLAND06', 'ISLAND10', 'ISLAND15']:
                        ax = self.its2_axes_list[ax_count]
                        patches_list = []
                        ind = 0
                        colour_list = []

                        # for each set of location, site and spp, we basically want to get a list of the samples
                        # that meet the set criteria, we then want to plot samples according to the ordered_sample_list
                        # order which will be in IDs. As such we will have to convert the sample_name in the info_df
                        # to a sample ID using the smp_name_to_smp_id_dict.

                        # get sample_names that fit the requirements
                        sample_names_of_set = self.parent.all_samples_info_df.loc[
                            (self.parent.all_samples_info_df['location'] == location) &
                            (self.parent.all_samples_info_df['site'] == site) &
                            (self.parent.all_samples_info_df['spp_water'] == spp)
                            ].index.values.tolist()

                        if spp == 'PORITES':
                            if 'CO0001674' in sample_names_of_set:
                                sample_names_of_set.remove('CO0001674')
                                sample_names_of_set.remove('CO0001669')

                        # convert these to sample IDs
                        # The sample names in symportal are actually the full file names version rather than
                        # the shorter versions in the info_df. As such we should we will have to do a conversion here
                        smple_ids_of_set = []
                        for smp_name in sample_names_of_set:
                            smple_ids_of_set.append(self.smp_name_to_smp_id_dict_short[smp_name])

                        # now we want to plot in the order of the ordered_sample_list
                        ordered_smple_ids_of_set = [smpl_id for smpl_id in self.ordered_sample_list if
                                                    smpl_id in smple_ids_of_set]

                        num_smp_in_this_subplot = len(ordered_smple_ids_of_set)
                        x_tick_label_list = []

                        for smple_id_to_plot in ordered_smple_ids_of_set:
                            # General plotting
                            sys.stdout.write('\rPlotting sample: {}'.format(smple_id_to_plot))
                            x_tick_label_list.append(self.smp_id_to_smp_name_dict[smple_id_to_plot].split('_')[0])
                            # for each sample we will start at 0 for the y and then add the height of each bar to this

                            # PLOT DIVs
                            self._plot_div_over_type_its2(colour_list, ind, patches_list, smple_id_to_plot)

                            # PLOT type
                            self._plot_type_under_div_its2(colour_list, ind, patches_list, smple_id_to_plot)
                            ind += 1

                        self._paint_rect_to_axes_div_and_type_its2(
                            ax=ax, colour_list=colour_list, num_smp_in_this_subplot=num_smp_in_this_subplot,
                            patches_list=patches_list, max_num_smpls_in_subplot=10)

                        ax_count += 1

        def _plot_div_over_type_its2(self, colour_list, ind, patches_list, smple_id_to_plot):
            bottom_div = 0
            # for each sequence, create a rect patch
            # the rect will be 1 in width and centered about the ind value.
            for seq in list(self.sp_output_df_div):
                # class matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
                rel_abund_div = self.sp_output_df_div.loc[smple_id_to_plot, seq]
                if rel_abund_div > 0:
                    patches_list.append(
                        Rectangle((ind - 0.5, bottom_div), 1, rel_abund_div, color=self.colour_dict_div[seq]))
                    # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
                    colour_list.append(self.colour_dict_div[seq])
                    bottom_div += rel_abund_div

        def _plot_type_under_div_its2(self, colour_list, ind, patches_list, smple_id_to_plot):
            # the idea of the type is to put it as a reflection below the y=0 line
            # as such we should just want to make everything negative
            bottom_type = 0
            # for each sequence, create a rect patch
            # the rect will be 1 in width and centered about the ind value.
            # we want to plot the rects so that they add to 1. As such we want to divide
            # each value by the total for that sample.
            tot_for_sample = self.sp_output_df_type.loc[smple_id_to_plot].sum()
            for its2_profile in list(self.sp_output_df_type):
                rel_abund = self.sp_output_df_type.loc[smple_id_to_plot, its2_profile]
                if rel_abund > 0:
                    depth = -0.2 * (rel_abund / tot_for_sample)
                    patches_list.append(
                        Rectangle((ind - 0.5, bottom_type), 1, depth,
                                  color=self.colour_dict_type[its2_profile]))
                    # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
                    colour_list.append(self.colour_dict_type[its2_profile])
                    bottom_type += depth

        def _add_labels_its2(self):

            self.its2_axes_list[0].set_title('ISLAND06', fontsize='large', fontweight='bold')
            self.its2_axes_list[1].set_title('ISLAND10', fontsize='large', fontweight='bold')
            self.its2_axes_list[2].set_title('ISLAND15', fontsize='large', fontweight='bold')

            self.its2_axes_list[9].set_title('ISLAND06', fontsize='large', fontweight='bold')
            self.its2_axes_list[10].set_title('ISLAND10', fontsize='large', fontweight='bold')
            self.its2_axes_list[11].set_title('ISLAND15', fontsize='large', fontweight='bold')

            self.its2_axes_list[18].set_title('ISLAND06', fontsize='large', fontweight='bold')
            self.its2_axes_list[19].set_title('ISLAND10', fontsize='large', fontweight='bold')
            self.its2_axes_list[20].set_title('ISLAND15', fontsize='large', fontweight='bold')

            self.its2_axes_list[0].set_ylabel('SITE 1', fontsize='large', fontweight='bold')
            self.its2_axes_list[3].set_ylabel('SITE 2', fontsize='large', fontweight='bold')
            self.its2_axes_list[6].set_ylabel('SITE 3', fontsize='large', fontweight='bold')

        def _get_div_colour_dict_and_ordered_list_of_seqs(self):
            colour_palette_div = get_colour_list()
            grey_palette_div = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']
            # get a list of the sequences in order of their abundance and use this list to create the colour dict
            # the abundances can be got by simply summing up the columns making sure to ommit the last columns
            abundance_dict = {}
            for col in list(self.sp_output_df_div):
                abundance_dict[col] = sum(self.sp_output_df_div[col])
            # get the names of the sequences sorted according to their totalled abundance
            self.ordered_list_of_seqs = [x[0] for x in sorted(abundance_dict.items(), key=lambda x: x[1], reverse=True)]
            # create the colour dictionary that will be used for plotting by assigning a colour from the colour_palette
            # to the most abundant seqs first and after that cycle through the grey_pallette assigning colours
            # If we aer only going to have a legend that is cols x rows as shown below, then we should only use
            # that many colours in the plotting.

            colour_dict_div = {}
            for i in range(len(self.ordered_list_of_seqs)):
                if i < self.num_leg_cells_seqs:
                    colour_dict_div[self.ordered_list_of_seqs[i]] = colour_palette_div[i]
                else:
                    grey_index = i % len(grey_palette_div)
                    colour_dict_div[self.ordered_list_of_seqs[i]] = grey_palette_div[grey_index]
            self.colour_dict_div = colour_dict_div

        def _process_type_df(self):
            try:
                sp_output_df_type = pd.read_csv(self.parent.symportal_profile_output_relative_path, sep='\t', lineterminator='\n',
                                                skiprows=[0, 1, 2, 3, 5],
                                                header=None)
            except FileNotFoundError:
                sp_output_df_type = pd.read_csv(self.parent.symportal_profile_output_relative_path_input_dir, sep='\t',
                                                lineterminator='\n',
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
            index_to_drop_from = None
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
            self.sp_output_df_type = sp_output_df_type.set_index(keys='sample_id', drop=True).astype('float')


            # we will use the col headers as the its2 type profile order for plotting but we
            # we should colour according to the abundance of the its2 type profiles
            # as we don't want to run out of colours by the time we get to profiles that are very abundant.
            # The sorted_type_prof_names_by_local_abund object has the names of the its2 type profile in order of abundance
            # we will use the index order as the order of samples to plot
            # create the colour dictionary that will be used for plotting by assigning a colour from the colour_palette
            # to the most abundant seqs first and after that cycle through the grey_pallette assigning colours
            self.sorted_type_prof_names_by_local_abund = [a[0] for a in
                                                     sorted(type_profile_to_abund_tup_list, key=lambda x: x[1],
                                                            reverse=True)]

        def _process_div_df(self):
            try:
                sp_output_df = pd.read_csv(self.parent.symportal_seq_output_relative_path, sep='\t', lineterminator='\n', header=0,
                                           index_col=0)
            except FileNotFoundError:
                sp_output_df = pd.read_csv(self.parent.symportal_seq_output_relative_path_input_dir, sep='\t',
                                           lineterminator='\n', header=0,
                                           index_col=0)
            # In order to be able to drop the DIV row at the end and the meta information rows, we should
            # drop all rows that are after the DIV column. We will pass in an index value to the .drop
            # that is called here. To do this we need to work out which index we are working with
            meta_index_to_cut_from = None
            index_values_as_list = sp_output_df.index.values.tolist()
            for i in range(-1, -(len(index_values_as_list)), -1):
                if index_values_as_list[i].startswith('DIV'):
                    # then this is the index (in negative notation) that we need to cut from
                    meta_index_to_cut_from = i
                    break
            sp_output_df = sp_output_df.iloc[:meta_index_to_cut_from]

            # create sample id to sample name dict
            self.smp_id_to_smp_name_dict = {ID: '_'.join(nm.split('_')[:3]) for ID, nm in
                                       zip(sp_output_df.index.values.tolist(),
                                           sp_output_df['sample_name'].values.tolist())}
            self.smp_name_to_smp_id_dict = {'_'.join(nm.split('_')[:3]): ID for ID, nm in
                                       zip(sp_output_df.index.values.tolist(),
                                           sp_output_df['sample_name'].values.tolist())}

            # now lets drop the QC columns from the SP output df and also drop the clade summation columns
            # we will be left with just clumns for each one of the sequences found in the samples
            sp_output_df.drop(
                columns=['sample_name', 'noName Clade A', 'noName Clade B', 'noName Clade C', 'noName Clade D',
                         'noName Clade E', 'noName Clade F', 'noName Clade G', 'noName Clade H',
                         'noName Clade I', 'raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                         'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                         'post_taxa_id_absolute_non_symbiodinium_seqs',
                         'post_taxa_id_unique_non_symbiodinium_seqs',
                         'size_screening_violation_absolute', 'size_screening_violation_unique',
                         'post_med_absolute', 'post_med_unique'
                         ], inplace=True)
            self.sp_output_df_div = sp_output_df.astype('float')

        def _paint_rect_to_axes_div_and_type_its2(
                self, ax, colour_list, num_smp_in_this_subplot, patches_list,max_num_smpls_in_subplot=10):
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
            ax.set_ylim(-0.2, 1)
            # ax.set_xticks(range(num_smp_in_this_subplot))
            # ax.set_xticklabels(x_tick_label_list, rotation='vertical', fontsize=6)

            self._remove_axes_but_allow_labels(ax)

            # as well as getting rid of the top and right axis splines
            # I'd also like to restrict the bottom spine to where there are samples plotted but also
            # maintain the width of the samples
            # I think the easiest way to do this is to hack a bit by setting the x axis spines to invisible
            # and then drawing on a line at y = 0 between the smallest and largest ind (+- 0.5)
            # ax.spines['bottom'].set_visible(False)
            ax.add_line(Line2D((0 - 0.5, num_smp_in_this_subplot - 0.5), (0, 0), linewidth=2, color='black'))

        def _vert_leg_axis_with_its2(self, legend_ax):
            legend_ax.set_ylim(0, 1)
            legend_ax.set_xlim(0, 1)
            legend_ax.invert_yaxis()
            number_of_icons = 3
            island_list = ['ISLAND06', 'ISLAND10', 'ISLAND15']

            icon_list = []
            # first populate the island icons
            # for site in site_list:
            #     icon_list.append((colour_dict[site], marker_dict['ISLAND06']))
            # then populate the site icons
            for island in island_list:
                icon_list.append(('#808080', self.marker_dict[island]))
            # lets assume that the axis is divided into 20 spaces for icons
            max_number_icons = 10
            # the  icon position should be mid way so max_number_icons
            # first icon position
            first_icon_position = int((max_number_icons - number_of_icons) / 2)
            pos_counter = first_icon_position
            for i in range(len(icon_list)):
                y_val_for_icon_and_text = (1 / max_number_icons) * pos_counter
                x_val_for_icon = 0.4
                x_val_for_text = x_val_for_icon + 0.4
                legend_ax.scatter(x=x_val_for_icon, y=y_val_for_icon_and_text, c=icon_list[i][0], marker=icon_list[i][1],
                                  s=100)

                if int(i / 3) == 0:
                    legend_ax.text(s=island_list[i], x=x_val_for_text, y=y_val_for_icon_and_text)

                pos_counter += 1

    class PlotPCoASpp(Generic_PCOA_DIST_Methods, GenericPlottingMethods):

        def __init__(self, parent, distance_method):
            EighteenSAnalysis.GenericPlottingMethods.__init__(self)
            EighteenSAnalysis.Generic_PCOA_DIST_Methods.__init__(self, parent=parent)
            self.colour_dict = {'SITE01': '#C0C0C0', 'SITE02': '#808080', 'SITE03': '#000000'}
            self.marker_dict = {'ISLAND06': '^', 'ISLAND10': 'o', 'ISLAND15': 's'}
            self.distance_method = distance_method

        def plot(self, is_three_d=False):
            """
            This is the code for producing the 18s pcoa with either the 3rd and 4th PC underneath or a 3d graph.
            """

            # For each species, get the pcoa df
            if self.distance_method == 'braycurtis':
                pcoa_df_dict = self._generate_bray_curtis_distance_and_pcoa_spp()
            elif self.distance_method == 'unifrac':
                pcoa_df_dict = self._generate_unifrac_distance_and_pcoa_spp()

            # setup figure
            spp_list = ['Porites', 'Pocillopora', 'Millepora']
            axarr = []
            fig = plt.figure(figsize=(18, 10))

            gs = plt.GridSpec(3, 6, figure=fig, width_ratios=[1, 0.2, 1, 0.2, 1, 0.5], height_ratios=[1, 0.2, 1])
            for j in range(0, 3, 2):
                for i in range(0, 5, 2):
                    if is_three_d:
                        if j == 1:
                            axarr.append(fig.add_subplot(gs[j, i], projection='3d'))
                        else:
                            axarr.append(plt.subplot(gs[j, i]))
                    else:
                        axarr.append(plt.subplot(gs[j, i]))

            legend_ax = plt.subplot(gs[:, 5])
            self._remove_axes_but_allow_labels(legend_ax)
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

                colour_list = [self.colour_dict[self.parent.coral_info_df_for_figures.loc[smp, 'site']] for smp in samples_of_spp]

                # shapes can designate sites


                marker_list = [self.marker_dict[self.parent.coral_info_df_for_figures.loc[smp, 'island']] for smp in samples_of_spp]

                # plot the points
                for x_val, y_val, col, mark in zip(x_values, y_values, colour_list, marker_list):
                    ax.scatter(x_val, y_val, c=col, marker=mark)

                # add axes labels
                ax.set_xlabel('PC1; explained = {}'.format('%.3f' % df_of_spp['PC1'][-1]))
                ax.set_ylabel('PC2; explained = {}'.format('%.3f' % df_of_spp['PC2'][-1]))
                # set axis title
                ax.set_title('{}'.format(spp))
                if is_three_d:
                    # then lets plot the 3d equivalent below the 3d figs
                    # plot the points
                    for x_val, y_val, z_val, col, mark in zip(x_values, y_values, z_values, colour_list, marker_list):
                        ax_second.scatter(x_val, y_val, z_val, c=col, marker=mark)

                    # add axes labels
                    ax_second.set_xlabel('PC1; explained = {}'.format('%.3f' % df_of_spp['PC1'][-1]))
                    ax_second.set_ylabel('PC2; explained = {}'.format('%.3f' % df_of_spp['PC2'][-1]))
                    ax_second.set_zlabel('PC3; explained = {}'.format('%.3f' % df_of_spp['PC3'][-1]))


                else:
                    # else lets just plot out the 3rd and 4th PCs below
                    # plot the points
                    for z_val, pc4_val, col, mark in zip(z_values, pc4_values, colour_list, marker_list):
                        ax_second.scatter(z_val, pc4_val, c=col, marker=mark)

                    # add axes labels
                    ax_second.set_xlabel('PC3; explained = {}'.format('%.3f' % df_of_spp['PC3'][-1]))
                    ax_second.set_ylabel('PC4; explained = {}'.format('%.3f' % df_of_spp['PC4'][-1]))

            # here we should put together the legend axis
            self._vert_leg_axis(self.colour_dict, legend_ax, self.marker_dict)

            fig.show()
            if not is_three_d:
                plt.savefig(os.path.join(self.parent.figure_output_dir, f'spp_{self.distance_method}_pcoa_with_pc3_pc4.png'), dpi=1200)
                plt.savefig(os.path.join(self.parent.figure_output_dir, f'spp_{self.distance_method}_pcoa_with_pc3_pc4.svg'))

    class PlotPCoASppIsland(Generic_PCOA_DIST_Methods, GenericPlottingMethods):

        def __init__(self, parent, distance_method):
            EighteenSAnalysis.GenericPlottingMethods.__init__(self)
            EighteenSAnalysis.Generic_PCOA_DIST_Methods.__init__(self, parent=parent)
            self.coral_info_df_for_figures = self.parent.coral_info_df_for_figures
            self.cache_dir = self.parent.cache_dir
            self.figure_output_dir = self.parent.figure_output_dir
            self.distance_method = distance_method

        def plot(self):
            """
            This is the code for producing a pcoa per island per species.
            """

            # For each species get the pcoa df
            if self.distance_method == 'braycurtis':
                pcoa_df_dict = self._generate_bray_curtis_distance_and_pcoa_spp_and_island()
            elif self.distance_method == 'unifrac':
                pcoa_df_dict = self._generate_unifrac_distance_and_pcoa_spp_and_island()

            # setup figure
            spp_list = ['PORITES', 'POCILLOPORA', 'MILLEPORA']

            island_list = ['ISLAND06', 'ISLAND10', 'ISLAND15']

            axarr = []

            fig = plt.figure(figsize=(18, 10))

            gs = plt.GridSpec(5, 7, figure=fig, width_ratios=[1, 0.2, 1, 0.2, 1, 0.2, 0.5],
                              height_ratios=[1, 0.2, 1, 0.2, 1])
            for j in range(0, 5, 2):
                for i in range(0, 5, 2):
                    axarr.append(plt.subplot(gs[j, i]))

            legend_ax = plt.subplot(gs[:, 6])
            self._remove_axes_but_allow_labels(legend_ax)
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
                    colour_dict = {'SITE01': '#C0C0C0', 'SITE02': '#808080', 'SITE03': '#000000'}
                    colour_list = [colour_dict[self.coral_info_df_for_figures.loc[smp, 'site']] for smp in samples_of_spp]

                    # shapes can designate sites

                    marker_dict = {'ISLAND06': '^', 'ISLAND10': 'o', 'ISLAND15': 's'}
                    marker_list = [marker_dict[self.coral_info_df_for_figures.loc[smp, 'island']] for smp in samples_of_spp]

                    # plot the points
                    for x_val, y_val, col, mark in zip(x_values, y_values, colour_list, marker_list):
                        ax.scatter(x_val, y_val, c=col, marker=mark)

                    # add axes labels
                    ax.set_xlabel('PC1; explained = {}'.format('%.3f' % df_of_spp_island['PC1'][-1]))
                    ax.set_ylabel('PC2; explained = {}'.format('%.3f' % df_of_spp_island['PC2'][-1]))
                    # set axis title
                    if ax_count in [0, 1, 2]:
                        ax.set_title('{}'.format(island), fontsize='large', fontweight='bold')
                    if ax_count in [2, 5, 8]:
                        ax2 = ax.twinx()
                        ax2.yaxis.set_ticklabels([])
                        ax2.yaxis.set_ticks_position('none')
                        ax2.set_ylabel(spp, fontsize='large', fontweight='bold')

                    ax_count += 1

            # here we should put together the legend axis
            self._vert_leg_axis(colour_dict, legend_ax, marker_dict)

            fig.show()
            plt.savefig(os.path.join(self.figure_output_dir, f'spp_island_{self.distance_method}_pcoa.png'), dpi=1200)
            plt.savefig(os.path.join(self.figure_output_dir, f'spp_island_{self.distance_method}_pcoa.svg'))

        def _generate_unifrac_distance_and_pcoa_spp_and_island(self):
            if os.path.isfile(os.path.join(self.parent.cache_dir, 'minor_div_abundance_dict.p')):
                minor_div_abundance_dict = pickle.load(
                    open(os.path.join(self.parent.cache_dir, 'minor_div_abundance_dict.p'), 'rb'))
            else:
                minor_div_abundance_dict = self._generate_minor_div_abundance_dict_from_scratch()

            self._create_df_from_minor_div_dict(minor_div_abundance_dict)

            columns, seq_fasta_list = self._make_all_seq_fasta_and_hash_names()

            self._set_df_cols_as_hashes(columns)

            self._align_seqs(seq_fasta_list)

            self._make_and_root_tree()

            spp_island_unifrac_pcoa_df_dict = {}
            for spp in ['PORITES', 'POCILLOPORA', 'MILLEPORA']:
                for island in ['ISLAND06', 'ISLAND10', 'ISLAND15']:
                    if os.path.isfile(os.path.join(self.cache_dir, f'spp_island_unifrac_pcoa_df_{spp}_{island}.p')):
                        spp_island_unifrac_pcoa_df = pickle.load(
                            open(os.path.join(self.cache_dir, f'spp_island_unifrac_pcoa_df_{spp}_{island}.p'), 'rb'))
                    else:
                        spp_island_unifrac_pcoa_df = self._make_spp_island_unifrac_pcoa_df_from_scratch(spp, island)
                    spp_island_unifrac_pcoa_df_dict[f'{spp}_{island}'] = spp_island_unifrac_pcoa_df
            return spp_island_unifrac_pcoa_df_dict

        def _generate_bray_curtis_distance_and_pcoa_spp_and_island(self):
            """Read in the minor div dataframe which should have normalised abundances in them
            # For each sample we have a fasta that we can read in which has the normalised (to 1000) sequences
            # For feeding into med. We can use a default dict to collect the sequences and abundances from this fairly
            # simply.
            # This is likely best done on for each sample outside of the pairwise comparison to save on redoing the same
            # collection of the sequences.
            """

            if os.path.isfile(os.path.join(self.cache_dir, 'minor_div_abundance_dict.p')):
                minor_div_abundance_dict = pickle.load(
                    open(os.path.join(self.cache_dir, 'minor_div_abundance_dict.p'), 'rb'))

            else:
                minor_div_abundance_dict = self._generate_minor_div_abundance_dict_from_scratch()

            # For each of the spp.
            spp_island_pcoa_df_dict = {}
            for spp in ['PORITES', 'POCILLOPORA', 'MILLEPORA']:
                for island in ['ISLAND06', 'ISLAND10', 'ISLAND15']:
                    if os.path.isfile(os.path.join(self.cache_dir, f'spp_island_pcoa_df_{spp}_{island}.p')):
                        spp_island_pcoa_df = pickle.load(
                            open(os.path.join(self.cache_dir, f'spp_island_pcoa_df_{spp}_{island}.p'), 'rb'))
                    else:
                        spp_island_pcoa_df = self._make_spp_island_braycurtis_pcoa_df_from_scratch(island,
                                                                                                   minor_div_abundance_dict,
                                                                                                   spp)
                    spp_island_pcoa_df_dict['{}_{}'.format(spp, island)] = spp_island_pcoa_df
            return spp_island_pcoa_df_dict

        def _make_spp_island_unifrac_pcoa_df_from_scratch(self, spp, island):
            # Get a list of the samples that we should be working with
            sample_names_of_spp_island = self.coral_info_df_for_figures.loc[
                (self.coral_info_df_for_figures['genus'] == spp.upper()) &
                (self.coral_info_df_for_figures['island'] == island.upper())].index.values.tolist()
            # Remove the two porites species from this that seem to be total outliers
            if spp == 'PORITES':
                if 'CO0001674' in sample_names_of_spp_island:
                    sample_names_of_spp_island.remove('CO0001674')
                    sample_names_of_spp_island.remove('CO0001669')

            # This is a subset of the main df that contains only the samples of the species in question
            spp_island_df = self.abundance_df.loc[sample_names_of_spp_island]
            spp_island_df = spp_island_df.loc[:, (spp_island_df != 0).any(axis=0)]

            # perform unifrac
            print('Performing unifrac calculations')
            wu = beta_diversity(
                metric='weighted_unifrac', counts=spp_island_df.to_numpy(),
                ids=[str(_) for _ in list(spp_island_df.index)],
                tree=self.rooted_tree, otu_ids=[str(_) for _ in list(spp_island_df.columns)])

            # compute the pcoa
            pcoa_output = pcoa(wu.data)
            pcoa_output.samples['sample'] = sample_names_of_spp_island
            spp_unifrac_pcoa_df = pcoa_output.samples.set_index('sample')
            # now add the variance explained as a final row to the renamed_dataframe
            spp_unifrac_pcoa_df = spp_unifrac_pcoa_df.append(
                pcoa_output.proportion_explained.rename('proportion_explained'))
            pickle.dump(spp_unifrac_pcoa_df,
                        open(os.path.join(self.parent.cache_dir, f'spp_island_unifrac_pcoa_df_{spp}_{island}.p'), 'wb'))
            return spp_unifrac_pcoa_df

        def _make_spp_island_braycurtis_pcoa_df_from_scratch(self, island, minor_div_abundance_dict, spp):
            # Get a list of the samples that we should be working with
            sample_names_of_spp_island = self.coral_info_df_for_figures.loc[
                (self.coral_info_df_for_figures['genus'] == spp.upper()) &
                (self.coral_info_df_for_figures['island'] == island.upper())].index.values.tolist()
            # Remove the two porites species from this that seem to be total outliers
            if spp == 'PORITES':
                if 'CO0001674' in sample_names_of_spp_island:
                    sample_names_of_spp_island.remove('CO0001674')
                    sample_names_of_spp_island.remove('CO0001669')
            spp_island_distance_dict = self._get_spp_island_distance_dict(
                minor_div_abundance_dict, sample_names_of_spp_island, spp, island)
            distance_out_file = self._make_and_output_distance_file_spp_island(
                sample_names_of_spp_island, spp_island_distance_dict, spp, island)
            # Feed this into the generate_PCoA_coords method
            spp_island_pcoa_df = self._generate_PCoA_coords(distance_out_file, spp)
            pickle.dump(spp_island_pcoa_df,
                        open(os.path.join(self.cache_dir, f'spp_island_pcoa_df_{spp}_{island}.pickle'), 'wb'))
            return spp_island_pcoa_df

        def _make_and_output_distance_file_spp_island(self, sample_names_of_spp, spp_island_distance_dict, spp, island):
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
                distance_out_file.append(
                    '\t'.join([str(distance_item) for distance_item in temp_sample_dist_string]))
            # from here we can hopefully rely on the rest of the methods as they already are. The .dist file should be
            # written out
            dist_out_path = os.path.join \
                (self.parent.dist_output_dir, f'bray_curtis_within_spp_island_sample_distances_{spp}_{island}.dist')
            with open(dist_out_path, 'w') as f:
                for line in distance_out_file:
                    f.write('{}\n'.format(line))
            return distance_out_file

        def _get_spp_island_distance_dict(self, minor_div_abundance_dict, sample_names_of_spp, spp, island):
            if os.path.isfile(os.path.join(self.cache_dir, f'spp_island_distance_dict_{spp}_{island}.p')):
                spp_island_distance_dict = pickle.load(
                    open(os.path.join(self.cache_dir, f'spp_island_distance_dict_{spp}_{island}.p'), 'rb'))

            else:
                spp_island_distance_dict = self._generate_spp_island_distance_dict_from_scratch(
                    minor_div_abundance_dict, sample_names_of_spp, spp, island)
            return spp_island_distance_dict

        def _generate_spp_island_distance_dict_from_scratch(self, minor_div_abundance_dict, sample_names_of_spp, spp, island):
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
            pickle.dump(
                spp_island_distance_dict,
                open(os.path.join(self.cache_dir, f'spp_island_distance_dict_{spp}_{island}.p'), 'wb'))
            return spp_island_distance_dict

    class TaxaAnnotation:

        """Methods associated with doing the blasts to generate taxa annotations for sequences in samples"""
        def __init__(self, parent):
            self.parent = parent

        def annotate(self, numProc=20):
            """The sequence_QC method has taken care of the basic mothur qc for us. The next step will be to run a blast
            analysis for each of the samples that we have
            In each directory we should have a stability.trim.contigs.good.unique.abund.pcr.fasta paired with a
            stability.trim.contigs.good.abund.pcr.names. We should work with these pairs and produce a dictionary
            that will give sequence name to taxonomic id.
            We are only interested in the coral samples.
            As such we should filter these out when doing the MPing.
            """

            ncbircFile = []
            ncbircFile.extend(["[BLAST]", "BLASTDB={}".format(self.parent.blast_nt_db_path)])

            # write out the ncbircFile
            with open(os.path.join(self.parent.root_dir, '.ncbirc'), 'w') as f:
                for line in ncbircFile:
                    f.write('{}\n'.format(line))

            input_q_for_blast = Queue()

            for ind in self.parent.all_samples_info_df.index.values.tolist():
                if 'CORAL' in self.parent.all_samples_info_df.loc[ind, 'fastq_fwd_file_path']:
                    input_q_for_blast.put(self.parent.all_samples_info_df.loc[ind, 'fastq_fwd_file_path'])

            for n in range(numProc):
                input_q_for_blast.put('STOP')

            # generate the dictionaries used for taxa_designation
            # send a copy to each of the processes
            node_dict, name_dict = self._generate_taxa_designation_from_staxid_dicts()

            all_procs = []
            for n in range(numProc):
                p = Process(target=self._blast_worker, args=(input_q_for_blast, node_dict, name_dict))
                all_procs.append(p)
                p.start()

            for p in all_procs:
                p.join()

            return

        def _generate_taxa_designation_from_staxid_dicts(self,
                taxdump_dir='/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/ntdbdownload/taxdump'):
            # read in the .nodes file. This file tells us which tax level the node is and which node is the parent level
            with open('{}/nodes.dmp'.format(taxdump_dir), 'r') as f:
                node_file = [line.rstrip() for line in f]
            # now make a dict from this where key is the tax id and the value is a tup where 0 = parent 1 = tax level
            node_dict = {line.split('\t|\t')[0]: (line.split('\t|\t')[1], line.split('\t|\t')[2]) for line in node_file}

            # next read in the names file. This file hold the name of the node.
            with open('{}/names.dmp'.format(taxdump_dir), 'r') as f:
                name_file = [line.rstrip() for line in f]

            # now make a dict from the names file where the key is the staxid and the value is the name
            name_dict = {line.split('\t|\t')[0]: line.split('\t|\t')[1] for line in name_file if
                         line.split('\t|\t')[3].replace('\t|', '') == 'scientific name'}

            return node_dict, name_dict

        def _get_taxa_designation_from_staxid(self, staxid, tax_level_list, node_dict=None, name_dict=None,
                                             taxdump_dir='/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/ntdbdownload/taxdump'):
            # I have set this up so that you can either feed in the already made node_dict and name_dict or
            # you they can be generated by this method. If you are running this method many times it will make
            # sense to run the above generate_taxa_designation_from_staxid_dicts methods
            # to generate the dicts to feed into this
            # This will take in a list as its argument and find the levels that are listed in the list
            list_to_return = [False for i in tax_level_list]
            if node_dict == None or name_dict == None:
                node_dict, name_dict = self._generate_taxa_designation_from_staxid_dicts()

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

        def _blast_worker(self, input_q, node_dict, name_dict):
            # this is the mp worker. We will aim to go to each of the directories and perform a blast and
            # and produce a dictionary that will be sequence id to taxonomic designation

            for path in iter(input_q.get, 'STOP'):
                sys.stdout.write('\nSample {}\n'.format(path.split('/')[-1]))

                sample_dir = '/'.join(path.split('/')[:-1])
                # we only need to be doing this sample if it hasn't already been done.
                if os.path.isfile('{}/sample_tax_dict.pickle'.format(sample_dir)) and os.path.isfile(
                        '{}/coral_dict.pickle'.format(sample_dir)) and os.path.isfile(
                        '{}/symbiodiniaceae_dict.pickle'.format(sample_dir)):
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
                            genus_level, family_level, order_level = self._get_taxa_designation_from_staxid(
                                staxid=comp_list[6], tax_level_list=['genus', 'family', 'order'],
                                node_dict=node_dict, name_dict=name_dict)
                        except:
                            # there are some tax_ids that we don't seem to be able to find
                            # for the time being we will just continue over this seqeunce
                            continue
                        if False in [genus_level, family_level, order_level]:
                            continue

                        if family_level == 'Symbiodiniaceae':
                            symbiodiniaceae_count += 1
                            symbiodiniaceae_genus_dict[genus_level] += 1

                        elif order_level == 'Scleractinia' or order_level == 'Anthoathecata':
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
                        sample_tax_dict[blast_key] = \
                        [a[0] for a in sorted(other_genus_count_dict.items(), key=lambda x: x[1], reverse=True)][0]
                    elif coral_count > max(max_other_taxa_count, symbiodiniaceae_count):
                        # then this is a scleratinian and we should store the order
                        sample_tax_dict[blast_key] = 'Scleractinia'
                        # for the scleractinian dictionary we should associate to the most abundant genus
                        most_abundant_genus_sclerac = \
                        [a[0] for a in sorted(coral_genus_dict.items(), key=lambda x: x[1], reverse=True)][0]
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

    class SeqStackedBarPlotter(GenericPlottingMethods):
        """This method produced stacked bar charts. It can produce different charts depending on the plot_type.
        'full' means all of the sequences including the maj
        'low' means without the maj
        'med' means without the maj and having been through med
        'qc_taxa_rel_abund' means plot the relative abundance of the post qc taxa categories
        'qc_absolute means' plot the absolute abundance of the post-qc sequences (all tax categories)"""

        def __init__(self, parent, plot_type):
            EighteenSAnalysis.GenericPlottingMethods.__init__(self)
            self.plot_type = plot_type
            self.parent = parent
            self.ax_list, self.fig = self._setup_axes()
            self.sample_abundance_df = None
            # This df is only used for MED plotting
            self.minor_intra_med_abund_df = None
            self.colour_dict = None

        def _setup_axes(self):
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
                                                                              subplot_spec=gs[
                                                                                  row_ind + increaser, col_ind])
                    grid_spec_subplot_list.append(temp_grid_spec_subplot)
                    for i in range(3):
                        # NB this might be a 2d array, lets see.
                        ax = plt.Subplot(fig, temp_grid_spec_subplot[i])
                        ax_list.append(ax)
                        fig.add_subplot(ax)
                # now put in the spacer row
                if increaser < 2:
                    ax_space = plt.subplot(gs[row_ind + 1 + increaser, :])
                    self._remove_axes_but_allow_labels(ax_space)
                    increaser += 1
            return ax_list, fig


        def plot(self):
            if self.plot_type == 'full' or self.plot_type == 'low':
                self.sample_abundance_df = self._generate_seq_abundance_df()
                # bob = sample_abundance_df.loc['CO0001677'].idxmax()
                self.colour_dict = self._generate_colour_dict(is_med=False, is_qc=False)
            elif self.plot_type == 'med':
                self.sample_abundance_df = self._generate_seq_abundance_df()
                self.minor_intra_med_abund_df = self._create_med_node_sample_abundance_df_for_minor_intras()
                self.colour_dict = self._generate_colour_dict(is_med=True, is_qc=False)
            elif 'qc' in self.plot_type:
                self.sample_abundance_df = self._generate_post_qc_taxa_df()
                self.colour_dict = self._generate_colour_dict(is_med=True, is_qc=True)

            if self.plot_type == 'low':
                self._plot_data_axes_18s(minor_DIV=True)
            elif self.plot_type == 'full' or self.plot_type == 'med':
                self._plot_data_axes_18s(minor_DIV=False)
            elif 'qc' in self.plot_type:
                self._plot_data_axes_18s(minor_DIV=False,
                                   qc=True)

            self._add_labels(self.ax_list)

            if self.plot_type == 'full':
                plt.savefig(os.path.join(self.parent.figure_output_dir, 'raw_seqs_abund_stacked.png'))
                plt.savefig(os.path.join(self.parent.figure_output_dir, 'raw_seqs_abund_stacked.svg'))
            elif self.plot_type == 'low':
                plt.savefig(os.path.join(self.parent.figure_output_dir, 'raw_seqs_abund_minor_div_only_stacked.png'))
                plt.savefig(os.path.join(self.parent.figure_output_dir, 'raw_seqs_abund_minor_div_only_stacked.svg'))
            elif self.plot_type == 'med':
                plt.savefig(os.path.join(self.parent.figure_output_dir, 'raw_seqs_abund_minor_div_only_stacked_with_med.png'))
                plt.savefig(os.path.join(self.parent.figure_output_dir, 'raw_seqs_abund_minor_div_only_stacked_with_med.svg'))
            elif self.plot_type == 'qc_taxa_rel_abund':
                plt.savefig(os.path.join(self.parent.figure_output_dir, 'post_qc_taxa_rel_abund.png'))
                plt.savefig(os.path.join(self.parent.figure_output_dir, 'post_qc_taxa_rel_abund.svg'))
            elif self.plot_type == 'qc_absolute':
                plt.savefig(os.path.join(self.parent.figure_output_dir, 'post_qc_absolute.png'))
                plt.savefig(os.path.join(self.parent.figure_output_dir, 'post_qc_absolute.svg'))

        def _create_med_node_sample_abundance_df_for_minor_intras(self):
            """
            # This code will create MED node profiles for each of the samples, disregarding the most abundant sequence
            # NB after creating the MED sequences it was not much different to the raw sequences. This is likely becauase
            # all of the sequences were found at such low and even abundances. I think, moving forwards we should just stick with
            # working with the raw sequencs rather than the MED nodes.

            # We should go sample by sample using the abundance_df and write out a fasta for the MED
            # this will be a subsample of 1000 sequences of the minor DIVs, which will be run with a dynamic
            # m value on MED.
            # we will then read in the MED and create a dict for each sample that is MED node sequence to relabund
            # we will then come out of MP and go serial to collect all of the MED nodes and get a cumulative abundance
            # across all samples. We will order this and use the MED nodes as columns for the new minor_intras_abund_df
            # we will then go back through all of the samples in serial and populate the MED df
            # at this point we can pickle this df out and we can put it directly into plotting code of below
            # to create the same minor intra figure that we had before but with having done the MED.
            """

            if os.path.isfile(os.path.join(self.parent.cache_dir, 'minor_intra_med_abund_df.p')):
                minor_intra_med_abund_df = pickle.load(
                    open(os.path.join(self.parent.cache_dir, 'minor_intra_med_abund_df.p'), 'rb'))
            else:
                # MP sample by sample to do the MED
                input_q = Queue()

                # put in a tup of series. The first series being the info, the second being the current raw seq abundance df series
                for ind in self.parent.coral_info_df_for_figures.index.values.tolist():
                    input_q.put((self.parent.coral_info_df_for_figures.loc[ind], self.sample_abundance_df.loc[ind]))

                numProc = 20
                for n in range(numProc):
                    input_q.put('STOP')

                all_procs = []
                for n in range(numProc):
                    p = Process(target=self._med_worker, args=(input_q,))
                    all_procs.append(p)
                    p.start()

                for p in all_procs:
                    p.join()

                # here we have pickled out dictionaries of the MED nodes and their abundnaaces in each of the sample direcotires
                # these abundances are relative abundances
                # now have a master default dict and collect the abundances of the med nodes
                master_abund_med_dict = defaultdict(float)
                for ind in self.parent.coral_info_df_for_figures.index.values.tolist():
                    print('adding med nodes to master dict for {}'.format(ind))
                    sample_dir = self.parent.coral_info_df_for_figures.loc[ind, 'sample_dir']

                    # read in the pickled abundance dict
                    sample_med_abundance_dict = pickle.load(
                        open('{}/sample_node_dict_rel_abund.pickle'.format(sample_dir), 'rb'))
                    for seq_key, seq_rel_abund in sample_med_abundance_dict.items():
                        master_abund_med_dict[seq_key] += seq_rel_abund

                # here we have the master dict updated with all of the med node sequences
                # sort and get list of med_node_seqs
                sorted_med_node_seqs_tups = sorted(master_abund_med_dict.items(), key=lambda x: x[1], reverse=True)
                med_node_seqs_sorted = [a[0] for a in sorted_med_node_seqs_tups]

                # now create a df that will hold all of the sample abundances
                minor_intra_med_abund_df = pd.DataFrame(columns=med_node_seqs_sorted)

                # now go back through the samples again and populate the df by making a temp series
                for ind in self.parent.coral_info_df_for_figures.index.values.tolist():
                    print('populating minor_intra_med_abund_df with {}'.format(ind))
                    sample_dir = self.parent.coral_info_df_for_figures.loc[ind, 'sample_dir']
                    sample_med_abundance_dict = pickle.load(
                        open('{}/sample_node_dict_rel_abund.pickle'.format(sample_dir), 'rb'))
                    sample_data_in_order = []
                    for seq in med_node_seqs_sorted:
                        if seq in sample_med_abundance_dict.keys():
                            sample_data_in_order.append(sample_med_abundance_dict[seq])
                        else:
                            sample_data_in_order.append(float(0))

                    minor_intra_med_abund_df = minor_intra_med_abund_df.append(
                        pd.Series(name=ind, data=sample_data_in_order, index=med_node_seqs_sorted))

                # here we have the minor_intra_med_abund_df populated and we can pickle it out
                pickle.dump(minor_intra_med_abund_df,
                            open(os.path.join(self.parent.cache_dir, 'minor_intra_med_abund_df.p'), 'wb'))
            return minor_intra_med_abund_df

        def _med_worker(self, input_q):
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
                pathToFile = None
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
                    sample_node_dict[node_file[i + 1]] = int(node_file[i].split(':')[1])

                # now do the rel abunds again and then pickle out
                tot = sum(sample_node_dict.values())
                sample_node_dict_rel_abund = {k: v / tot for k, v in sample_node_dict.items()}

                pickle.dump(sample_node_dict_rel_abund,
                            open('{}/sample_node_dict_rel_abund.pickle'.format(sample_dir), 'wb'))

            return

        def _generate_post_qc_taxa_df(self):

            if os.path.isfile(os.path.join(self.parent.cache_dir, 'tax_absolute_abund_df.p')):
                tax_absolute_abund_df = pickle.load(open(os.path.join(self.parent.cache_dir, 'tax_absolute_abund_df.p'), 'rb'))
            else:
                input_q = Queue()

                for ind in self.parent.figure_output_dir.index.values.tolist():
                    input_q.put((ind, self.parent.figure_output_dir.loc[ind, 'sample_dir']))

                numProc = 20
                for n in range(numProc):
                    input_q.put('STOP')

                all_procs = []
                for n in range(numProc):
                    p = Process(target=self._qc_abundance_worker, args=(input_q,))
                    all_procs.append(p)
                    p.start()

                for p in all_procs:
                    p.join()

                # at this point we have the sample_tax_category_dict items for each sample pickled out
                # we can now go through each of the sample_dir again and use them to populate the master df
                # which we will then return
                tax_columns = ['Porites', 'Millepora', 'Pocillopora', 'other_coral', 'Symbiodiniaceae', 'other_taxa']
                tax_absolute_abund_df = pd.DataFrame(columns=tax_columns)
                for ind in self.parent.figure_output_dir.index.values.tolist():
                    sample_dir = self.parent.figure_output_dir.loc[ind, 'sample_dir']

                    # load the sample_tax_category_dict
                    sample_tax_category_dict = pickle.load(
                        open('{}/sample_tax_category_dict.pickle'.format(sample_dir), 'rb'))

                    sample_data = []
                    for cat in tax_columns:
                        if cat in sample_tax_category_dict.keys():
                            sample_data.append(sample_tax_category_dict[cat])
                        else:
                            sample_data.append(0)
                    tax_absolute_abund_df = tax_absolute_abund_df.append(
                        pd.Series(data=sample_data, index=tax_columns, name=ind))

                # here we have the tax_absolute_abund_df populated
                # now pickle it out
                pickle.dump(tax_absolute_abund_df, open(os.path.join(self.parent.cache_dir, 'tax_absolute_abund_df.p'), 'wb'))
            return tax_absolute_abund_df

        def _qc_abundance_worker(self, input_q):
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
                pickle.dump(sample_tax_category_dict,
                            open('{}/sample_tax_category_dict.pickle'.format(sample_dir), 'wb'))

        def _plot_data_axes_18s(self, minor_DIV=False, qc=False):
            ax_count = 0
            for site in ['SITE01', 'SITE02', 'SITE03']:
                for location in ['ISLAND06', 'ISLAND10', 'ISLAND15']:
                    for spp in ['PORITES', 'POCILLOPORA', 'MILLEPORA']:
                        ax = self.ax_list[ax_count]
                        patches_list = []
                        ind = 0
                        colour_list = []

                        # for each set of location, site and spp, we basically want to get a list of the samples
                        # that meet the set criteria, we then want to plot samples according to the ordered_sample_list
                        # order which will be in IDs. As such we will have to convert the sample_name in the info_df
                        # to a sample ID using the smp_name_to_smp_id_dict.

                        # get sample_names that fit the requirements
                        sample_names_of_set = self.parent.coral_info_df_for_figures.loc[
                            (self.parent.coral_info_df_for_figures['island'] == location) &
                            (self.parent.coral_info_df_for_figures['site'] == site) &
                            (self.parent.coral_info_df_for_figures['genus'] == spp)
                            ].index.values.tolist()

                        # get the above sample_names_of_set in order of the sample_order list
                        ordered_sample_names_of_set = []
                        for samp_name in self.parent.sample_order:
                            if samp_name in sample_names_of_set:
                                ordered_sample_names_of_set.append(samp_name)

                        # there may be some samples in the sample_names_of_set that aren't in the sample_order
                        # add these here
                        for sample_name in sample_names_of_set:
                            if sample_name not in self.parent.sample_order:
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
                                self._plot_div_over_type_minor_div_only(colour_list, ind, patches_list,
                                                                  smple_id_to_plot)
                            elif qc:
                                self._plot_div_over_type_qc_abund(
                                    colour_list=colour_list, ind=ind, patches_list=patches_list,
                                    smple_id_to_plot=smple_id_to_plot)
                            else:
                                self._plot_div_over_type_18s(colour_list, ind, patches_list, smple_id_to_plot)

                            ind += 1
                        if qc:
                            self._paint_rect_to_axes_div_and_type_18s(
                                ax=ax, colour_list=colour_list, num_smp_in_this_subplot=num_smp_in_this_subplot,
                                patches_list=patches_list, x_tick_label_list=x_tick_label_list,
                                max_num_smpls_in_subplot=10, ax_count=ax_count)
                        else:
                            self._paint_rect_to_axes_div_and_type_18s(
                                ax=ax, colour_list=colour_list, num_smp_in_this_subplot=num_smp_in_this_subplot,
                                patches_list=patches_list, x_tick_label_list=x_tick_label_list,
                                max_num_smpls_in_subplot=10)

                        ax_count += 1

        def _paint_rect_to_axes_div_and_type_18s(
                self, ax, colour_list, num_smp_in_this_subplot, patches_list, x_tick_label_list=None,
                max_num_smpls_in_subplot=10, ax_count=None):
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
            if self.plot_type != 'qc_absolute':
                ax.set_ylim(0, 1)
            ax.set_xticks(range(num_smp_in_this_subplot))
            ax.set_xticklabels(x_tick_label_list, rotation='vertical', fontsize=6)

            if self.plot_type == 'qc_absolute':
                ax.set_ylim(0, 2000000)
                ax.set_yscale('symlog')
                if ax_count % 9 != 0:
                    self._remove_axes_but_allow_labels(ax, x_tick_label_list)
                else:
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)

            else:
                self._remove_axes_but_allow_labels(ax, x_tick_label_list)

            # as well as getting rid of the top and right axis splines
            # I'd also like to restrict the bottom spine to where there are samples plotted but also
            # maintain the width of the samples
            # I think the easiest way to do this is to hack a bit by setting the x axis spines to invisible
            # and then drawing on a line at y = 0 between the smallest and largest ind (+- 0.5)
            # ax.spines['bottom'].set_visible(False)
            if self.plot_type != 'qc_absolute':
                ax.add_line(Line2D((0 - 0.5, num_smp_in_this_subplot - 0.5), (0, 0), linewidth=2, color='black'))
            elif self.plot_type == 'qc_absolute':
                ax.add_line(Line2D((0 - 0.5, num_smp_in_this_subplot - 0.5), (0.1, 0.1), linewidth=2, color='black'))

        def _add_labels(self, ax_list):
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

        def _plot_div_over_type_18s(self, colour_list, ind, patches_list, smple_id_to_plot):
            bottom_div = 0
            # for each sequence, create a rect patch
            # the rect will be 1 in width and centered about the ind value.

            sample_series = self.sample_abundance_df.loc[smple_id_to_plot]
            # https://pandas.pydata.org/pandas-docs/stable/generated/pandas.Series.nonzero.html
            non_zero_series = sample_series[sample_series.nonzero()[0]]
            non_zero_series_index_list = non_zero_series.index.values.tolist()
            print('\nplotting {} intra seqs'.format(len(non_zero_series_index_list)))
            for seq in non_zero_series_index_list:
                # class matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
                rel_abund_div = self.sample_abundance_df.loc[smple_id_to_plot, seq]
                if rel_abund_div > 0:
                    patches_list.append(Rectangle((ind - 0.5, bottom_div), 1, rel_abund_div, color=self.colour_dict[seq]))
                    # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
                    colour_list.append(self.colour_dict[seq])
                    bottom_div += rel_abund_div

        def _plot_div_over_type_qc_abund(self, colour_list, ind, patches_list, smple_id_to_plot):

            bottom_div = 0
            # for each sequence, create a rect patch
            # the rect will be 1 in width and centered about the ind value.

            sample_series = self.sample_abundance_df.loc[smple_id_to_plot]
            # https://pandas.pydata.org/pandas-docs/stable/generated/pandas.Series.nonzero.html
            non_zero_series = sample_series[sample_series.nonzero()[0]]

            if self.plot_type == 'qc_absolute':
                for tax_category in non_zero_series.index.values.tolist():
                    # class matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
                    abs_abund = self.sample_abundance_df.loc[smple_id_to_plot, tax_category]
                    patches_list.append(Rectangle((ind - 0.5, bottom_div), 1, abs_abund, color='#808080'))
                    colour_list.append('#808080')
                    bottom_div += abs_abund

            elif self.plot_type == 'qc_taxa_rel_abund':
                cat_abund_list = []
                cat_name_list = []
                for tax_category in non_zero_series.index.values.tolist():
                    cat_abs_abund = self.sample_abundance_df.loc[smple_id_to_plot, tax_category]
                    cat_abund_list.append(cat_abs_abund)
                    cat_name_list.append(tax_category)
                tot = sum(cat_abund_list)
                rel_abunds_list = [abund / tot for abund in cat_abund_list]
                for i, cat_name in enumerate(cat_name_list):
                    patches_list.append(
                        Rectangle((ind - 0.5, bottom_div), 1, rel_abunds_list[i], color=self.colour_dict[cat_name]))
                    colour_list.append(self.colour_dict[cat_name])
                    bottom_div += rel_abunds_list[i]

        def _plot_div_over_type_minor_div_only(self, colour_list, ind, patches_list, smple_id_to_plot):
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
            sample_series = self.sample_abundance_df.loc[smple_id_to_plot]
            # https://pandas.pydata.org/pandas-docs/stable/generated/pandas.Series.nonzero.html
            non_zero_series = sample_series[sample_series.nonzero()[0]]

            for seq in non_zero_series.index.values.tolist()[1:]:
                rel_abund_div = self.sample_abundance_df.loc[smple_id_to_plot, seq]
                if rel_abund_div > 0 and rel_abund_div < 0.5:
                    seq_list_in_order.append(seq)
                    relabund_list_in_order.append(rel_abund_div)
            # now re normalise the abundances
            tot = sum(relabund_list_in_order)
            for i in range(len(seq_list_in_order)):
                normalised_seq_abund = relabund_list_in_order[i] / tot
                patches_list.append(Rectangle((ind - 0.5, bottom_div), 1, normalised_seq_abund,
                                              color=self.colour_dict[seq_list_in_order[i]]))
                # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
                colour_list.append(self.colour_dict[seq_list_in_order[i]])
                bottom_div += normalised_seq_abund

        def _generate_colour_dict(self, is_med=False, is_qc=False):
            # the purpose of this is to return a seuence to colour dictionary that we will pickle out to maintain
            # continuity
            # make the colour list and the grey list
            # then get the sequence order from the sample_abundance_df
            # then associate to the colours until we are out of colours
            colour_dict = {}
            if is_qc:
                return {'Porites': '#FFFF00', 'Pocillopora': '#87CEFA', 'Millepora': '#FF6347',
                        'other_coral': '#C0C0C0', 'Symbiodiniaceae': '#00FF00', 'other_taxa': '#696969'}

            if is_med:
                colour_list = get_colour_list()[3:]
            else:
                colour_list = get_colour_list()
            grey_palette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']
            for i, seq_name in enumerate(list(self.sample_abundance_df)):
                if i < len(colour_list):
                    colour_dict[seq_name] = colour_list[i]
                else:
                    colour_dict[seq_name] = grey_palette[i % 6]
            return colour_dict

        def _generate_seq_abundance_df(self, numProc=20):
            """ The purpose of this will be to get a dataframe that contains the
            sequence abundances for each of the samples
            This will only contain the sequence abundances for the coral sequences
            so it will not contain any zooxs seqs or non-scleractinian seqs.

            One thing that will be tough is getting a global collection of the sequences that have been found.
            I think the best way to do this is to MP this and return a dictionary for each sample
            where the key will be the sample name and the value will be a dictionary which will be
            key of actual sequences and relabund of the seq once we have all of these we can get a global diversity
            of sequence and we can also get a global abundance of the seuqences so that we can plot the
            most abundant sequenecs first.
            """

            if os.path.isfile(os.path.join(self.parent.cache_dir, 'sample_abundance_df.p')):
                sample_abundance_df = pd.read_pickle(os.path.join(self.parent.cache_dir, 'sample_abundance_df.p'))
            else:
                input_q = Queue()

                for ind in self.parent.coral_info_df_for_figures.index.values.tolist():
                    input_q.put((ind, self.parent.coral_info_df_for_figures.loc[ind, 'sample_dir']))

                for n in range(numProc):
                    input_q.put('STOP')

                all_procs = []
                for n in range(numProc):
                    p = Process(target=self._abundance_worker, args=(input_q,))
                    all_procs.append(p)
                    p.start()

                for p in all_procs:
                    p.join()

                # At this stage we have pickled out abundance dictionaries for each of the samples in their
                # respective directories

                # Now we go through each of samples and check to see if each of the seqs is found in the master seq dict
                # or if it is a subset or superset of any of the strings. This has the potential to be very slow
                # so we could maybe chunk this and do it MP. But for the time being lets try to do it serial.

                # First we will go through each of the samples and
                # generate a set of sequence for the master_abund_dict
                # that will be the largest superset sequences of all the sequences found in the samples.

                # Once we've done that, then we will go back through the sequences and see which
                # sequences the sample's sequences should be associated to.
                # We should pickle out along the way to save us time if we have to re-run

                master_abund_dict = self._get_master_abund_dict()

                # Now we once again go sample by sample, seq by seq and check all of the seqs in each sample and
                # identify which of the master abund dict sequences the sequence in question should be associated to.
                # We then use that representative sequence to populate the dataframe

                header_seq_s = self._make_sample_abundance_df_columns_headers(master_abund_dict)
                sample_abundance_df = self._make_empty_abundance_df(header_seq_s)

                # For every coral sample
                for ind in self.parent.coral_info_df_for_figures.index.values.tolist():

                    print(f'Populating df with {ind}')
                    sample_dir = self.parent.coral_info_df_for_figures.loc[ind, 'sample_dir']
                    if os.path.isfile(os.path.join(sample_dir, 'sample_abundance_series.pkl')):
                        sample_abundance_series = pd.read_pickle(os.path.join(sample_dir, 'sample_abundance_series.pkl'))
                    else:
                        sample_abundance_series = self._generate_sample_abundance_series_from_scratch(
                            header_seq_s, ind, master_abund_dict, sample_dir)
                    sample_abundance_df = sample_abundance_df.append(sample_abundance_series)
                # at this point we have fully populated the sample_abundance_df and we should pickle it out
                sample_abundance_df.to_pickle(os.path.join(self.parent.cache_dir, 'sample_abundance_df.p'))

            return sample_abundance_df

        def _generate_sample_abundance_series_from_scratch(self, header_seq_s, ind, master_abund_dict, sample_dir):
            # create an empty series that will be full of float 0s to start with
            sample_abundance_series = pd.Series(name=ind, data=[float(0) for i in header_seq_s],
                                                index=header_seq_s)
            sample_abund_dict = pickle.load(open(os.path.join(sample_dir, 'seq_abund_dict.pickle'), 'rb'))
            # For every sequence in the sample
            for sample_key_sequence, v_sample in sample_abund_dict.items():
                # list that holds the sequences that are either a subset or superset of the sample sequence
                sub_or_super_match_list = []
                exact_match = False
                # For every sequence in the master abundance dictionary
                for master_key_sequence, v_master in master_abund_dict.items():
                    if sample_key_sequence == master_key_sequence:
                        # then we have an exact match and this is the sequence that
                        # the abundance should be associated to.
                        # else continue to work with the match_list etc.
                        sample_abundance_series[master_key_sequence] += v_sample
                        exact_match = True
                        break
                    if sample_key_sequence in master_key_sequence or master_key_sequence in sample_key_sequence:
                        sub_or_super_match_list.append(master_key_sequence)
                if exact_match:
                    continue
                # here we have not acheived an exact match and we should work with the list
                # here we have checked each of the master sequences
                if not sub_or_super_match_list:
                    # this should never happen
                    continue
                elif len(sub_or_super_match_list) == 1:
                    # then we have a single match and this is the column of the df that we should associate the
                    # abundance to
                    sample_abundance_series[sub_or_super_match_list[0]] += v_sample

                elif len(sub_or_super_match_list) > 1:
                    # then we have a bit of a problem here and we need to work out what to do from looking at the
                    # situation
                    # this is a particular situation and has happened when there is a sequence that is shorter than the
                    # majority and one of the sequences has a differntiating snp in the difference in sequence length.
                    # for the time being I will attribute the abundance of this sequence to the most abundant match
                    sorted_tups_temp = sorted([(seq, master_abund_dict[seq]) for seq in sub_or_super_match_list],
                                              key=lambda x: x[1],
                                              reverse=True)
                    most_abund_seq = [a[0] for a in sorted_tups_temp][0]
                    sample_abundance_series[most_abund_seq] += v_sample
            # here we have the series for the sample populated.
            # now add this to the dataframe
            # pickle out the series to save time
            sample_abundance_series.to_pickle('{}/sample_abundance_series.pkl'.format(sample_dir))
            return sample_abundance_series

        def _make_empty_abundance_df(self, header_seq_s):
            sample_abundance_df = pd.DataFrame(columns=header_seq_s)
            return sample_abundance_df

        def _make_sample_abundance_df_columns_headers(self, master_abund_dict):
            sorted_tups = sorted(list(master_abund_dict.items()), key=lambda x: x[1], reverse=True)
            header_seq_s = [a[0] for a in sorted_tups]
            return header_seq_s

        def _get_master_abund_dict(self):
            if os.path.isfile(os.path.join(self.parent.cache_dir, 'master_abund_dict.p')):
                master_abund_dict = pickle.load(
                    open(os.path.join(self.parent.cache_dir, 'master_abund_dict.p'), 'rb'))
            else:
                master_abund_dict = self._make_master_abund_dict()
            return master_abund_dict

        def _make_master_abund_dict(self):
            master_abund_dict = defaultdict(float)
            # For every coral sample
            tot = len(self.parent.coral_info_df_for_figures.index.values.tolist())
            count = 0
            for ind in self.parent.coral_info_df_for_figures.index.values.tolist():
                count += 1
                print(f'Counting {ind}: {count} out of {tot} samples')
                sample_dir = self.parent.coral_info_df_for_figures.loc[ind, 'sample_dir']
                sample_abund_dict = pickle.load(open(os.path.join(sample_dir, 'seq_abund_dict.pickle'), 'rb'))
                # For every sequence in the coral
                for sample_key_sequence, v_sample in sample_abund_dict.items():
                    # For each of the sequences, see if it is found within one of the
                    # master_abund_dicts sequences.
                    # Check to see that it is not found in more than one of the sequences too.
                    # If it is then this is a problem but there may be a simple solution,
                    # i.e. consolidating the two sequences.

                    # list that holds the sequences that are either a subset or superset of the sample sequence
                    sub_or_super_match_list = []
                    # For every sequence in the master dictionary
                    for master_key_sequence, v_master in master_abund_dict.items():
                        if sample_key_sequence in master_key_sequence or master_key_sequence in sample_key_sequence:
                            sub_or_super_match_list.append(master_key_sequence)

                    # here we have checked each of the master sequences
                    if not sub_or_super_match_list:
                        # if the match list is empty then this is a new sequence and we can add it to the master_abund_dict
                        master_abund_dict[sample_key_sequence] += v_sample
                    elif len(sub_or_super_match_list) == 1:
                        self._update_representative_master_sequence_abundance(
                            master_abund_dict, sample_key_sequence, sub_or_super_match_list, v_sample)
                    elif len(sub_or_super_match_list) > 1:
                        # Then we have a bit of a problem here and we need to work out what to do.
                        # From looking at the situation
                        # this is a particular situation and has happened when there is a sequence that
                        # is shorter than the majority and one of the sequences has a differntiating snp
                        # in the difference in sequence length.
                        # For the time being I will attribute the abundance
                        # of this sequence to the most abundant match
                        most_abund_seq = [a[0] for a in
                                          sorted([(seq, master_abund_dict[seq]) for seq in sub_or_super_match_list],
                                                 key=lambda x: x[1], reverse=True)][0]
                        master_abund_dict[most_abund_seq] += v_sample
                        continue
            # at this point we have the master_abund_dict populated and we should pickle it out
            pickle.dump(master_abund_dict, open(os.path.join(self.parent.cache_dir, 'master_abund_dict.p'), 'wb'))
            return master_abund_dict

        def _update_representative_master_sequence_abundance(self, master_abund_dict, sample_key_sequence,
                                                             sub_or_super_match_list, v_sample):
            # Then we have a single match.
            if self._if_sample_seq_is_subset_of_master_seq(
                    master_seq=sub_or_super_match_list[0], sample_seq=sample_key_sequence):
                self._attribute_sample_seq_abund_to_master_seq_representative(
                    master_abund_dict, sub_or_super_match_list, v_sample)
            elif self._if_master_seq_is_subset_of_sample_seq(
                    master_seq=sub_or_super_match_list[0], sample_seq=sample_key_sequence):
                self._make_sample_seq_representative_and_update_abund(
                    master_abund_dict, sample_key_sequence, sub_or_super_match_list, v_sample)

        def _attribute_sample_seq_abund_to_master_seq_representative(self, master_abund_dict, sub_or_super_match_list,
                                                                     v_sample):
            # if the sample is a subset of one of the master_abund_dicts then its
            # abundance should be attributed to this sequence
            master_abund_dict[sub_or_super_match_list[0]] += v_sample

        def _make_sample_seq_representative_and_update_abund(self, master_abund_dict, sample_key_sequence,
                                                             sub_or_super_match_list, v_sample):
            # else if the sample sequences contains the master sequence, then the master sequence
            # should be updated to the seq in question, the current abundance should be associated
            # to this new sequence and the abundance of the seq in q should also be attributed (added)
            temp_abundance = master_abund_dict[sub_or_super_match_list[0]]
            del master_abund_dict[sub_or_super_match_list[0]]
            master_abund_dict[sample_key_sequence] = temp_abundance
            master_abund_dict[sample_key_sequence] += v_sample

        def _if_master_seq_is_subset_of_sample_seq(self, master_seq, sample_seq):
            return master_seq in sample_seq

        def _if_sample_seq_is_subset_of_master_seq(self, master_seq, sample_seq):
            return sample_seq in master_seq

        def _abundance_worker(self, input_q):
            """This worker will be run in a multiprocessing framework and will work on a sample
            by sample basis to create a dictionary per sample. This dictionary will
            be the nucleotide sequence as key and the relative abundance as the value. These
            dictionaries will then be collected outside of the multiprocessing framework
            so that overall abundances across samples can be calculated."""
            for sample_name, sample_dir in iter(input_q.get, 'STOP'):
                # check to see if the sample's abundance dictionary has already been created
                if os.path.isfile(os.path.join(sample_dir, 'seq_abund_dict.pickle')):
                    continue

                # I want this dictionary to have key as the raw sequence and the relative abundance
                # as the value
                sample_dict = defaultdict(int)

                sys.stdout.write('\nSample {}\n'.format(sample_dir))

                # for each sample we will only need two files
                # we will need the name file inorder to create an abundance dict
                # then we will need the coral_genus_only.fasta file so that we know which sequences we should be counting

                # read in the name file
                with open('{}/stability.trim.contigs.good.abund.pcr.names'.format(sample_dir), 'r') as f:
                    name_file = [line.rstrip() for line in f]

                abund_dict = {line.split('\t')[0]: len(line.split('\t')[1].split(',')) for line in name_file}

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
                normalised_sample_dict = {k: v / tot_seqs for k, v in sample_dict.items()}

                # here we have the samples abundance dictionary populated. Now pickle out
                pickle.dump(normalised_sample_dict, open(os.path.join(sample_dir, 'seq_abund_dict.pickle'), 'wb'))

    class SequenceQC():
        """Mothur qc of the sequences"""

        def __init__(self, parent):
            self.parent = parent

        def do_seq_qc(self, numProc=20):
            """this function will be responsible for creating the processed name and fasta pairs from the fastq files
            # that can be found in the info df
            """

            # lets make an iput queue that is going to be each of the samples
            mothur_qc_input_queue = Queue()

            # for each sample in the info_df I will add a tuple that is a pair of the fwd and rev directories to the fastq.gz
            for ind in self.parent.all_samples_info_df.index.values.tolist():
                mothur_qc_input_queue.put(
                    (self.parent.all_samples_info_df.loc[ind, 'fastq_fwd_file_path'],
                     self.parent.all_samples_info_df.loc[ind, 'fastq_rev_file_path']))

            for n in range(numProc):
                mothur_qc_input_queue.put('STOP')

            all_procs = []
            for n in range(numProc):
                p = Process(target=self._mothur_worker, args=(mothur_qc_input_queue,))
                all_procs.append(p)
                p.start()

            for p in all_procs:
                p.join()

            return

        def _mothur_worker(self, input_q):

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


    def _generate_info_df_for_samples(self):
        """Generate two dataframes, the first that contains information from all samples
        The second will contain information from only the coral samples and will be used in the plotting.
        This information will be gathered by parsing the directory
        structure that the fastq sequencing files are housed in
        """
        if os.path.isfile(os.path.join(self.cache_dir, 'all_samples_info_df.p')):
            self.all_samples_info_df = pickle.load(open(os.path.join(self.cache_dir, 'all_samples_info_df.p'), 'rb'))
            self.coral_info_df_for_figures = pickle.load(
                open(os.path.join(self.cache_dir, 'coral_info_df_for_figures.p'), 'rb'))
        else:
            columns_for_df = ['sample_name', 'fastq_fwd_file_path', 'fastq_rev_file_path', 'coral_plankton',
                              'spp_water', 'location', 'site', 'size_fraction']
            self.all_samples_info_df = pd.DataFrame(columns=columns_for_df)

            # now lets parse through the directories using them to get some of the information facts above
            self._parse_data_dir_structure_to_infer_sample_info_and_populate_df()

            pickle.dump(self.all_samples_info_df, open(os.path.join(self.cache_dir,'all_samples_info_df.p'), 'wb'))

            self._generate_coral_info_df_for_figures_from_all_sample_info_df()

            pickle.dump(self.coral_info_df_for_figures,
                        open(os.path.join(self.cache_dir, 'coral_info_df_for_figures.p'), 'wb'))



    def _generate_coral_info_df_for_figures_from_all_sample_info_df(self):
        if os.path.isfile(os.path.join(self.cache_dir, 'coral_info_df_for_figures.p')):
            self.coral_info_df_for_figures = pickle.load(
                open(os.path.join(self.cache_dir, 'coral_info_df_for_figures.p'), 'rb'))
        else:
            self.coral_info_df_for_figures = pd.DataFrame(columns=['island', 'site', 'genus', 'individual', 'sample_dir'])
            for ind in self.all_samples_info_df.index.values.tolist():
                if 'CORAL' in self.all_samples_info_df.loc[ind, 'fastq_fwd_file_path']:
                    fastq_string = self.all_samples_info_df.loc[ind, 'fastq_fwd_file_path']
                    components = fastq_string.split('/')
                    smp_site = components[-6]
                    smp_island = components[-7]
                    smp_individual = components[-3]
                    smp_genus = components[-4]
                    smp_dir = '/'.join(components[:-1])
                    self.coral_info_df_for_figures = self.coral_info_df_for_figures.append(
                        pd.Series(name=ind, data=[smp_island, smp_site, smp_genus, smp_individual, smp_dir],
                                  index=['island', 'site', 'genus', 'individual', 'sample_dir']))


    def _parse_data_dir_structure_to_infer_sample_info_and_populate_df(self):
        """Parse through the directory strucutre that holds the fastq files to get
        the sample infomation and populate the info_df """
        for location in os.listdir(self.data_root_dir):
            print(f'Parsing directory strucutre for {location}')
            info_dir_parser = self.InfoDirectoryParser(parent=self, location=location)
            info_dir_parser.parse()

    class InfoDirectoryParser:
        def __init__(self, parent, location):
            self.parent = parent
            self.location = location
            self.current_dir = None
            self.sample_type = None
            self.site = None
            self.spp_water = None
            self.size_fraction = None
            self.individual = None

        def parse(self):
            if 'ISLAND' in self.location:
                self._parse_island_directory()

            elif 'OA' in self.location:
                self._parse_oa_directory()

        def _parse_oa_directory(self,):
            self.current_dir = os.path.join(self.parent.data_root_dir, self.location, 'PLANKTON', 'SURFACE')
            for size_fraction_indi in os.listdir(self.current_dir):
                self.size_fraction = size_fraction_indi
                self.current_dir = os.path.join(self.current_dir, size_fraction_indi)
                self.spp_water = 'PLANKTON'
                self.sample_type = 'OA'
                self.site = 'OA'
                self._populate_info_df_sample_row()
                self._set_current_dir_up_one_dir()

        def _set_current_dir_up_one_dir(self):
            self.current_dir = '/'.join(self.current_dir.split('/')[:-1])

        def _parse_island_directory(self):
            self.current_dir = os.path.join(self.parent.data_root_dir, self.location)
            for site in os.listdir(self.current_dir):
                self.site = site
                self.current_dir = os.path.join(self.current_dir, self.site)
                for sample_type in os.listdir(self.current_dir):
                    self.sample_type = sample_type
                    self.current_dir = os.path.join(self.current_dir, sample_type)

                    if sample_type == 'CORAL':
                        self.size_fraction = 'coral'
                        self._parse_coral_directory()

                    elif sample_type == 'PLANKTON':
                        self._parse_plankton_directory()
                self._set_current_dir_up_one_dir()
            self._set_current_dir_up_one_dir()

        def _parse_plankton_directory(self):
            for water_type in os.listdir(self.current_dir):
                self.spp_water = water_type
                if water_type == 'CSW':
                    self._parse_csw_directory()
                elif water_type == 'SURFACE':
                    self._parse_surface_directory()
            self._set_current_dir_up_one_dir()

        def _parse_surface_directory(self):
            # then this is a SURFACE sample and there are no individuals
            self.size_fraction='S320'
            self.current_dir = os.path.join(self.current_dir, self.spp_water, self.size_fraction)
            # collect the information we need
            self._populate_info_df_sample_row()
            self._set_current_dir_up_one_dir()
            self._set_current_dir_up_one_dir()

        def _parse_csw_directory(self):
            self.current_dir = os.path.join(self.current_dir, self.spp_water)
            for individual in os.listdir(self.current_dir):
                self.current_dir = os.path.join(self.current_dir, individual)
                for size_fraction_indi in os.listdir(self.current_dir):
                    self.size_fraction = size_fraction_indi
                    self.current_dir = os.path.join(self.current_dir, self.size_fraction)
                    # now we are in the directory that contains the actual paired fastq.gz files for a
                    # given water sample
                    self._populate_info_df_sample_row()
                    self._set_current_dir_up_one_dir()
                self._set_current_dir_up_one_dir()
            self._set_current_dir_up_one_dir()

        def _parse_coral_directory(self):
            for species in os.listdir(self.current_dir):
                self.spp_water = species
                self.current_dir = os.path.join(self.current_dir, self.spp_water)
                for individual in os.listdir(self.current_dir):
                    self.individual = individual
                    self.current_dir = os.path.join(self.current_dir, self.individual, 'CS4L')
                    # now we are in the directory that contains the actual paired fastq.gz files for a
                    # given coral individual
                    # collect the information we need
                    self._populate_info_df_sample_row()
                    self._set_current_dir_up_one_dir()
                    self._set_current_dir_up_one_dir()
                self._set_current_dir_up_one_dir()
            self._set_current_dir_up_one_dir()

        def _populate_info_df_sample_row(self):

            sample_name = [file for file in os.listdir(self.current_dir) if 'fastq.gz' in file][0].split('_')[0]

            fwd_path, rev_path = self._get_fwd_and_rev_fastq_paths()

            temp_dict = {
                'sample_name':sample_name ,'fastq_fwd_file_path': fwd_path, 'fastq_rev_file_path': rev_path,
                'coral_plankton': self.sample_type, 'spp_water': self.spp_water, 'location': self.location,
                'site': self.site, 'size_fraction': self.size_fraction
            }

            temp_series = pd.Series(data=temp_dict, name=sample_name)
            self.parent.all_samples_info_df = self.parent.all_samples_info_df.append(temp_series)

        def _get_fwd_and_rev_fastq_paths(self):
            fastq_files = [file for file in os.listdir(self.current_dir) if 'fastq.gz' in file]
            fastq_files = self._marge_fastq_files_if_more_than_two_exist(fastq_files)
            fwd_path, rev_path = self._infer_fastq_paths(fastq_files)
            return fwd_path, rev_path

        def _infer_fastq_paths(self, files):
            fwd_found = False
            rev_found = False
            fwd_path = None
            rev_path = None
            for file_name in files:
                if 'R1' in file_name:
                    fwd_path = os.path.join(self.current_dir, file_name)
                    fwd_found = True
                elif 'R2' in file_name:
                    rev_path = os.path.join(self.current_dir, file_name)
                    rev_found = True
            # make sure that both the fwd and rev paths have been identified
            if not fwd_found or not rev_found:
                print('fwd or rev read not found')
                sys.exit(1)
            return fwd_path, rev_path

        def _marge_fastq_files_if_more_than_two_exist(self, files):
            if len(files) != 2:
                print('merging fastqs in {}'.format(self.current_dir))
                # this method will merge the fastqs together so that 4 turn into 2.
                self._merge_fastqs(files)

                # finally re-read in the files
                files = os.listdir(self.current_dir)
                if len(files) != 2:
                    print('we still gotta problem {}'.format(self.current_dir))
                    sys.exit(0)
            return files

        def _merge_fastqs(self, files):
            # If we get here then there have been multiple sequencing runs for the same sample
            # we will aim to simply merge the fastqs together insitu
            # first get the name that we want to use (this is the one containing BG
            one_found = False
            two_found = False
            R1_name = None
            R2_name = None
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
                print('couldnt find the right files {}'.format(self.current_dir))
            # second unzip all of the files
            for file_name in files:
                subprocess.run(['gunzip', '{}/{}'.format(self.current_dir, file_name)])
            # now we have all of the files unzipped
            # now look at each of the files and add them to a master R1 and R2.
            # it doesn't matter whether we add the R1 or R2 first. So long as we add the pair at the same time
            un_files = os.listdir(self.current_dir)
            master_fastq_R1 = []
            master_fastq_R2 = []
            for un_file in un_files:
                # we should go into this if twice once for each pair
                if 'R1' in un_file:
                    # then we can read this and its pair
                    rone_path = '{}/{}'.format(self.current_dir, un_file)
                    with open(rone_path, 'r') as f:
                        rone_file = [line.rstrip() for line in f]
                    master_fastq_R1.extend(rone_file)

                    # now do the same for the corresponing R2 file
                    rtwo_path = '{}/{}'.format(self.current_dir, un_file.replace('R1', 'R2'))
                    with open(rtwo_path, 'r') as f:
                        rtwo_file = [line.rstrip() for line in f]
                    master_fastq_R2.extend(rtwo_file)

                    # now delte the files
                    os.remove(rone_path)
                    os.remove(rtwo_path)
            # here we have a master file for the R1 and R2. We can now write it out to the name that we have
            rone_path_to_write = '{}/{}'.format(self.current_dir, R1_name.replace('.gz', ''))
            with open(rone_path_to_write, 'w') as f:
                for line in master_fastq_R1:
                    f.write('{}\n'.format(line))
            rtwo_path_to_write = '{}/{}'.format(self.current_dir, R2_name.replace('.gz', ''))
            with open(rtwo_path_to_write, 'w') as f:
                for line in master_fastq_R2:
                    f.write('{}\n'.format(line))
            # now we simply need to recompress the files
            un_files = os.listdir(self.current_dir)
            for un_file in un_files:
                subprocess.run(['gzip', '{}/{}'.format(self.current_dir, un_file)])

    class NetworkStuff:
        """
        # TODO this is still work in progress. I have just shoved all of the methods that were related to this
        in progress work in to this class as a holder. It still needs to be properly refactored.
        I want to make splits tree networks of the 18s sequences. To do this i will need to get the seuqences in better
        shape so that I can align them. Currently they are of very different lengths. I will go back into each of the
        sample folders and I will align the latest fasta file and then do cropping using 90% cutoff or something.
        this means that I will get rid of gaps at the beginning and end of the alignment if the column position has gaps for
        90% of the sequences.
        On these cropped sequences I will then re-unique and then finally realign the sequences. These sequences will then
        be ready for making networks from. It may also be a good idea to work with these sequences for all of the work we have
        done up until now.

        """
        def __init__(self, parent):
            self.parent = parent

        def prepare_sequences_for_networking(self):

            # for each coral sample do the alignment and cropping and re-alignment of sequences first
            for ind in self.parent.coral_info_df_for_figures.index.values.tolist():
                sample_dir = self.parent.coral_info_df_for_figures.loc[ind, 'sample_dir']

                if os.path.isfile('{}/coral_aligned_fasta_for_networks.fasta'.format(sample_dir)) and os.path.isfile(
                        '{}/coral_fasta_aligned_and_cropped.names'.format(sample_dir)):
                    # then this sample has already been completed
                    continue

                print('Processing {}'.format(ind))
                sample_genus = self.parent.coral_info_df_for_figures.loc[ind, 'genus']

                current_fasta_file_path = '{}/stability.trim.contigs.good.unique.abund.pcr.fasta'.format(sample_dir)
                with open(current_fasta_file_path, 'r') as f:
                    current_fasta_file = [line.rstrip() for line in f]
                current_fasta_dict = {current_fasta_file[i][1:].split('\t')[0]: current_fasta_file[i + 1] for i in
                                      range(0, len(current_fasta_file), 2)}

                current_names_file_path = '{}/stability.trim.contigs.good.abund.pcr.names'.format(sample_dir)
                with open(current_names_file_path, 'r') as f:
                    current_names_file = [line.rstrip() for line in f]
                current_names_dict = {line.split('\t')[0]: line for line in current_names_file}

                # this fasta currently contains all of the sequences including non-genus specific seqs
                # we are only interested in the genus specific seqs so lets pull these out using the
                # coral_dict.pickle file
                coral_seq_dict = pickle.load(open('{}/coral_dict.pickle'.format(sample_dir), 'rb'))

                fasta_for_alignment = []
                names_for_alignment = []
                for seq_key, coral_genus in coral_seq_dict.items():
                    if coral_genus.upper() == sample_genus.upper():
                        fasta_for_alignment.extend(['>{}'.format(seq_key), current_fasta_dict[seq_key]])
                        names_for_alignment.append(current_names_dict[seq_key])

                # here we have a fasta file and a names file pair that are just the coral genus in question.
                # we should now write these out and then align them. Then do the cropping
                path_to_fasta_file_to_align = '{}/coral_fasta_to_align_and_crop.fasta'.format(sample_dir)
                path_to_names_file_to_align = '{}/coral_names.names'.format(sample_dir)

                # write out the fasta
                with open(path_to_fasta_file_to_align, 'w') as f:
                    for line in fasta_for_alignment:
                        f.write('{}\n'.format(line))

                # write out the .names file
                with open(path_to_names_file_to_align, 'w') as f:
                    for line in names_for_alignment:
                        f.write('{}\n'.format(line))

                aligned_fasta_path = '{}/coral_fasta_aligned_to_crop.fasta'.format(sample_dir)

                self.align_fasta(input_fasta_path=path_to_fasta_file_to_align, output_fasta_path=aligned_fasta_path)

                # this method takes an input path and an output path
                path_to_fasta_cropped = '{}/coral_fasta_aligned_and_cropped.fasta'.format(sample_dir)

                # at this point we have the original fasta aligned and cropped.
                cropping.crop_fasta(input_path=aligned_fasta_path, output_path=path_to_fasta_cropped, cutoff=0.9)

                # we now need to remove the gaps from the sequences else we end up with a strange '.' character in our fasta
                # file after the seqs.unique
                # read in the aligned fasta file

                self.remove_gaps_from_alignment(input_fasta_alignment_path=path_to_fasta_cropped)

                # now we re run mothur to do a unique on the fasta a names pair
                mBatchFile = [
                    r'set.dir(input={})'.format(sample_dir),
                    r'set.dir(output={})'.format(sample_dir),
                    r'unique.seqs(fasta={}, name={})'.format(path_to_fasta_cropped, path_to_names_file_to_align),
                    r'summary.seqs(fasta={0}/coral_fasta_aligned_and_cropped.unique.fasta, name={0}/coral_fasta_aligned_and_cropped.names)'.format(
                        sample_dir)
                ]

                mBatchFile_path = '{}/mBatchFile_two'.format(sample_dir)

                # write out batch file
                with open(mBatchFile_path, 'w') as f:
                    for line in mBatchFile:
                        f.write('{}\n'.format(line))

                # run the mothur processing
                # subprocess.run(['mothur', r'{0}'.format(mBatchFile_path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                subprocess.run(['mothur', r'{}'.format(mBatchFile_path)])

                # at this point we should have a new .fasta file and a new .names file
                # we will then want to align these
                # this seems to have worked really well. We are down to 149 sequences

                # now we need to re-align these sequencs to get an alignment to work with
                input_fasta_for_alignment = '{}/coral_fasta_aligned_and_cropped.unique.fasta'.format(sample_dir)
                output_fasta_aligned = '{}/coral_aligned_fasta_for_networks.fasta'.format(sample_dir)
                self.align_fasta(input_fasta_path=input_fasta_for_alignment, output_fasta_path=output_fasta_aligned)
                # at this point we have the .names file which is coral_fasta_aligned_and_cropped.names and the aligned
                # fasta that is coral_aligned_fasta_for_networks.fasta
            apples = 'asdf'

        def remove_gaps_from_alignment(self, input_fasta_alignment_path):
            with open(input_fasta_alignment_path, 'r') as f:
                fasta_to_remove_gaps = [line.rstrip() for line in f]
            fasta_without_gaps = []
            for i in range(len(fasta_to_remove_gaps)):
                if i % 2 == 1:
                    fasta_without_gaps.append(fasta_to_remove_gaps[i].replace('-', ''))
                else:
                    fasta_without_gaps.append(fasta_to_remove_gaps[i])
            # now write out the fasta without gaps
            with open(input_fasta_alignment_path, 'w') as f:
                for line in fasta_without_gaps:
                    f.write('{}\n'.format(line))

        def align_fasta(self, input_fasta_path, output_fasta_path):
            # now perform the alignment with MAFFT
            mafft = local["mafft-linsi"]
            out_file = input_fasta_path.replace('.fasta', '_aligned.fasta')
            # now run mafft including the redirect
            (mafft['--thread', -1, input_fasta_path] > out_file)()
            # read in the interleaved aligned fasta
            with open(out_file, 'r') as f:
                aligned_fasta_interleaved = [line.rstrip() for line in f]
            # make a serial fasta from the interleaved fasta
            aligned_fasta = self.convert_interleaved_to_sequencial_fasta_two(aligned_fasta_interleaved)
            # write out the fasta to be cropped
            with open(output_fasta_path, 'w') as f:
                for line in aligned_fasta:
                    f.write('{}\n'.format(line))

        def convert_interleaved_to_sequencial_fasta_two(self, fasta_in):
            fasta_out = []
            for i in range(len(fasta_in)):
                if fasta_in[i].startswith('>'):
                    if fasta_out:
                        # if the fasta is not empty then this is not the first
                        fasta_out.append(temp_seq_str)
                    # else then this is the first sequence and there is no need to add the seq.
                    temp_seq_str = ''
                    fasta_out.append(fasta_in[i])
                else:
                    temp_seq_str = temp_seq_str + fasta_in[i]
            # finally we need to add in the last sequence
            fasta_out.append(temp_seq_str)
            return fasta_out


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

eighteen_s_analysis = EighteenSAnalysis()
# eighteen_s_analysis.do_qc()
# eighteen_s_analysis.do_taxa_annotations()
# eighteen_s_analysis.plot_seq_stacked_bar_plots(plot_type='full')
# eighteen_s_analysis.plot_seq_stacked_bar_plots(plot_type='low')
# eighteen_s_analysis.plot_seq_stacked_bar_plots(plot_type='qc_taxa_rel_abund')
# eighteen_s_analysis.plot_seq_stacked_bar_plots(plot_type='qc_absolute')
# eighteen_s_analysis.plot_seq_stacked_bar_plots(plot_type='med')
# eighteen_s_analysis.plot_pcoa_spp(distance_method='braycurtis')
# eighteen_s_analysis.plot_pcoa_spp_island(distance_method='braycurtis')
eighteen_s_analysis.plot_pcoa_spp_18s_its2(distance_method='braycurtis')
# TODO draw up some networks

