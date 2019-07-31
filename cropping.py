#!/usr/bin/env python3.6
import pandas as pd
import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser(description='simple script for corpping the gaps on the begninning and end of alignments', epilog='For support email: benjamin.hume@kaust.edu.sa')
    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument('-i', '--input', help='path to aligned fasta to crop', required=True)
    required_named.add_argument('-o', '--output', help='path to write out the cropped and aligned fasta to', required=True)

    parser.add_argument('-c', '--cutoff', help='The minimum proportion of sequences that have a gap at any given column '
                                               'in order for the column to be cut from the alignment. Once column has '
                                               'been reached that has a proportion of -s lower than this, '
                                               'the cropping of that end will stop. Default = 1.0. e.g. 0.8', required=False, default=1, type=float)

    args = parser.parse_args()

    cropper = Cropper(input_path=args.input, output_path=args.output, cutoff=args.cutoff)
    cropper.crop()

class Cropper:
    def __init__(self, input_path=None, output_path=None, aligned_fasta_as_list = None, cutoff=1.0, return_as='list'):
        # if working with input path
        if input_path:
            self.fasta_as_pandas_df = self._fasta_to_pandas_df(self._read_file_to_list(input_path))
        elif aligned_fasta_as_list:
            self.fasta_as_pandas_df = self._fasta_to_pandas_df(aligned_fasta_as_list)
        else:
            raise RuntimeError('Please provide either the path to the aligned fasta or the aligned fasta as a list')

        # if output path then we will write out the cropped fasta, else we will return as a list
        if output_path:
            self.return_as = 'fasta_file'
            self.output_path = output_path
        else:
            self.return_as = return_as

        self.cutoff = cutoff

        self.cropped_fasta_as_df = None
        self.corpped_fasta_as_list = None

    def crop(self):

        # crop the df
        self.cropped_fasta_as_df = self._crop_aligned_fasta_df()
        # convert the df back to list object ready to write out
        if self.return_as == 'list':
            return self._pandas_df_to_fasta()
        elif self.return_as == 'fasta_file':
            self.cropped_fasta_list = self._pandas_df_to_fasta()
            # now write out the corpped fasta list
            self._write_list_as_file_to_path()

    def _read_file_to_list(self, filename):
        with open(filename, mode='r') as f:
             return [line.rstrip() for line in f]

    def _write_list_as_file_to_path(self):
        with open(self.output_path, 'w') as f:
            for line in self.cropped_fasta_list:
                f.write('{}\n'.format(line))

    def _fasta_to_pandas_df(self, fasta_as_list):
        temp_df = pd.DataFrame([list(line) for line in fasta_as_list if not line.startswith('>')])
        seq_names = [line[1:] for line in fasta_as_list if line.startswith('>')]
        temp_df.index=seq_names
        return temp_df

    def _crop_aligned_fasta_df(self):
        columns_to_drop = []
        for i in list(self.fasta_as_pandas_df):
            # if there is a gap in the column at the beginning
            #add functionality to incorporate the cutoff.
            # this cutoff will mean that the proportion of the gaps in a column must be greater than this for the
            # column to be removed from the df. Once a column has been reached that has a proportion lower than this value
            # cropping will be complete for that end.
            list_of_col_characters = self.fasta_as_pandas_df.loc[:, i].values.tolist()
            count_of_gaps = list_of_col_characters.count('-') + list_of_col_characters.count('*')
            total_characters = len(list_of_col_characters)
            prop_of_gaps = count_of_gaps / total_characters
            if prop_of_gaps >= self.cutoff:
                columns_to_drop.append(i)
            else:
                break

        for i in reversed(list(self.fasta_as_pandas_df)):
            # if there is a gap in the column at the end
            list_of_col_characters = self.fasta_as_pandas_df.loc[:, i].values.tolist()
            count_of_gaps = list_of_col_characters.count('-') + list_of_col_characters.count('*')
            total_characters = len(list_of_col_characters)
            prop_of_gaps = count_of_gaps / total_characters
            if prop_of_gaps >= self.cutoff:
                columns_to_drop.append(i)
            else:
                break

        # get a list that is the columns indices that we want to keep
        col_to_keep = [col_index for col_index in list(self.fasta_as_pandas_df) if col_index not in columns_to_drop]

        # drop the gap columns
        return self.fasta_as_pandas_df[col_to_keep]

    def _pandas_df_to_fasta(self):
        temp_fasta = []
        for ind in self.cropped_fasta_as_df.index.tolist():
            temp_fasta.extend(['>{}'.format(ind), ''.join(list(self.cropped_fasta_as_df.loc[ind]))])
        return temp_fasta

if __name__ == "__main__":
    main()