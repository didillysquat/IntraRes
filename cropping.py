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

    input_path = args.input
    output_path = args.output
    cutoff = args.cutoff

    crop_fasta(input_path, output_path, cutoff)


def crop_fasta(input_path, output_path, cutoff=1.0):
    # read in the path to the input fasta as list object
    try:
        fasta_to_crop = read_file_to_list(input_path)
    except:
        sys.exit('error reading in the aligned fasta')

    # input_path = '{}/test_fastas/tester_aligned.fasta'.format(os.getcwd())
    # output_path = '{}/test_fastas/tester_aligned_cropped.fasta'.format(os.getcwd())
    # try:
    #     fasta_to_crop = read_file_to_list(input_path)
    # except:
    #     sys.exit('error reading in the aligned fasta')
    # create a df from the fasta
    fasta_as_pandas_df = fasta_to_pandas_df(fasta_as_list=fasta_to_crop)
    # crop the df
    cropped_fasta_as_df = crop_aligned_fasta_df(aligned_fasta_as_pandas_df_to_crop=fasta_as_pandas_df, cutoff=cutoff)
    # convert the df back to list object ready to write out
    cropped_fasta_list = pandas_df_to_fasta(cropped_fasta_df=cropped_fasta_as_df)
    # now write out the corpped fasta list
    write_list_as_file_to_path(cropped_fasta_list, output_path)


def read_file_to_list(filename):
    with open(filename, mode='r') as f:
         return [line.rstrip() for line in f]

def write_list_as_file_to_path(file, path):
    with open(path, 'w') as f:
        for line in file:
            f.write('{}\n'.format(line))

def fasta_to_pandas_df(fasta_as_list):
    temp_df = pd.DataFrame([list(line) for line in fasta_as_list if not line.startswith('>')])
    seq_names = [line[1:] for line in fasta_as_list if line.startswith('>')]
    temp_df.index=seq_names
    return temp_df

def crop_aligned_fasta_df(aligned_fasta_as_pandas_df_to_crop, cutoff):
    columns_to_drop = []
    for i in list(aligned_fasta_as_pandas_df_to_crop):
        # if there is a gap in the column at the beginning
        #add functionality to incorporate the cutoff.
        # this cutoff will mean that the proportion of the gaps in a column must be greater than this for the
        # column to be removed from the df. Once a column has been reached that has a proportion lower than this value
        # cropping will be complete for that end.
        list_of_col_characters = aligned_fasta_as_pandas_df_to_crop.loc[:, i].values.tolist()
        count_of_gaps = list_of_col_characters.count('-') + list_of_col_characters.count('*')
        total_characters = len(list_of_col_characters)
        prop_of_gaps = count_of_gaps / total_characters
        if prop_of_gaps >= cutoff:
            columns_to_drop.append(i)
        else:
            break
    for i in reversed(list(aligned_fasta_as_pandas_df_to_crop)):
        # if there is a gap in the column at the end
        list_of_col_characters = aligned_fasta_as_pandas_df_to_crop.loc[:, i].values.tolist()
        count_of_gaps = list_of_col_characters.count('-') + list_of_col_characters.count('*')
        total_characters = len(list_of_col_characters)
        prop_of_gaps = count_of_gaps / total_characters
        if prop_of_gaps >= cutoff:
            columns_to_drop.append(i)
        else:
            break

    # get a list that is the columns indices that we want to keep
    col_to_keep = [col_index for col_index in list(aligned_fasta_as_pandas_df_to_crop) if col_index not in columns_to_drop]

    # drop the gap columns
    return aligned_fasta_as_pandas_df_to_crop[col_to_keep]

def pandas_df_to_fasta(cropped_fasta_df):
    temp_fasta = []
    for ind in cropped_fasta_df.index.tolist():
        temp_fasta.extend(['>{}'.format(ind), ''.join(list(cropped_fasta_df.loc[ind]))])
    return temp_fasta

if __name__ == "__main__":
    main()