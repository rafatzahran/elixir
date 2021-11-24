# python3.6
# Given two types of text files of genomic information, using the SEGMENT format and the FUNCTION format.

import pandas as pd
import numpy as np
from argparse import ArgumentParser
import os.path
# ---------------------------------------------------------------
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error('The file %s does not exist!' % arg)
    elif not(arg.endswith('.s') or arg.endswith('.f')):
        parser.error('The file %s not accepted. File format must be in SEGMENT format \".s\" or FUNCTION format \".f\".' % arg)
    else:
        return arg


parser = ArgumentParser(description='Development task for Elixir.no')
parser.add_argument('-f1', dest='file_name1', required=True,
                    help='Input file, it can be SEGMENT format or FUNCTION format', metavar='FILE',
                    type=lambda x: is_valid_file(parser, x))
parser.add_argument('-f2', dest='file_name2', required=True,
                    help='Input file, it can be SEGMENT format or FUNCTION format', metavar='FILE',
                    type=lambda x: is_valid_file(parser, x))                   

args = parser.parse_args()
# ---------------------------------------------------------------
def files_with_function_format(args):
    return args.file_name1.endswith('.f') and args.file_name2.endswith('.f')

def files_with_segment_format(args):
    return args.file_name1.endswith('.s') and args.file_name2.endswith('.s')

def files_with_segment_and_function_format(args):
    return (args.file_name1.endswith('.s') and args.file_name2.endswith('.f')) or \
        (args.file_name1.endswith('.s') and args.file_name2.endswith('.f'))
# ---------------------------------------------------------------
"""
2 FUNCTION files: calculate the sample Pearson correlation coefficient of the two number lists.
"""
def get_sample_pearson_correlation_coefficient(args):
    file1 = open(args.file_name1, 'r')
    df_file1 = pd.read_csv(file1, header=None, names=["X"])

    file2 = open(args.file_name2, 'r')
    df_file2 = pd.read_csv(file2, header=None, names=["Y"])

    frames = [df_file1, df_file2]
    df = pd.concat(frames, axis=1)

    df['X_MINUS_MEAN'] = df['X'] - df['X'].mean()
    df['Y_MINUS_MEAN'] = df['Y'] - df['Y'].mean()

    df['(X_MINUS_MEAN)(Y_MINUS_MEAN)'] = df['X_MINUS_MEAN'] * df['Y_MINUS_MEAN']

    df['X_MINUS_MEAN_SQUARE'] = df['X_MINUS_MEAN'].pow(2)
    df['Y_MINUS_MEAN_SQUARE'] = df['Y_MINUS_MEAN'].pow(2)

    total1 = df['(X_MINUS_MEAN)(Y_MINUS_MEAN)'].sum()
    total2_1 = np.sqrt(df['X_MINUS_MEAN_SQUARE'].sum())
    total2_2 = np.sqrt(df['Y_MINUS_MEAN_SQUARE'].sum())

    return total1 / (total2_1 * total2_2)
# ---------------------------------------------------------------
"""
2 SEGMENT files: calculate the overlap (in number of positions) of the regions from file X.s with regions from file Y.s.
"""
def calculate_regions_overlap(args):
    file1 = open(args.file_name1, 'r')
    df_file1 = pd.read_csv(file1, header=None, delimiter=r"\s+", names=["start", "end"])

    file2 = open(args.file_name2, 'r')
    df_file2 = pd.read_csv(file2, header=None, delimiter=r"\s+", names=["start", "end"])

    segment_a_positions = set()
    for start, end in zip(df_file1['start'], df_file1['end']):
        segment_a_positions.update(range(start,end))

    segment_b_positions = set()
    for start, end in zip(df_file2['start'], df_file2['end']):
        segment_b_positions.update(range(start,end))

    return set.intersection(segment_a_positions, segment_b_positions)
# ---------------------------------------------------------------
"""
1 SEGMENT and one FUNCTION file: The mean of the numbers in the FUNCTION file 
whose positions are covered by the regions in the SEGMENT file.
"""
def get_mean_of_FUNCTION_covered_by_SEGMENT(args):
    file1 = open(args.file_name1, 'r')
    file2 = open(args.file_name2, 'r')
    df_file1 = pd.read_csv(file1, header=None, delimiter=r"\s+", names=["start", "end"]) if args.file_name1.endswith('.s') else \
        pd.read_csv(file2, header=None, delimiter=r"\s+", names=["start", "end"])
    df_file2 = pd.read_csv(file1, header=None, names=["Y"]) if args.file_name1.endswith('.f') else \
        pd.read_csv(file2, header=None, names=["Y"])

    segment_positions = set()
    for start, end in zip(df_file1['start'], df_file1['end']):
        segment_positions.update(range(start,end))
    return df_file2.loc[segment_positions,['Y']].mean()[0]
# ---------------------------------------------------------------
def main():
    if files_with_function_format(args):
        total = get_sample_pearson_correlation_coefficient(args)
        print(total)
    elif files_with_segment_format(args):
        segment_intersection = calculate_regions_overlap(args)
        #print(segment_intersection)
        print('Number of overlaped positions= ', len(segment_intersection))
    else:
        mean = get_mean_of_FUNCTION_covered_by_SEGMENT(args)
        print(mean)
# ---------------------------------------------------------------
main()