#!/usr/bin/env python2

"""
Python 2.7 code for performing a Benjamini-Hochberg procedure to correct for
multiple-testing false discovery rate.
Required input:
    File containing SRA (or other experiment ID numbers) and their p-values
Optional input:
    alpha value - desired significance level for discovery of strandedness/
        rejection of null hypothesis (unstrandedness)
    Path to output directory to store list of experiments declared "stranded"
"""

import argparse
from datetime import datetime
import glob
import logging
from operator import itemgetter
import os


def get_sorted_pvals(file_list):
    '''Returns a list of SRA accession numbers and their p-values from file.
    Input the path to a file with SRA numbers and p-values.
    Returns a list of data from the file, sorted from lowest to highest
    p-value.
    '''
    list = []
    ticker = 0
    for file in file_list:
        with open(file, 'r') as pval_file:
            for line in pval_file:
                if not line == '\n':
                    line_list = line.split(',')
                    line_list[1] = float(line_list[1].strip('\n'))
                    # line_list[1] = float(line_li)
                    if line_list not in list:
                        list.append(line_list)
                    # list.append(line_list)
                    else:
                        ticker += 1
    list.sort(key=itemgetter(1))
    logging.info('The total number of experiments tested is:')
    logging.info(len(list))
    logging.info('The sorted pvalues are:')
    logging.info(list)
    return list


def BH_correction_procedure(pvals, alpha):
    '''Performs Benjamini-Hochberg multiple testing correction procedure.
     Input list of SRA numbers and their p-values, sorted from smallest to
     largest p-value, and the acceptable significance level for declaring a
     discovery (alpha).
     Returns whether or not a discovery can be found, and if discovery exists,
     the list number that represents the last (highest p-value) discovery.
     '''
    m = len(pvals)
    for k_rev, expt in enumerate(list(reversed(pvals))):
        k = m - k_rev
        p_value = float(expt[1])
        if p_value <= alpha * k / m:
            return True, k
        elif k == 1:
            return False, k


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Determine sample '
                                                 'strandedness.')
    parser.add_argument('--pvaluedir', '-p', required=True, help='Folder with '
                        'SRA accession numbers and their p-values; generated'
                        'by strandedness.py.')
    parser.add_argument('--alpha', '-a', type=float, default=0.05,
                        help='Acceptable false discovery rate for strandedness'
                        ' call; must be between 0 and 1.')
    parser.add_argument('--outputpath', '-o', default='./', help='Give path '
                        'for output files: junction reads and p-values.')
parser.add_argument('--log-level', '-l', choices=['DEBUG', 'INFO', 'ERROR'
                                                  'WARNING', 'CRITICAL'],
                    default='INFO', help='Choose what logging mode to run')

args = parser.parse_args()
alpha = float(args.alpha)
pval_file_path = args.pvaluedir
out_path = args.outputpath
log_mode = args.log_level

name_tag = "BH_analyzer"
now = str(datetime.now())
log_file = os.path.join(out_path, '{}_{}_log.txt'.format(name_tag, now))
logging.basicConfig(filename=log_file, level=log_mode)

random_path = os.path.join(pval_file_path, '*/output/random*')
weighted_path = os.path.join(pval_file_path, '*/output/weighted*')
random_files = glob.glob(random_path)
weighted_files = glob.glob(weighted_path)

weighted_pvals = get_sorted_pvals(weighted_files)
random_pvals = get_sorted_pvals(random_files)

weighted_disc, weighted_large_k = BH_correction_procedure(weighted_pvals,
                                                          alpha)
random_disc, random_large_k = BH_correction_procedure(random_pvals, alpha)

logging.info('For weighted pvals:')
if weighted_disc:
    stranded_expts_w = weighted_pvals[:weighted_large_k]
    logging.info('BH complete, largest k is {}'.format(weighted_large_k))
    logging.info('Weighted stranded experiments are these:')
    logging.info(stranded_expts_w)
    weighted_outpath = os.path.join(out_path, 'weighted_stranded_SRAs.txt')
    with open(weighted_outpath, 'w') as weighted_file:
        for set in stranded_expts_w:
            weighted_file.write('{},{}\n'.format(set[0], set[1]))
else:
    logging.info('BH complete, all weighted experiments were unstranded.')

logging.info('For random pvals:')
if random_disc:
    stranded_expts_r = random_pvals[:random_large_k]
    logging.info('BH complete, largest k is {}'.format(random_large_k))
    logging.info('Random stranded experiments are these:')
    logging.info(stranded_expts_r)
    random_outpath = os.path.join(out_path, 'random_stranded_SRAs.txt')
    with open(random_outpath, 'w') as random_file:
        for set in stranded_expts_r:
            random_file.write('{},{}\n'.format(set[0], set[1]))
else:
    logging.info('BH complete, all random experiments are unstranded.')
