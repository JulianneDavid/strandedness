#!/usr/bin/env python2

"""
Python 2.7 code for collecting p-values for strandedness of RNA-seq experiments

Required input:
    .csv File from SRA with SRA accession numbers to check.
    Reference genome for hisat2 alignment

Optional input:
    Path to fastq-dump
    Path to hisat2
    Path to directory to store output files
    Number of "useful" reads required after running.
    Multiplier: how many reads are believed to be necessary to download to
        obtain the desired number of junction reads.
    Max attmpts: how many times to try fastq-dump read download, and hisat2
        read alignment, before moving on to the next set of downloads or the
        next SRA accession number.
    Logging level: INFO is the only option supported right now.

Improvements to be made:
    - checking for paired-end experiment - is just the "PAIRED" tag OK?
"""

import argparse
import csv
from datetime import datetime
import gzip
from itertools import groupby
import logging
import mmh3
from operator import itemgetter
import os
import random
import re
import subprocess as sp
from scipy import stats
from time import sleep


class TimeOut(Exception): pass


def get_hisat_input(required, multiplier, total, fastq_path, acc, output,
                    pairedtag, fail_file):
    """Samples & downloads fastq reads, and prepares them for a hisat2 -12 run.

    Input:
        required: the target number of "useful"/junction reads desired (int)
        multiplier: a multiplier to correct for the fact that most reads are
            not junction reads, and some will not be downloaded due to low
            quality scores (int)
        total: the total number of spots in the SRA experiment available to
            sample (int)
        fastq_path: the path to fastq-dump (string)
        acc: the SRA accession number for the current experiment (string)
        output: the output path for writing out (string)
        paired_tag: whether the experiment is paired or not, only necessary for
            processing multiple reads per spot (true/false)
        fail_file: the file to which to write any failed accession numbers for
            re-running later.

    A list of unique random numbers between 1 and the total number of spots is
    generated, then the reads at those spots are downloaded with fastq-dump.
    Each downloaded read is then formated appropriately for hisat2 with the -12
    option, which requires the following form: one read or read pair per line,
    tab separated as follows: read name, read 1 alignment, read 1 quality
    scores, then (if applicable) read 2 alignment and read 2 quality scores.

    Writes a file, SRAx_spots.txt, with the sampled spots used to get data.

    Returns all of the sampled reads formatted for batch alignment by hisat2,
    and the list of random spots (to be saved for later reference).
    """
    max_time = 300
    dl_start = datetime.now()
    spot_path = os.path.join(output, '{}_spots.txt'.format(acc))
    with open(spot_path, 'w') as spot_file:
        required_spots = required * multiplier
        num_bins = 100
#         num_bins = 10
        bin_spots = required_spots/num_bins
        bin_size = total/num_bins
        bin = 1
        read_format = []
        while bin <= num_bins:
            bin_start = (bin - 1) * bin_size + 1
            bin_stop = bin * bin_size
            bin += 1
            start_spot = random.randint(bin_start, bin_stop - bin_spots)
            spot_file.write('{}\n'.format(start_spot))
            spot = str(start_spot)
            stop_spot = str(start_spot + bin_spots)
            bin_time = 0
            bin_start = datetime.now()
            delay = 1
            attempt = 1
            back_off = 3
            while bin_time <= max_time:
                try:
                    fastq = sp.check_output(['{}'.format(fastq_path), '-I',
                                             '-B', '-W', '-E', '--split-spot',
                                             '--skip-technical', '-N', spot,
                                             '-X', stop_spot, '-Z', acc])
                    break
                except:
                    logging.info('acc {} failed: attempt {}'.format(acc,
                                                                    attempt))
                    delay = delay * back_off
                    logging.info('waiting {} sec before retry'.format(delay))
                    sleep(delay)
                    attempt += 1
                    bin_time = (datetime.now() - bin_start).total_seconds()
            else:
                logging.info('Accession number {} download timed out.  Moving '
                             'on to the next accession number.\n'.format(acc))
                with open(fail_file, 'w') as failed:
                    failed.write('{}\n'.format(acc))
                raise TimeOut('This SRA accession number download has timed '
                              'out.  See the standard output for the '
                              'fastq-dump error messages.  Moving on to the '
                              'next accession number.')

            lines = fastq.split('\n')
            format_lines = []
            last_line = 4 * (pairedtag + 1)
            for i, line in enumerate(lines, 1):
                if i % 2 == 0:
                    format_lines += [line]
                elif (i - 1) % last_line == 0:
                    format_lines += [line]
                if i % last_line == 0:
                    read_format.extend(['\t'.join(format_lines) + '\n'])
                    format_lines = []
        read_input = ''.join(read_format)
        hisat_formatted_input = read_input
        dl_end = datetime.now()
        time_difference = dl_end - dl_start
        elapsed_time = time_difference.total_seconds()
        logging.info('total download time was {} seconds'.format(elapsed_time))
        return hisat_formatted_input


def align_reads(hisat2_path, reference_genome, outpath, acc, reads):
    """Returns SAM format reads aligned by hisat2.
    Input:
        hisat2_path: the path to hisat2 (string)
        reference_genome: the path to the reference genome to be used for
            alignment (string)
        reads: fastq reads correctly formatted for the hisat2 -12 option
            (string, one read per line tab separated, name seq qual if
            single-end or name seq qual seq qual if paired-end)
        temp_dir: dir in which to put temporary fastq

    Returns a list of aligned reads in SAM format.
    """
    align_start = datetime.now()
    reads_path = os.path.join(outpath, '{}_reads.sam.gz'.format(acc))
    align_command = ('set -exo pipefail; {h2} --no-head --12 - -x {ref} | gzip'
                     ' > {file}').format(h2=hisat2_path, ref=reference_genome,
                                         file=reads_path)
    align_process = sp.Popen(align_command, stdin=sp.PIPE, stderr=sp.PIPE,
                             shell=True, executable='/bin/bash')
    align_process.communicate(input=reads)

    align_end = datetime.now()
    time_difference = align_end - align_start
    elapsed_time = time_difference.total_seconds()
    logging.info('total alignment time was {}'.format(elapsed_time))
    return reads_path, align_process.returncode


def filter_alignments(alignments, paired_tag):
    """Returns only alignments with highest alignment scores.

    Input:
        alignment_list: list of all returned alignments for one read or pair.
        paired_tag: whether the experiment is paired end or single end.

    Checks each alignment for AS:i: tag, then filters the list by highest
    possible AS:i: value for the group of alignments.

    Then, if the experiment is paired-end: chooses randomly whether to pick
    first reads or second reads, then checks the SAM flag & 64 for first reads.
    (This assumes that the "PAIRED" entry in the SRA csv is correct, and that
    all reads here will be paired-end, and that therefore, any read that does
    not have flag & 64 must have flag & 128, i.e. be a second read.)

    Returns a list of alignments with the highest alignment scores, either
    first or second read only if paired end.
    """
    # results = re.findall('(^.*AS:i:([+-]?\d+).*$)', '\n'.join(alignment_list),
    #                      flags=re.M)
    # if not results:
    #     return 0
    # max_score = max(results, key=itemgetter(1))[1]
    # alignments = filter(lambda item: item[1] == max_score, results)
    # primary_alignments = [item[0] for item in alignments]
    # return primary_alignments
    alignment_scores = []
    for entry in alignments:
        for tag in entry.split('\t'):
            if tag[:5] == 'AS:i:':
                alignment_scores.extend([int(tag[5:])])
    if not alignment_scores:
        return 0
    max_alignment_score = max(alignment_scores)
    primary_alignments = [alignments[i] for i, alignment_score
                          in enumerate(alignment_scores)
                          if alignment_score == max_alignment_score]

    if paired_tag:
        first_reads = random.getrandbits(1)
        single_end = []
        for entry in primary_alignments:
            read = entry.split('\t')
            flag = int(read[1])
            if not flag & 1:
                break
            if first_reads == (flag & 64 == 64):
                single_end.extend([entry])
        prepared_alignments = single_end
    else:
        prepared_alignments = primary_alignments
    return prepared_alignments


def read_sense(SAM_flag, plus_or_minus):
    """Checks a read's SAM flag and XS:A: tag, and determines its "direction".

    Input the SAM flag (int) and XS:A:? tag (string) from the aligned read.

    We have two states, arbitrarily called "sense" and "antisense," indicating
    whether all the first/second reads align with a gene or with its reverse
    complement, or whether this is random.

    Return the read's "sense"ness - if bit = 1, the read is "sense", otherwise
    it is "antisense."
    """
    fwd_gene = plus_or_minus[0] == 'XS:A:+'
    paired = SAM_flag & 1 == 1
    rev_read = SAM_flag & 16 == 16
    first_read = SAM_flag & 64 == 64
    return (fwd_gene != rev_read) == (first_read or not paired)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Determine sample '
                                                 'strandedness.')
    parser.add_argument('--sra-file', '-s', required=True, help='File with '
                        'SRA accession numbers to be downloaded and '
                        'checked for strandedness.')
    parser.add_argument('--ref-genome', '-r', required=True, help='Path to '
                        'reference genome to use for the aligner.')
    parser.add_argument('--fastq-dump-path', '-f', default='fastq-dump',
                        help='specify the path for fastq-dump')
    parser.add_argument('--hisat-path', '-p', default='hisat2',
                        help='specify the path for hisat2')
    parser.add_argument('--output-path', '-o', default='./', help='give path '
                        'for output files: sampled spots, aligned junction '
                        'reads, and SRA numbers with their p-values.')
    parser.add_argument('--required-reads', '-n', type=int, default=100,
                        help='give the target number of useful reads.')
    parser.add_argument('--multiplier', '-m', type=int, default=10, help='a '
                        'multiplier for generating the number of reads to '
                        'download from fastq-dump, to account for quality '
                        'filtering and for not all reads being useful.')
    parser.add_argument('--max-attempts', '-a', type=int, default=3, help='the'
                        ' number of times to attempt hisat2 alignment one one '
                        'SRA accession number after alingment failure before '
                        'continuing.')
    parser.add_argument('--log-level', '-l', choices=['DEBUG', 'INFO', 'ERROR'
                                                      'WARNING', 'CRITICAL'],
                        default='INFO', help='choose what logging mode to run')

    args = parser.parse_args()
    sra_file = args.sra_file
    ref_genome = args.ref_genome
    fastq_dump = args.fastq_dump_path
    hisat2 = args.hisat_path
    out_path = args.output_path
    required_reads = args.required_reads
    read_multiplier = args.multiplier
    max_attempts = args.max_attempts
    log_mode = args.log_level

    minimum_spots = 10000000
    useful = True
    name_tag = os.path.basename(sra_file).split('.')[0]
    now = str(datetime.now())
    log_file = os.path.join(out_path, '{}_{}_log.txt'.format(name_tag, now))
    logging.basicConfig(filename=log_file, level=log_mode)
    fail_file = os.path.join(out_path, 'failed_expts_{}.txt'.format(name_tag))
    pv_rand = os.path.join(out_path, 'random_pvals_{}.txt'.format(name_tag))
    pv_weigh = os.path.join(out_path, 'weighted_pvals_{}.txt'.format(name_tag))
    with open(sra_file) as sra_array, \
         open(pv_rand, 'w', 0) as pval_rand_file,\
         open (pv_weigh, 'w', 0) as pval_weigh_file:
        sra_array.next()
        csv_reader = csv.reader(sra_array)
        for experiment in csv_reader:
            sra_acc, num_spots, paired = (experiment[0], int(experiment[3]),
                                          experiment[15] == 'PAIRED')

            # to filter out single cell reads
            if num_spots < minimum_spots:
                logging.info('SRA {} skipped due to too few spots\n'
                             ''.format(sra_acc))
                continue

            seed = mmh3.hash(sra_acc)
            random.seed(seed)
            logging.info('the random seed is {}'.format(seed))

            try:
                hisat_input = get_hisat_input(required_reads, read_multiplier,
                                              num_spots, fastq_dump, sra_acc,
                                              out_path, paired, fail_file)
            except TimeOut:
                continue

            attempt = 0
            while attempt < max_attempts:
                reads_path, failure = align_reads(hisat2, ref_genome, out_path,
                                                  sra_acc, hisat_input)
                attempt += 1
                if not failure:
                    break
            else:
                logging.info('alignment failed for SRA {}'.format(sra_acc))
                continue

            random_sense = 0
            random_checked = 0
            weighted_sense = 0
            weighted_checked = 0
            with gzip.open('{}'.format(reads_path), 'r') as aligned_reads:
                sort_by_names = lambda x: x.split('\t')[0]
                for key, group in groupby(aligned_reads, sort_by_names):
                    alignments = list(group)
                    primary_alignments = filter_alignments(alignments, paired)
                    if not primary_alignments:
                        continue

                    num_primaries = float(len(primary_alignments))
                    weight = 1 / num_primaries
                    random.shuffle(primary_alignments)
                    one_read = not useful
                    for entry in primary_alignments:
                        read = entry.split('\t')
                        flag = int(read[1])
                        XS_A = re.findall('XS:A:[+-]', entry)
                        if XS_A:
                            weighted_sense += read_sense(flag, XS_A) * weight
                            weighted_checked += weight
                            if not one_read == useful:
                                random_sense += read_sense(flag, XS_A)
                                random_checked += 1
                                one_read = useful

            # stats.binom_test is a 2-sided & symmetrical cdf calculation:
            # same result whether sense or antisense is used
            r = 0.5
            random_anti = random_checked - random_sense
            random_p = stats.binom_test(random_sense, random_checked, r)
            weighted_anti = weighted_checked - weighted_sense
            weighted_p = stats.binom_test(weighted_sense, weighted_checked, r)

            logging.info('The SRA accession number is {}'.format(sra_acc))
            logging.info('pval file record is {}'.format(name_tag))
            logging.info('Drawing one random primary alignment, there are:')
            logging.info('{} sense reads.'.format(random_sense))
            logging.info('{} antisense reads.'.format(random_anti))
            logging.info('{} junction reads.'.format(random_checked))
            logging.info('The random p-value is {}'.format(random_p))
            logging.info('Looking at all primary alignments, there are:')
            logging.info('{} weighted sense reads.'.format(weighted_sense))
            logging.info('{} weighted antisense reads.'.format(weighted_anti))
            logging.info('The weighted p-value is {}.\n'.format(weighted_p))

            print >>pval_rand_file, '{},{}'.format(sra_acc, random_p)
            print >>pval_weigh_file, '{},{}'.format(sra_acc, weighted_p)
            # pval_file.write('{},{}\n'.format(sra_acc, p_value))
