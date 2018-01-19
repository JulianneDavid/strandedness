#!/usr/bin/env python2

"""
Python 2.7 code to make a stranded fastq file randomly, artificially unstranded

Required input:
    .fastq file: stranded

Optional input:
    .fastq file of paired end reads, if the experiment is paired end.
    output directory: where the artificially unstranded fastq should be stored.

"""

import argparse
import mmh3
import os
import random
import string


def reverse_complement(dna_sequence):
    complement = string.maketrans("ATCG", "TAGC")
    rev_comp = dna_sequence.translate(complement)[::-1]
    return rev_comp


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Determine sample '
                                                 'strandedness.')
    parser.add_argument('--fastq-file', '-f', required=True, help='Stranded '
                        'RNA-seq experiment .fastq file to artificially '
                        'unstrand.')
    parser.add_argument('--paired-fastq', '-p', default='', help='If the '
                        'experiment is stranded, give the file containing the '
                        'paired end reads here.')
    parser.add_argument('--output-directory', '-o', default='./',
                        help='Desired path for storing output FASTA files.')

    args = parser.parse_args()
    fastq_file_1 = args.fastq_file
    fastq_file_2 = args.paired_fastq
    out_path = args.output_directory

    name_tag = os.path.basename(fastq_file_1).split('.')[0]
    seed = mmh3.hash(name_tag)
    random.seed(seed)
    line_number = 0
    coinflip = random.randint(0, 1)

    if fastq_file_2:
        shuffled_fastq1 = os.path.join(out_path,
                                       '{}_shuffled_1.fastq'.format(name_tag))
        shuffled_fastq2 = os.path.join(out_path,
                                       '{}_shuffled_2.fastq'.format(name_tag))
        with open(fastq_file_1) as fastq1, \
             open(fastq_file_2) as fastq2,\
             open(shuffled_fastq1, 'w') as out_fastq_1, \
             open(shuffled_fastq2, 'w') as out_fastq_2:
            for line1, line2 in zip(fastq1, fastq2):
                line_number += 1
                line1 = line1.strip('\n')
                line2 = line2.strip('\n')
                if line_number == 5:
                    line_number = 1
                    coinflip = random.randint(0, 1)
                if line_number == 2 and coinflip:
                    line1 = reverse_complement(line1)
                    line2 = reverse_complement(line2)
                if line_number == 4 and coinflip:
                    line1 = line1[::-1]
                    line2 = line2[::-1]
                print >>out_fastq_1, '{}'.format(line1)
                print >>out_fastq_2, '{}'.format(line2)

    else:
        shuffled_fastq = os.path.join(out_path,
                                      '{}_shuffled.fastq'.format(name_tag))
        with open(fastq_file_1) as fastq, \
             open(shuffled_fastq, 'w') as out_fastq:
            for line in fastq:
                line_number += 1
                line = line.strip('\n')
                if line_number == 5:
                    line_number = 1
                    coinflip = random.randint(0, 1)
                if line_number == 2 and coinflip:
                    line = reverse_complement(line)
                if line_number == 4 and coinflip:
                    line = line[::-1]
                print >>out_fastq, '{}'.format(line)
