#!/usr/bin/env python2

"""
Python 2.7 code to make a stranded fastq file randomly, artificially unstranded

Required input:
    .fastq file: stranded

Optional input:
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
    parser.add_argument('--output-directory', '-o', default='./',
                        help='Desired path for storing output FASTA files.')

    args = parser.parse_args()
    fastq_file = args.fastq_file
    out_path = args.output_directory

    name_tag = os.path.basename(fastq_file).split('.')[0]
    shuffled_fastq = os.path.join(out_path,
                                  '{}_shuffled.fastq'.format(name_tag))
    seed = mmh3.hash(name_tag)
    random.seed(seed)
    line_number = 0
    coinflip = random.randint(0, 1)
    with open(fastq_file) as fastq, open(shuffled_fastq, 'w') as out_fastq:
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
