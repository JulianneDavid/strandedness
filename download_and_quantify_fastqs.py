#!/usr/bin/env python2

"""
Python 2.7 code for processing randomly sampled "stranded"-called experiments.

Required input:
    .txt file of called "stranded" experiment accession numbers generated by
        "exacloud_p_val_analyzer.py"

Optional input:


This code will randomly sample the stranded experiments, use fastq-dump to
collect the fastq files, artificially unstrand them, and then perform
quantification with Salmon.
"""


import argparse
from datetime import datetime
import logging
import mmh3
import os
import random
import string
import subprocess as sp
from time import sleep


class TimeOut(Exception): pass


def sample_stranded_experiments(sra_containing_file, sample_size):
    """Randomly samples accession numbers from a previously generated list
    Input:
        sra_containing_file: generated by p_val_analyzer
        sample_size: the number of sra numbers to randomly collect

    Returns a list of random SRA accession numbers previously classified as
    stranded.
    """
    name_tag = os.path.basename(sra_containing_file)
    seed = mmh3.hash(name_tag)
    random.seed(seed)
    with open(sra_containing_file) as sra_file:
        lines = random.sample(sra_file.readlines(), sample_size)

    accession_numbers = []
    for line in lines:
        accession_numbers.append(line.split(',')[0])
    return accession_numbers


def collect_fastq_files(fastq_path, accs, fail_file, success_file):
    """Downloads fastq reads and returns their directory

    Input:
        fastq_path: the path to fastq-dump (string)
        accs: list of SRA accession numbers to download (list of strings)
        output: the output path for writing out (string)
        fail_file: the file to which to write any failed accession numbers for
            re-running later. (string)

    Returns the path where the new fastq files are stored.
    """
    existing_fastqs = []
    with open(success_file, 'r') as success:
        for line in success:
            existing_fastqs.append(line.strip('\n'))
            
    # downloaded_accs 
    max_time = 3000
    dl_start = datetime.now()
    fastq_output = os.path.join(out_path, 'original_fastq_files')
    back_off = 3
    for acc in accs:
        logging.info('current acc is {}'.format(acc))
        if acc in existing_fastqs or acc in ['SRR5575952', 'SRR2960573']:
            logging.info('skipping acc {}'.format(acc))
            continue
        logging.info('\ncollecting fastq {}'.format(acc))
        bin_time = 0
        bin_start = datetime.now()
        delay = 1
        attempt = 1
        while bin_time <= max_time:
            try:
                fastq_stdout = sp.check_output(['{}'.format(fastq_path), '-I',
                                                '-B', '-W', '-E',
                                                '--split-files',
                                                '--skip-technical', '-O',
                                                fastq_output, acc])
                logging.info('fastq-dump was successful; std output was:')
                logging.info(fastq_stdout)
                with open(success_file, 'a') as success:
                    success.write('{}\n'.format(acc))
                break
            except:
                logging.info('acc {} failed: attempt {}'.format(acc, attempt))
                delay = delay * back_off
                logging.info('waiting {} sec before retry'.format(delay))
                sleep(delay)
                attempt += 1
                bin_time = (datetime.now() - bin_start).total_seconds()
        else:
            logging.info('Accession number {} download timed out.  Moving '
                         'on to the next accession number.\n'.format(acc))
            with open(fail_file, 'a') as failed:
                failed.write('{}\n'.format(acc))
            raise TimeOut('This SRA accession number download has timed '
                          'out.  See the standard output for the '
                          'fastq-dump error messages.  Moving on to the '
                          'next accession number.\n')

    dl_end = datetime.now()
    time_difference = dl_end - dl_start
    elapsed_time = time_difference.total_seconds()
    logging.info('total download time was {} seconds\n'.format(elapsed_time))
    return fastq_output


def reverse_complement(dna_sequence):
    """Creates and returns the reverse complement of a nucleotide sequence

    Input:
        dna_sequence: the nucleotide sequence to reverse complement

    Returns: the sequence's reverse complement
    """
    complement = string.maketrans("ATCG", "TAGC")
    rev_comp = dna_sequence.translate(complement)[::-1]
    return rev_comp


def artificially_unstrand(accession, fastq_path):
    """Makes a stranded fastq file randomly, artificially "unstranded."

    Required input:
        accession: SRA accession number for generating the fastq file names
        path for fastq file storage

    This function generates the path names for the fastq files, artificially, 
    randomly "unstrands" the reads, and writes the new fastq output to the
    fastq storage directory.

    Output:
        Tag labeling the experiment as paired or single-end.
    """
    # name_tag = os.path.basename(fastq_file_1).split('.')[0].split('_')[0]
    seed = mmh3.hash(accession)
    random.seed(seed)
    line_number = 0
    coinflip = random.randint(0, 1)

    fastq_file_1 = os.path.join(fastq_path, '{}_1.fastq'.format(accession))
    fastq_file_2 = os.path.join(fastq_path, '{}_2.fastq'.format(accession))

    if os.path.isfile(fastq_file_2):
        paired = True
        shuffled_fastq1 = os.path.join(fastq_path,
                                       '{}_shuffled_1.fastq'.format(accession))
        shuffled_fastq2 = os.path.join(fastq_path   ,
                                       '{}_shuffled_2.fastq'.format(accession))
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
        paired = False
        shuffled_fastq = os.path.join(fastq_path,
                                      '{}_shuffled_1.fastq'.format(accession))
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
    return paired


def call_salmon_quantification(salmon_path, salmon_ind, outpath, acc, paired):
    """Runs salmon in quantification mode
    Input:
        salmon_path: the path to salmon (string)
        salmon_ind: the path to the reference genome to be used for
            alignment (string)
        outpath: where we will store the new salmon quantification information
        acc: the base name of the fastq files to quantify
        paired: whether or not the experiment is paired end

    Runs Salmon quantification, and stores results in a directory specific to
    the accession number of the experiment.
    """
    logging.info('starting quantification for {} reads'.format(acc))
    quant_start = datetime.now()
    orig_out = os.path.join(outpath, 'salmon_output_{}_original'.format(acc))
    shuf_out = os.path.join(outpath, 'salmon_output_{}_unstranded'.format(acc))

    if paired:
        fq_orig_1 = os.path.join(fastq_path, '{}_1.fastq'.format(acc))
        fq_orig_2 = os.path.join(fastq_path, '{}_2.fastq'.format(acc))
        salmon_command_orig = ('set -exo pipefail; {sal} quant -i {ref} -l A '
                               '-1 {f} -2 {s} -o {out} --posBias --gcBias'
                               ).format(sal=salmon_path, ref=salmon_ind,
                                        f=fq_orig_1, s=fq_orig_2,
                                        out=orig_out)

        fq_shuf_1 = os.path.join(fastq_path, '{}_shuffled_1.fastq'.format(acc))
        fq_shuf_2 = os.path.join(fastq_path, '{}_shuffled_2.fastq'.format(acc))
        salmon_command_shuf = ('set -exo pipefail; {sal} quant -i {ref} -l A '
                               '-1 {f} -2 {s} -o {out} --posBias --gcBias'
                               ).format(sal=salmon_path, ref=salmon_ind,
                                        f=fq_shuf_1, s=fq_shuf_2,
                                        out=shuf_out)
    else:
        fq_orig_1 = os.path.join(fastq_path, '{}_1.fastq'.format(acc))
        salmon_command_orig = ('set -exo pipefail; {sal} quant -i {ref} -l A '
                               '-r {r} -o {out} --posBias --gcBias'
                               ).format(sal=salmon_path, ref=salmon_ind,
                                        r=fq_orig_1, out=orig_out)
        fq_shuf_1 = os.path.join(fastq_path, '{}_shuffled_1.fastq'.format(acc))
        salmon_command_shuf = ('set -exo pipefail; {sal} quant -i {ref} -l A '
                               '-r {r} -o {out} --posBias --gcBias'
                               ).format(sal=salmon_path, ref=salmon_ind,
                                        r=fq_shuf_1, out=shuf_out)

    salmon_process_orig = sp.Popen(salmon_command_orig, stdin=sp.PIPE,
                                   stderr=sp.PIPE, shell=True,
                                   executable='/bin/bash')
    salmon_process_orig.communicate()

    salmon_process_shuf = sp.Popen(salmon_command_shuf, stdin=sp.PIPE,
                                   stderr=sp.PIPE, shell=True,
                                   executable='/bin/bash')
    salmon_process_shuf.communicate()

    quant_end = datetime.now()
    time_difference = quant_end - quant_start
    elapsed_time = time_difference.total_seconds()
    logging.info('quantification for {} reads complete'.format(acc))
    logging.info('output stored in {} and {}'.format(orig_out, shuf_out))
    logging.info('salmon quantification time was {}\n'.format(elapsed_time))
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Determine sample '
                                                 'strandedness.')
    parser.add_argument('--stranded_list', '-s', required=True, help='File '
                        'with SRA accession numbers to be sampled, downloaded,'
                        ' and quantified')
    parser.add_argument('--fastq-dump-path', '-f', default='fastq-dump',
                        help='specify the path for fastq-dump')
    parser.add_argument('--salmon-path', '-q', default='salmon',
                        help='specify the path for salmon')
    parser.add_argument('--salmon-index', '-i', required=True, help='Provide '
                        'the path to the index to be used for salmon '
                        'quantification.')
    parser.add_argument('--sample-size', '-ss', type=int, default=100, 
                        help='Give the number of fastq files to download.')
    parser.add_argument('--output-path', '-o', default='./', help='give path '
                        'for output files: shuffled and unshuffled fastq '
                        'files and salmon quantification output.')
    parser.add_argument('--log-level', '-l', choices=['DEBUG', 'INFO', 'ERROR'
                                                      'WARNING', 'CRITICAL'],
                        default='INFO', help='choose what logging mode to run')
    parser.add_argument('--downloaded-sras', '-d', help='Give the path to file'
                        ' with already-downloaded sra accession numbers.')

    args = parser.parse_args()
    sra_file = args.stranded_list
    fastq_dump = args.fastq_dump_path
    salmon = args.salmon_path
    salmon_index = args.salmon_index
    sample_size = args.sample_size
    out_path = args.output_path
    log_mode = args.log_level
    success_file = args.downloaded_sras

    name_tag = os.path.basename(sra_file).split('.')[0]
    now = str(datetime.now())
    log_file = os.path.join(out_path, '{}_{}_log.txt'.format(name_tag, now))
    logging.basicConfig(filename=log_file, level=log_mode)
    fail_file = os.path.join(out_path, 'failed_expts_{}.txt'.format(name_tag))
    if success_file is None:
        success_file = os.path.join(out_path, 
                                    'downloaded_sras_{}.txt'.format(name_tag))
    pv_rand = os.path.join(out_path, 'random_pvals_{}.txt'.format(name_tag))
    pv_weigh = os.path.join(out_path, 'weighted_pvals_{}.txt'.format(name_tag))

    accession_numbers = sample_stranded_experiments(sra_file, sample_size)
    print('starting collection of fastq files')
    fastq_path = collect_fastq_files(fastq_dump, accession_numbers, fail_file,
                                     success_file)

    for acc in accession_numbers:
        paired = artificially_unstrand(acc, fastq_path)
        call_salmon_quantification(salmon, salmon_index, out_path, acc, paired)
