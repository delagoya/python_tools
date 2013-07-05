#!/usr/bin/env python
from __future__ import print_function, division
from glob import glob
import os
import re
import argparse
import logging
import csv
import locale

def get_arguments():
    '''Parses the CLI arguments'''
    args = argparse.ArgumentParser()
    args.add_argument(
        "glob_pattern",
        type=str,
        help="A glob pattern like you would give to 'ls' to find a set of RUM mapping_stats.txt files. Example: '**/Sample_*/mapping_stats.txt'"
    )
    args.add_argument(
        "sample_name_regex",
        type=str,
        help="A regular expression pattern to pull the sample name from the STAR result directories' names. Example '(Sample_\w+)\/'"
    )
    args.add_argument(
        '--output','-o',
        default="mapping_stat_report.csv",
        type=str,
        help="An output CSV file name. Default: 'mapping_stat_report.csv'"
    )
    args.add_argument(
        '--debug', '-d',
        action='store_true',
        help="Print debugging information"
    )
    args.add_argument(
        '--verbose', '-v',
        action='store_true',
        help="Print verbose information. Basically this will print the mapping stats table to your screen."
    )
    return args.parse_args()

def setup_logging(args):
    '''Sets up normal and verbose logging'''
    logging.basicConfig(level=logging.ERROR,
                        format='%(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

def format_percentage(num):
    '''Formats a percentage for printing'''
    logging.debug("Formatting as %: " + str(num))
    return ("%0.2f" % (num * 100))  + "%"

def format_integer(num):
    '''Formats an integer for printing'''
    logging.debug("Formatting as Integer: " + str(num))
    return locale.format('%d', num, grouping=True)

def get_sample_name(pathname, args):
    logging.debug("Grabbing sample name from basedir: " + pathname)
    sn  = None
    if os.path.exists(pathname):
        m = re.compile(args.sample_name_regex).search(pathname)
        if m:
            sn = m.group(1)
        else:
            log.error("No sample name found from {0} using {1}".format(pathname,args.sample_name_regex))
    logging.debug("Grabbed sample_name: " + sn)
    return sn

def main():
    '''Run routine to pull mapping statistics from RUM mapping_stats.txt files.
    '''
    args = get_arguments()
    logging.debug("ARGS: {0}".format(args))
    locale.setlocale(locale.LC_ALL,'en_US')
    setup_logging(args)
    report = csv.writer(open(args.output,'wb'))

    # grab the list of mapping stats files from the glob pattern
    mapping_stats_files =  glob(args.glob_pattern)
    mapping_stats_files.sort()

    min_fragments_num = float('inf')
    max_fragments_num = 0
    min_aligned_num = float('inf')
    max_aligned_num = 0
    min_uniq_num = float('inf')
    max_uniq_num = 0
    min_nu_num = float('inf')
    max_nu_num = 0
    min_aligned_percentage = float('inf')
    max_aligned_percentage = 0
    min_uniq_percentage = float('inf')
    max_uniq_percentage = 0
    min_nu_percentage = float('inf')
    max_nu_percentage = 0

    headers = [
        'SampleName',
        "Num_fragments",
        "Num_uniq_aligned",
        "Num_NU_aligned",
        "Num_ALL_aligned",
        "%_aligned_uniq",
        "%_aligned_NU",
        "%_aligned_ALL"
    ]

    report.writerow(headers)
    logging.info("\t".join(headers))
    logging.info("-" * 100)

    for i,ms_fn in enumerate(mapping_stats_files):
        logging.debug("Grabbing mapping stats from: " + ms_fn)
        # get the sample name
        sample_name = get_sample_name(ms_fn, args)

        # now read in the mapping_stats and pull the information
        data = open(ms_fn,'rU').readlines()

        # grab the total reads and percentage mapped for this data set
        total_fragments_num = re.search(r'Number of input reads.+\|\s+(\d+)',data[5]).group(1)
        total_fragments_num = int(total_fragments_num)
        uniq_num = int(re.search(r'Uniquely mapped reads number\s+\|\s+(\d+)', data[8]).group(1))
        nu_num = int(re.search(r'Number of reads mapped to multiple loci\s+\|\s+(\d+)', data[23]).group(1))
        nu_num += int(re.search(r'Number of reads mapped to too many loci\s+\|\s+(\d+)', data[25]).group(1))
        total_aligned_num = uniq_num + nu_num

        total_aligned_percentage = total_aligned_num / total_fragments_num
        uniq_percentage = uniq_num / total_fragments_num
        nu_percentage = nu_num / total_fragments_num

        logging.debug(sample_name + ": total reads = " + format_integer(total_fragments_num))
        logging.debug(sample_name + ": total aligned num = " + format_integer(total_aligned_num))
        logging.debug(sample_name + ": total aligned % = " + format_percentage(total_aligned_percentage))

        # grab the number of aligned reads for chrM and contigs
        # Adjust for Min & Max numbers
        min_fragments_num = min(min_fragments_num, total_fragments_num)
        max_fragments_num = max(max_fragments_num, total_fragments_num)


        min_aligned_num = min(min_aligned_num, total_aligned_num)
        max_aligned_num = max(max_aligned_num, total_aligned_num)
        min_aligned_percentage = min(min_aligned_percentage, total_aligned_percentage)
        max_aligned_percentage = max(max_aligned_percentage, total_aligned_percentage)


        min_uniq_num = min(min_uniq_num, uniq_num)
        max_uniq_num = max(max_uniq_num, uniq_num)
        min_uniq_percentage = min(min_uniq_percentage, uniq_percentage)
        max_uniq_percentage = max(max_uniq_percentage, uniq_percentage)

        min_nu_num = min(min_nu_num, nu_num)
        max_nu_num = max(max_nu_num, nu_num)
        min_nu_percentage = min(min_nu_percentage, nu_percentage)
        max_nu_percentage = max(max_nu_percentage, nu_percentage)

        # now print out the information
        report.writerow([sample_name,
            total_fragments_num,
            total_aligned_num,
            uniq_num,
            nu_num,
            total_aligned_num,
            uniq_percentage,
            nu_percentage,
            total_aligned_percentage,
            ])
        logging.info("\t".join([sample_name,
            format_integer(total_fragments_num),
            format_integer(uniq_num),
            format_integer(nu_num),
            format_integer(total_aligned_num),
            format_percentage(uniq_percentage),
            format_percentage(nu_percentage),
            format_percentage(total_aligned_percentage)]))

    report.writerow(['Minimums',
        min_fragments_num,
        min_uniq_num,
        min_nu_num,
        min_aligned_num,
        min_uniq_percentage,
        min_nu_percentage,
        min_aligned_percentage])

    report.writerow(['Maximums',
        max_fragments_num,
        max_uniq_num,
        max_nu_num,
        max_aligned_num,
        max_uniq_percentage,
        max_nu_percentage,
        max_aligned_percentage])

    logging.info("-" * 100 )
    logging.info("\t".join(['Minimums',
        format_integer(min_fragments_num),
        format_integer(min_uniq_num),
        format_integer(min_nu_num),
        format_integer(min_aligned_num),
        format_percentage(min_uniq_percentage),
        format_percentage(min_nu_percentage),
        format_percentage(min_aligned_percentage)]))
    logging.info("\t".join(['Maximums',
        format_integer(max_fragments_num),
        format_integer(max_uniq_num),
        format_integer(max_nu_num),
        format_integer(max_aligned_num),
        format_percentage(max_uniq_percentage),
        format_percentage(max_nu_percentage),
        format_percentage(max_aligned_percentage)]))

if __name__ == '__main__':
    main()
