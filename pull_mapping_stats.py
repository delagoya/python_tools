#!/usr/bin/env python
from __future__ import print_function, division
import glob
import os, re, sys
import argparse
import logging
import csv
import locale
import prettytable
from prettytable import from_csv

def get_arguments():
    '''Parses the CLI arguments'''
    args = argparse.ArgumentParser()
    args.add_argument(
        "GLOB_PATTERN",
        type=str,
        help="A glob pattern like you would give to 'ls' to find a set of RUM mapping_stats.txt files. Example: '**/Sample_*/mapping_stats.txt'"
    )
    args.add_argument(
        "SAMPLE_REGEX",
        type=str,
        help="A regular expression pattern to pull the sample number from the mapping stats file path"
    )

    args.add_argument(
        '--output','-o',
        default="mapping_stat_report.csv",
        type=str,
        help="An CSV report output file path. Default: './mapping_stat_report.csv'"
    )

    args.add_argument(
        '--debug', '-d',
        action='store_true',
        help="Print debugging information"
    )

    args.add_argument(
        '--verbose', '-v',
        action='store_true',
        help="Print verbose information. Basically this will print a nicely formatted mapping stats table to your screen."
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

def get_sample_name(bdir,sample_regex=None):
    logging.debug("Grabbing sample name from basedir: " + bdir)
    sn  = None
    if os.path.exists(os.path.join(bdir,'rum_job_report.txt')):
        fn = os.path.join(bdir,'rum_job_report.txt')
        logging.debug("Grabbing sample name from file: " + fn )
        data = open(fn,'rU').readlines()
        sn= re.findall(r'Job name : (\S+)', data[9])[0]
    elif os.path.exists(os.path.join(bdir,'rum.log_master')):
        fn = os.path.join(bdir,'rum.log_master')
        logging.debug("Grabbing sample name from file: " + fn )
        data = open(fn,'rU').readlines()
        sn =  re.findall(r'name: (\S+)', data[6])[0]
    else:
        # use the basedir regex to get the sample name
        rgx = re.compile(sample_regex)
        sn = rgx.match(bdir).group(1)
    logging.debug("Grabbed sample_name: " + sn)
    return sn

def main():
    '''Run routine to pull mapping statistics from RUM mapping_stats.txt files.
    '''
    args = get_arguments()
    locale.setlocale(locale.LC_ALL,'en_US')
    setup_logging(args)
    report_file = open(args.output,'w')
    report = csv.writer(report_file)

    mapping_stats_files = []
    # grab the list of mapping stats files from the glob pattern
    for fn in glob.glob(args.GLOB_PATTERN):
        mapping_stats_files.append(fn)

    mapping_stats_files.sort()

    min_reads_num = float('inf')
    max_reads_num = 0
    min_aligned_num = float('inf')
    max_aligned_num = 0
    min_aligned_percentage = float('inf')
    max_aligned_percentage = 0
    min_chrm_percentage    = float('inf')
    max_chrm_percentage    = 0
    min_contig_percentage  = float('inf')
    max_contig_percentage  = 0

    # need to add percent of unique and the percent of non-unique mappers
    headers = ['SampleName',
        "Num_read_pairs",
        "Num_uniq_aligned",
        "Num_NU_aligned",
        "Num_ALL_aligned",
        "Num_chrM_aligned",
        "Num_contigs_aligned",
        "%_aligned_uniq",
        "%_aligned_NU",
        "%_aligned_ALL",
        "%_aligned_chrM",
        "%_aligned_contigs"]
    report.writerow(headers)

    for i,ms_fn in enumerate(mapping_stats_files):
        logging.debug("Grabbing mapping stats from: " + ms_fn)
        # ms_fn = os.path.realpath(ms_fn)

        # get the sample name
        base_dir = os.path.dirname(ms_fn)
        sample_name = get_sample_name(base_dir,args.SAMPLE_REGEX)

        # now read in the mapping_stats and pull the information
        data = open(ms_fn,'rU').readlines()

        # grab the total reads and percentage mapped for this data set
        total_reads_num = int(re.findall(r'Number of read pairs\: (\S+)',data[0])[0].replace(',',''))
        # uniq mapped
        int_percent_regex = re.compile(".+\:\s(\S+)\s\((\S+)\%\)$")
        m =  int_percent_regex.match(data[11])
        uniq_num = int(m.group(1).replace(",",''))
        uniq_percent = float(m.group(2))/100
        # NU mapped
        m =  int_percent_regex.match(data[18])
        nu_num = int(m.group(1).replace(",",''))
        nu_percent = float(m.group(2))/100
        # all mapped
        m =  int_percent_regex.match(data[26])
        all_num = int(m.group(1).replace(",",''))
        all_percent = float(m.group(2))/100

        # Grab number of uniq mapped reads per chromosome
        chrM = [0,0]
        chr_contigs =[0,0];
        index = 0
        for l in data[22:-1]:
            logging.debug(l)
            if not re.match(r'^chr',l):
                continue
            if re.match('Non-Uniquely',l):
                index = 1
            ch, n = l.split()
            n = int(n)
            if ch == 'chrM':
                chrM[index] = n
            if ch.find("_") > 0:
                chr_contigs[index] += n

        chrm_reads_percentage =  sum(chrM) / total_reads_num
        contig_reads_percentage = sum(chr_contigs) / total_reads_num

        # Adjust for Min & Max numbers
        min_reads_num = min(min_reads_num,total_reads_num)
        max_reads_num =  max(max_reads_num,total_reads_num)
        min_aligned_num =  min(min_aligned_num,all_num)
        max_aligned_num =  max(max_aligned_num,all_num)
        min_aligned_percentage = min(min_aligned_percentage,all_percent)
        max_aligned_percentage = max(max_aligned_percentage,all_percent)
        min_chrm_percentage = min(min_chrm_percentage,chrm_reads_percentage)
        max_chrm_percentage = max(max_chrm_percentage, chrm_reads_percentage)
        min_contig_percentage = min(min_contig_percentage, contig_reads_percentage)
        max_contig_percentage = max(max_contig_percentage, contig_reads_percentage)

        # now print out the information
        report.writerow([
            sample_name,
            total_reads_num,
            uniq_num,
            nu_num,
            all_num,
            sum(chrM),
            sum(chr_contigs),
            uniq_percent,
            nu_percent,
            all_percent,
            chrm_reads_percentage,
            contig_reads_percentage
        ])

    report.writerow(['Minimums',
        min_reads_num,
        min_aligned_percentage])
    report.writerow(['Maximums',
        max_reads_num,
        max_aligned_percentage])

    report_file.close()

    # if verbose, read in report and write to stdout the prettytable
    if args.verbose:
        with open(args.output) as report:
            reader = csv.reader(report)
            pt = prettytable.PrettyTable(headers)
            for row in reader:
                if len(row) < len(headers):
                    break
                logging.debug(row)
                pt.add_row(row)
            print(pt)

if __name__ == '__main__':
    main()
