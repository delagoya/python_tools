#!/usr/bin/env python
from __future__ import division, print_function
import argparse
import logging
import os
import sys
import re
from math import ceil
import numpy as np
import pysam

"""Filters a SAM/BAM file, sorted in read ID order, to three SAM/BAM files of unique, non-unique, and unmapped reads.
"""

def get_arguments():
    '''Parses the CLI arguments'''
    args = argparse.ArgumentParser(description="Filters a SAM/BAM alignment file (sorted by read IDs) to canonical chromosomal alignments. Optionally outputs rejected reads to another file.")
    args.add_argument(
        "input",
        type=str,
        help="Input SAM/BAM file for filtering."
    )
    args.add_argument(
        'output_prefix',
        type=str,
        help="Output SAM/BAM filename prefix. Two files will be produced. <prefix>.uniq.sam and <prefix>.nuniq.sam"
    )
    args.add_argument(
        '-B', '--output-bam',
        action='store_true',
        help="Output files are BAM format, instead of SAM. E.b. <prefix>.uniq.bam, <prefix>.nuniq.bam, and <prefix>.unmapped.bam"
    )
    args.add_argument(
        '--verbose', '-v',
        action='store_true',
        help="Print verbose information"
    )
    args.add_argument(
        '--debug', '-d',
        action='store_true',
        help="Print debugging information"
    )
    return args.parse_args()

def setup_logging(args):
    '''Sets up normal and verbose logging'''
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    # filename=args.output_prefix + ".log",
    # filemode='w')
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

def get_valid_chrs(header):
    """Returns a Dictionary of the valid target sequence entries from the Sam file headers
    """
    canonical_chr = re.compile(r'chr[0-9XY]{1,2}$')
    valid_chrs = {}
    for i,e in enumerate(header['SQ']):
        if canonical_chr.match(e['SN']):
            valid_chrs[i] = True
    logging.debug("Valid Chromosomes: [{0}]".format(valid_chrs))
    return valid_chrs

def get_next_alignments(samfile,last_entry):
    """Returns a set of mapping entries for a read from a SAM file that map to one of the standard chromosomes.
    """
    entries = []
    try:
        if last_entry == None:
            last_entry = samfile.next()
        entries.append(last_entry)
        while True:
            tmp_entry = samfile.next()
            # go on to the next entry if not on a standard CHR
            if last_entry.qname == tmp_entry.qname:
                entries.append(tmp_entry)
                last_entry = tmp_entry
            else:
                last_entry = tmp_entry
                break
    except StopIteration, e:
        logging.info("Reached end of SAM file " + samfile.filename)
        entries.append(last_entry)
        last_entry = None
    return (entries,last_entry)

def write_entries(entries,samfile,stats=None):
    for entry in entries:
        # determine is mate is valid
        samfile.write(entry)
        if stats:
            stats[entry.rname - 1] += 1

def write_stats(total_count, stats, stats_log, sam):
    stats_log.write("Total Fragments\t{0}\n".format(total_count))
    for i,c in enumerate(stats):
        stats_log.write("{chr}\t{cnt}\n".format(chr=sam.getrname(i),
            cnt=c
        ))

def main():
    args = get_arguments()
    setup_logging(args)
    logging.debug(args)

    # define output filenames
    pf = args.output_prefix
    if args.output_bam:
        o_ext = 'bam'
    else:
        o_ext = 'sam'

    # handle read and write modes, taking into account BAM files
    read_mode = 'r'
    write_mode = 'wh'
    if re.search(r'\.bam$',args.input):
        read_mode = 'rb'
    if args.output_bam:
        write_mode = 'wb'

    # Input/output/rejected SAM/BAM files
    src =  pysam.Samfile(args.input,read_mode)

    uniq = pysam.Samfile(".".join((pf,'uniq',o_ext)), write_mode, template=src)
    uniq_stats_log = open(args.output_prefix + ".uniq.mapping_stats.txt",'w')

    nuniq = pysam.Samfile(".".join((pf,'nuniq',o_ext)), write_mode, template=src)
    nuniq_stats_log = open(args.output_prefix + ".nuniq.mapping_stats.txt",'w')

    unmapped = pysam.Samfile(".".join((pf,'unmappes',o_ext)),write_mode, template=src)

    # Total number of reads output
    total_uniq = total_unmapped = total_nu = 0
    uniq_stats = [0] * len(src.header['SQ'])
    nuniq_stats = [0] * len(src.header['SQ'])

    # Classic Python iteration using PySAM's fetch() function is
    # **much** slower than using the loop below.
    last_entry = None
    count = 0
    while True:
        entries, last_entry = get_next_alignments(src,last_entry)
        if not last_entry:
            break
        e = entries[0]
        if e.is_unmapped and e.mate_is_unmapped:
            total_unmapped += 1
            # may be prudent to write these out to a SAM file as well
            write_entries(entries,unmapped)
        else:
            map_count = [x for x in e.tags if re.match(r'IH|NH',x[0])][0][1]

            if map_count > 1:
                write_entries(entries,nuniq,nuniq_stats)
                total_nu += 1
            else:
                write_entries(entries,uniq,uniq_stats)
                total_uniq += 1
        count += 1

    # write out the stats
    write_stats(total_uniq, uniq_stats,uniq_stats_log, uniq)
    write_stats(total_nu, nuniq_stats, nuniq_stats_log, nuniq)

    src.close()
    uniq.close()
    nuniq.close()
    unmapped.close()

    uniq_stats_log.close()
    nuniq_stats_log.close()

    logging.info("Number of unique entries: %d\n" % total_uniq)
    logging.info("Number of non-unique entries: %d\n" % total_nu)
    logging.info("Number of unmapped entries: %d\n" % total_unmapped)

if __name__ == '__main__':
    main()
