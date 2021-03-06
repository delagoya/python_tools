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

"""Filters a SAM/BAM file, sorted in read ID order, to sequences that only come from the canonical chromosomes. E.g. not the contigs or mitochondial chromosomes. Likely only to work for human and mouse.
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
        'output',
        type=str,
        help="Output SAM/BAM filename."
    )
    args.add_argument(
        '-r','--rejected',
        type=str,
        help="Optional output SAM/BAM filename for rejected reads."
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

# A regex for grabbing canonical chromosomes from human/mouse
_RE_CONANICAL_CHR  = re.compile(r'chr[0-9XY]{1,2}$')

def setup_logging(args):
    '''Sets up normal and verbose logging'''
    logging.basicConfig(level=logging.ERROR,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

def get_valid_chrs(headers):
    """Returns a Dictionary of the valid target sequence entries from the Sam file headers
    """
    valid_chrs = {}
    for i,e in enumerate(headers['SQ']):
        if _RE_CONANICAL_CHR.match(e['SN']):
            valid_chrs[i] = True
    logging.debug("Valid Chromosomes: [{0}]".format(valid_chrs))
    return valid_chrs

def main():
    args = get_arguments()
    setup_logging(args)
    logging.debug(args)

    # handle STDIN/STDOUT streams
    if not args.input:
        args.input = '-'
    logging.debug("I=%s ; O=%s" % (args.input, args.output))

    # When input filename given, handle missing output filename
    if not args.output:
        if args.input == '-':
            args.output = '-'
        else:
            d,f = os.path.split(args.input)
            f = re.sub(r"\.([sb]am)$",r'.filtered.\1',f)
            args.output = os.path.join(d,f)

    # handle read and write modes, taking into account BAM files
    read_mode = 'r'
    write_mode = 'wh'
    if args.input_bam or re.search(r'\.bam$',args.input):
        read_mode = 'rb'
    if args.output_bam or re.search(r'\.bam$',args.output):
        write_mode = 'wb'

    # Input/output/rejected SAM/BAM files
    src =  pysam.Samfile(args.input,read_mode)
    target = pysam.Samfile(args.output,write_mode,template=src)
    rejected = None
    if args.rejected:
        rejected =  pysam.Samfile(args.rejected,write_mode,template=src)

    # some handy vars for bookkeeping and workflow
    valid_chrs = get_valid_chrs(src.header)

    # Total number of reads output
    total_output_tally = 0
    rejected_output_tally = 0
    totally_rejected_tally = 0
    last_entry = None
    try:
        while True:
            entry = src.next()
            if not last_entry:
                last_entry = entry
            valid_written = False
            reject_written = False
            if valid_chrs.has_key(entry.rname):
                target.write(entry)
                valid_written = True
            else:
                entry_rejected = True
                if rejected:
                    rejected.write(entry)
                    reject_written = True

            if last_entry.qname != entry.qname:
                if valid_written:
                    total_output_tally += 1
                elif reject_written:
                    rejected_output_tally += 1
                    totally_rejected_tally += 1
                elif entry_rejected and not valid_written and not reject_written:
                    totally_rejected_tally += 1
            last_entry = entry
    except StopIteration,e:
        pass

    src.close()
    target.close()
    sys.stderr.write("Number of valid entries: %d\n" % total_output_tally)
    sys.stderr.write("Number of totally rejected entries: %d\n" % totally_rejected_tally)
    if args.rejected:
        rejected.close()
        sys.stderr.write("Number of output rejected entries: %d\n" % rejected_output_tally)

if __name__ == '__main__':
    main()
