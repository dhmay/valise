#!/usr/bin/env python
"""
merge two fasta files, keeping unique sequences. Keeps the first sequence from the first file.
Does not handle ambiguity characters in any interesting way.
This is boneheaded and bloated, doesn't try to save memory at all.
"""

import argparse
import logging
from datetime import datetime
from pyvalise.io import fasta_io

__author__ = "Damon May"
__copyright__ = "Copyright (c) 2015 Damon May"
__license__ = ""
__version__ = ""

logger = logging.getLogger(__name__)


def declare_gather_args():
    """
    Declare all arguments, parse them, and return the args dict.
    Does no validation beyond the implicit validation done by argparse.
    return: a dict mapping arg names to values
    """

    # declare args
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fasta1', type=argparse.FileType('r'),
                        help='input fasta file 1')
    parser.add_argument('fasta2', type=argparse.FileType('r'),
                        help='input fasta file 2')
    parser.add_argument('--out', required=True, type=argparse.FileType('w'),
                        help='output file')
    parser.add_argument('--debug', action="store_true", help='Enable debug logging')
    return parser.parse_args()


def main():
    args = declare_gather_args()
    # logging
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s: %(message)s")
    if args.debug:
        logger.setLevel(logging.DEBUG)
        # any module-specific debugging goes below

    script_start_time = datetime.now()
    logger.debug("Start time: %s" % script_start_time)

    seq_name_map = {}

    n_withinfirst_dups = 0
    n_seqs_processed = 0
    n_written = 0
    for record in fasta_io.load_fasta_records(args.fasta1):
        n_seqs_processed += 1
        seq = str(record.seq)
        if seq in seq_name_map:
            n_withinfirst_dups += 1
        else:
            seq_name_map[seq] = record.name
            fasta_io.write_record_to_fasta(record, args.out)
            n_written += 1
    print("File 1 processed. %d duplicates from %d total sequences" % (n_withinfirst_dups, n_seqs_processed))
    n_written_firstfile = n_written
    n_secondfile_dups = 0
    for record in fasta_io.load_fasta_records(args.fasta2):
        n_seqs_processed += 1
        seq = str(record.seq)
        if seq in seq_name_map:
            n_secondfile_dups += 1
        else:
            seq_name_map[seq] = record.name
            fasta_io.write_record_to_fasta(record, args.out)
            n_written += 1
    print("Done. %d duplicates in second file, %d new seqs written" % (n_secondfile_dups, n_written - n_written_firstfile))
    print("Total sequences processed: %d" % n_seqs_processed)
    print("Total sequences written: %d" % n_written)

    print("Done.")




    logger.debug("End time: %s. Elapsed time: %s" % (datetime.now(), datetime.now() - script_start_time))


main()
