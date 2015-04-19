#!/usr/bin/env python

"""
FASTA reading/writing
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import logging


log = logging.getLogger(__name__)


def load_fasta_records(fasta_file):
    """load a list of all records in a fasta"""
    result = list()

    return SeqIO.parse(fasta_file, "fasta")


def load_fasta_name_record_map(fasta_file):
    """load the records from the fasta and assemble a dict mapping names to
    records"""
    result = dict()
    for record in load_fasta_records(fasta_file):
        result[record.name] = record
    return result


def load_fasta_names(fasta_file):
    """Load a list of the names of the records in a fasta file"""
    return [record.name for record in load_fasta_records(fasta_file)]

def load_fasta_sequence_record_map(fasta_file):
    """load the records from the fasta and assemble a dict mapping sequences to
    records"""
    result = dict()
    for record in load_fasta_records(fasta_file):
        result[str(record.seq)] = record
    return result


def load_fasta_name_sequence_map(fasta_file):
    """load the records from the fasta and assemble a dict mapping names to
    record sequences"""
    result = dict()
    for record in load_fasta_records(fasta_file):
        result[record.name] = str(record.seq)
    return result


def write_records_to_fasta(records, fasta_file):
    """write a list of entries to a fasta file"""
    SeqIO.write(records, fasta_file, "fasta")


def build_record_fasta_text(record):
    return ">%s %s\n%s" % (record.name, record.description, record.sequence)


def make_dna_seq_record(id, description, seq):
    return SeqRecord(Seq(seq, IUPAC.ambiguous_dna), id=id, description=description)

