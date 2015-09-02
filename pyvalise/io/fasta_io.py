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


def write_dna_name_seq_to_fasta(name, seq, fasta_file, description=None):
    """
    create a single record for a name and a sequence; write it to a fasta file.
    Should be able to call this multiple times on a single file
    :param name:
    :param seq:
    :param fasta_file:
    :return:
    """
    write_record_to_fasta(make_dna_seq_record(name, description, seq), fasta_file)


def write_protein_name_seq_to_fasta(name, seq, fasta_file, description=""):
    """
    create a single record for a name and a sequence; write it to a fasta file.
    Should be able to call this multiple times on a single file
    :param name:
    :param seq:
    :param fasta_file:
    :return:
    """
    write_record_to_fasta(make_protein_seq_record(name, description, seq), fasta_file)


def write_record_to_fasta(record, fasta_file):
    """write a single record to a fasta file. Should be able to call this
    multiple times on a single file
    """
    SeqIO.write([record], fasta_file, "fasta")


def build_record_fasta_text(record):
    """
    This is only to be used if other options are difficult, for some reason.
    Doesn't split seq across lines
    :param record:
    :return:
    """
    return ">%s %s\n%s" % (record.name, record.description, record.sequence)


def make_seq_record(id, description, seq, alphabet):
    return SeqRecord(Seq(seq, alphabet), id=id, description=description)


def make_dna_seq_record(id, description, seq):
    return make_seq_record(id, description, seq, IUPAC.ambiguous_dna)


def make_protein_seq_record(id, description, seq):
    return make_seq_record(id, description, seq, IUPAC.protein)

