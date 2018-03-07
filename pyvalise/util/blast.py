#!/usr/bin/env python
"""
parse BLAST results
"""

import logging
import csv

__author__ = "Damon May"
__copyright__ = "Copyright (c) 2015 Damon May"
__license__ = ""
__version__ = ""

logger = logging.getLogger(__name__)


def load_blast_protein_proteins_evalues_map(blastfile, max_e=1000, baitproteins_tokeep=None,
                                            trim_accessions=True):
    """

    :param blastfile:
    :param max_e:
    :param baitproteins_tokeep:  a list of bait protein IDs. Only keep results that are on this list
    :param trim_accessions: if accession numbers contain two ||'s, trim down to the value in the middle
    :return:
    """
    result = {}
    for line in blastfile:
        # this is embarrassing... fault tolerance for screwed-up files
        try:
            blast_hit = parse_blast_line(line)
        except ValueError:
            logger.debug("Bad blast line:  %s" % line)
            continue
        if blast_hit.evalue > max_e:
            continue
        if baitproteins_tokeep and blast_hit.query_protein not in baitproteins_tokeep:
            continue
        if blast_hit.query_protein not in result:
            result[blast_hit.query_protein] = []
        if trim_accessions and blast_hit.hit_protein.count('|') == 2:
            blast_hit.hit_protein = blast_hit.hit_protein.split('|')[1]
        result[blast_hit.query_protein].append((blast_hit.hit_protein, blast_hit.evalue))
    return result


def load_blast_protein_proteins_map(blastfile, max_e=1000, baitproteins_tokeep=None):
    """

    :param blastfile:
    :param max_e:
    :param baitproteins_tokeep:  a list of bait protein IDs. Only keep results that are on this list
    :return:
    """
    result = {}
    for line in blastfile:
        # this is embarrassing... fault tolerance for screwed-up files
        try:
            blast_hit = parse_blast_line(line)
        except ValueError:
            logger.debug("Bad blast line:  %s" % line)
            continue
        if blast_hit.evalue > max_e:
            continue
        if baitproteins_tokeep and blast_hit.query_protein not in baitproteins_tokeep:
            continue
        if blast_hit.query_protein not in result:
            result[blast_hit.query_protein] = []
        result[blast_hit.query_protein].append(blast_hit.hit_protein)
    return result


def iterparse_blast_lines(blastfile):
    """
    Yield bait, hit, evalue
    :param file:
    :return:
    """
    for line in blastfile:
        try:
            blast_hit = parse_blast_line(line)
            yield blast_hit
        except ValueError as err:
            logger.debug("Bad blast line %s" % line)

def parse_blast_line(line):
    """
    Returns bait, hit, evalue, taxon_id
    :param line:
    :return:
    """
    chunks = line.split('\t')
    if len(chunks) < 11:
        raise ValueError("BLAST result line too short (%d fields)\n%s" % (len(chunks), line))
    #bait is 0. match is 1. e-value is 10
    # chunks are:
    # qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
    try:
        qseqid = chunks[0]
        sseqid = chunks[1]
        pident = float(chunks[2])
        length = int(chunks[3])
        mismatch = int(chunks[4])
        gapopen = int(chunks[5])
        qstart = int(chunks[6])
        qend = int(chunks[7])
        sstart = int(chunks[8])
        send = int(chunks[9])
        evalue = float(chunks[10])
        bitscore = float(chunks[11])
    except ValueError:
        raise ValueError("BLAST result %s won't parse." % line)
    hit_taxon_id = None
    if len(chunks) > 12:
        if chunks[12]:
            hit_taxon_id = int(chunks[12])
    return BlastHit(qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send,
                    evalue, bitscore, hit_taxon_id)


class BlastHit(object):
    def __init__(self, query_protein, hit_protein, percent_id, length, mismatch, gapopen, qstart, qend, sstart, send,
                 evalue, bitscore, hit_taxon_id = None):
        self.query_protein = query_protein
        self.hit_protein = hit_protein
        self.percent_id = percent_id
        self.length = length
        self.mismatch = mismatch
        self.gapopen = gapopen
        self.qstart = qstart
        self.qend = qend
        self.sstart = sstart
        self.send = send
        self.evalue = evalue
        self.bitscore = bitscore
        self.hit_taxon_id = hit_taxon_id

