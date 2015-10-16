#!/usr/bin/env python

"""Client for getting information from UniProt"""

import logging
import urllib2
import time
import sys
import requests
import StringIO
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

# requests gets really annoying otherwise
logging.getLogger("requests").setLevel(logging.WARNING)


__author__ = "Damon May"
__copyright__ = "Copyright (c) 2015 Damon May"
__license__ = ""
__version__ = ""

logger = logging.getLogger(__name__)

# types of secondary structure features
SECONDARY_STRUCTURE_TYPES = ["helix", "turn", "strand"]

UNIPROT_SERVER = "http://uniprot.org/uniprot"

# batch size for queries. I know for sure 1800 is too big!
DEFAULT_BATCH_SIZE = 200


def make_node_name(raw_node_name):
    return "{%s}%s" % (UNIPROT_SERVER, raw_node_name)


class UniProtClient(object):

    def __init__(self, server=UNIPROT_SERVER, reqs_per_sec=20):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def fetch_entries(self, seqids, batch_size=DEFAULT_BATCH_SIZE):
        """
        fetch UniprotEntry objects for a bunch of sequences, divided up into batches
        """
        seqid_batches = []
        for i in xrange(0, len(seqids), batch_size):
            seqid_batches.append(seqids[i:i+batch_size])
        logger.debug("Splitting %d entries into %d batches of size <= %d" %
                     (len(seqids), len(seqid_batches), batch_size))
        entries = []
        for seqid_batch in seqid_batches:
            entries.extend(self.fetch_entries_onebatch(seqid_batch))
        return entries


    def fetch_entries_onebatch(self, seqids):
        """
        fetch UniprotEntry objects for a bunch of sequences, in one batch
        """
        primary_seqids = [s[:6] for s in seqids]
        logger.debug("Fetching metadata for %d Uniprot IDs from http://uniprot.org ...\n" % len(primary_seqids))
        r = requests.post(
            'http://www.uniprot.org/batch/',
            files={'file': StringIO.StringIO(' '.join(primary_seqids))},
            params={'format': 'xml',
                    'columns': 'id,reviewed',
                    'compress': 'no'
                    })
        while 'Retry-After' in r.headers:
            t = int(r.headers['Retry-After'])
            logger.debug('Waiting %d\n' % t)
            time.sleep(t)
            r = requests.get(r.url)
        cache_txt = r.text

        root = ET.fromstring(cache_txt)
        entries = []
        for entry_xml in root.findall("{" + UNIPROT_SERVER + "}entry"):
            entries.append(UniprotEntry(entry_xml))
        return entries

    def fetch_seqid_entry_map(self, seqids):
        entries = self.fetch_entries(seqids)
        result = {}
        for entry in entries:
            result[entry.accession] = entry
        return result

    def fetch_result(self, uniprot_id):
        hdrs = {}

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        content = None
        try:
            url = self.server + "/" + uniprot_id + ".xml"
            logger.debug("Hitting URL: %s" % url)
            request = urllib2.Request(url, headers=hdrs)
            response = urllib2.urlopen(request)
            content = response.read()
            self.req_count += 1

        except urllib2.HTTPError, e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    content = self.fetch_result(uniprot_id)
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.code} Reason: {1}\n'.format(uniprot_id, e))

        return content

    def get_uniprot_sec_structure_and_length(self, uniprot_id):
        """
        Pull secondary structure elements and length from the UniProt entry. Also return length,
        to prevent having to do a second roundtrip to calculate percent of sequence in secondary structure
        elements.
        :param uniprot_id:
        :return: start, end
        """

        feature_type_stretches_map = {}
        for feature_type in SECONDARY_STRUCTURE_TYPES:
            feature_type_stretches_map[feature_type] = []

        logger.debug("get_uniprot_sec_structure, id=%s" % uniprot_id)
        xml_metadata = self.fetch_result(uniprot_id)
        if not xml_metadata:
            return None, None
        root = ET.fromstring(xml_metadata)

        has_pdb = False
        entry = root.find("{" + UNIPROT_SERVER + "}entry")
        seq_length = int(entry.find("{" + UNIPROT_SERVER + "}sequence").get("length"))
        for dbref_elem in entry.findall("{" + UNIPROT_SERVER + "}dbReference"):
            if dbref_elem.get("type") == "PDB":
                has_pdb = True
                break
        if not has_pdb:
            return None, seq_length
        for feature in entry.findall("{" + UNIPROT_SERVER + "}feature"):
            feature_type = feature.get("type")
            if feature_type in SECONDARY_STRUCTURE_TYPES:
                location = feature.find("{" + UNIPROT_SERVER + "}location")
                begin_1based = int(location.find("{" + UNIPROT_SERVER + "}begin").get("position"))
                end_1based = int(location.find("{" + UNIPROT_SERVER + "}end").get("position"))
                feature_type_stretches_map[feature_type].append((begin_1based, end_1based))

        return feature_type_stretches_map, seq_length


class UniprotEntry:
    def __init__(self, xml_entry):
        self.accession = xml_entry.find(make_node_name('accession')).text
        self.name = xml_entry.find(make_node_name('name')).text
        self.ncbi_taxonomy_id = None
        xml_organism = xml_entry.find(make_node_name('organism'))
        if xml_organism:
            for dbref_elem in xml_organism.findall(make_node_name("dbReference")):
                if dbref_elem.get("type") == "NCBI Taxonomy":
                    self.ncbi_taxonomy_id = int(dbref_elem.get('id'))
                    break


