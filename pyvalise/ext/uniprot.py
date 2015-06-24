#!/usr/bin/env python

"""Client for getting information from UniProt"""

import logging
import urllib
import urllib2
import time
import sys
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET


__author__ = "Damon May"
__copyright__ = "Copyright (c) 2015 Damon May"
__license__ = ""
__version__ = ""

logger = logging.getLogger(__name__)

# types of secondary structure features
SECONDARY_STRUCTURE_TYPES = ["helix", "turn", "strand"]

UNIPROT_SERVER = "http://uniprot.org/uniprot"

class UniProtClient(object):

    def __init__(self, server=UNIPROT_SERVER, reqs_per_sec=20):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

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

