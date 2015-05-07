#!/usr/bin/env python

"""REST client for Ensembl"""

import logging
import urllib
import urllib2
import time
import sys
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom


__author__ = "Damon May"
__copyright__ = "Copyright (c) 2015 Damon May"
__license__ = ""
__version__ = ""

logger = logging.getLogger(__name__)

class EnsemblRestClient(object):
    '''copied from:
    https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client
    '''
    def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=20):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def do_rest_action(self, endpoint, hdrs=None, params=None):
        if hdrs is None:
            hdrs = {}

        if params:
            endpoint += '?' + urllib.urlencode(params)

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        content = None
        try:
            request = urllib2.Request(self.server + endpoint, headers=hdrs)
            response = urllib2.urlopen(request)
            content = response.read()
            self.req_count += 1

        except urllib2.HTTPError, e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    content = self.do_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.code} Reason: {1}\n'.format(endpoint, e))

        return content

    def get_transcript_seq(self, transcript_id):
        logger.debug("get_transcript, id=%s" % transcript_id)
        transcript_seq = self.do_rest_action(
            '/sequence/id/{0}'.format(transcript_id),
            params={'content-type': 'text/plain'}
        )
        return transcript_seq

    def get_transcript_cds(self, transcript_id):
        logger.debug("get_transcript, id=%s" % transcript_id)
        transcript_seq = self.do_rest_action(
            '/sequence/id/{0}'.format(transcript_id),
            params={'content-type': 'text/plain',
                    'type': 'cds'}
        )
        return transcript_seq

    def get_transcript_startend(self, transcript_id):
        """

        :param transcript_id:
        :return: start, end
        """
        logger.debug("get_transcript_metadata, id=%s" % transcript_id)
        xml_metadata = self.do_rest_action(
            '/lookup/id/{0}'.format(transcript_id),
            params={'content-type': 'text/xml',
                    'expand': '0'}
        )
        root = ET.fromstring(xml_metadata)

        data_elem = root.find("data")
        start = int(data_elem.get("start"))
        end = int(data_elem.get("end"))
        return start, end

    def get_transcript_cds_startends(self, transcript_id):
        """

        :param transcript_id:
        :return:
        """
        logger.debug("get_transcript_cds_startends, id=%s" % transcript_id)
        xml_metadata = self.do_rest_action(
            '/overlap/id/{0}'.format(transcript_id),
            params={'feature': 'cds',
                    'content-type': 'text/xml'}

        )
        root = ET.fromstring(xml_metadata)

        cds_starts_ends = []
        for data_elem in root.findall("data"):
            if data_elem.get("Parent") != transcript_id:
                continue
            cds_starts_ends.append((int(data_elem.get("start")), int(data_elem.get("end"))))
        return cds_starts_ends

    def get_transcript_startend_strand(self, transcript_id):
        """
        get transcript start and end coordinates and strand
        :param transcript_id:
        :return: transcript_start, transcript_end, cds_start, cds_end
        """
        logger.debug("get_transcript_metadata, id=%s" % transcript_id)
        xml_metadata = self.do_rest_action(
            '/lookup/id/{0}'.format(transcript_id),
            params={'content-type': 'text/xml',
                    'expand': '1'}
        )
        root = ET.fromstring(xml_metadata)

        data_elem = root.find("data")
        start = int(data_elem.get("start"))
        end = int(data_elem.get("end"))
        strand = int(data_elem.get("strand"))
        return start, end, strand

    def translate_position_relative(self, abs_start, abs_end, position, strand):
        """
        Given an absolute start and end position, and a strand, give the position
        in relative terms from the effective start
        :param abs_start:
        :param abs_end:
        :param position:
        :param strand: 1 or -1
        :return:
        """
        assert(strand == 1 or strand == -1)
        assert(abs_start <= position <= abs_end)
        if strand == 1:
            return position - abs_start
        return abs_end - position

    def get_transcript_with_relcdspos(self, transcript_id):
        """
        :param transcript_id:
        :return: transcript sequence, and 0-based relative positions of cds start and end
        """
        logger.debug("get_transcript_with_relcdspos 0")
        transcript_seq = self.get_transcript_seq(transcript_id)
        transcript_startpos, transcript_endpos, strand = self.get_transcript_startend_strand(transcript_id)
        logger.debug("get_transcript_with_relcdspos 1")
        cds_abs_starts_ends = self.get_transcript_cds_startends(transcript_id)
        logger.debug("get_transcript_with_relcdspos 2")
        if not cds_abs_starts_ends:
            raise Exception('no coding region')
        logger.debug("get_transcript_with_relcdspos 3")

        if len(transcript_seq) != abs(transcript_startpos - transcript_endpos) + 1:
            raise Exception('Transcript seq length %d does not match start and end %d, %d' %
                            (len(transcript_seq), transcript_startpos, transcript_endpos))
        logger.debug("get_transcript_with_relcdspos 4")
        cds_starts_ends_rel = []
        for cds_abs_start, cds_abs_end in cds_abs_starts_ends:
            logger.debug("get_transcript_with_relcdspos cds loop")
            if strand == -1:
                rel_start = transcript_endpos - cds_abs_end
                rel_end = transcript_endpos - cds_abs_start
            else:
                rel_start = cds_abs_start - transcript_startpos
                rel_end = cds_abs_end - transcript_startpos
            cds_starts_ends_rel.append((rel_start, rel_end))
        cds_starts_ends_rel.sort(key=lambda x: x[0])

        return transcript_seq, cds_starts_ends_rel

    def get_transcript_cds_5prime_3prime_seqs(self, transcript_id):
        """
        Cover method to get a transcript's CDS, 3' UTR and 5' UTR
        :param transcript_id:
        :return: CDS, 5' UTR, 3' UTR
        """
        transcript_seq, cds_starts_ends_rel = self.get_transcript_with_relcdspos(transcript_id)

        cds_start = cds_starts_ends_rel[0][0]
        cds_end = cds_starts_ends_rel[len(cds_starts_ends_rel) - 1][1]

        utr_5prime = transcript_seq[0:cds_start]
        utr_3prime = transcript_seq[cds_end+1:]
        cds = ''
        for cds_chunk in cds_starts_ends_rel:
            cds = cds + transcript_seq[cds_chunk[0]: cds_chunk[1] + 1]
        return cds, utr_5prime, utr_3prime





