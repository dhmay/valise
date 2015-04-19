#!/usr/bin/env python

"""REST client for Ensembl"""

import logging
import urllib
import urllib2
import time
import sys


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

    def perform_rest_action(self, endpoint, hdrs=None, params=None):
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
                    content = self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.code} Reason: {1}\n'.format(endpoint, e))

        return content

    def get_masked_transcript(self, transcript_id):
        logger.debug("get_masked_transcript, id=%s" % transcript_id)
        masked_transcript = self.perform_rest_action(
            '/sequence/id/{0}'.format(transcript_id),
            params={'content-type': 'text/plain',
                    'mask_feature': '1'}
        )
        return masked_transcript

    def get_transcript_5_cdswithintrons_3(self, transcript_id):
        '''
        :param transcript_id:
        :return: a tuple with (5' UTR, cds with introns, 3' UTR), all uppercase
        '''
        masked_transcript = self.get_masked_transcript(transcript_id)
        cds_start_idx = -1
        for i in xrange(0, len(masked_transcript)):
            if masked_transcript[i].islower():
                cds_start_idx = i
                break
        assert(cds_start_idx >= 0)
        cds_last_idx = -1
        for i in xrange(len(masked_transcript)-1, 0, -1):
            if masked_transcript[i].islower():
                cds_last_idx = i
                break
        assert(cds_last_idx > cds_start_idx)
        utr_5prime = masked_transcript[0:cds_start_idx].upper()
        cds_with_introns = masked_transcript[cds_start_idx:cds_last_idx + 1].upper()
        utr_3prime = ''
        if cds_last_idx < len(masked_transcript):
            utr_3prime = masked_transcript[cds_last_idx + 1:].upper()
        return utr_5prime, cds_with_introns, utr_3prime



