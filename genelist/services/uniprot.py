#!/usr/bin/env python
# encoding: utf-8

import sys
import requests
import requests_cache
import time
import xmltodict

from ..utils import cleanup_description

class Uniprot(object):
    """Basic interface to the public UniProt API.

    Args:
        response_format (str): format for response (xml, json, etc.)
    """

    def __init__(self, response_format='application/xml'):
        super(Uniprot, self).__init__()
        self.base_url = 'http://www.uniprot.org/uniprot'
        self.format = response_format

        requests_cache.install_cache('uniprot_cache', backend='sqlite', expire_after=8460000)

    def get(self, handler):
        """Compose url and universal headers for any request handler.

        Args:
            handler (str): API entry point

        Returns (dict): a parsed xml object
        """
        url = "%s/%s" % (self.base_url, handler)
        headers = {
            'Accept': self.format
        }

        res = requests.get(url, headers=headers)

        if not res.from_cache:
            time.sleep(0.25) # wait for 250ms, play nice

        return xmltodict.parse(res.text)

    def fetch_description(self, uniprot_id):
        """Fetches the description of a uniprot id.

        Args:
            uniprot_id (str): a UniProt ID

        Returns (str): descriptive text

        """
        data = self.get("%s.xml" % uniprot_id)
        try:
            info = data['uniprot']['entry']['protein']['recommendedName']['fullName']
        except (KeyError, IndexError):
            return None

        return cleanup_description(info)
