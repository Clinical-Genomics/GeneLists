#!/usr/bin/env python
# encoding: utf-8

import sys
import requests
import requests_cache
import time

class Genenames(object):
    """Basic interface to the public genenames API.

    Args:
        response_format (str): format for response (xml, json, etc.)
    """

    def __init__(self, response_format='application/json'):
        super(Genenames, self).__init__()
        self.base_url = 'http://rest.genenames.org/'
        self.format = response_format

        requests_cache.install_cache('genenames_cache', backend='sqlite', expire_after=8460000)

    def get(self, handler):
        """Compose url and universal headers for any request handler.

        Args:
            handler (str): API entry point

        Returns (json): a parsed json object
        """
        url = "%s/%s" % (self.base_url, handler)
        headers = {
            'Accept': self.format
        }

        res = requests.get(url, headers=headers)

        if not res.from_cache:
            time.sleep(0.25) # wait for 250ms, play nice

        return res.json()

    def fetch_symbol(self, hgnc_symbol, key):
        """Fetches the HGNC symbol and information related to it.

        Args:
            hgnc_symbol (str): an HGNC symbol.
            key (str): e.g. alias_symbol, uniprotids, refseqids

        Returns (dict): result set.

        """
        data = self.get("fetch/symbol/%s" % hgnc_symbol)
        try:
            info = data['response']['docs'][0][key]
        except (KeyError, IndexError):
            return None

        return info

    def omim(self, hgnc_symbol):
        """Fetches the OMIM morbid id.

        Args:
            hgnc_symbol (str): an HGNC symbol.

        Returns (list): OMIM morbid.

        """
        return self.fetch_symbol(hgnc_symbol, 'omim_id')

    def aliases(self, hgnc_symbol):
        """Fetches the HGNC aliases

        Args:
            hgnc_symbol (str): an HGNC symbol.

        Returns (list): a list of HGNC symbols

        """
        return self.fetch_symbol(hgnc_symbol, 'alias_symbol')

    def uniprot(self, hgnc_symbol):
        """Fetches the HGNC aliases

        Args:
            hgnc_symbol (str): an HGNC symbol.

        Returns (list): a list of HGNC symbols

        """
        return self.fetch_symbol(hgnc_symbol, 'uniprot_ids')

    def refseq(self, hgnc_symbol):
        """Fetches the HGNC aliases

        Args:
            hgnc_symbol (str): an HGNC symbol.

        Returns (list): a list of HGNC symbols

        """
        return self.fetch_symbol(hgnc_symbol, 'refseq_accession')

    def _parse_official(self, data, omim_morbid=None):
        """Parses the official HGNC symbol out of a JSON result set
        fetched from genenames.org.
        On missing OMIM morbid number, the first HGNC symbol is returned

        Args:
            data (json): result set from querying rest.genenames.org/fetch/symbol/%s
            omim_morbid (str, opt): an option omim morbid number

        Returns (str): the official HGNC identifier

        """
        if omim_morbid == None:
            return data['response']['docs'][0]['symbol']

        for symbols in data['response']['docs']:
            if str(symbols['omim_id'][0]) == str(omim_morbid):
                return symbols['symbol']

        return None

    def official(self, hgnc_symbol, omim_morbid=None):
        """Fetches the HGNC official symbol for this HGNC symbol. On multiple matches,
        it is better to also provide the OMIM morbid number.

        Args:
            hgnc_symbol (str): an HGNC symbol.
            omim_morbid (str, opt): the associated omim morbid number

        Returns (str): the official HGNC symbol

        """

        data = self.get("fetch/symbol/%s" % hgnc_symbol)

        official_symbol = None
        try:
            official_symbol = self._parse_official(data, omim_morbid)
        except (KeyError, IndexError) as e:
            pass

        if official_symbol == None:
            # ok, no results found, maybe try its previous symbol?
            data = self.get("fetch/prev_symbol/%s" % hgnc_symbol)

            try:
                official_symbol = self._parse_official(data, omim_morbid)
            except (KeyError, IndexError) as e:
                return None

        return official_symbol

if __name__ == '__main__':
    main(sys.argv[1:])
