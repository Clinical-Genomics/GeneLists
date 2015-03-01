# -*- coding: utf-8 -*-
from __future__ import absolute_import, unicode_literals
from datetime import datetime

import time
import requests

def format_entry(json_entry):
  """Extract interesting information from a single OMIM entry."""
  # extract nested titles section
  titles = json_entry.get('titles', {})

  if 'includedTitles' in titles:
    other_entities = titles['includedTitles'].split(';;\n')
  else:
    other_entities = []

  # extract "geneMap"
  gene_map = json_entry.get('geneMap', {})

  # extract phenotypes for the gene
  phenotypes_raw = (item['phenotypeMap']
                    for item in gene_map.get('phenotypeMapList', []))
  phenotypes = [{
    'phenotype_mim_number': phenotype.get('phenotypeMimNumber'),
    'mim_number': phenotype.get('mimNumber'),
    'phenotype': phenotype.get('phenotype'),
    'mapping_key': phenotype.get('phenotypeMappingKey'),
    'inheritance': phenotype.get('phenotypeInheritance')
  } for phenotype in phenotypes_raw]

  data = {
    'prefix': json_entry.get('prefix'),
    'mim_number': json_entry.get('mimNumber'),
    'status': json_entry.get('status'),
    'other_entities': other_entities,
    'gene_symbol': gene_map.get('geneSymbols'),
    'gene_name': gene_map.get('geneName'),
    'gene_location': gene_map.get('computedCytoLocation'),
    'phenotypes': phenotypes,
  }

  if 'epochCreated' in json_entry:
    data['created_at'] = datetime.fromtimestamp(json_entry.get('epochCreated'))
    data['updated_at'] = datetime.fromtimestamp(json_entry.get('epochUpdated'))

  return data


class OMIM(object):

  """Basic interface to the public OMIM API.

  Can be lazy loaded like other Flask extensions by calling ``init_app``
  post initialization.

  Args:
    app (Flask): Flask instance
    response_format (str): format for response (xml, json, etc.)
  """

  def __init__(self, api_key, response_format='json'):
    super(OMIM, self).__init__()
    self.base_url = 'http://api.europe.omim.org/api'
    self.format = response_format
    self.api_key = api_key

  def base(self, handler):
    """Compose url and universal params for any request handler.

    Args:
      handler (str): API entry point

    Returns:
      str, dict: URL entry point and params as a ``dict``
    """
    url = "%s/%s" % (self.base_url, handler)
    params = {
      'apiKey': self.api_key,
      'format': self.format
    }

    return url, params

  def parse_phenotypic_disease_model(self, phenotypes):
    """Compose Phenotypic_disease_model entry based on the inheritance models of a gene.
    Such entry might look like: 614096>AR|615889>AR

    Args:
        phenotypes (list of dicts): each dict represents a phenotype. This is the 'phenotypes' key in the data returned from OMIM

    Returns: TODO

    """
    TERMS_MAPPER = {
      'Autosomal recessive': 'AR',
      'Autosomal dominant': 'AD',
      'X-linked dominant': 'XD',
      'X-linked recessive': 'XR',
    }

    TERMS_BLACKLIST = [
      'Isolated cases',
    ]

    phenotypic_disease_model = []
    models = set()
    for phenotype in phenotypes:
      if phenotype['inheritance'] is not None:
        models.update([model.strip('? ') for model in phenotype['inheritance'].split(';')])
        models = models.difference(TERMS_BLACKLIST) # remove blacklisted terms
        models = set([TERMS_MAPPER.get(model_human, model_human) for model_human in models]) # rename them if possible

        phenotypic_disease_model.append('%s>%s' % (phenotype['phenotype_mim_number'], '/'.join(models)))
      else:
        if (phenotype['phenotype_mim_number'] is not None):
          phenotypic_disease_model.append(str(phenotype['phenotype_mim_number']))

    if len(phenotypic_disease_model) == 0:
      return None
    return '|'.join(phenotypic_disease_model)

  def gene(self, hgnc_symbol):
    entry = self.search_gene(hgnc_symbol)

    return format_entry(entry)

  def search_gene(self, hgnc_symbol, include=('geneMap', 'dates')):
    """Search for MIM number for a HGNC approved symbol.

    Args:
      hgnc_symbol (str): HGNC approved symbol
      include (list, optional): additional sections to include

    Returns:
      dict: first matching entry in the results
    """
    url, params = self.base('entry/search')

    params['search'] = "approved_gene_symbol:%s" % hgnc_symbol
    params['include'] = include

    res = False
    sleep = 0
    retry = True # Execute the first
    while retry or res.status_code == requests.codes.conflict:
        time.sleep(sleep)

        try:
            retry = False
            res = requests.get(url, params=params)
        except Exception:
            retry = True
        #except TypeError:
        #    retry = True
        #except ProtocolError:
        #    retry = True

        # when we get trottled, give it a sec
        if sleep > 1000: # if sleeping for 15mins, reset
            sleep = 1
        else:
            sleep *= 2

    data = res.json()

    entries = data['omim']['searchResponse']['entryList']

    if entries:
      return entries[0]['entry']
    else:
      return {}

  def clinical_synopsis(self, mim, include=('clinicalSynopsis',),
                        exclude=None):
    """Get data from ``clinicalSynopsis`` handler.

    Args:
      mim (str): OMIM Phenotype MIM number
      include (list, optional): included sections in response
      exclude (list, optional): exclude sections in response

    Returns:
      response: response object from ``requests``
    """
    url, params = self.base('clinicalSynopsis')
    params['mimNumber'] = mim
    params['include'] = include

    if exclude:
      params['exclude'] = exclude

    return requests.get(url, params=params)

  def entry(self, mim, includes=None, include_all=False):
    """Get data from ``entry`` handler.

    Args:
      mim (str or list): OMIM Phenotype MIM number(s)
      includes (str or list, optional): included sections in response
      include_all (bool, optional): same as ``includes='all'``

    Returns:
      response: response object from ``requests``
    """
    url, params = self.base('entry')

    # add MIM number(s)
    params['mimNumber'] = mim

    # add include section(s)
    if include_all:
      params['include'] = 'all'
    else:
      params['include'] = includes

    # send request
    return requests.get(url, params=params)

if __name__ == '__main__':
    omim = OMIM(api_key='<fill in key>')
    print(omim.gene('AP1S1'))
