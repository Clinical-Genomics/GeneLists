# -*- coding: utf-8 -*-
from __future__ import absolute_import, unicode_literals
from datetime import datetime

#import json
import time
import requests
import requests_cache

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
    'mim_number': json_entry.get('mimNumber', False),
    'status': json_entry.get('status'),
    'other_entities': other_entities,
    'gene_symbol': gene_map.get('geneSymbols'),
    'gene_name': gene_map.get('geneName'),
    'gene_location': gene_map.get('computedCytoLocation', False),
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
    api_key (str): the OMIM API key http://www.omim.org/api
    response_format (str): format for response (xml, json, etc.)
  """

  def __init__(self, api_key, response_format='json'):
    super(OMIM, self).__init__()
    self.base_url = 'http://api.omim.org/api'
    self.format = response_format
    self.api_key = api_key

    requests_cache.install_cache('omim_cache', backend='sqlite', expire_after=8460000)

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

  def parse_phenotypic_disease_models_ext(self, phenotypes, chromosome='', phenotype_number=False):
    """Extract and clean up the inheritance models returned from an OMIM query.

    Args:
        phenotypes (list of dicts): each dict represents a phenotype. This is the 'phenotypes' key in the data returned from OMIM

    Returns
        dict of sorted list: { phenotype_mim_number: [ inheritance_model, inheritance_model, ...] }

    """
    TERMS_MAPPER = {
      'Autosomal recessive': 'AR',
      'Autosomal dominant': 'AD',
      'X-linked dominant': 'XD',
      'X-linked recessive': 'XR',
    }

    TERMS_X = [
      'X-linked dominant',
      'X-linked recessive'
    ]

    TERMS_AUTOSOMAL = [
      'Autosomal recessive',
      'Autosomal dominant'
    ]

    TERMS_BLACKLIST = [
      'Isolated cases',
      'Mitochondrial'
    ]

    phenotypic_disease_model = {}
    for phenotype in phenotypes:

      # only use certain phenotype number
      if phenotype_number and str(phenotype_number) != str(phenotype['phenotype_mim_number']):
          continue

      models = set()
      if phenotype['phenotype_mim_number'] not in phenotypic_disease_model:
        phenotypic_disease_model[ phenotype['phenotype_mim_number'] ] = []

      if phenotype['inheritance'] is not None:
        for model in phenotype['inheritance'].split(';'):
          model = model.strip(' ')

          # add a model
          models.update([ model ])

        models = models.difference(TERMS_BLACKLIST) # remove blacklisted terms
        # remove models that don't belong on this chromosome
        models = models.difference(TERMS_AUTOSOMAL) if chromosome.upper() == 'X' else models.difference(TERMS_X)
        models = set([TERMS_MAPPER.get(model_human, model_human) for model_human in models]) # rename them if possible

      phenotypic_disease_model[ phenotype['phenotype_mim_number'] ].append({
        'models': sorted(list(models)) if len(models) else [],
        'description': phenotype['phenotype']
      })

    return phenotypic_disease_model

  def parse_phenotypic_disease_models(self, phenotypes, chromosome=''):
    """Extract and clean up the inheritance models returned from an OMIM query.

    Args:
        phenotypes (list of dicts): each dict represents a phenotype. This is the 'phenotypes' key in the data returned from OMIM

    Returns
        dict of sorted list: { phenotype_mim_number: [ inheritance_model, inheritance_model, ...] }

    """
    TERMS_MAPPER = {
      'Autosomal recessive': 'AR',
      'Autosomal dominant': 'AD',
      'X-linked dominant': 'XD',
      'X-linked recessive': 'XR',
    }

    TERMS_X = [
      'X-linked dominant',
      'X-linked recessive'
    ]

    TERMS_AUTOSOMAL = [
      'Autosomal recessive',
      'Autosomal dominant'
    ]

    TERMS_BLACKLIST = [
      'Isolated cases',
      'Mitochondrial'
    ]

    phenotypic_disease_model = {}
    for phenotype in phenotypes:
      models = set()
      if phenotype['inheritance'] is not None:
        for model in phenotype['inheritance'].split(';'):
          model = model.strip(' ')

          # exclude unconfirmed phenotypes
          if model.startswith('?'): continue

          # exclude unconfirmed inheritance models
          if phenotype['phenotype'].startswith('?'): continue

          # exclude susceptibility to phenotypes
          if phenotype['phenotype'].startswith('{'): continue

          # add a model
          models.update([ model ])

        models = models.difference(TERMS_BLACKLIST) # remove blacklisted terms
        # remove models that don't belong on this chromosome
        models = models.difference(TERMS_AUTOSOMAL) if chromosome.upper() == 'X' else models.difference(TERMS_X)
        models = set([TERMS_MAPPER.get(model_human, model_human) for model_human in models]) # rename them if possible

      phenotypic_disease_model[ phenotype['phenotype_mim_number'] ] = sorted(list(models)) if len(models) else None

    return phenotypic_disease_model

  def parse_phenotypic_descriptions(self, phenotypes, chromosome):
      descriptions = {}
      for phenotype in phenotypes:
          descriptions[phenotype['phenotype_mim_number']] = phenotype['phenotype']

      return descriptions

  def gene(self, hgnc_symbol=None, mim_number=None):
    entries = self.search_gene(hgnc_symbol=hgnc_symbol, mim_number=mim_number)

    # don't check further if we don't have anything
    if not entries:
      return format_entry({})

    for entry in entries:
      if 'geneMap' in entry['entry']:
        if 'phenotypeMapList' in entry['entry']['geneMap']:
          return format_entry(entry['entry'])

    # no phenotypes found, return something
    return format_entry(entries[0]['entry'])

  def search_gene(self, hgnc_symbol=None, mim_number=None, include=('geneMap', 'dates')):
    """Search for MIM number for a HGNC approved symbol.
    hgnc_symbol or mim_number need to be provided. hgnc_symbol takes precedence.

    Args:
      hgnc_symbol (str, opt): HGNC approved symbol.
      mim_number (str, opt): the omim morbid number.
      include (list, optional): additional sections to include

    Returns:
      list: of dicts. Each dict contains the response for one entry.
    """
    url, params = self.base('entry/search')

    #params['search'] = "%s" % hgnc_symbol # leaving out approved_gene_symbol to get a match on aliases
    if mim_number:
        params['search'] = "number:%s" % mim_number
    else:
        params['search'] = "approved_gene_symbol:%s" % hgnc_symbol
    params['include'] = include

    res = False
    sleep = 0
    retry = True # Execute the first
    while retry or res.status_code == requests.codes.conflict:
        try:
          retry = False
          res = requests.get(url, params=params)
          if not res.from_cache:
              time.sleep(0.25) # wait for 250ms as according to OMIM specs
        except TypeError:
            retry = True
        except ProtocolError:
            retry = True

        # when we get trottled, give it a sec
        if sleep > 1000: # if sleeping for 15mins, reset
            sleep = 1
        else:
            sleep *= 2

    data = res.json()

    entries = data['omim']['searchResponse']['entryList']
#    import json
#    print(json.dumps(entries, sort_keys=True, indent=4, separators=(',', ': ')))

    if entries:
      return entries
    else:
      return []

  def search_gene_raw(self, hgnc_symbol, include=('geneMap', 'dates')):
    """Search for MIM number for a HGNC approved symbol.

    Args:
      hgnc_symbol (str): HGNC approved symbol
      include (list, optional): additional sections to include

    Returns: raw json object.
    """
    url, params = self.base('entry/search')

    params['search'] = "approved_gene_symbol:%s" % hgnc_symbol
    params['include'] = include

    res = False
    sleep = 0
    retry = True # Execute the first
    while retry or res.status_code == requests.codes.conflict:
        try:
            retry = False
            res = requests.get(url, params=params)
            if not res.from_cache:
                time.sleep(0.25) # wait for 250ms as according to OMIM specs
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

    return res.text

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
