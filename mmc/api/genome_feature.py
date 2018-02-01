# -*- coding: utf_8 -*-
import json
import re
import requests

API_SNPS = 'http://churchill-lab.jax.org/ensimpl_snps/api/snps'
API_VERSIONS = 'http://churchill-lab.jax.org/ensimpl_snps/api/versions'
API_GENES = 'http://churchill-lab.jax.org/ensimpl/api/genes'
REGEX_RS_NUMBER = re.compile("RS[0-9]{1,}", re.IGNORECASE)
#REGEX_ENSEMBL_MOUSE_ID = re.compile("ENSMUS([EGTP])[0-9]{11}", re.IGNORECASE)
#REGEX_ENSEMBL_HUMAN_ID = re.compile("ENS([EGTP])[0-9]{11}", re.IGNORECASE)


def get_genomic_locations(ids, version, species):
    snps = []
    genes = []

    for _ in ids:
        if REGEX_RS_NUMBER.match(_):
            snps.append(_)
        else:
            # assume it's a gene
            genes.append(_)

    print('snps=', snps)
    print('genes=', genes)

    if len(genes) > 0:
        genes = get_genes(genes, version, species)

    if len(snps) > 0:
        snps = get_snps(snps, version, species)

    print('genes=', genes)
    print('snps=', snps)

    return {
        'genes': genes if len(genes) > 0 else None,
        'snps': snps if len(snps) > 0 else None
    }



def get_snps(ids, version, species):

    try:
        url = API_SNPS
        post_data = {
            'ids': ids,
            'version': version,
            'species': species
        }

        print('url=', url)
        print('post_data=', post_data)
        r = requests.post(url, data=post_data)
        data = json.loads(r.text)

        return data
    except Exception as e:
        print(str(e))
        return None


def get_genes(ids, version, species):
    try:
        url = API_GENES
        post_data = {
            'id': ids,
            'version': version,
            'species': species
        }
        print('url=', url)
        print('post_data=', post_data)
        r = requests.post(url, data=post_data)
        data = json.loads(r.text)

        return data
    except Exception as e:
        print(str(e))
        return None