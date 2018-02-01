# -*- coding: utf_8 -*-
from flask import (
    Blueprint,
    redirect,
    jsonify,
g,
    request,
    url_for,
    Response,
json,


    render_template)
from datetime import datetime
import time
from mmc.utils import format_time
import json
from werkzeug.datastructures import Headers
from flask import current_app
from mmc.api.converter import MouseMapFactory
from mmc.api.converter import MouseMapFactory
from mmc.api.genome_feature import get_genes
from mmc.api.genome_feature import get_snps
from werkzeug.exceptions import BadRequest
from collections import OrderedDict

from mmc import utils

api = Blueprint('api', __name__, template_folder='templates', url_prefix='/api')


MMC = None


def fix_chrom(chrom):
    chrom = str(chrom).upper()
    for prefix in ('CHR', 'CHROM', 'CHROMSOME'):
        chrom = chrom.replace(prefix, '')
        break

    if chrom == 'X':
        return '20'
    elif chrom == 'Y':
        return '21'

    return chrom


def split_url(url):
    """Splits the given URL into a tuple of (protocol, host, uri)"""
    proto, rest = url.split(':', 1)
    rest = rest[2:].split('/', 1)
    host, uri = (rest[0], rest[1]) if len(rest) == 2 else (rest[0], "")
    return proto, host, uri


def parse_url(url):
    """Parses out Referer info indicating the request is from a previously proxied page.
    For example, if:
        Referer: http://localhost:8080/proxy/google.com/search?q=foo
    then the result is:
        ("google.com", "search?q=foo")
    """

    proto, host, uri = split_url(url)

    print('url=', url)
    print('proto=', proto)
    print('host=', host)
    print('uri=', uri)

    if uri.find("/") < 0:
        return None

    rest = uri.split("/", 2)[2]

    return {'proto': proto, 'host': host, 'uri': uri, 'url': rest}


@api.record_once
def record_params(setup_state):
    print(setup_state.app.config)


@api.route('/units', methods=['GET', 'POST'])
def api_units():
    if 'MMF' not in g:
        g.MMF = MouseMapFactory(current_app.config['MMC_CONFIG'])

    mmf = g.MMF
    conversion_maps = mmf.conversion_map

    display = []
    for idx, elem in enumerate(current_app.config['CONVERSION_UNITS']):
        display.append(elem)

    units = {}
    for curr_from_unit, entry in sorted(conversion_maps.items()):
        print(curr_from_unit, list(sorted(entry.keys())))
        units[curr_from_unit] = list(sorted(entry.keys()))

    # genes and snps
    units['genes'] = sorted(units[mmf.config['genomeFeatureIDFromUnit']] + ['mm10'])
    units['snps'] = sorted(units[mmf.config['genomeFeatureIDFromUnit']] + ['mm10'])

    return jsonify({
        'display': display,
        'conversion': units
    })


@api.route('/convert', methods=['POST'])
def api_convert():
    if 'MMF' not in g:
        g.MMF = MouseMapFactory(current_app.config['MMC_CONFIG'])

    mmf = g.MMF
    json_data = request.get_json()

    if 'includeErrors' in json_data:
        include_errors = utils.str2bool(json_data['includeErrors'])
    else:
        include_errors = True

    if 'includeInput' in json_data:
        include_input = utils.str2bool(json_data['includeInput'])
    else:
        include_input = True

    if 'includeExtra' in json_data:
        include_extra = utils.str2bool(json_data['includeExtra'])
    else:
        include_extra = True

    if 'fromUnit' in json_data:
        from_unit = json_data['fromUnit']
    else:
        from_unit = 'genomeFeatureID'

    to_unit = json_data['toUnit']
    data = json_data['dataToConvert']

    print('includeErrors: ', include_errors)
    print('includeInput: ', include_input)
    print('includeExtra: ', include_extra)
    print('fromUnit: ', from_unit)
    print('toUnit: ', to_unit)
    print('dataToConvert: ', data)

    converted = []
    errors = []

    if from_unit == 'snps':
        from_unit = mmf.config['genomeFeatureIDFromUnit']

        print('from_unit now = ', from_unit)

        # make a map of al the snps
        ordered_data = OrderedDict()
        for d in data:
            ordered_data[d[0]] = 1

        snps = get_snps(data, 91, 'Mm')

        if from_unit == to_unit:
            # special case of just getting the locations and returning
            for snp in snps['snps']:
                if include_input:
                    ordered_data[snp[2]] = [snp[2], snp[0], snp[1]]
                else:
                    ordered_data[snp[2]] = [snp[0], snp[1]]
        else:
            the_map = mmf.get_map(from_unit, to_unit)

            for snp in snps['snps']:
                try:
                    # original data
                    original_chrom = snp[0]
                    original_value = snp[1]

                    # massaged data
                    fixed_chrom = fix_chrom(original_chrom)
                    fixed_value = float(original_value)

                    # converted value
                    converted_value = the_map.convert(fixed_chrom, fixed_value)

                    if include_input:
                        ordered_data[snp[2]] = [snp[2], fixed_chrom, converted_value]
                    else:
                        ordered_data[snp[2]] = [fixed_chrom, converted_value]
                except Exception as e:
                    current_app.logger.error(e)
                    if include_errors:
                        if include_input:
                            ordered_data[snp[2]] = [snp[2], original_chrom, original_value]
                        else:
                            ordered_data[snp[2]] = [None, None]
                    errors.append([original_chrom, original_value, str(e)])

        for unknown in snps['unknown']:
            if include_errors:
                if include_input:
                    ordered_data[unknown] = [unknown, None, None]
                else:
                    ordered_data[unknown] = [None, None]
            errors.append([unknown, 'Unknown SNP'])

        to_return = {'results': list(ordered_data.values())}
    elif from_unit == 'genes':
        from_unit = mmf.config['genomeFeatureIDFromUnit']

        print('from_unit now = ', from_unit)

        # make a map of al the genes
        ordered_data = OrderedDict()
        for d in data:
            ordered_data[d[0]] = 1

        genes = get_genes(data, 91, 'Mm')
        print(genes)

        if from_unit == to_unit:
            # special case of just getting the locations and returning
            for gene_id, gene in genes['genes'].items():
                if gene:
                    if include_input:
                        ordered_data[gene_id] = [gene_id, gene['chromosome'], gene['start'], gene['end']]
                    else:
                        ordered_data[gene_id] = [gene['chromosome'], gene['start'], gene['end']]
                else:
                    if include_input:
                        ordered_data[gene_id] = [gene_id, None, None, None]
                    else:
                        ordered_data[gene_id] = [None, None, None]
        else:
            the_map = mmf.get_map(from_unit, to_unit)

            for gene_id, gene in genes['genes'].items():
                if gene:
                    try:
                        # original data
                        original_chrom = gene['chromosome']
                        original_value_1 = gene['start']
                        original_value_2 = gene['end']

                        # massaged data
                        fixed_chrom = fix_chrom(original_chrom)
                        fixed_value_1 = int(original_value_1)
                        fixed_value_2 = int(original_value_2)

                        # converted value
                        converted_value_1 = the_map.convert(fixed_chrom,
                                                            fixed_value_1)

                        converted_value_2 = the_map.convert(fixed_chrom,
                                                            fixed_value_2)

                        if include_input:
                            ordered_data[gene_id] = [gene_id, fixed_chrom, converted_value_1, converted_value_2]
                        else:
                            ordered_data[gene_id] = [fixed_chrom, converted_value_1, converted_value_2]

                    except Exception as e:
                        current_app.logger.error(e)
                        if include_errors:
                            if include_input:
                                ordered_data[gene_id] = [gene_id, None, None, None]
                            else:
                                ordered_data[gene_id] = [None, None, None]
                        errors.append([gene_id, str(e)])
                else:
                    if include_input:
                        ordered_data[gene_id] = [gene_id, None, None, None]
                    else:
                        ordered_data[gene_id] = [None, None, None]
                    errors.append([gene_id, 'Unknown gene'])

        to_return = {'results': list(ordered_data.values())}
    else:
        the_map = mmf.get_map(from_unit, to_unit)

        for to_convert in data:
            try:
                if len(to_convert) < 2:
                    raise ValueError('{} does not contain at least 2 values'.format(str(to_convert)))

                # original data
                original_chrom = to_convert[0]
                original_value_1 = to_convert[1]
                original_value_2 = None

                if len(to_convert) > 2:
                    original_value_2 = to_convert[2]

                # massaged data
                fixed_chrom = fix_chrom(original_chrom)
                fixed_value_1 = float(original_value_1)

                input_list = [original_chrom, fixed_value_1, None]

                if original_value_2:
                    fixed_value_2 = float(original_value_2)
                    input_list[2] = fixed_value_2

                # converted value
                converted_value_1 = the_map.convert(fixed_chrom, fixed_value_1)

                output_list = [fixed_chrom, converted_value_1, None]

                if original_value_2:
                    converted_value_2 = the_map.convert(fixed_chrom,
                                                        fixed_value_2)

                    output_list[2] = converted_value_2

                if include_input:
                    converted.append(input_list + output_list)
                else:
                    converted.append(output_list)
            except Exception as e:
                current_app.logger.error(e)
                if include_errors:
                    if include_input:
                        output_list = [None] * 6
                        for i, x in enumerate(to_convert):
                            output_list[i] = x

                        converted.append(output_list)
                    else:
                        output_list = [None] * 3
                        converted.append(output_list)

                error_list = [None] * 3
                for i, x in enumerate(to_convert):
                    error_list[i] = x

                error_list.append(str(e))

                errors.append(error_list)

        # construct return object
        to_return = {'results': converted}

    if include_extra:
        to_return['extra'] = {'errors': errors}

    return jsonify(to_return)


@api.route('/save', methods=['POST'])
def api_filesave():
    '''
     * The whole purpose of this function is to echo the conversion results
     * back to the client so that they get a "save as" dialog. There are client
     * side ways to do this but unfortunately they're not yet suported by all
     * browsers
     * @param conversionResults text of the conversion results
     * @return  the same text that was passed in
     */

    :return:
    '''
    results = request.form.get('textarea-to')

    file_name = 'mmc-results.txt'


    headers = Headers()
    headers.add('Cache-Control', 'must-revalidate')
    headers.add('Pragma', 'must-revalidate')
    headers.add('Content-type', 'text/plain')
    headers.add('Content-Disposition', 'attachment; filename={}'.format(file_name))

    def generate():
        for row in results.split('\n'):
            yield row.strip() + '\n'

    return Response(generate(), mimetype='text/plain', headers=headers)


'''
@api.route('/<path:path>')
def page_not_found(path):
    return jsonify({'error': 'not found',
                    'path': path})

'''