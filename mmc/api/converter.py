import json
import os
import sys
#from collections import OrderedDict

from mmc import utils

GENO_COORD_MAP = {}


class GenoCoordMapNormalizer():
    def __init__(self, from_coord_map, to_coord_map):
        self.from_coord_map = from_coord_map
        self.to_coord_map = to_coord_map

    def validate_normalize(self):
        if self.from_coord_map.keys() != self.to_coord_map.keys():
            raise ValueError(
                'chromosomes between from coordinates and to coordinates do not match up.')

        for chrom in self.from_coord_map:
            from_chr_coords = self.from_coord_map[chrom]
            to_chr_coords = self.to_coord_map[chrom]

            if len(from_chr_coords) == 0:
                raise ValueError(
                    'coordinate list for {} is empty'.format(chrom))

            if len(from_chr_coords) != len(to_chr_coords):
                raise ValueError(
                    'coordinate lists for {} are not equal size'.format(chrom))

            if not utils.is_sorted(from_chr_coords):
                raise ValueError(
                    'from coordinates must be given in sorted order')

            if not utils.is_sorted(to_chr_coords):
                raise ValueError('to coordinates must be given in sorted order')

            # we will fix expanses of identical coordinate values as long as
            # they don't overlap between the two maps
            any_same_runs = False
            prev_from_same = False
            prev_to_same = False

            for i in range(1, len(from_chr_coords)):
                curr_from_same = (from_chr_coords[i] == from_chr_coords[i - 1])
                curr_to_same = (to_chr_coords[i] == to_chr_coords[i - 1])

                if curr_from_same or curr_to_same:
                    any_same_runs = True

                if (curr_from_same or prev_from_same) and (
                        curr_to_same or prev_to_same):
                    raise ValueError(
                        'found overlapping identical coordinate runs')

            if any_same_runs:
                self.normalize_identical_runs(from_chr_coords, to_chr_coords);
                self.normalize_identical_runs(to_chr_coords, from_chr_coords);

                self.from_coord_map[chrom] = from_chr_coords
                self.to_coord_map[chrom] = to_chr_coords


    def collapse_identical(self, start, count, coords_from, coords_to):
        min_to = coords_to[start]
        max_to = coords_to[start + count - 1]
        new_val_to = (min_to + max_to) / 2.0

        del coords_from[start + 1:start + count]
        del coords_to[start + 1:start + count]

        coords_to[start] = new_val_to

    def normalize_identical_runs(self, chords_from, chords_to):
        size = len(chords_from)
        ident_count = 0
        prev_val = chords_from[-1]

        for i in range(size - 2, -1, -1):
            curr_val = chords_from[i]

            if curr_val == prev_val:
                ident_count += 1
            else:
                if ident_count >= 1:
                    self.collapse_identical(i + 1, ident_count + 1, chords_from,
                                            chords_to)
                ident_count = 0

            prev_val = curr_val

        if ident_count >= 1:
            self.collapse_identical(0, ident_count + 1, chords_from, chords_to)


class DirectGenoCoordMap():
    def __init__(self, from_coord_map, to_coord_map):
        self.from_coord_map = from_coord_map
        self.to_coord_map = to_coord_map

    def extrapolate_from_to(self, chrom, from_value, from_coord_map,
                            to_coord_map):
        if chrom not in from_coord_map or chrom not in to_coord_map:
            raise ValueError('chrom {} is not valid'.format(chrom))

        from_points = from_coord_map[chrom]
        to_points = to_coord_map[chrom]

        from_start = from_points[0]
        from_end = from_points[-1]
        from_total_dist = from_end - from_start

        to_start = to_points[0]
        to_end = to_points[-1]
        to_total_dist = to_end - to_start

        to_units_per_from_unit = to_total_dist / from_total_dist

        from_val_offset = from_value - from_start
        to_val_offset = from_val_offset * to_units_per_from_unit

        return to_start + to_val_offset

    def convert_from_to(self, chrom, from_value, from_coord_map, to_coord_map):
        if chrom not in from_coord_map or chrom not in to_coord_map:
            raise ValueError('chrom {} is not valid'.format(chrom))

        from_points = from_coord_map[chrom]
        to_points = to_coord_map[chrom]

        index = utils.binary_search(from_points, from_value)

        if index >= 0:
            return to_points[index]

        # convert the index from the negative values which is documented
        # as (-(insertion point) - 1) in the binary search function
        index = -(index + 1)

        if index == 0 or index == len(from_points):
            return self.extrapolate_from_to(chrom, from_value, from_coord_map,
                                            to_coord_map)

        high_from_value = from_points[index]
        low_from_value = from_points[index - 1]
        within_range_from_diff = high_from_value - low_from_value

        high_to_value = to_points[index]
        low_to_value = to_points[index - 1]
        within_range_to_diff = high_to_value - low_to_value

        to_units_per_from_unit = within_range_to_diff / within_range_from_diff

        return low_to_value + to_units_per_from_unit * (
                from_value - low_from_value)

    def get_inverse(self):
        return DirectGenoCoordMap(self.to_coord_map, self.from_coord_map)

    def convert(self, chrom, from_value):
        return self.convert_from_to(chrom, from_value, self.from_coord_map,
                                    self.to_coord_map)


class GenoCoordParser():
    def __init__(self):
        pass

    def parse_genome_coord_file(self, file_name):
        #from_map = OrderedDict()
        #to_map = OrderedDict()
        from_map = {}
        to_map = {}

        i = 0

        with open(file_name, 'r') as fd:
            for line in fd:
                i += 1

                elems = line.strip().split()
                chrom = elems[0].strip()
                from_value = float(elems[1].strip())
                to_value = float(elems[2].strip())

                if chrom not in from_map:
                    from_map[chrom] = []
                    to_map[chrom] = []

                from_map[chrom].append(from_value)
                to_map[chrom].append(to_value)

        # validate and normalize
        normalizer = GenoCoordMapNormalizer(from_map, to_map)
        normalizer.validate_normalize()

        return DirectGenoCoordMap(normalizer.from_coord_map,
                                  normalizer.to_coord_map)


class IndirectGenoCoordMap():
    def __init__(self, from_map, to_map):
        self.from_map = from_map
        self.to_map = to_map

    def convert(self, chrom, from_value):
        return self.to_map.convert(chrom,
                                   self.from_map.convert(chrom, from_value))

    def get_inverse(self):
        return IndirectGenoCoordMap(self.to_map.get_inverse(),
                                    self.from_map.get_inverse())


class MouseMapFactory():
    def __init__(self, config_file):
        self.config_file = config_file
        #self.conversion_map = OrderedDict()
        self.conversion_map = {}

        with open(self.config_file, 'r') as fd:
            self.config = json.load(fd)

        for resource in self.config['conversionFiles']:
            file_to_parse = os.path.join('resources', resource['resourcePath'])

            # print('File to parse: {}'.format(file_to_parse))
            gcm_parser = GenoCoordParser()
            curr_map = gcm_parser.parse_genome_coord_file(file_to_parse)

            self.put_geno_coord_map(resource['fromUnits'],
                                    resource['toUnits'],
                                    curr_map)

            self.put_geno_coord_map(resource['toUnits'],
                                    resource['fromUnits'],
                                    curr_map.get_inverse())

        for ic in self.config['indirectConversions']:
            connect_unit = ic['connectingUnit']
            for from_unit in ic['connectedUnits']:
                for to_unit in ic['connectedUnits']:
                    if (from_unit != to_unit) and (
                    not self.has_map(from_unit, to_unit)):
                        from_map = self.conversion_map[from_unit][connect_unit]
                        to_map = self.conversion_map[connect_unit][to_unit]
                        self.put_geno_coord_map(from_unit, to_unit,
                                                IndirectGenoCoordMap(from_map,
                                                                     to_map))

    def put_geno_coord_map(self, from_coord, to_coord, gcm):
        if from_coord not in self.conversion_map:
            #from_map = OrderedDict()
            from_map = {}
            from_map[to_coord] = gcm
            self.conversion_map[from_coord] = from_map
        else:
            from_map = self.conversion_map[from_coord]
            from_map[to_coord] = gcm
            self.conversion_map[from_coord] = from_map

    def has_map(self, from_coord, to_coord):
        return (from_coord in self.conversion_map) and (
                to_coord in self.conversion_map[from_coord])

    def get_map(self, from_coord, to_coord):
        print('get_map({}, {})'.format(from_coord, to_coord))
        return self.conversion_map[from_coord][to_coord]

    '''
if __name__ == '__main__':
    mmf = init('resources/mmc-conf.json')
    conversion_maps = mmf.conversion_map
    genome_feature_id_to_units = mmf.config['genomeFeatureIDToUnits']

    # for curr_from_units, entry in conversion_maps.items():
    #    print(curr_from_units, entry.keys())

    from_coord = sys.argv[1]
    to_coord = sys.argv[2]
    value = sys.argv[3]

    print(from_coord, to_coord, value)

    # print(from_coord)
    # print(to_coord)
    # print(value)

    the_map = mmf.get_map(from_coord, to_coord)

    # print('map=', the_map)
    print('map-from', the_map.from_coord_map['1'][0:10])
    print('map-to', the_map.to_coord_map['1'][0:10])

    elems = value.split()
    new_value = the_map.convert(elems[0], float(elems[1]))
    print(new_value)
    '''



