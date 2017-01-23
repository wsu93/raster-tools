# (c) Nelen & Schuurmans, see LICENSE.rst.
# -*- coding: utf-8 -*-
"""
Fills nodata in a large raster dataset.

The main principle is to aggregate the raster per quad of pixels and repeat
until only one pixel is left, and then zooming back in and smoothing at each
zoom step.

However, things get complicated since a typical dataset does not fit in
memory.  Therefore the dataset is processed in tiles and per tile the edge
voids are remembered for later merging and storing. If all edge voids for a
tile are resolved, it is written a second time. It is up to the user to
combine the result tiles into a single void dataset.
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from __future__ import division

from os.path import dirname, exists, join
import argparse
import itertools
import os
from hashlib import md5
from struct import pack

from scipy import ndimage
import numpy as np

from raster_tools import datasets
# from raster_tools import datasources
from raster_tools import groups
from raster_tools import utils

from osgeo import gdal
from osgeo import ogr

#
#
#
class BaseExtractor(object):
    """ Extract objects from a void mask. """

    def __init__(self, mask):
        self.label, total = ndimage.label(mask)
        self.objects = ndimage.find_objects(self.label)

    def grow(self, slices):
        """
        Return slices, limited.

        :param slices: 2-tuple of slice objects

        The result is a 2-tuple of grown slice objects. Limited is a boolean
        indicating whether growth was limited by one or more array edges.
        """
        h, w = self.label.shape
        y1, y2, x1, x2 = (slices[0].start, slices[0].stop,
                          slices[1].start, slices[1].stop)

        Y1 = max(0, y1 - 1)
        Y2 = min(h, y2 + 1)
        X1 = max(0, x1 - 1)
        X2 = min(w, x2 + 1)

        slices = (slice(Y1, Y2), slice(X1, X2))
        limited = (y1 == Y1 or y2 == Y2 or x1 == X1 or x2 == X2)

        return slices, limited

    def get_objects(self, objects):
        """ Return slices, void, edge. """
        for slices, label in objects:
            void = self.label[slices] == label
            edge = ndimage.binary_dilation(void) - void
            yield slices, void, edge

    def get_interior_objects(self):
        return self.get_objects(self.interior_objects)

    def get_exterior_objects(self):
        return self.get_objects(self.exterior_objects)


class MultiExtractor(BaseExtractor):
    def __init__(self, *args, **kwargs):
        super(MultiExtractor, self).__init__(*args, **kwargs)
        # split interior and exterior objects in lists of (slices,
        # label)-tuples
        self.interior_objects = []
        self.exterior_objects = []
        for label, slices in enumerate(self.objects, 1):
            slices, limited = self.grow(slices)
            if limited:
                self.exterior_objects.append((slices, label))
            else:
                self.interior_objects.append((slices, label))


class SingleExtractor(BaseExtractor):
    def __init__(self, *args, **kwargs):
        super(SingleExtractor, self).__init__(*args, **kwargs)
        # find the main object in the list of objects
        height, width = self.label.shape
        slices = slice(1, height - 1), slice(1, width - 1)
        label = self.objects.index(slices)
        slices = self.grow(slices)[0]
        self.interior_objects = [(slices, label)]
        self.exterior_objects = []


class Filler(object):
    """
    Fill voids, one by one. The voids have to be filled one by one because a
    recursive aggregation producedure is followed. The edges of one void may
    interfere on some aggregation level with the edges of another void, thereby
    affecting both fillings.
    """
    def __init__(self, source_path, target_path, func, single):
        """
        If edges_path is given, write edge void geometries to a shapefile.
        Otherwise process only a single void spanning the whole geometry.
        """
        self.Extractor = SingleExtractor if single else MultiExtractor
        self.source = groups.Group(gdal.Open(source_path))
        self.target_path = target_path
        self.single = single
        self.func = func

    def fill(self, feature):
        """
        Write filled raster and edge shapes to target directory.
        """
        # determine path
        name = feature[str('unit')]
        root = join(self.target_path, name[:3], name)
        path = root + '.tif'
        if exists(path):
            return

        # retrieve data
        geometry = feature.geometry()
        geo_transform = self.source.geo_transform.shifted(geometry)
        values = self.source.read(geometry)
        no_data_value = self.source.no_data_value
        mask = (values == no_data_value)
        if mask.all():
            return

        extractor = self.Extractor(mask)

        # part one - filling of interior voids
        func = self.func
        result = np.full_like(values, no_data_value)
        for slices, void, edge in extractor.get_interior_objects():
            break
            work = values[slices].copy()
            work[~edge] = no_data_value
            filled = fill(func=func, values=work, no_data_value=no_data_value)
            result[slices][void] = filled[void]

        # create directory
        try:
            os.makedirs(dirname(path))
        except OSError:
            pass  # no problem

        # write tiff
        kwargs = {
            'geo_transform': geo_transform,
            'projection': self.source.projection,
            'no_data_value': no_data_value.item(),
        }
        array = result[np.newaxis]
        with datasets.Dataset(array, **kwargs) as dataset:
            TIF.CreateCopy(path, dataset, options=OPTIONS)

        if self.single:
            return

        # part two - keeping track of exterior voids
        data_source = MEM.CreateDataSource(str(''))
        layer = data_source.CreateLayer(
            str(name),
            srs=geometry.GetSpatialReference(),
        )
        field_defn = ogr.FieldDefn(str('unit'), ogr.OFTString)
        layer.CreateField(field_defn)
        kwargs = {
            'no_data_value': 0,
            'projection': self.source.projection,
        }
        for slices, void, edge in extractor.get_exterior_objects():
            origin = slices[0].start, slices[1].start
            kwargs['geo_transform'] = geo_transform.rebased(origin)
            # array = np.uint8(void | edge)[np.newaxis]
            array = np.uint8(void)[np.newaxis]
            with datasets.Dataset(array, **kwargs) as dataset:
                band = dataset.GetRasterBand(1)
                gdal.Polygonize(band, band, layer, 0)

        # write shape
        SHP.CreateDataSource(root + '.shp').CopyLayer(layer, name)
#
#
#

TIF = gdal.GetDriverByName(str('GTiff'))
MEM = ogr.GetDriverByName(str('Memory'))
SHP = ogr.GetDriverByName(str('ESRI Shapefile'))

SIZE = 600


OPTIONS = ['compress=deflate', 'tiled=yes']

KERNEL = np.array([[0.0625, 0.1250,  0.0625],
                   [0.1250, 0.2500,  0.1250],
                   [0.0625, 0.1250,  0.0625]])


def zoom(values):
    """ Return zoomed array. """
    return values.repeat(2, axis=0).repeat(2, axis=1)


def grow(slices):
    """ Return 2-tuple of grown slice objects. """
    s1, s2 = slices
    return slice(s1.start - 1, s1.stop + 1), slice(s2.start - 1, s2.stop + 1)


def smooth(values):
    """ Two-step uniform for symmetric smoothing. """
    return ndimage.correlate(values, KERNEL)


def block(slices):
    """ Return 2-tuple of superblock slice objects. """
    s1, s2 = slices
    i1, i2 = s1.start, s1.stop - 1
    j1, j2 = s2.start, s2.stop - 1
    d = 1
    while i1 // d != i2 // d or j1 // d != j2 // d:
        d *= 2
    print(d)
    print(d, i1, i2, j1, j2)
    return (slice(i1 // d * d, d * (1 + i1 // d)),
            slice(j1 // d * d, d * (1 + j1 // d)))


from math import ceil, log
def block1(slices):
    s1, s2 = slices
    i1, i2 = s1.start, s1.stop - 1
    j1, j2 = s2.start, s2.stop - 1
    h, w = i2 - i1, j2 - j1
    d = max(2 ** int(ceil(log(h, 2))),
            2 ** int(ceil(log(w, 2))))
    if i1 // d != (i2 - 1) // d or j1 // d != (j2 - 1) // d:
        d *= 2
    return (slice(i1 // d * d, d * (1 + i1 // d)),
            slice(j1 // d * d, d * (1 + j1 // d)))




def fill(func, values, no_data_value):
    """
    Fill must return a filled array. It does so by aggregating, requesting a
    fill for that (so this is a recursive thing), and zooming back. After
    zooming back, it smooths the filled values and returns.
    """
    mask = values == no_data_value
    if not mask.any():
        # this should end the recursion
        return values

    # aggregate
    aggregated = utils.aggregate_uneven(func=func,
                                        values=values,
                                        no_data_value=no_data_value)

    filled = fill(func=func, **aggregated)
    zoomed = zoom(filled)[:values.shape[0], :values.shape[1]]
    return np.where(mask, smooth(zoomed), values)






class TileManager(object):

    def __init__(self, shape_path, source_path):
        # remember source dataset and derivatives
        source = groups.Group(gdal.Open(source_path))
        self.no_data_value = source.no_data_value.item()
        self.geo_transform = source.geo_transform
        self.projection = source.projection
        self.source = source

        # find indices into source from defining shape
        shape = ogr.Open(shape_path)
        layer = shape[0]
        feature = layer[0]
        geometry = feature.geometry()
        bounds = self.geo_transform.get_indices(geometry, inflate=True)

        # have to fix: get_indices currently swaps x and y
        self.bounds = bounds[1], bounds[0], bounds[3], bounds[2]

    def __iter__(self):
        """ Yield tile objects with data load capability. """
        i1, j1, i2, j2 = self.bounds
        i_range = range(i1, i2, SIZE)
        j_range = range(j1, j2, SIZE)
        total = len(i_range) * len(j_range)
        count = itertools.count(1)
        for i in range(i1, i2, SIZE):
            for j in range(j1, j2, SIZE):
                bounds = i, j, min(i + SIZE, i2), min(j + SIZE, j2)
                yield self._get_tile(bounds)
                gdal.TermProgress_nocb(count.next() / total)

    def _get_tile(self, bounds):
        tile = Tile(
            manager=self,
            bounds=bounds,
            projection=self.projection,
            no_data_value=self.no_data_value,
            geo_transform=self.geo_transform.rebased(origin=bounds[:2]),
        )
        return tile

    def get_above_of(self, tile):
        i1, j1, i2, j2 = tile.bounds
        if i1 != self.bounds[0]:
            bounds = i1 - SIZE, j1, i1, j2
            return self._get_tile(bounds)

    def get_left_of(self, tile):
        i1, j1, i2, j2 = tile.bounds
        if j1 != self.bounds[1]:
            bounds = i1, j1 - SIZE, i2, j1
            return self._get_tile(bounds)

    def get_data_for(self, tile):
        # have to fix: source.read() currently swaps x and y
        i1, j1, i2, j2 = tile.bounds
        bounds = j1, i1, j2, i2
        return self.source.read(bounds)


class Tile(object):
    def __init__(self, manager, bounds,
                 projection, no_data_value, geo_transform):

        self.manager = manager
        self.bounds = bounds
        self.projection = projection
        self.no_data_value = no_data_value
        self.geo_transform = geo_transform
        self.hash = md5(pack('4q', *bounds)).hexdigest()

    def get_left(self):
        return self.manager.get_left_of(self)

    def get_above(self):
        return self.manager.get_above_of(self)

    def get_data(self):
        return self.manager.get_data_for(self)

    def __repr__(self):
        i1, j1, i2, j2 = self.bounds
        return str((self.bounds, i2 - i1, j2 - j1))


def process(tile):
    """
    Returns array with interior voids filled.

    Modifies the tile object by setting sparse exterior void edges as
    attributes on the tile.

    # here's a way to restart at value changes
    a = np.array([1,1,1,1,1,1,2,2,2,2,5,5,5,5,5])
    b = np.ones_like(a)
    c = ndimage.sum(b, a, np.arange(-1, 6)).cumsum()
    np.arange(a.size) - c[a]

    # and here is the algo to relate edge voids:
    result = np.logical_and(label_a, label_b)
    join on the bases of this result
    """
    # load and label data
    data = tile.get_data()
    mask = np.equal(data, tile.no_data_value)
    label, total = ndimage.label(mask)
    objects = ndimage.find_objects(label)
    objects = [(s1.start, s1.stop, s2.start, s2.stop) for s1, s2 in objects]
    objects = np.array(objects)

    # identify edge voids from slices
    top, bottom, left, right = objects.transpose()
    height, width = label.shape
    inside = ~np.logical_or.reduce([
        top == 0, bottom == height, left == 0, right == width,
    ])

    # interior = objects[~outside] + [[-1, 1, -1, 1]]  # grown

    # initial conditions
    squares = (objects[inside] - [[0, 1, 0, 1]])
    tracker = np.arange(1, 1 + total)[inside]
    size = 2

    # loop superpixel sizes
    while len(tracker):
        regions = squares // size
        top, bottom, left, right = regions.transpose()
        s1 = (top == bottom) & (left == right)
        # loop until all voids within the same superpixel are filled
        subtracker = tracker[s1]
        subregions = regions[s1]
        while len(subtracker):
            # now tracker[s1] contains the label numbers
            # and squares[s1] contains the superpixel indices
            # s1 is a boolean index into inside objects that fit in quads
            # s2 is a boolean index selecting unique objects per superpixel
            s2 = np.zeros(len(subtracker), dtype='b1')
            s2[np.unique(
                SIZE * subregions[:, 0] + subregions[:, 2],
                return_index=True,
            )[1]] = True

            # 
            # what to do whith these? construct edge array with data and mask
            # array with labels. Both can be just returned to check the result
            import ipdb
            ipdb.set_trace()


            # remove processed voids from
            subtracker = subtracker[~s2]
            subregions = subregions[~s2]
            break
        # drop s1 selection from tracker
        tracker = tracker[~s1]
        # squares? regions?
        size *=2
        break

           



    # edges
    # top, bottom, left, right
    tile.voids = (
    )

    # keep edges that refer to voids via labels
    tile.edges = (
        label[:, 0].copy(),
        label[0].copy(),
        label[:, -1].copy(),
        label[-1].copy(),
    )

    # determine the maximum aggregation factor per void

    
    # determine batch? Just add objects to the batch as long as the object does
    # not fit in the


def fillnodata(shape_path, source_path, target_path):
    """
    Determine extent and create a tile manager.
    """
    os.mkdir(target_path)

    tile_manager = TileManager(shape_path=shape_path, source_path=source_path)
    for tile in tile_manager:
        array = process(tile)[np.newaxis]  # modifies tile, too
        kwargs = {
            'projection': tile.projection,
            'geo_transform': tile.geo_transform,
            'no_data_value': tile.no_data_value,
        }
        path = join(target_path, tile.hash + '-1.tif')
        with datasets.Dataset(array, **kwargs) as dataset:
            TIF.CreateCopy(path, dataset, options=OPTIONS)


def get_parser():
    """ Return argument parser. """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        'shape_path',
        metavar='SHAPE',
        help='source OGR dataset indicating data selection to fill.'
    )
    parser.add_argument(
        'source_path',
        metavar='SOURCE',
        help='source GDAL raster dataset with voids'
    )
    parser.add_argument(
        'target_path',
        metavar='OUTPUT',
        help='output folder for rasters containing the fillings',
    )
    return parser


def main():
    """ Call fillnodata with args from parser. """
    kwargs = vars(get_parser().parse_args())
    fillnodata(**kwargs)
