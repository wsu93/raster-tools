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

from os.path import join
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

TIF = gdal.GetDriverByName(str('GTiff'))
# MEM = ogr.GetDriverByName(str('Memory'))
# SHP = ogr.GetDriverByName(str('ESRI Shapefile'))

UP = 0
LEFT = 1
DOWN = 2
RIGHT = 3

SIZE = 256

OPTIONS = ['COMPRESS=DEFLATE', 'TILED=YES']

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
    Returns filling array for interior voids.

    Modifies the tile object by setting sparse exterior void edges as
    attributes on the tile.

    # and here is the algo to relate edge voids:
    result = np.logical_and(label_a, label_b)
    join on the bases of this result
    """
    # load and label data
    no_data_value = tile.no_data_value
    data = tile.get_data()
    mask = np.equal(data, no_data_value)
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

    # define object bboxes as having inclusive ENDS as opposed to slices
    bboxes = (objects[inside] + [[-1, 0, -1, 0]])

    # have an array that tracks the labels that still have to be processed
    tracker = np.arange(1, 1 + total)[inside]

    # start at 4 since a single pixel including its edges won't fit size 2
    size = 4

    # loop superpixel sizes
    result = np.full_like(data, no_data_value)
    while len(tracker):
        # floor dividing by size results the index into quads of current size
        regions = bboxes // size
        top, bottom, left, right = regions.transpose()
        # s1 is a boolean index selecting objects that fit in the quads
        s1 = (top == bottom) & (left == right)

        # have an array that tracks the subset of labels for the current size
        subtracker = tracker[s1]
        subregions = regions[s1]
        # loop until all voids within the same superpixel are filled
        while len(subtracker):
            # now for the current quad size, subtracker contains the label
            # numbers and subregions contains the indices into a the quad grid

            # s2 is a boolean index that select unique indices
            s2 = np.zeros(len(subtracker), dtype='b1')
            s2[np.unique(
                SIZE * subregions[:, 0] + subregions[:, 2],
                return_index=True,
            )[1]] = True

            # using the subtracker
            void = np.in1d(label, subtracker[s2]).reshape(label.shape)
            edge = ndimage.binary_dilation(void) - void
            values = np.where(edge, data, no_data_value).astype('f4')
            filled = fill(
                func='min', values=values, no_data_value=no_data_value
            )
            result[void] = filled[void]

            # remove processed voids from tracker and regions
            subtracker = subtracker[~s2]
            subregions = subregions[~s2]

        # drop s1 selection from tracker, because it is now processed
        tracker = tracker[~s1]
        bboxes = bboxes[~s1]

        # increase the size for the next loop
        size *= 2

    # write edge info
    # top, bottom, left, right
    tile.edges = (
        label[:, 0].copy(),
        label[0].copy(),
        label[:, -1].copy(),
        label[-1].copy(),
    )

    # keep edges that refer to voids via labels
    tile.voids = []
    for e in tile.edges:
        side = {}
        for l in np.unique(e):
            if l == 0:
                continue
            void = label == l
            edge = ndimage.binary_dilation(void) - void
            index = edge.nonzero()
            (x, y), z = index, data[index]
        try:
            side[l] = np.vstack([x, y, z])
        except UnboundLocalError:
            pass
    tile.voids.append(side)

    return result


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

def fillnodata(source_path, **kwargs):
    from scipy import ndimage
    f1 = np.array([[0, 1, 0],
                   [1, 1, 0],
                   [0, 0, 0]], dtype='b1')
    f2 = np.array([[0, 0, 0],
                   [0, 1, 1],
                   [0, 1, 0]], dtype='b1')

    l = np.array([[0, 0, 0, 0, 0],
                  [0, 1, 1, 1, 0],
                  [0, 0, 1, 1, 0],
                  [0, 2, 0, 1, 0],
                  [0, 0, 0, 0, 0]], 'u1')

    e1 = ndimage.grey_dilation(l, footprint=f1) - l
    e2 = ndimage.grey_dilation(l, footprint=f2) - l

    n1 = e1.nonzero()
    n2 = e2.nonzero()

    i, j = np.mgrid[0:5, 0:5]

    
    m = e1[n1].tolist() + e2[n2].tolist()
    u = i[n1].tolist() + i[n2].tolist()
    v = j[n1].tolist() + j[n2].tolist()

    import ipdb
    ipdb.set_trace() 
    




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
    # import cProfile
    # pr = cProfile.Profile()
    # pr.enable()

    kwargs = vars(get_parser().parse_args())
    fillnodata(**kwargs)

    # pr.disable()
    # pr.print_stats(sort='time')
