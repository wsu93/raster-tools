# (c) Nelen & Schuurmans, see LICENSE.rst.
# -*- coding: utf-8 -*-
"""
Fills nodata in a raster, partitioning the work by geometries supplied via
shapefile, writing the output in a separate rasters.

The main principle is to aggregate the raster per quad of pixels and repeat
until only one pixel is left, and then zooming back in and smoothing at each
zoom step.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from __future__ import division

from os.path import dirname, exists, join
import argparse
import os

from scipy import ndimage
import numpy as np

from raster_tools import datasets
from raster_tools import datasources
from raster_tools import groups
from raster_tools import utils

from osgeo import gdal
from osgeo import ogr

TIF = gdal.GetDriverByName(str('GTiff'))
MEM = ogr.GetDriverByName(str('Memory'))
SHP = ogr.GetDriverByName(str('ESRI Shapefile'))

OPTIONS = ['compress=deflate', 'tiled=yes']

KERNEL = np.array([[0.0625, 0.1250,  0.0625],
                   [0.1250, 0.2500,  0.1250],
                   [0.0625, 0.1250,  0.0625]])


def zoom(values):
    """ Return zoomed array. """
    return values.repeat(2, axis=0).repeat(2, axis=1)


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


def fillnodata(feature_path, part, **kwargs):
    """
    """
    # select some or all polygons
    features = datasources.PartialDataSource(feature_path)
    if part is not None:
        features = features.select(part)

    filler = Filler(**kwargs)

    for feature in features:
        filler.fill(feature)


def get_parser():
    """ Return argument parser. """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        'feature_path',
        metavar='FEATURES',
        help='shapefile with geometries and names of output tiles',
    )
    parser.add_argument(
        'source_path',
        metavar='RASTER',
        help='source GDAL raster dataset with voids'
    )
    parser.add_argument(
        'target_path',
        metavar='OUTPUT',
        help='output folder for rasters containing the fillings',
    )
    parser.add_argument(
        '-s', '--single',
        action='store_true',
        help='Only process one matching void per feature - implies 2nd pass',
    )
    parser.add_argument(
        '-a', '--aggregation',
        dest='func',
        default='min',
        choices=('min', 'mean', 'max'),
        help='statistic to use for quad aggregations',
    )
    parser.add_argument(
        '-p', '--part',
        help='partial processing source, e.g. "2/3"',
    )
    return parser


def main():
    """ Call fillnodata with args from parser. """
    kwargs = vars(get_parser().parse_args())
    fillnodata(**kwargs)
