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
    
    INTERIOR = 1
    EXTERIOR = 2

    def __init__(self):
        self.values = values
        self.no_data_value = no_data_value
        
        mask = (values == no_data_value)
        if mask.all():
            raise ValueError('All values are no data.')

        self.label, self.total = ndimage.label(mask)
        objects = ndimage.find_objects(self.label)
    
    def grow(self, slices):
        """ Grow slices by one, but do not exceed shape dims. """
        if s1 == 0
        return tuple(slice(
            max(0, s.start - 1),
            min(l, s.stop + 1))
            for s, l in zip(slices, self.shape))



class MultiExtractor(object):

    def __init__(self, values, no_data_value):
        super(SingleExtractor, self).__init__()

        # split interior and exterior objects
        self.internal_objects = []
        self.external_objects = []

        for count, slices in enumerate(objects, 1):
            x1, x2, y1, y2 = (
                slices[1].start,
                slices[1].stop,
                slices[0].start,
                slices[0].stop,
            )
            if x1 == 0 or y1 == 0 or x2 == width or y2 == height:
                # write including edge. Is that possible?
                continue
            import ipdb
            ipdb.set_trace() 
            slices = grower.grow(slices)

            # determine edges
            # create work array with only edges
            # fill void and write into result array

    def get_external_objects(self)
        # if self.single:
            # slices = slice(1, height - 1), slice(1, width -1)
            # index = objects.index(slices) + 1
            # void = (label == index)

        yield {
            'id':
            'slices':


    def get_internal_objects(self)
        yield 

    def get_single_object(self)


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
        # source group
        self.source = groups.Group(gdal.Open(source_path))
        self.target_path = target_path

    def fill(self, feature):
        """
        Write filled raster and edge shapes to target directory.
        """
        # target path
        name = feature[str('unit')]
        root = join(self.target_path, name[:3], name)
        path = root + 'tif'
        
        if exists(path):
            return

        # retrieve data
        values = self.source.read(feature.geometry())
        no_data_value = self.source.no_data_value
        try:
            extractor = Extractor(values=values, no_data_value=no_data_value)
        except ValueError:
            return  # all values are no data
        
        result = np.full_like(values, no_data_value)

        # create directory
        try:
            os.makedirs(dirname(path))
        except OSError:
            pass  # no problem
        # write tiff
        # kwargs = {'projection': source_dataset.GetProjection(),
                  # 'geo_transform': source_dataset.GetGeoTransform()}
        # kwargs['array'] = result['values'][np.newaxis]
        # kwargs['no_data_value'] = result['no_data_value'].item()
        # with datasets.Dataset(**kwargs) as dataset:
            # GTIF.CreateCopy(path, dataset, options=OPTIONS)

        # write shape


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
