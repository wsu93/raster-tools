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

import argparse

import numpy as np
from scipy import ndimage

from raster_tools import datasets
from raster_tools import groups
from raster_tools import utils

from raster_tools import gdal

GTIF = gdal.GetDriverByName(str('gtiff'))
SHAPE = ogr.GetDriverByName(str('esri shapefile'))
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


def fill(values, no_data_value):
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
    aggregated_shape = values.shape[0] / 2, values.shape[1] / 2
    aggregated = utils.aggregate_uneven(func='mean',
                                        values=values,
                                        no_data_value=no_data_value)

    filled = fill(ceiling=ceiling, **aggregated)
    zoomed = zoom(filled)[:values.shape[0], :values.shape[1]]
    return np.where(mask, smooth(zoomed), values)


class Filler(object):
    def __init__(self, source_path, target_path, single=True):
        """ 
        If edges_path is given, write edge void geometries to a shapefile.
        Otherwise process only a single void spanning the whole geometry.
        """
        # source group
        self.source = groups.Group(gdal.Open(source_path))
        self.target_path = target_path

    def fill(self, feature):
        """
        Return dictionary with data and no data value of filling.
        Write filled raster and edge shapes to target directory.
        """
        # target path
        name = feature[str('name')]
        root_path = os.path.join(self.target_path, name[:3], name)
        tif_path = root_path + '.tif'
        shp_path = root_path + '.shp'
        
        if os.path.exists(tif_path):
            continue

        # create directory
        try:
            os.makedirs(os.path.dirname(path))
        except OSError:
            pass  # no problem

        values = self.source.read(geometry)
        no_data_value = self.source.no_data_value
        result = fill(values=values, no_data_value=no_data_value)
        result[values != no_data_value] = no_data_value
        return {'values': result, 'no_data_value': no_data_value}


def fillnodata(source_path, target_path, ceiling_path):
    """
    Fill a single raster.
    """
    source_dataset = gdal.Open(source_path)
    geometry = utils.get_geometry(source_dataset)

    filler = Filler(source_path=source_path, ceiling_path=ceiling_path)
    result = filler.fill(geometry)

    kwargs = {'projection': source_dataset.GetProjection(),
              'geo_transform': source_dataset.GetGeoTransform()}
    kwargs['array'] = result['values'][np.newaxis]
    kwargs['no_data_value'] = result['no_data_value'].item()

    with datasets.Dataset(**kwargs) as dataset:
        GTIF.CreateCopy(target_path, dataset, options=OPTIONS)

def fillnodata(feature_path, source_path, target_path, single, part):
    """
    """
    # select some or all polygons
    features = datasources.PartialDataSource(features_path)
    if part is not None:
        features = features.select(part)

    filler = Filler(source_path=source_path)

    for feature in features:
        # target path
        name = feature[str('name')]
        path = os.path.join(output_path, name[:2], '{}.tif'.format(name))
        if os.path.exists(path):
            continue

        # create directory
        try:
            os.makedirs(os.path.dirname(path))
        except OSError:
            pass  # no problem

        # geometries
        inner_geometry = feature.geometry()
        outer_geometry = inner_geometry.Buffer(32, 1)

        # geo transforms
        geo_transform = filler.source.geo_transform
        inner_geo_transform = geo_transform.shifted(inner_geometry)
        outer_geo_transform = geo_transform.shifted(outer_geometry)

        # fill
        result = filler.fill(outer_geometry)

        # cut out
        slices = outer_geo_transform.get_slices(inner_geometry)
        values = result['values'][slices]
        no_data_value = result['no_data_value']
        if np.equal(values, no_data_value).all():
            continue

        # save
        options = ['compress=deflate', 'tiled=yes']
        kwargs = {'projection': filler.source.projection,
                  'geo_transform': inner_geo_transform,
                  'no_data_value': no_data_value.item()}

        with datasets.Dataset(values[np.newaxis], **kwargs) as dataset:
            GTIF.CreateCopy(path, dataset, options=options)


def get_parser():
    """ Return argument parser. """
    parser = argparse.ArgumentParser(
        description=__doc__
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
        'output_path',
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
        help='statistic to use per pixel quad prior to smoothing',
    )
    parser.add_argument(
        '-p', '--part',
        help='partial processing source, for example "2/3"',
    )
    return parser


def main():
    """ Call fillnodata with args from parser. """
    kwargs = vars(get_parser().parse_args())
    fillnodata(**kwargs)
