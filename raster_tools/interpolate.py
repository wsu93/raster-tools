# -*- coding: utf-8 -*-
""" TODO Docstring. """

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from __future__ import division

import argparse
import logging
import sys

logger = logging.getLogger(__name__)

import numpy as np
from osgeo import gdal
from osgeo import gdal_array
from osgeo import ogr
from osgeo import osr
from scipy import ndimage

from raster_tools import datasets
from raster_tools import utils

gdal.UseExceptions()
ogr.UseExceptions()
osr.UseExceptions()


def get_parser():
    """ Return argument parser. """
    parser = argparse.ArgumentParser(
        description=__doc__
    )
    parser.add_argument(
        'index_path',
        metavar='INDEX',
    )
    parser.add_argument(
        'mask_path',
        metavar='MASK',
    )
    parser.add_argument(
        'raster_path',
        metavar='RASTER',
    )
    parser.add_argument(
        'output_path',
        metavar='OUTPUT',
    )
    return parser


class Interpolator(object):
    def __init__(self, mask_layer, output_path, raster_dataset, margin=200):
        self.margin = margin
        self.mask_layer = mask_layer
        self.output_path = output_path
        self.raster_dataset = raster_dataset

        self.geometry = utils.get_geometry(raster_dataset)
        self.geo_transform = utils.GeoTransform(
            raster_dataset.GetGeoTransform(),
        )

        # no data value
        band = raster_dataset.GetRasterBand(1)
        data_type = band.DataType
        no_data_value = band.GetNoDataValue()
        self.no_data_value = gdal_array.flip_code(data_type)(no_data_value)

    def get_arrays(self, geometry):
        """ Meta is the three class array. """
        # read the data
        window = self.geo_transform.get_window(geometry)
        data = self.raster_dataset.ReadAsArray(**window)
        meta = np.where(np.equal(data, self.no_data_value), 0, 1).astype('u1')

        # rasterize the water
        kwargs = {'geo_transform': self.geo_transform.shifted(geometry)}
        with datasets.Dataset(meta[np.newaxis], **kwargs) as dataset:
            gdal.RasterizeLayer(dataset,
                                [1],
                                self.mask_layer,
                                burn_values=[2])

        return data, meta

    def interpolate(self, index_feature):
        print(index_feature.items())
        inner_geometry = index_feature.geometry()
        outer_geometry = (inner_geometry
                          .Buffer(self.margin, 1)
                          .Intersection(self.geometry))

        inner_geo_transform = self.geo_transform.shifted(inner_geometry)
        outer_geo_transform = self.geo_transform.shifted(outer_geometry)
        slices = outer_geo_transform.get_slices(inner_geometry)

        data, meta = self.get_arrays(inner_geometry)

        label = ndimage.label(np.equal(meta, 0))[0]
        print(meta[75:85, 165:175])
        print(label[75:85, 165:175])
        #from pylab import imshow, show
        #imshow(meta)
        #show()


def command(index_path, mask_path, raster_path, output_path):
    """
    - use label to find the edge of islands of nodata
    - interpolate or take the min from the part of the edges that have data
    - write to output according to index
    """
    index_data_source = ogr.Open(index_path)
    index_layer = index_data_source[0]

    mask_data_source = ogr.Open(mask_path)
    mask_layer = mask_data_source[0]

    raster_dataset = gdal.Open(raster_path)

    interpolator = Interpolator(mask_layer=mask_layer,
                                output_path=output_path,
                                raster_dataset=raster_dataset)

    for index_feature in index_layer:
        print(index_feature.geometry())
        interpolator.interpolate(index_feature)
        exit()
    return 0


def main():
    """ Call command with args from parser. """
    logging.basicConfig(stream=sys.stderr,
                        level=logging.DEBUG,
                        format='%(message)s')
    try:
        return command(**vars(get_parser().parse_args()))
    except SystemExit:
        raise  # argparse does this
    except:
        logger.exception('An exception has occurred.')
