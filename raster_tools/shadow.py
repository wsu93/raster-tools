# -*- coding: utf-8 -*-
# (c) Nelen & Schuurmans, see LICENSE.rst.
"""
Calculate shadows.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from __future__ import division

import argparse
import itertools
import logging
import math
import os
import sys

import numpy as np
from scipy import ndimage

from raster_tools import gdal
from raster_tools import datasets
from raster_tools import utils
from raster_tools import groups

logger = logging.getLogger(__name__)
driver = gdal.GetDriverByName(str('gtiff'))


class Shadower(object):
    """
    2. Do shifts and subtractions, compare to original and mark where
       shifted > original. And take out of index.
    3. Stop when nothing changes anymore.
    """
    def __init__(self, raster_path, output_path):

        # put the input raster(s) in a group
        if os.path.isdir(raster_path):
            datasets = [gdal.Open(os.path.join(raster_path, path))
                        for path in sorted(os.listdir(raster_path))]
        else:
            datasets = [gdal.Open(raster_path)]
        self.group = groups.Group(*datasets)

        # TODO Check this is june 12, 15:00
        azimuth = 236
        elevation = 50
        slope = math.tan(math.radians(elevation))
        pixel = self.group.geo_transform[1]

        dx = -math.sin(math.radians(azimuth))
        dy = math.cos(math.radians(azimuth))

        # calculate shift and corresponding elevation change
        self.ds = 1 / max(dx, dy)          # pixels
        self.dz = self.ds * slope * pixel  # meters
        self.dx = dx * self.ds             # pixels
        self.dy = dy * self.ds             # pixels

        # calculate margin for input data
        mz = 367  # gerbrandy tower, in meters
        ms = mz / slope / pixel                                     # pixels
        self.mx = int(math.copysign(math.ceil(abs(dx * ms)), -dx))  # pixels
        self.my = int(math.copysign(math.ceil(abs(dy * ms)), -dy))  # pixels

        self.output_path = output_path

    def get_bounds_and_slices(self, geometry):
        """
        Return the window into the source raster that includes the
        required margin, and the slices that return from that window
        the part corresponding to geometry.
        """
        x1, y1, x2, y2 = self.group.geo_transform.get_indices(geometry)

        # bounds
        bounds = (min(x1, x1 + self.mx),
                  min(y1, y1 + self.my),
                  max(x2, x2 + self.mx),
                  max(y2, y2 + self.my))

        # slices
        w, h = x2 - x1, y2 - y1
        slices = (slice(None, h) if self.my > 0 else slice(-h, None),
                  slice(None, w) if self.mx > 0 else slice(-w, None))

        return bounds, slices

    def shadow(self, feature):
        geometry = feature.geometry()
        bounds, slices = self.get_bounds_and_slices(geometry)

        # prepare
        array1 = self.group.read(bounds)
        array2 = np.empty_like(array1)
        view1 = array1[slices]
        view2 = array2[slices]
        target = np.zeros_like(view1, dtype='b1')

        # calculate shadow
        for iteration in itertools.count(1):
            shift = iteration * self.dy, iteration * self.dx
            ndimage.shift(array1, shift, array2, order=0)
            view2 -= iteration * self.dz
            index = np.logical_and(~target, view2 > view1)
            print(index.sum())
            if not index.any():
                break
            target[index] = True
        target = target.astype('u1')

        # target path
        leaf_number = feature[b'BLADNR']
        path = os.path.join(self.output_path,
                            leaf_number[:3],
                            '{}.tif'.format(leaf_number))
        if os.path.exists(path):
            logger.debug('Target already exists.')
            return

        # create directory
        try:
            os.makedirs(os.path.dirname(path))
        except OSError:
            pass  # no problem

        kwargs = {
            'no_data_value': 0,
            'projection': self.group.projection,
            'geo_transform': self.group.geo_transform.shifted(geometry),
        }
        options = [
            'tiled=yes',
            'compress=deflate',
        ]
        with datasets.Dataset(target[np.newaxis], **kwargs) as dataset:
            driver.CreateCopy(path, dataset, options=options)


def shadow(index_path, raster_path, output_path, part):
    """
    """
    index = utils.PartialDataSource(index_path)
    if part is not None:
        index = index.select(part)

    shadower = Shadower(output_path=output_path,
                        raster_path=raster_path)

    for feature in index:
        shadower.shadow(feature)
    return 0


def get_parser():
    """ Return argument parser. """
    parser = argparse.ArgumentParser(
        description=__doc__
    )
    parser.add_argument(
        'index_path',
        metavar='INDEX',
        help='shapefile with geometries and names of output tiles',
    )
    parser.add_argument(
        'raster_path',
        metavar='RASTER',
        help='source GDAL raster dataset with voids'
    )
    parser.add_argument(
        'output_path',
        metavar='OUTPUT',
        help='target folder',
    )
    parser.add_argument(
        '-p', '--part',
        help='partial processing source, for example "2/3"',
    )
    parser.add_argument('-v', '--verbose', action='store_true')
    return parser


def main():
    """ Call shadow with args from parser. """
    # logging
    kwargs = vars(get_parser().parse_args())
    if kwargs.pop('verbose'):
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(stream=sys.stderr, level=level, format='%(message)s')

    # run or fail
    try:
        shadow(**kwargs)
        return 0
    except:
        logger.exception('An exception has occurred.')
        return 1
