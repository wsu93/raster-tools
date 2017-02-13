# (c) Nelen & Schuurmans, see LICENSE.rst.
# -*- coding: utf-8 -*-
"""
Label large raster dataset and store it in HDF5 format.
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from __future__ import division

import argparse
import os

from osgeo import gdal
from scipy import ndimage
import numpy as np
import h5py

from raster_tools import groups

BATCH_SIZE = 4 * 256**3
COMPRESSION = 'lzf'


def fillnodata(raster_path, label_path):
    # raster init
    if os.path.isdir(raster_path):
        raster_datasets = [gdal.Open(os.path.join(raster_path, path))
                           for path in sorted(os.listdir(raster_path))]
    else:
        raster_datasets = [gdal.Open(raster_path)]
    group = groups.Group(*raster_datasets)
    width = group.width
    height = group.height
    no_data_value = group.no_data_value

    # determine strip height from the batchsize and group dimensions
    line_size = width * group.dtype.itemsize
    strip_height = BATCH_SIZE // line_size
    strip_height = 100

    offsets = xrange(0, height, strip_height)

    # create hdf5 dataset in rooftile layout
    shape = len(offsets), height, width
    chunks = 1, 256, 256
    dtype = 'u4'

    label = h5py.File(label_path)
    rooftile = label.create_dataset(
        'rooftile',
        dtype=dtype,
        shape=shape,
        chunks=chunks,
        compression=COMPRESSION,
    )

    # label in strips
    x1 = 0
    x2 = width
    for count, offset in enumerate(offsets):
        y1 = offset
        if y1 + strip_height < height:
            y2 = y1 + strip_height + 1
        else:
            y2 = height
        values = group.read((x1, y1, x2, y2))
        mask = values == no_data_value
        label['rooftile'][count, y1:y2, x1:x2] = ndimage.label(mask)
        
    label.close()


def get_parser():
    """ Return argument parser. """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        'raster_path',
        metavar='SOURCE',
        help='Source directory with GDAL raster datasets to process'
    )
    parser.add_argument(
        'label_path',
        metavar='OUTPUT',
        help='Outputfile containing the labels.',
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
