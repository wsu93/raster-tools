# (c) Nelen & Schuurmans, see LICENSE.rst.
# -*- coding: utf-8 -*-
"""
Aggregate based on labels.
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from __future__ import division

import argparse

from os import mkdir, listdir
from os.path import isdir, join

from osgeo import gdal
from scipy import ndimage
from scipy import sparse
from scipy.sparse.csgraph import connected_components

import numpy as np
import h5py

from raster_tools import groups
from raster_tools import utils

DTYPE = 'u8'
COMPRESSION = 'lzf'
CHUNKS = 256, 256

FOOTPRINT1 = np.array([[0, 1, 0],
                       [1, 1, 0],
                       [0, 0, 0]], dtype='b1')

FOOTPRINT2 = np.array([[0, 0, 0],
                       [0, 1, 1],
                       [0, 1, 0]], dtype='b1')

def aggregate(raster_path, label_path, aggregate_path):
    # raster init
    if isdir(raster_path):
        raster_datasets = [gdal.Open(join(raster_path, path))
                           for path in sorted(listdir(raster_path))]
    else:
        raster_datasets = [gdal.Open(raster_path)]

    # label init
    h5_label = h5py.File(join(label_path, 'final.h5'), 'r')
    ds_label = h5_label['final']
    depth = ds_label.attrs['maximum'] + 1

    # convenient variables
    group = groups.Group(*raster_datasets)
    width = group.width
    height = group.height
    no_data_value = group.no_data_value

    # open and create hdf5 datasets

    mkdir(aggregate_path)
    ds_kwargs = {
        'dtype': group.dtype,
        'chunks': (1,) + CHUNKS,
        'compression': COMPRESSION,
        'shape': (depth, height, width),
    }
    h5_aggregate = h5py.File(join(aggregate_path, 'aggregate.h5'))
    ds_aggregate = h5_aggregate.create_dataset('aggregate', **ds_kwargs)

    # label in chunks
    partitioner = utils.Partitioner(size=(height, width), chunks=CHUNKS)
    for i1, j1, i2, j2 in partitioner:

        # determine width, height and padding
        ni, nj = i2 - i1, j2 - j1
        pi1 = 0 if i1 == 0 else -1
        pj1 = 0 if j1 == 0 else -1
        pi2 = 1 if i2 < height else 0
        pj2 = 1 if j2 < width else 0

        slices = slice(i1 + pi1, i2 + pi2), slice(j1 + pj1, j2 + pj2)

        # # determine labels and apply the offset
        values = group[slices]
        label = ds_label[slices]
        
        mask = values == no_data_value
        
        edge1 = ndimage.grey_dilation(label, footprint=FOOTPRINT1) - label
        edge2 = ndimage.grey_dilation(label, footprint=FOOTPRINT2) - label

        nonzero1 = edge1.nonzero()
        nonzero2 = edge2.nonzero()

        depth1 = edge1[nonzero1]
        depth2 = edge2[nonzero2]

        depths = np.union1d(depth1, depth2)
        if depths.size < 10:
            continue

        convert = np.empty(int(depths.max() + 1), dtype='i8')
        convert[depths] = np.arange(depths.size)

        deep_index1 = convert[depth1], nonzero1[0], nonzero1[1]
        deep_index2 = convert[depth2], nonzero2[0], nonzero2[1]

        work = ds_aggregate[depths, slices[0], slices[1]]
        work[deep_index1] = values[nonzero1]
        work[deep_index2] = values[nonzero2]
        ds_aggregate[depths, slices[0], slices[1]] = work

        import ipdb
        ipdb.set_trace() 


    h5_label.close()
    h5_aggregate.close()


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
        metavar='LABEL',
        help='Inputfile containing the labels.',
    )
    parser.add_argument(
        'aggregate_path',
        metavar='LABEL',
        help='Outputfile containing the aggregations.',
    )
    return parser


def main():
    """ Call aggregate with args from parser. """
    kwargs = vars(get_parser().parse_args())
    aggregate(**kwargs)
