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

CHUNKS = 1250, 1000
DTYPE = 'u8'
COMPRESSION = 'lzf'


def label(raster_path, label_path):
    # raster init
    if isdir(raster_path):
        raster_datasets = [gdal.Open(join(raster_path, path))
                           for path in sorted(listdir(raster_path))]
    else:
        raster_datasets = [gdal.Open(raster_path)]

    # convenient variables
    group = groups.Group(*raster_datasets)
    width = group.width
    height = group.height
    no_data_value = group.no_data_value

    # create hdf5 datasets
    shape = height, width
    mkdir(label_path)
    ds_kwargs = {
        'compression': COMPRESSION,
        'chunks': CHUNKS,
        'dtype': DTYPE,
        'shape': shape,
    }
    h5_label = h5py.File(join(label_path, 'label.h5'))
    h5_south = h5py.File(join(label_path, 'south.h5'))
    h5_east = h5py.File(join(label_path, 'east.h5'))
    h5_final = h5py.File(join(label_path, 'final.h5'))
    ds_label = h5_label.create_dataset('label', **ds_kwargs)
    ds_south = h5_south.create_dataset('south', **ds_kwargs)
    ds_east = h5_east.create_dataset('east', **ds_kwargs)
    ds_final = h5_final.create_dataset('final', **ds_kwargs)

    # label in chunks
    offset = 0
    partitioner = utils.Partitioner(size=shape, chunks=CHUNKS)
    for i1, j1, i2, j2 in partitioner:

        # determine width, height and padding
        ni, nj = i2 - i1, j2 - j1
        pi = 1 if i2 < height else 0
        pj = 1 if j2 < width else 0

        # determine labels and apply the offset
        values = group[i1:i2 + pi, j1:j2 + pj]
        mask = values == no_data_value
        label, count = ndimage.label(mask)
        label[mask] += offset

        # store offset labels in the label dataset
        ds_label[i1:i2, j1:j2] = label[:ni, :nj]

        # store the south and east edges as well
        ds_south[i2:i2 + pi, j1:j2] = label[ni:, :nj]
        ds_east[i1:i2, j2:j2 + pj] = label[:ni, nj:]

        # increment the offset to prevent duplicate void numbers
        offset += count

    # determine connections linking labels to south and east strips
    links = set()
    for i in range(CHUNKS[0], height, CHUNKS[0]):
        links.update(filter(any, zip(ds_label[i, :], ds_south[i, :])))
    for j in range(CHUNKS[1], width, CHUNKS[1]):
        links.update(filter(any, zip(ds_label[:, j], ds_east[:, j])))

    h5_south.close()
    h5_east.close()

    # build a sparse matrix using the links
    csr_indices = zip(*links)
    csr_values = np.zeros(len(csr_indices[0]))
    csr_kwargs = {
        'dtype': DTYPE,
        'shape': (offset + 1, offset + 1),
    }
    csr = sparse.csr_matrix((csr_values, csr_indices), **csr_kwargs)

    # do components analysis and build a conversion array
    components = connected_components(csr, directed=False)[1]
    voids = np.arange(offset + 1)
    minima = ndimage.minimum(voids, components, voids)
    convert = minima[components]

    # apply the conversion to a final dataset
    partitioner = utils.Partitioner(size=shape, chunks=CHUNKS)
    for i1, j1, i2, j2 in partitioner:
        ds_final[i1:i2, j1:j2] = convert[ds_label[i1:i2, j1:j2]]
    ds_final.attrs['maximum'] = offset

    h5_label.close()
    h5_final.close()


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
    """ Call label with args from parser. """
    kwargs = vars(get_parser().parse_args())
    label(**kwargs)
