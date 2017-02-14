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

CHUNKS = 256, 256
COMPRESSION = 'lzf'


def partition(size, chunks=CHUNKS):
    """
    Return generator of (i1, j1, i2, j2) index tuples.

    :param size: Size of the array that is processed
    :param chunks: Maximum size of the chunks. Chunks at lower right
        edges may be smaller.

    :type size: 2-tuple
    :type chunks: 2-tuple

    The returned offset and shape can be used to process an array in
    chunks.
    """

    # offsets of the chunks into the raster
    i1, j1 = np.ogrid[0:size[0]:chunks[0], 0:size[1]:chunks[1]]

    # stops per chunk, either offset plus chunksize, or size if at edge
    i2, j2 = (np.minimum(i1 + chunks[0], size[0]),
              np.minimum(j1 + chunks[1], size[1]))

    # use broadcasting to turn this into a generator of indices
    for lower, upper in zip(np.broadcast(i1, j1), np.broadcast(i2, j2)):
        yield lower + upper


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


    # create hdf5 dataset, one for labels and one for strips
    shape = height, width
    kwargs = {
        'compression': COMPRESSION,
        'chunks': CHUNKS,
        'shape': shape,
        'dtype': 'u4',
    }
    h5_file = h5py.File(label_path)
    h5_label = h5_file.create_dataset('label', **kwargs)
    h5_strip = h5_file.create_dataset('strip', **kwargs)

    # label in chunks
    slices = []
    offset = 0

    for i1, j1, i2, j2 in partition(shape):
        
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
        h5_label[i1:i2, j1:j2] = label[:ni, :nj]

        # store it in the strips, too
        h5_strip[i2:i2 + pi, j1:j2] = label[ni:, :nj]  # south
        h5_strip[i1:i2, j2:j2 + pj] = label[:ni, nj:]  # east

        # increment the offset to prevent duplicate void numbers
        offset += count

    # determine connections
    label = np.concatenate([h5_label[256::256, :].flatten(),
                            h5_label[:, 256::256].flatten()])
    strip = np.concatenate([h5_strip[256::256, :].flatten(),
                            h5_strip[:, 256::256].flatten()])
    links = np.logical_and(label, strip)
    label, strip = label[links], strip[links]

    # now what.
    return

    from pylab import show, imshow
    imshow(h5_strip, interpolation='none')
    show()
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
