# -*- coding: utf-8 -*-

# (c) Nelen & Schuurmans, see LICENSE.rst.
"""
Vectorize flow.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from __future__ import division

import argparse
import os

import numpy as np

from raster_tools import utils

from raster_tools import gdal
from raster_tools import ogr
from raster_tools import osr

SHAPE = ogr.GetDriverByName(str('esri shapefile'))
COURSES = np.array([(64, 128, 1),
                    (32,   0, 2),
                    (16,   8, 4)], 'u1')

INDICES = COURSES.nonzero()
NUMBERS = COURSES[INDICES][np.newaxis, ...]
OFFSETS = np.array(INDICES).transpose() - 1


def get_traveled(courses):
    """ Return indices when travelling along courses. """
    # turn indices into points array
    height, width = courses.shape
    indices = (np.arange(height).repeat(width),
               np.tile(np.arange(width), height))
    points = np.array(indices).transpose()

    # determine direction and apply offset
    encode = courses[indices][:, np.newaxis]     # which codes
    select = np.bool8(encode & NUMBERS)          # which courses
    target = points + OFFSETS[select.argmax(1)]  # apply offsets

    return tuple(target.transpose())             # return tuple


def vectorize(direction, accumulation):
    """
    Vectorize flow.

    Key principle is the flow array, that relates the source cell A to
    the target cell B as B = flow[A].

    A special value equal to size indicates the stream leaving the flow.
    """
    # construct a mapping array for the flow
    size = direction.size
    height, width = direction.shape
    traveled = get_traveled(direction)

    # construct the flow array
    flow = np.empty(size + 1, dtype='i8')
    flow[-1] = size
    flow[:size] = np.where(np.logical_or.reduce([
        direction.ravel() == 0,    # undefined cells
        traveled[0] < 0,        # flow-off to the top
        traveled[0] >= height,  # ... bottom
        traveled[1] < 0,        # ... left
        traveled[1] >= width,   # ... right
    ]), size, traveled[0] * width + traveled[1])

    for klass in 1, 2, 3, 4, 5:

        # select points that match klass
        points = (accumulation.ravel() >= klass).nonzero()[0]

        # determine sources, merges and sinks
        sources = points[np.logical_and(flow[points] != size,
                         np.in1d(points, flow[points], invert=True))]
        merges = (np.bincount(flow[points]) > 1).nonzero()[0][:-1]
        sinks = points[flow[points] == size]

        # determine starts and stops
        starts = np.unique(np.concatenate([sources, merges]))
        stops = set(np.concatenate([merges, sinks]).tolist())  # native set

        # travel them and yield per section
        for x in starts:
            if x in sinks:
                continue
            l = [x]
            while True:
                x = flow[x]
                l.append(x)
                if x in stops:
                    break
            a = np.array(l)
            yield klass, (a // width + 0.5, a % width + 0.5)


class Vectorizer(object):
    def __init__(self, direction_path, accumulation_path, output_path):
        # paths and source data
        self.direction_dataset = gdal.Open(direction_path)
        self.accumulation_dataset = gdal.Open(accumulation_path)
        self.output_path = output_path

        # geospatial reference
        geo_transform = self.direction_dataset.GetGeoTransform()
        self.geo_transform = utils.GeoTransform(geo_transform)
        self.projection = self.direction_dataset.GetProjection()

    def vectorize(self, index_feature):
        # target path
        name = index_feature[str('bladnr')]
        path = os.path.join(self.output_path, name[:3], '{}'.format(name))
        if os.path.exists(path):
            return

        # create directory
        try:
            os.makedirs(os.path.dirname(path))
        except OSError:
            pass  # no problem

        index_geometry = index_feature.geometry().Buffer(0)
        geo_transform = self.geo_transform.shifted(index_geometry)

        # data
        window = self.geo_transform.get_window(index_geometry)
        direction = self.direction_dataset.ReadAsArray(**window)
        accumulation = self.accumulation_dataset.ReadAsArray(**window)

        # processing
        data_source = SHAPE.CreateDataSource(str(path))
        layer_sr = osr.SpatialReference(self.projection)
        layer_name = str(os.path.basename(path))
        layer = data_source.CreateLayer(layer_name, layer_sr)
        layer.CreateField(ogr.FieldDefn(str('class'), ogr.OFTInteger))
        layer_defn = layer.GetLayerDefn()
        generator = vectorize(direction=direction, accumulation=accumulation)
        for klass, indices in generator:
            feature = ogr.Feature(layer_defn)
            points = geo_transform.get_coordinates(indices)
            feature[str('class')] = klass
            geometry = ogr.Geometry(ogr.wkbLineString)
            for p in zip(*points):
                geometry.AddPoint_2D(*p)
            feature.SetGeometry(geometry)
            layer.CreateFeature(feature)


def flow_vec(index_path, part, **kwargs):
    """
    """
    # select some or all polygons
    index = utils.PartialDataSource(index_path)
    if part is not None:
        index = index.select(part)

    vectorizer = Vectorizer(**kwargs)

    for feature in index:
        vectorizer.vectorize(feature)
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
        'direction_path',
        metavar='DIRECTION',
        help='GDAL direction raster dataset',
    )
    parser.add_argument(
        'accumulation_path',
        metavar='ACCUMULATION',
        help='GDAL accumulation raster dataset',
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
    return parser


def main():
    """ Call aggregate with args from parser. """
    kwargs = vars(get_parser().parse_args())
    flow_vec(**kwargs)