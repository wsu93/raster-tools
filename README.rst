raster-tools
==========================================

Dependencies
------------
python-gdal
python-psycopg2
python-scipy


Rasterizing landuse tables
--------------------------
For rasterization of landuse tables from a postgres datasource a special
wrapper command is available at bin/rasterize-landuse, use --help for args.


Creating a seamless large-scale void-filled raster (deprecated)
---------------------------------------------------------------
1. Make index with power-of-two dimensions for some extent [reindex] 
2. Aggregate base dataset to that index for 32x32 pixels [aggregate]
3. Merge result of 2. into single raster [gdal_merge.py]
4. Fill result of 3. with fillnodata algorithm [fillnodata]
5. Combine with result of 3. to single filled dem [gdalwarp]
6. Fill base data using result of 5. as
   'ceiling' and result of 1. as index [mfillnodata]
7. Cut result back into desired tiling using [retile]

Procedure for filling completely filling internal nodata
--------------------------------------------------------
Rebuild the current file to not use an internal index of optional tiling size
Add the batch speedup to the extractor (yield multiple non interfering voids)
Create SparseEdge object:
    - incremental growth
    - mergeable
Create SparseEdgeAggregation
    - aggregate on init from sparse edge
    - rasterize method

# relate current tile edges to previous tile edges how? A hash in a dictionary
of the coordinates?

process a tile
swiftly smooth internals - write the tile to disk
create edgevoid sparse objects for all exterior voids
on next related tile, see if void can be closed from data
otherwise, it can be merged.
remember voids by tile link
remember voids by tile

2 big questions remain.
sparse merge how?
    - missing data points from adjacent tile can be completed
    - now voids can be merged from still missing points:
    - loop current tile missing edge pixels => pop corresponding voids from
    - something with labels and an & operation.
aggregated rasterize how? slice from sparse arrays!




Define edge voids as unresolvable
Per void keep track of
- edge data
- yet unknown edges to compare to upcoming tiles
- tiles related to this void

Expanded void fills must be properly aligned to tile, how? Easy, the origin
will stay the same and can be related to the rectangular grid


Creating streamlines
--------------------

flow-fil index raster cover output/f
flow-dir index output/f/all.vrt cover output/d              # fill depressions
flow-acc index output/d/all.vrt output/a                    # accumulate
flow-vec index output/d/all.vrt output/a/all.vrt output/v   # makes features
flow-rst index output/v output/r                            # features to tifs


Post-nensskel setup TODO
------------------------

Here are some instructions on what to do after you've created the project with
nensskel.

- Add a new jenkins job at
  http://buildbot.lizardsystem.nl/jenkins/view/djangoapps/newJob or
  http://buildbot.lizardsystem.nl/jenkins/view/libraries/newJob . Job name
  should be "raster-tools", make the project a copy of the existing "lizard-wms"
  project (for django apps) or "nensskel" (for libraries). On the next page,
  change the "github project" to ``https://github.com/nens/raster-tools/`` and
  "repository url" fields to ``git@github.com:nens/raster-tools.git`` (you might
  need to replace "nens" with "lizardsystem"). The rest of the settings should
  be OK.

- The project is prepared to be translated with Lizard's
  `Transifex <http://translations.lizard.net/>`_ server. For details about
  pushing translation files to and fetching translation files from the
  Transifex server, see the ``nens/translations`` `documentation
  <https://github.com/nens/translations/blob/master/README.rst>`_.
