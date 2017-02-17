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
Available: fill-label. Script to make an almost arbitrary size array with void
labels.


creating streamlines
--------------------

flow-fil index raster cover output/f
flow-dir index output/f/all.vrt cover output/d              # fill depressions
flow-acc index output/d/all.vrt output/a                    # accumulate
flow-vec index output/d/all.vrt output/a/all.vrt output/v   # makes features
flow-rst index output/v output/r                            # features to tifs


post-nensskel setup todo
------------------------

here are some instructions on what to do after you've created the project with
nensskel.

- add a new jenkins job at
  http://buildbot.lizardsystem.nl/jenkins/view/djangoapps/newjob or
  http://buildbot.lizardsystem.nl/jenkins/view/libraries/newjob . job name
  should be "raster-tools", make the project a copy of the existing "lizard-wms"
  project (for django apps) or "nensskel" (for libraries). on the next page,
  change the "github project" to ``https://github.com/nens/raster-tools/`` and
  "repository url" fields to ``git@github.com:nens/raster-tools.git`` (you might
  need to replace "nens" with "lizardsystem"). the rest of the settings should
  be ok.

- the project is prepared to be translated with lizard's
  `transifex <http://translations.lizard.net/>`_ server. for details about
  pushing translation files to and fetching translation files from the
  transifex server, see the ``nens/translations`` `documentation
  <https://github.com/nens/translations/blob/master/readme.rst>`_.
