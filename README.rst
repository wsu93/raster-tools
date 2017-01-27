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


resolving multiple entries point to same sparsevoid

Objects
~~~~~~~
TileManager
    - read data
    - label data
    - create sparse void pieces using grey_dilation, register the voids
    - per void register links
    - register dangling (that is, edge) voids in the dangledict as a list of
      absolute indices

Tile
    - data
    - label
    - offset
    - origin

Void
    - tile, label list, so in each tile know which labels to materialize 
    - indices list
    - value list
    - aggregated indices and values lists
    - resolve dictionary with links to adjacent tiles and indices
    - unknown: dictionary with connected tiles

        pop backward link to us
        merge (auto removes 2-way links, adds tiles, registers)

    def merge(self, void):
        update
            links
            indices
            data
            counter
            there may be data or voids on the other side
            void points can be merged
            mutual data must be added, but there will be an edge that is not
            to be used twice
            

    def resolve():
        while links:
            find a linked void from linked tile

    def rasterize(self, array, origin):
        from points selectively paste and smooth at increasing resolutions
        decrease tile counter

Another idea. Process data in tiles. Keep a connection table relating
voids per tile. How to not miss over the edge data? Well, take 1px buffer
around each tile,

labelrelator



Procedure
~~~~~~~~~
    voids dictionary linking tile,void to sparsevoid objects


    pick a tile. now all voids are in void dictionary.
    find objects.
    per object:
        get void from dictionary
        void.resolve()
        void.rasterize()
        remove if void rasterize counter done.
        transfer using the label and the object into the result array


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
