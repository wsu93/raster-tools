from setuptools import setup

version = '0.2.dev0'

long_description = '\n\n'.join([
    open('README.rst').read(),
    open('CREDITS.rst').read(),
    open('CHANGES.rst').read(),
    ])

install_requires = [
    'gdal',
    'h5py',
    'psycopg2',
    'scipy',
    'setuptools',
    'unipath',
    'numpy >= 1.8.2',
    ],

tests_require = [
    'nose',
    'coverage',
    ]

setup(name='raster-tools',
      version=version,
      description="Tools for processing of elevation and landuse raster data.",
      long_description=long_description,
      # Get strings from http://www.python.org/pypi?%3Aaction=list_classifiers
      classifiers=[],
      keywords=[],
      author='Arjan Verkerk',
      author_email='arjan.verkerk@nelen-schuurmans.nl',
      url='',
      license='GPL',
      packages=['raster_tools'],
      zip_safe=False,
      install_requires=install_requires,
      tests_require=tests_require,
      extras_require={'test': tests_require},
      entry_points={
          'console_scripts': [
              # ahn2
              'ahn2-aig2tif       = raster_tools.ahn2.aig2tif:main',
              'ahn2-zip2tif       = raster_tools.ahn2.zip2tif:main',
              'ahn2-constant      = raster_tools.ahn2.constant:main',
              # # ahn3
              'ahn3-download      = raster_tools.ahn3.download:main',
              # srtm
              'srtm-make-index    = raster_tools.srtm.make_index:main',
              'srtm-fix-nodata    = raster_tools.srtm.fix_nodata:main',
              'srtm-organize      = raster_tools.srtm.organize:main',
              # flow analysis
              'flow-fil           = raster_tools.flow.flow_fil:main',
              'flow-dir           = raster_tools.flow.flow_dir:main',
              'flow-acc           = raster_tools.flow.flow_acc:main',
              'flow-vec           = raster_tools.flow.flow_vec:main',
              'flow-rst           = raster_tools.flow.flow_rst:main',
              # void filling
              'fill-label         = raster_tools.fill.label:main',
              # modification
              'mfillnodata        = raster_tools.mfillnodata:main',
              'fillnodata         = raster_tools.fillnodata:main',
              'hillshade          = raster_tools.hillshade:main',
              'aggregate          = raster_tools.aggregate:main',
              'rebase             = raster_tools.rebase:main',
              'shadow             = raster_tools.shadow:main',
              # organization
              'retile             = raster_tools.retile:main',
              'reindex            = raster_tools.reindex:main',
              # rasterization
              'bag2tif            = raster_tools.bag2tif:main',
              'rasterize-landuse  = raster_tools.rasterize_landuse:main',
              'rasterize          = raster_tools.rasterize:main',
              # modeling
              'extract            = raster_tools.extract:main',
              'line-up            = raster_tools.line_up:main',
              'vselect            = raster_tools.vselect:main',
              # pointclouds
              'txt2tif            = raster_tools.txt2tif:main',
              'pol2laz            = raster_tools.pol2laz:main',
              'roof               = raster_tools.roof:main',
              # analysis
              'zonal              = raster_tools.zonal:main',
              'rgb-zonal          = raster_tools.rgb_zonal:main',
              'green-factor       = raster_tools.green_factor:main',
          ]},
      )
