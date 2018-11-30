# -*- coding: utf-8 -*-
import sys

try:
    from osgeo import gdal
    from osgeo import osr
except:
    import gdal

    print 'You are using "old gen" bindings. gdal2tiles needs "new gen" bindings.'
    sys.exit(1)
import os

try:
    from PIL import Image
    import numpy
    import osgeo.gdal_array as gdalarray
except:
    pass
__version__ = '$Id: gdal2tiles.py 27349 2014-05-16 18:58:51Z rouault $'

resampling_list = (
    'average',
    'near',
    'bilinear',
    'cubic',
    'cubicspline',
    'lanczos',
    'antialias',
)
profile_list = ('mercator', 'geodetic', 'raster')  # ,'zoomify')
webviewer_list = ('all', 'google', 'openlayers', 'none')

import math
import time
import numpy as np
from array import array

MAXZOOMLEVEL = 32


class GlobalMercator(object):
    """
    TMS Global Mercator Profile
    ---------------------------

  Functions necessary for generation of tiles in Spherical Mercator projection,
  EPSG:900913 (EPSG:gOOglE, Google Maps Global Mercator), EPSG:3785, OSGEO:41001.

  Such tiles are compatible with Google Maps, Bing Maps, Yahoo Maps,
  UK Ordnance Survey OpenSpace API, ...
  and you can overlay them on top of base maps of those web mapping applications.

    Pixel and tile coordinates are in TMS notation (origin [0,0] in bottom-left).

    What coordinate conversions do we need for TMS Global Mercator tiles::

         LatLon      <->       Meters      <->     Pixels    <->       Tile

     WGS84 coordinates   Spherical Mercator  Pixels in pyramid  Tiles in pyramid
         lat/lon            XY in metres     XY pixels Z zoom      XYZ from TMS
        EPSG:4326           EPSG:900913
         .----.              ---------               --                TMS
        /      \     <->     |       |     <->     /----/    <->      Google
        \      /             |       |           /--------/          QuadTree
         -----               ---------         /------------/
       KML, public         WebMapService         Web Clients      TileMapService

    What is the coordinate extent of Earth in EPSG:900913?

      [-20037508.342789244, -20037508.342789244, 20037508.342789244, 20037508.342789244]
      Constant 20037508.342789244 comes from the circumference of the Earth in meters,
      which is 40 thousand kilometers, the coordinate origin is in the middle of extent.
      In fact you can calculate the constant as: 2 * math.pi * 6378137 / 2.0
      $ echo 180 85 | gdaltransform -s_srs EPSG:4326 -t_srs EPSG:900913
      Polar areas with abs(latitude) bigger then 85.05112878 are clipped off.

    What are zoom level constants (pixels/meter) for pyramid with EPSG:900913?

      whole region is on top of pyramid (zoom=0) covered by 256x256 pixels tile,
      every lower zoom level resolution is always divided by two
      initialResolution = 20037508.342789244 * 2 / 256 = 156543.03392804062

    What is the difference between TMS and Google Maps/QuadTree tile name convention?

      The tile raster itself is the same (equal extent, projection, pixel size),
      there is just different identification of the same raster tile.
      Tiles in TMS are counted from [0,0] in the bottom-left corner, id is XYZ.
      Google placed the origin [0,0] to the top-left corner, reference is XYZ.
      Microsoft is referencing tiles by a QuadTree name, defined on the website:
      http://msdn2.microsoft.com/en-us/library/bb259689.aspx

    The lat/lon coordinates are using WGS84 datum, yeh?

      Yes, all lat/lon we are mentioning should use WGS84 Geodetic Datum.
      Well, the web clients like Google Maps are projecting those coordinates by
      Spherical Mercator, so in fact lat/lon coordinates on sphere are treated as if
      the were on the WGS84 ellipsoid.

      From MSDN documentation:
      To simplify the calculations, we use the spherical form of projection, not
      the ellipsoidal form. Since the projection is used only for map display,
      and not for displaying numeric coordinates, we don't need the extra precision
      of an ellipsoidal projection. The spherical projection causes approximately
      0.33 percent scale distortion in the Y direction, which is not visually noticable.

    How do I create a raster in EPSG:900913 and convert coordinates with PROJ.4?

      You can use standard GIS tools like gdalwarp, cs2cs or gdaltransform.
      All of the tools supports -t_srs 'epsg:900913'.

      For other GIS programs check the exact definition of the projection:
      More info at http://spatialreference.org/ref/user/google-projection/
      The same projection is degined as EPSG:3785. WKT definition is in the official
      EPSG database.

      Proj4 Text:
        +proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0
        +k=1.0 +units=m +nadgrids=@null +no_defs

      Human readable WKT format of EPGS:900913:
         PROJCS["Google Maps Global Mercator",
             GEOGCS["WGS 84",
                 DATUM["WGS_1984",
                     SPHEROID["WGS 84",6378137,298.257223563,
                         AUTHORITY["EPSG","7030"]],
                     AUTHORITY["EPSG","6326"]],
                 PRIMEM["Greenwich",0],
                 UNIT["degree",0.0174532925199433],
                 AUTHORITY["EPSG","4326"]],
             PROJECTION["Mercator_1SP"],
             PARAMETER["central_meridian",0],
             PARAMETER["scale_factor",1],
             PARAMETER["false_easting",0],
             PARAMETER["false_northing",0],
             UNIT["metre",1,
                 AUTHORITY["EPSG","9001"]]]
    """

    def __init__(self, tileSize=256):
        '''Initialize the TMS Global Mercator pyramid'''

        self.tileSize = tileSize
        self.initialResolution = 2 * math.pi * 6378137 / self.tileSize

        # 156543.03392804062 for tileSize 256 pixels

        self.originShift = 2 * math.pi * 6378137 / 2.0

        # 20037508.342789244

    def LatLonToMeters(self, lat, lon):
        '''Converts given lat/lon in WGS84 Datum to XY in Spherical Mercator EPSG:900913'''

        mx = lon * self.originShift / 180.0
        my = math.log(math.tan((90 + lat) * math.pi / 360.0)) \
             / (math.pi / 180.0)

        my = my * self.originShift / 180.0
        return (mx, my)

    def MetersToLatLon(self, mx, my):
        '''Converts XY point from Spherical Mercator EPSG:900913 to lat/lon in WGS84 Datum'''

        lon = mx / self.originShift * 180.0
        lat = my / self.originShift * 180.0

        lat = 180 / math.pi * (2 * math.atan(math.exp(lat * math.pi
                                                      / 180.0)) - math.pi / 2.0)
        return (lat, lon)

    def PixelsToMeters(
            self,
            px,
            py,
            zoom,
    ):
        '''Converts pixel coordinates in given zoom level of pyramid to EPSG:900913'''

        res = self.Resolution(zoom)
        mx = px * res - self.originShift
        my = py * res - self.originShift
        return (mx, my)

    def MetersToPixels(
            self,
            mx,
            my,
            zoom,
    ):
        '''Converts EPSG:900913 to pyramid pixel coordinates in given zoom level'''

        res = self.Resolution(zoom)
        px = (mx + self.originShift) / res
        py = (my + self.originShift) / res
        return (px, py)

    def PixelsToTile(self, px, py):
        '''Returns a tile covering region in given pixel coordinates'''

        tx = int(math.ceil(px / float(self.tileSize)) - 1)
        ty = int(math.ceil(py / float(self.tileSize)) - 1)
        return (tx, ty)

    def PixelsToRaster(
            self,
            px,
            py,
            zoom,
    ):
        '''Move the origin of pixel coordinates to top-left corner'''

        mapSize = self.tileSize << zoom
        return (px, mapSize - py)

    def MetersToTile(
            self,
            mx,
            my,
            zoom,
    ):
        '''Returns tile for given mercator coordinates'''

        (px, py) = self.MetersToPixels(mx, my, zoom)
        return self.PixelsToTile(px, py)

    def TileBounds(
            self,
            tx,
            ty,
            zoom,
    ):
        '''Returns bounds of the given tile in EPSG:900913 coordinates'''

        (minx, miny) = self.PixelsToMeters(tx * self.tileSize, ty
                                           * self.tileSize, zoom)
        (maxx, maxy) = self.PixelsToMeters((tx + 1) * self.tileSize,
                                           (ty + 1) * self.tileSize, zoom)
        return (minx, miny, maxx, maxy)

    def TileLatLonBounds(
            self,
            tx,
            ty,
            zoom,
    ):
        '''Returns bounds of the given tile in latutude/longitude using WGS84 datum'''

        bounds = self.TileBounds(tx, ty, zoom)
        (minLat, minLon) = self.MetersToLatLon(bounds[0], bounds[1])
        (maxLat, maxLon) = self.MetersToLatLon(bounds[2], bounds[3])

        return (minLat, minLon, maxLat, maxLon)

    def Resolution(self, zoom):
        '''Resolution (meters/pixel) for given zoom level (measured at Equator)'''

        # return (2 * math.pi * 6378137) / (self.tileSize * 2**zoom)

        return self.initialResolution / 2 ** zoom

    def ZoomForPixelSize(self, pixelSize):
        '''Maximal scaledown zoom of the pyramid closest to the pixelSize.'''

        for i in range(MAXZOOMLEVEL):
            if pixelSize > self.Resolution(i):
                if i != 0:
                    return i - 1
                else:
                    return 0  # We don't want to scale up

    def GoogleTile(
            self,
            tx,
            ty,
            zoom,
    ):
        '''Converts TMS tile coordinates to Google Tile coordinates'''

        # coordinate origin is moved from bottom-left to top-left corner of the extent

        return (tx, 2 ** zoom - 1 - ty)

    def QuadTree(
            self,
            tx,
            ty,
            zoom,
    ):
        '''Converts TMS tile coordinates to Microsoft QuadTree'''

        quadKey = ''
        ty = 2 ** zoom - 1 - ty
        for i in range(zoom, 0, -1):
            digit = 0
            mask = 1 << i - 1
            if tx & mask != 0:
                digit += 1
            if ty & mask != 0:
                digit += 2
            quadKey += str(digit)

        return quadKey


# ---------------------

class GDAL2Tiles(object):
    def process(self):
        """The main processing function, runs all the main steps of processing"""

        # Opening and preprocessing of the input file

        self.open_input()
        self.generate_metadata()
        self.generate_base_tiles()

    def open_input(self):
        """Initialization of the input raster, reprojection if necessary"""

        gdal.AllRegister()

        # Initialize necessary GDAL drivers
        self.out_drv = gdal.GetDriverByName(self.tiledriver)
        self.mem_drv = gdal.GetDriverByName('MEM')

        if not self.out_drv:
            raise Exception("The '%s' driver was not found, is it available in this GDAL build?"
                            , self.tiledriver)
        if not self.mem_drv:
            raise Exception("The 'MEM' driver was not found, is it available in this GDAL build?"
                            )
        # Open the input file
        if self.input:
            self.in_ds = gdal.Open(self.input, gdal.GA_ReadOnly)
            self.proj = self.in_ds.GetProjection()
        else:
            raise Exception('No input file was specified')
        if not self.in_ds:
            # Note: GDAL prints the ERROR message too

            self.error("It is not possible to open the input file '%s'."
                       % self.input)
        # Read metadata from the input file

        if self.in_ds.RasterCount == 0:
            self.error("Input file '%s' has no raster band"
                       % self.input)
        # Get NODATA value

        self.in_nodata = []
        for i in range(1, self.in_ds.RasterCount + 1):
            if self.in_ds.GetRasterBand(i).GetNoDataValue() != None:
                self.in_nodata.append(self.in_ds.GetRasterBand(i).GetNoDataValue())
        # Spatial Reference System of the input raster

        self.in_srs = None
        if self.options.s_srs:
            self.in_srs = osr.SpatialReference()
            self.in_srs.SetFromUserInput(self.options.s_srs)
            self.in_srs_wkt = self.in_srs.ExportToWkt()
        else:
            self.in_srs_wkt = self.in_ds.GetProjection()
            if not self.in_srs_wkt and self.in_ds.GetGCPCount() != 0:
                self.in_srs_wkt = self.in_ds.GetGCPProjection()
            if self.in_srs_wkt:
                self.in_srs = osr.SpatialReference()
                self.in_srs.ImportFromWkt(self.in_srs_wkt)
        # Spatial Reference System of tiles

        self.out_srs = osr.SpatialReference()

        if self.options.profile == 'mercator':
            self.out_srs.ImportFromEPSG(900913)
        elif self.options.profile == 'geodetic':
            self.out_srs.ImportFromEPSG(4326)
        else:
            self.out_srs = self.in_srs

        self.out_ds = gdal.AutoCreateWarpedVRT(self.in_ds,
                                               self.in_srs_wkt, self.out_srs.ExportToWkt())
        # Output Bounds - coordinates in the output SRS
        self.out_gt = self.out_ds.GetGeoTransform()

        self.ominx = self.out_gt[0]
        self.omaxx = self.out_gt[0] + self.out_ds.RasterXSize \
                                      * self.out_gt[1]
        self.omaxy = self.out_gt[3]
        self.ominy = self.out_gt[3] - self.out_ds.RasterYSize \
                                      * self.out_gt[1]

        if self.options.profile == 'mercator':
            self.mercator = GlobalMercator()  # from globalmaptiles.py

            # Function which generates SWNE in LatLong for given tile

            self.tileswne = self.mercator.TileLatLonBounds
            # Generate table with min max tile coordinates for all zoomlevels
            self.tminmax = list(range(0, 32))
            for tz in range(0, 32):
                (tminx, tminy) = self.mercator.MetersToTile(self.ominx,
                                                            self.ominy, tz)
                (tmaxx, tmaxy) = self.mercator.MetersToTile(self.omaxx,
                                                            self.omaxy, tz)

                # crop tiles extending world limits (+-180,+-90)

                (tminx, tminy) = (max(0, tminx), max(0, tminy))
                (tmaxx, tmaxy) = (min(2 ** tz - 1, tmaxx), min(2 ** tz - 1, tmaxy))
                self.tminmax[tz] = (tminx, tminy, tmaxx, tmaxy)
            # Get the minimal zoom level (map covers area equivalent to one tile)

            if self.tminz == None:
                self.tminz = \
                    self.mercator.ZoomForPixelSize(self.out_gt[1]
                                                   * max(self.out_ds.RasterXSize,
                                                         self.out_ds.RasterYSize) / float(self.tilesize))

                # Get the maximal zoom level (closest possible zoom level up on the resolution of raster)

            if self.tmaxz == None:
                self.tmaxz = \
                    self.mercator.ZoomForPixelSize(self.out_gt[1])

    def generate_metadata(self):
        """Generation of main metadata files and HTML viewers (metadata related to particular tiles are generated during the tile processing)."""
        if not os.path.exists(self.output):
            os.makedirs(self.output)
        if self.options.profile == 'mercator':
            (south, west) = self.mercator.MetersToLatLon(self.ominx,
                                                         self.ominy)
            (north, east) = self.mercator.MetersToLatLon(self.omaxx,
                                                         self.omaxy)
            (south, west) = (max(-85.05112878, south), max(-180.0,
                                                           west))
            (north, east) = (min(85.05112878, north), min(180.0, east))
            self.swne = (south, west, north, east)

    def generate_base_tiles(self):
        """Generation of the base tiles (the lowest in the pyramid) directly from the input raster"""

        print 'Generating Base Tiles:'

        if self.options.verbose:
            # mx, my = self.out_gt[0], self.out_gt[3] # OriginX, OriginY
            # px, py = self.mercator.MetersToPixels( mx, my, self.tmaxz)
            # print "Pixel coordinates:", px, py, (mx, my)

            print ''
            print 'Tiles generated from the max zoom level:'
            print '----------------------------------------'
            print ''

            # Set the bounds
        for mtz in xrange(self.tmaxz + 1):
            self.tmaxz = mtz
            if (mtz >= 7):
                self.tilesize = 128
                self.querysize = self.tilesize
            (tminx, tminy, tmaxx, tmaxy) = self.tminmax[self.tmaxz]

            querysize = self.querysize
            ds = self.out_ds

            tcount = (1 + abs(tmaxx - tminx)) * (1 + abs(tmaxy - tminy))
            ti = 0
            tz = self.tmaxz
            yrange = range(tmaxy, tminy - 1, -1)
            for ty in yrange:
                for tx in range(tminx, tmaxx + 1):
                    if self.stopped:
                        break
                    ti += 1
                    tilefilename = os.path.join(self.output, str(tz), str(tx),
                                                '%s.%s' % (pow(2, tz) - 1 - ty, self.tileext))
                    if self.options.verbose:
                        print (ti, '/', tcount, tilefilename)  # , "( TileMapService: z / x / y )"

                    if self.options.resume and os.path.exists(tilefilename):
                        if self.options.verbose:
                            print 'Tile generation skiped because of --resume'
                        else:
                            self.progressbar(ti / float(tcount))
                        continue

                    # Create directories for the tile
                    if not os.path.exists(os.path.dirname(tilefilename)):
                        os.makedirs(os.path.dirname(tilefilename))

                    if self.options.profile == 'mercator':
                        # Tile bounds in EPSG:900913

                        b = self.mercator.TileBounds(tx, ty, tz)
                        minLat, minLon = self.mercator.MetersToLatLon(b[0], b[1])
                        maxLat, maxLon = self.mercator.MetersToLatLon(b[2], b[3])
                        # Don't scale up by nearest neighbour, better change the querysize
                        # to the native resolution (and return smaller query tile) for scaling

                    if self.options.profile in ('mercator', 'geodetic'):
                        (rb, wb) = self.geo_query(ds, b[0], b[3], b[2],
                                                  b[1])
                        nativesize = wb[0] + wb[2]  # Pixel size in the raster covering query geo extent
                        if self.options.verbose:
                            print ('\tNative Extent (querysize',
                                   nativesize, '): ', rb, wb)
                        # Tile bounds in raster coordinates for ReadRaster query
                        (rb, wb) = self.geo_query(
                            ds,
                            b[0],
                            b[3],
                            b[2],
                            b[1],
                            querysize=querysize,
                        )
                        (rx, ry, rxsize, rysize) = rb
                        (wx, wy, wxsize, wysize) = wb
                        # -------------------------------------------------------------------------
                    data = ds.ReadRaster(
                        rx,
                        ry,
                        rxsize,
                        rysize,
                        wxsize,
                        wysize,
                        band_list=[1],
                    )

                    mydata = np.asarray(array("f", data))
                    newMydata = np.reshape(mydata, (wysize, wxsize))
                    driver = gdal.GetDriverByName("GTiff");  # 数据类型必须有，因为要计算需要多大内存空间
                    dataset = driver.Create(tilefilename, querysize, querysize, 1, gdal.GDT_Float32)
                    self.fill_value = np.finfo(np.float32).min
                    noData = np.full(shape=(querysize, querysize), fill_value=self.fill_value)
                    dataset.GetRasterBand(1).WriteArray(noData)  # 写入数组数据
                    dataset.GetRasterBand(1).WriteArray(newMydata, wx, wy)  # 写入数组数据
                    dataset.GetRasterBand(1).SetNoDataValue(np.double(self.fill_value))
                    (minLat, minLon, maxLat, maxLon) = self.mercator.TileLatLonBounds(tx, ty, tz)

                    geotrans = (minLon, (maxLon - minLon) / querysize, 0.0, maxLat, 0.0, -(maxLat - minLat) / querysize)
                    dataset.SetGeoTransform(geotrans)
                    dataset.SetProjection(self.proj)  # 写入投影
                    del dataset
                    del mydata
                    del newMydata
                    # dstile.WriteRaster(wx, wy, wxsize, wysize, data, band_list=[1], dtype="float32")
                    # print "max", np.max(dstile.ReadAsArray(0, 0, dstile.RasterYSize, dstile.RasterYSize))

    def geo_query(self, ds, ulx, uly, lrx, lry, querysize=0):
        """For given dataset and query in cartographic coordinates
        returns parameters for ReadRaster() in raster coordinates and
        x/y shifts (for border tiles). If the querysize is not given, the
        extent is returned in the native resolution of dataset ds."""
        geotran = ds.GetGeoTransform()
        rx = int((ulx - geotran[0]) / geotran[1] + 0.001)
        ry = int((uly - geotran[3]) / geotran[5] + 0.001)
        rxsize = int((lrx - ulx) / geotran[1] + 0.5)
        rysize = int((lry - uly) / geotran[5] + 0.5)
        if not querysize:
            (wxsize, wysize) = (rxsize, rysize)
        else:
            (wxsize, wysize) = (querysize, querysize)

        # Coordinates should not go out of the bounds of the raster

        wx = 0
        if rx < 0:
            rxshift = abs(rx)
            wx = int(wxsize * (float(rxshift) / rxsize))
            wxsize = wxsize - wx
            rxsize = rxsize - int(rxsize * (float(rxshift) / rxsize))
            rx = 0
        if rx + rxsize > ds.RasterXSize:
            wxsize = int(wxsize * (float(ds.RasterXSize - rx) / rxsize))
            rxsize = ds.RasterXSize - rx

        wy = 0
        if ry < 0:
            ryshift = abs(ry)
            wy = int(wysize * (float(ryshift) / rysize))
            wysize = wysize - wy
            rysize = rysize - int(rysize * (float(ryshift) / rysize))
            ry = 0
        if ry + rysize > ds.RasterYSize:
            wysize = int(wysize * (float(ds.RasterYSize - ry) / rysize))
            rysize = ds.RasterYSize - ry

        return ((rx, ry, rxsize, rysize), (wx, wy, wxsize, wysize))

        # -------------------------------------------------------------------------

    def __init__(self, arguments):
        self.input = None
        self.output = None
        self.stopped = False
        # Tile format

        self.tilesize = 256
        self.tiledriver = 'GTiff'
        self.tileext = 'tif'
        self.querysize = self.tilesize

        self.optparse_init()
        (self.options, self.args) = self.parser.parse_args(args=arguments)
        # print self.args[-1]
        # print self.args[-2]
        if os.path.isdir(self.args[-1]) or len(self.args) > 1 and not os.path.exists(self.args[-1]):
            self.output = self.args[-1]
            self.args = self.args[:-1]
        if not self.args:
            self.error('No input file specified')

        self.input = self.args[0]
        self.tminz = None
        self.tmaxz = None
        if self.options.zoom:
            minmax = self.options.zoom.split('-', 1)
            minmax.extend([''])
            (min, max) = minmax[:2]
            self.tminz = int(min)
            if max:
                self.tmaxz = int(max)
            else:
                self.tmaxz = int(min)

    # -------------------------------------------------------------------------

    def optparse_init(self):
        """Prepare the option parser for input (argv)"""

        from optparse import OptionParser, OptionGroup
        usage = 'Usage: %prog [options] input_file(s) [output]'
        p = OptionParser(usage, version='%prog ' + __version__)
        p.add_option(
            '-p',
            '--profile',
            dest='profile',
            type='choice',
            choices=profile_list,
            help="Tile cutting profile (%s) - default 'mercator' (Google Maps compatible)"
                 % ','.join(profile_list),
        )
        p.add_option(
            '-r',
            '--resampling',
            dest='resampling',
            type='choice',
            choices=resampling_list,
            help="Resampling method (%s) - default 'average'"
                 % ','.join(resampling_list),
        )
        p.add_option('-s', '--s_srs', dest='s_srs', metavar='SRS',
                     help='The spatial reference system used for the source input data'
                     )
        p.add_option('-z', '--zoom', dest='zoom',
                     help="Zoom levels to render (format:'2-5' or '10')."
                     )
        p.add_option('-e', '--resume', dest='resume',
                     action='store_true',
                     help='Resume mode. Generate only missing files.')
        p.add_option('-a', '--srcnodata', dest='srcnodata',
                     metavar='NODATA',
                     help='NODATA transparency value to assign to the input data'
                     )
        p.add_option('-d', '--tmscompatible', dest='tmscompatible',
                     action='store_true',
                     help='When using the geodetic profile, specifies the base resolution as 0.703125 or 2 tiles at zoom level 0.'
                     )
        p.add_option('-l', '--leaflet', action='store_true',
                     dest='leaflet',
                     help="Set 0,0 point to north. For use with 'leaflet'. Requires -p raster. "
                     )
        p.add_option('-v', '--verbose', action='store_true',
                     dest='verbose',
                     help='Print status messages to stdout')

        # KML options

        g = OptionGroup(p, 'KML (Google Earth) options',
                        'Options for generated Google Earth SuperOverlay metadata'
                        )
        g.add_option('-k', '--force-kml', dest='kml',
                     action='store_true',
                     help="Generate KML for Google Earth - default for 'geodetic' profile and 'raster' in EPSG:4326. For a dataset with different projection use with caution!"
                     )
        g.add_option('-n', '--no-kml', dest='kml', action='store_false'
                     ,
                     help='Avoid automatic generation of KML files for EPSG:4326'
                     )
        g.add_option('-u', '--url', dest='url',
                     help='URL address where the generated tiles are going to be published'
                     )
        p.add_option_group(g)

        # HTML options

        g = OptionGroup(p, 'Web viewer options',
                        'Options for generated HTML viewers a la Google Maps'
                        )
        g.add_option(
            '-w',
            '--webviewer',
            dest='webviewer',
            type='choice',
            choices=webviewer_list,
            help="Web viewer to generate (%s) - default 'all'"
                 % ','.join(webviewer_list),
        )
        g.add_option('-t', '--title', dest='title',
                     help='Title of the map')
        g.add_option('-c', '--copyright', dest='copyright',
                     help='Copyright for the map')
        g.add_option('-g', '--googlekey', dest='googlekey',
                     help='Google Maps API key from http://code.google.com/apis/maps/signup.html'
                     )

        (g.add_option('-b', '--bingkey', dest='bingkey',
                      help='Bing Maps API key from https://www.bingmapsportal.com/'
                      ),)
        p.add_option_group(g)

        # TODO: MapFile + TileIndexes per zoom level for efficient MapServer WMS
        # g = OptionGroup(p, "WMS MapServer metadata", "Options for generated mapfile and tileindexes for MapServer")
        # g.add_option("-i", "--tileindex", dest='wms', action="store_true"
        #                 help="Generate tileindex and mapfile for MapServer (WMS)")
        # p.add_option_group(g)

        p.set_defaults(
            verbose=False,
            profile='mercator',
            kml=False,
            url='',
            webviewer='all',
            copyright='',
            resampling='average',
            resume=False,
            googlekey='INSERT_YOUR_KEY_HERE',
            bingkey='INSERT_YOUR_KEY_HERE',
        )

        self.parser = p

        # -------------------------------------------------------------------------


        # 确定最好等级的区域


if __name__ == '__main__':
    startTime = time.time()
    argv = gdal.GeneralCmdLineProcessor(sys.argv)
    if argv:
        gdal2tiles = GDAL2Tiles(argv[1:])
        gdal2tiles.process()
    print "end time ", time.time() - startTime
