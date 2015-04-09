from __future__ import print_function

# Copyright (c) 2014, Roger Lew (rogerlew.gmail.com)
#
# The project described was supported by NSF award number IIA-1301792
# from the NSF Idaho EPSCoR Program and by the National Science Foundation.

"""
Data Handler for iSNOBAL IPW Files
"""
__version__ = "0.0.1"

from collections import namedtuple
from glob import glob
import os
import time
import warnings

import h5py
import numpy as np
from numpy.testing import assert_array_almost_equal

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, ListedColormap

from osgeo import gdal, gdalconst, ogr, osr

in_db__vars = tuple('I_lw T_a e_a u T_g S_n'.split())
out_em__vars = tuple('R_n H L_v_E G M delta_Q E_s '
                     'melt ro_predict cc_s'.split())
out_snow__vars = tuple('z_s rho m_s h2o T_s_0 T_s_l T_s z_s_l h2o_sat'.split())

_unpackindx  = lambda L: int(L.split()[2])
_unpackint   = lambda L: int(L.split('=')[1].strip())
_unpackfloat = lambda L: float(L.split('=')[1].strip())
_unpackstr   = lambda L: L.split('=')[1].strip()


class Band:
    """
    Represents a raster band of geospatial data
    """
    def __init__(self, nlines, nsamps):

        # Using classes instead of dicts makes things faster
        # even if though they are more verbose to implement

        # width, height
        self.nlines, self.nsamps = nlines, nsamps

        # basic_image
        self.name = None
        self.bytes = None
        self.fmt = None
        self.bits = None
        self.annot = None
        self.history = []

        # geo
        self.bline = None
        self.bsamp = None
        self.dline = None
        self.dsamp = None
        self.geounits = None
        self.coord_sys_ID = None
        self.geotransform = None

        # lq
        self.x0 = None
        self.xend = None
        self.y0 = None
        self.yend = None
        self.units = None
        self.transform = None

        self.data = None

    def _parse_geo(self, L0, L1, L2, L3, L4, L5):
        """
        sets attributes and builds GDAL ordered
        geotransform list
        """
        bline = self.bline = _unpackfloat(L0)
        bsamp = self.bsamp = _unpackfloat(L1)
        dline = self.dline = _unpackfloat(L2)
        dsamp = self.dsamp = _unpackfloat(L3)
        self.geounits = _unpackstr(L4)
        self.coord_sys_ID = _unpackstr(L5)
        self.geotransform = [bsamp - dsamp / 2.0,
                             dsamp, 0.0,
                             bline - dline / 2.0,
                             0.0, dline]

    def _parse_lq(self, L0, L1):
        """
        Pulls the scale values from line0 and line1 builds
        a function for transforming integer values to floats
        """
        [[x0, y0], [xend, yend]] = [L0.split(), L1.split()]
        x0, y0, xend, yend = map(float, [x0, y0, xend, yend])
        self.transform = lambda x: (yend - y0) * (x / xend) + y0
        self.x0, self.xend = x0, xend
        self.y0, self.yend = y0, yend

    def __str__(self):
        return '''\
    nlines (width): {0.nlines}
    nsamps (height): {0.nsamps}
    bytes: {0.bytes}
    transform: {0.x0} -> {0.xend}
               {0.y0} -> {0.yend}
    geo: {0.bline} {0.bsamp} {0.dline} {0.dsamp}
'''.format(self)


class IPW:
    """
    Represents a IPW file container
    """
    def __init__(self, fname, rescale=True, epsg=32611):
        """
        IPW(fname[, rescale=True])

        Parameters
        ----------
        fname : string
           path to the IPW container
           if "in.", "em.", or "snow." are in fname
           band names will be pulled from the cooresponding
           global varlist:
               in_db__vars
               out_em__vars
               out_snow__vars

        rescale : bool (default = True)
            Specifies whether data should be rescaled to singles
            based on the lq map or whether they should remain uint8/16

        epsg : int (default = 32611)
             North America Centric Cheat Sheet
            -----------------------------------------------
            UTM Zone 10 Northern Hemisphere (WGS 84)  32610
            UTM Zone 11 Northern Hemisphere (WGS 84)  32611
            UTM Zone 12 Northern Hemisphere (WGS 84)  32612
            UTM Zone 13 Northern Hemisphere (WGS 84)  32613
            UTM Zone 14 Northern Hemisphere (WGS 84)  32614
            UTM Zone 15 Northern Hemisphere (WGS 84)  32615
            UTM Zone 16 Northern Hemisphere (WGS 84)  32616
            UTM Zone 17 Northern Hemisphere (WGS 84)  32617
            UTM Zone 18 Northern Hemisphere (WGS 84)  32618
            UTM Zone 19 Northern Hemisphere (WGS 84)  32619

        """
        global in_db__vars, out_em__vars, out_snow__vars

        self.epsg = epsg    # this should just be stored as an attribute
                            # it produces alot of book-keeping otherwise

        # read the data to a list of lines
        fid = open(fname, 'rb')
        readline = fid.readline  # dots make things slow (and ugly)
        tell = fid.tell

        bands = []
        byteorder = None
        nlines = None
        nsamps = None
        nbands = None

        # size of the file we are reading in bytes
        st_size = os.fstat(fid.fileno()).st_size
        while 1:  # while 1 is faster than while True
            line = readline()

            # fail safe, haven't needed it, but if this is running
            # on a server not a bad idea to have it here
            if tell() == st_size:
                raise Exception('Unexpectedly reached EOF')

            if '!<header> basic_image_i' in line:
                byteorder = map(int, readline().split('=')[1].strip())
                nlines = _unpackint(readline())
                nsamps = _unpackint(readline())
                nbands = _unpackint(readline())

                # initialize the band instances in the bands list
                bands = [Band(nlines, nsamps) for j in xrange(nbands)]

            elif '!<header> basic_image' in line:
                indx = _unpackindx(line)
                bytes = bands[indx].bytes = _unpackint(readline())
                bands[indx].fmt = ('uint8', 'uint16')[bytes == 2]
                bands[indx].bits = _unpackint(readline())

                while line[0] != '!':
                    bands[indx].history.append(_unpackstr(line))

            elif '!<header> geo' in line:
                indx = _unpackindx(line)
                bands[indx]._parse_geo(*[readline() for i in xrange(6)])

            elif '!<header> lq' in line:
                indx = _unpackindx(line)
                line1 = fid.readline()
                if 'units' in line1:
                    bands[indx].units = _unpackstr(line1)

                    bands[indx]._parse_lq(_unpackstr(readline()),
                                          _unpackstr(readline()))
                else:
                    bands[indx]._parse_lq(_unpackstr(line1),
                                          _unpackstr(readline()))

            if '\f' in line:  # feed form character separates the
                break         # image header from the binary data

        # attempt to assign names to the bands
        assert nbands == len(bands)

        if 'in.' in fname:
            varlist = in_db__vars
        elif 'em.' in fname:
            varlist = out_em__vars
        elif 'snow.' in fname:
            varlist = out_snow__vars
        else:
            varlist = ['band%02i'%i for i in xrange(nbands)]

        assert len(varlist) >= nbands

        for b, name in zip(bands, varlist[:nbands]):
            b.name = name

        # Unpack the binary data using numpy.fromfile
        # because we have been reading line by line fid is at the
        # first data byte, we will read it all as one big chunk
        # and then put it into the appropriate bands
        #
        # np.types allow you to define heterogenous arrays of mixed
        # types and reference them with keys, this helps us out here
        dt = np.dtype([(b.name, b.fmt) for b in bands])

        bip = sum([b.bytes for b in bands])  # bytes-in-pixel
        required_bytes = bip * nlines * nsamps
        assert (st_size - tell()) >= required_bytes

        # this is way faster than looping with struct.unpack
        # struct.unpack also starts assuming there are pad bytes
        # when format strings with different types are supplied
        data = np.fromfile(fid, dt, count=nlines*nsamps)

        # Separate into bands
        data = data.reshape(nlines, nsamps)
        for b in bands:
            if rescale:
                b.data = np.array(b.transform(data[b.name]),
                                  dtype=np.float32)
            else:
                b.data = np.array(data[b.name],
                                  dtype=np.dtype(b.fmt))

        # clean things up
        self.fname = fname
        self.rescale = rescale
        self.name_dict = dict(zip(varlist, range(nbands)))
        self.bands = bands
        self.bip = bip
        self.byteorder = byteorder
        self.nlines = nlines
        self.nsamps = nsamps
        self.nbands = nbands

        fid.close()

    def __getitem__(self, key):
        return self.bands[self.name_dict[key]]

    def colorize(self, dst_fname, band, colormap, ymin=None, ymax=None,
                 drivername='Gtiff'):
        """
        colorize(band[, colormap][, ymin=None][, ymax=None])

        Build a colorized georeferenced raster of a band

        Parameters
        ----------
        band : int or string
            index of band, 1st band is "0"
            string name of raster band e.g. "ro_predicted"

        colormap : string
            name of a matplotlib.colors colormap

        ymin : None or float
            float specifies the min value for normalization
            None will use min value of data

        ymax : None or float
            float specifies the max value for normalization
            None will use max value of data
        """
        bands = self.bands
        nsamps, nlines = self.nsamps, self.nlines

        # find band
        try:
            band = bands[int(band)]
        except:
            band = self[band]

        # build normalize function
        if ymin is None:
            ymin = (0.0, band.y0)[self.rescale]

        if ymax is None:
            ymax = (2.0**band.bits - 1.0, band.yend)[self.rescale]

        assert ymax > ymin

        norm_func = Normalize(ymin, ymax, clip=True)

        # colorize band
        cm = plt.get_cmap(colormap)
        rgba = np.array(cm(norm_func(band.data)) * 255.0, dtype=np.uint8)

        # find geotransform
        for b in bands:
            gt0 = b.geotransform
            if gt0 is not None:
                break

        # initialize raster
        driver = GetDriverByName(drivername)
        ds = driver.Create(dst_fname, nsamps, nlines,
                           4, gdalconst.GDT_Byte)

        # set projection
        if epsg is not None:
            proj = osr.SpatialReference()
            status = proj.ImportFromEPSG(epsg)
            if status != 0:
                warnings.warn('Importing epsg %i return error code %i'
                              % (epsg, status))
            ds.SetProjection(proj.ExportToWkt())

        # set geotransform
        if gt0 is None:
            warnings.warn('Unable to find a geotransform')
        else:
            # set transform
            ds.SetGeoTransform(gt0)

        # write data
        for i in xrange(4):
            ds.GetRasterBand(i+1).WriteArray(rgba[:, :, i])

        ds = None  # Writes and closes file

    def translate(self, dst_fname, writebands=None,
                  drivername='Gtiff', multi=True):
        """
        translate(dst_dataset[, bands=None][, drivername='GTiff']
                  [, multi=True])

        translates the data to a georeferenced tif.

        Parameters
        ----------
        dst_fname : string
           path to destination file without extension.
           Assumes folder exists.

        writebands : None or iterable of integers
            Specifies which bands to write to file.
            Bands are written in the order specifed.

            If none, all bands will be written to file
            The first band is "0" (like iSNOBAL, not like GDAL)

        multi : bool (default True)
            True write each band to its own dataset
            False writes all the bands to a single dataset
        """
        if writebands is None:
            writebands = range(self.nbands)

        if multi:
            for i in writebands:
                self._translate(dst_fname + '.%02i'%i, [i], drivername)
        else:
            self._translate(dst_fname, writebands, drivername)

    def _translate(self, dst_fname, writebands=None, drivername='Gtiff'):
        epsg = self.epsg
        rescale = self.rescale
        bands = self.bands
        nbands = self.nbands
        nlines, nsamps = self.nlines, self.nsamps

        if writebands is None:
            writebands = range(nbands)

        num_wb = len(writebands)

        assert num_wb >= 1

        # The first band of the inputs doesn't have a
        # geotransform. I'm sure this is a feature and not a bug ;)
        #
        # Check to make sure the defined geotransforms are the same
        #
        # Haven't really found any cases where multiple bands have
        # different projections. Is this really a concern?
        gt_override = 0

        # search write bands for valid transform
        for i in writebands:
            gt0 = bands[i].geotransform
            if gt0 is not None:
                break

        if gt0 is None:
            # search all bands for valid transform
            for b in bands:
                gt0 = b.geotransform
                if gt0 is not None:
                    gt_override = 1
                    break
            if gt0 is None:
                raise Exception('No Projection Found')
            else:
                warnings.warn('Using Projection from another band')

        if not gt_override:
            for i in writebands:
                gt = bands[i].geotransform
                if gt is None:
                    continue
                assert_array_almost_equal(gt0, gt)

        # If the data hasn't been rescaled all bands are written as
        # Float 32. If the data has not been scaled the type is
        # Uint8 if all channels are Uint8 and Uint16 otherwise
        if rescale:
            gdal_type = gdalconst.GDT_Float32
        else:
            if all([bands[i].bytes == 1 for i in writebands]):
                gdal_type = gdalconst.GDT_Byte
            else:
                gdal_type = gdalconst.GDT_UInt16

        # initialize raster
        driver = gdal.GetDriverByName(drivername)
        ds = driver.Create(dst_fname + '.tif', nsamps, nlines,
                           num_wb, gdal_type)

        # set projection
        if epsg is not None:
            proj = osr.SpatialReference()
            status = proj.ImportFromEPSG(epsg)
            if status != 0:
                warnings.warn('Importing epsg %i return error code %i'
                              %(epsg, status))
            ds.SetProjection(proj.ExportToWkt())

        # set geotransform
        ds.SetGeoTransform(gt0)

        # write data
        j = 1
        for i in writebands:
            ds.GetRasterBand(j).WriteArray(bands[i].data)
            j += 1

        ds = None  # Writes and closes file

    def __str__(self):
        s = ['''\
IPW({0.fname})
---------------------------------------------------------
byteorder: {0.byteorder}
nlines: {0.nlines}
nsamps: {0.nsamps}
nbands: {0.nbands}
'''.format(self)]

        for i, b in enumerate(self.bands):
            s.append('\nBand %i\n'%i)
            s.append(str(b))

        return ''.join(s)


def _packgrp(root, grp, wc, varlist, nbands=None):
    fns = glob(wc)

    assert len(fns) > 0

    ipw0 = IPW(fns[0])
    shape = (ipw0.nlines * ipw0.nsamps,
             len(fns),
             (nbands, ipw0.nbands)[nbands is None])
    data = np.zeros(shape, np.float32)

    for i, fn in enumerate(fns):
        ipw = IPW(fn)
        for j, b in enumerate(ipw.bands):
            assert varlist[j] == b.name
            data[:, i, j] = b.data.flatten()

    root.create_group(grp)
    for i, key in enumerate(varlist):
        root[grp].create_dataset(key, data=data[:, :, i])


def packToHd5(in_path, out_path=None, fname=None):
    """
    packToHd5(in_path[, out_path][, fname=None])

    Packs input and output data into an hdf5 container

    Parameters
    in_path : string
        path to file containing input IPW files

    out_path : string
        path to file containing output IPW files

    """
    if fname is None:
        fname = 'insnobal_data.hd5'

    if not os.path.isdir(in_path):
        raise Exception('in_path should be a directory')

    if out_path is None:
        out_path = in_path

    if not os.path.isdir(out_path):
        raise Exception('out_path should be a directory')

    root = h5py.File(fname, 'w')

    # Process in_db files
    # Some of the input files have 5 bands, some have 6
    wc = os.path.join(in_path, 'in.*')
    _packgrp(root, 'in_db', wc, in_db__vars, nbands=6)

    # Process out_em files
    wc = os.path.join(out_path, 'em.*')
    _packgrp(root, 'out_em', wc, out_em__vars)

    # Process out_snow files
    wc = os.path.join(out_path, 'snow.*')
    _packgrp(root, 'out_snow', wc, out_snow__vars)

    root.close()
