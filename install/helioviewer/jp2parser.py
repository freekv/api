#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""Helioviewer.org JP2 Download Daemon (HVPull)
JPEG 2000 Image XML Box parser class
"""
import os
import sys
from xml.etree import cElementTree as ET
import numpy as np
from sunpy.io.jp2 import get_header
from sunpy.map import Map
from sunpy.util.xml import xml_to_dict
from sunpy.io.header import FileHeader
from glymur import Jp2k

from glymur.jp2box import UnknownBox, _BOX_WITH_ID
from glymur.codestream import Codestream

from sunpy.util.xml import xml_to_dict
from os import stat
from struct import pack, unpack
import codecs
import io

__HV_CONSTANT_RSUN__ = 959.644
__HV_CONSTANT_AU__ = 149597870700
def parsexml(fptr, offset, length):
    num_bytes = offset + length - fptr.tell()
    read_buffer = fptr.read(num_bytes)

    if sys.hexversion < 0x03000000 and codecs.BOM_UTF8 in read_buffer:
        msg = ('A BOM (byte order marker) was detected and '
               'removed from the XML contents in the box starting at byte '
               'offset {offset:d}.')
        msg = msg.format(offset=offset)
        warnings.warn(msg, UserWarning)
        read_buffer = read_buffer.replace(codecs.BOM_UTF8, b'')

    try:
        text = read_buffer.decode('utf-8')
    except UnicodeDecodeError as err:
        decl_start = read_buffer.find(b'<?xml')
        if decl_start <= -1:
            msg = ('A problem was encountered while parsing an XML box:'
                   '\n\n\t"{error}"\n\nNo XML was retrieved.')
            warnings.warn(msg.format(error=str(err)), UserWarning)
            return XMLBox(xml=None, length=length, offset=offset)

        text = read_buffer[decl_start:].decode('utf-8')

        msg = ('A UnicodeDecodeError was encountered parsing an XML box '
               'at byte position {offset:d} ({reason}), but the XML was '
               'still recovered.')
        msg = msg.format(offset=offset, reason=err.reason)
        warnings.warn(msg, UserWarning)

    text = text.rstrip(chr(0))
    return text

def hv_parse_this_box(fptr, box_id, start, num_bytes):
    try:
        parser = parsexml
    except KeyError:
        # We don't recognize the box ID, so create an UnknownBox and be
        # done with it.
        msg = 'Unrecognized box ({0}) encountered.'.format(box_id)
        warnings.warn(msg)
        return UnknownBox(box_id, offset=start, length=num_bytes)

    try:
        xmltxt = parser(fptr, start, num_bytes)
    except ValueError as err:
        xmltxt = None
        msg = ('Encountered an unrecoverable ValueError while parsing a {0} '
               'box at byte offset {1}.  The original error message was "{2}"')
        msg = msg.format(box_id.decode('utf-8'), start, str(err))
        warnings.warn(msg, UserWarning)
    return xmltxt

def find_xml(fptr, offset, length):

    fptr_read = fptr.read
    fptr_seek = fptr.seek
    fptr_tell = fptr.tell

    superbox = []

    if offset == 0:
        start = 0
    else:
        start = fptr_tell()

    while True:

        # Are we at the end of the superbox?
        if start >= offset + length:
            break

        read_buffer = fptr_read(8)
        if len(read_buffer) < 8:
            msg = 'Extra bytes at end of file ignored.'
            warnings.warn(msg)
            break

        (box_length, box_id) = unpack('>I4s', read_buffer)
        if box_length == 0:
            # The length of the box is presumed to last until the end of
            # the file.  Compute the effective length of the box.
            # num_bytes = os.path.getsize(fptr.name) - fptr.tell() + 8

            # !!! does not work if not top level box, unlikely to occur
            num_bytes = length - start # length - (start + 8) + 8
        elif box_length == 1:
            # The length of the box is in the XL field, a 64-bit value.
            read_buffer = fptr_read(8)
            num_bytes, = unpack('>Q', read_buffer)
        else:
            # The box_length value really is the length of the box!
            num_bytes = box_length
        if box_id == b"xml ":
            return hv_parse_this_box(fptr, box_id, start, num_bytes)

        if box_length == 0:
            # We're done, box lasted until the end of the file.
            break

        # Position to the start of the next box.
        start += num_bytes
        cur_pos = fptr_tell()

        if num_bytes > length:
            # Length of the current box goes past the end of the
            # enclosing superbox.
            msg = '{0} box has incorrect box length ({1})'
            msg = msg.format(box_id, num_bytes)
            warnings.warn(msg)
        elif cur_pos == start:
            # At the start of the next box, jump to it.
            continue
        elif cur_pos > start:
            # The box must be invalid somehow, as the file pointer is
            # positioned past the end of the box.
            msg = ('{0} box may be invalid, the file pointer is positioned '
                   '{1} bytes past the end of the box.')
            msg = msg.format(box_id, cur_pos - start)
            warnings.warn(msg)

        fptr_seek(start)


class JP2parser:
    _filepath = None
    _data = None
            
    def __init__(self, path):
        """Main application"""
        self._filepath = path
        
    def getData(self):
        """Create data object of JPEG 2000 image.
    
        Get image observatory, instrument, detector, measurement, date from image
        metadata and create an object.
        """
            
        imageData = Map(self.read_header_only_but_still_use_sunpy_map())
        image = dict()
        
        #Calculate sun position/size/scale
        dimensions              = self.getImageDimensions();
        refPixel                = self.getRefPixelCoords();
        imageScale              = self.getImagePlateScale();
        dsun                    = self.getDSun();
        layeringOrder           = self.getLayeringOrder();

        # Normalize image scale
        imageScale = imageScale * (dsun / __HV_CONSTANT_AU__);
        
        image['scale'] = imageScale
        image['width'] = dimensions['width']
        image['height'] = dimensions['height']
        image['refPixelX'] = refPixel['x']
        image['refPixelY'] = refPixel['y']
        image['layeringOrder'] = layeringOrder

        image['DSUN_OBS'] = self._data['DSUN_OBS'] if 'DSUN_OBS' in self._data else 'NULL' 
        image['SOLAR_R'] = self._data['SOLAR_R'] if 'SOLAR_R' in self._data else 'NULL'
        image['RADIUS'] = self._data['RADIUS'] if 'RADIUS' in self._data else 'NULL'
        image['NAXIS1'] = self._data['NAXIS1'] if 'NAXIS1' in self._data else 'NULL'
        image['NAXIS2'] = self._data['NAXIS2'] if 'NAXIS2' in self._data else 'NULL'
        image['CDELT1'] = self._data['CDELT1'] if 'CDELT1' in self._data else 'NULL'
        image['CDELT2'] = self._data['CDELT2'] if 'CDELT2' in self._data else 'NULL'
        image['CRVAL1'] = self._data['CRVAL1'] if 'CRVAL1' in self._data else 'NULL'
        image['CRVAL2'] = self._data['CRVAL2'] if 'CRVAL2' in self._data else 'NULL'
        image['CRPIX1'] = self._data['CRPIX1'] if 'CRPIX1' in self._data else 'NULL'
        image['CRPIX2'] = self._data['CRPIX2'] if 'CRPIX2' in self._data else 'NULL'
        image['XCEN'] = self._data['XCEN'] if 'XCEN' in self._data else 'NULL'
        image['YCEN'] = self._data['YCEN'] if 'YCEN' in self._data else 'NULL'
        image['CROTA1'] = self._data['CROTA1'] if 'CROTA1' in self._data else 'NULL'
        
        #Fix FITS NaN parameters
        for key, value in image.items():
            if self._is_string(value):
                if value.lower() == 'nan' or value.lower() == '-nan':
                    image[key] = 'NULL'
        
        image['nickname'] = imageData.nickname
        image['observatory'] = imageData.observatory.replace(" ","_")
        image['instrument'] = imageData.instrument.split(" ")[0]
        image['detector'] = imageData.detector
        measurement = str(imageData.measurement).replace(".0 Angstrom", "").replace(".0", "")
        # Convert Yohkoh measurements to be helioviewer compatible
        if image['observatory'] == "Yohkoh":
            if measurement == "AlMg":
                image['measurement'] = "AlMgMn"
            elif measurement == "Al01":
                image['measurement'] = "thin-Al"
            else:
                image['measurement'] = measurement
        elif image['observatory'] == "Hinode":
            image['filter1'] = measurement.split("-")[0].replace(" ", "_")
            image['filter2'] = measurement.split("-")[1].replace(" ", "_")
        else:
            image['measurement'] = measurement
        image['date'] = imageData.date
        image['filepath'] = self._filepath
        image['header'] = imageData.meta
    
        return image
    
    def read_header_only_but_still_use_sunpy_map(self):
        """
        Reads the header for a JPEG200 file and returns some dummy data.
        Why does this function exist?  SunPy map objects perform some important
        homogenization steps that we would like to take advantage of in
        Helioviewer.  The homogenization steps occur on the creation of the sunpy
        map object.  All SunPy maps have the same properties, some of which are
        useful for Helioviewer to use in order to ingest JPEG2000 data.  The SunPy
        map properties are based on the header information in the JPEG2000 file,
        which is a copy of the source FITS header (with some modifications in
        some cases - see JP2Gen).  So by using SunPy's maps, Helioviewer does not
        have to implement these homogenization steps.
    
        So what's the problem?  Why not use SunPy's JPEG2000 file reading
        capability?  Well let's explain. SunPy's JPEG2000 file reading reads
        both the file header and the image data.  The image data is then decoded
        ultimately creating a numpy array.  The decoding step is computationally
        expensive for the 4k by 4k images provided by AIA and HMI.  It takes long
        enough that the ingestion of AIA and HMI data would be severely impacted,
        possibly to the point that we would never catch up if we fell behind in
        ingesting the latest data.
    
        The solution is to not decode the image data, but to pass along only the
        minimal amount of information required to create the SunPy map.  This
        function implements this solution tactic, admittedly in an unsatisfying
        manner.  The actual image data is replaced by a 1 by 1 numpy array.  This
        is sufficient to create a SunPy map with the properties required by the
        Helioviewer Project.
    
        Parameters
        ----------
        filepath : `str`
            The file to be read.
    
        Returns
        -------
        pairs : `list`
            A (data, header) tuple
        """
        header = self.get_header()
    
        return [(np.zeros([1, 1]), header[0])]
    
    
    def get_header(self):
        """
        Reads the header from the file
    
        Parameters
        ----------
        filepath : `str`
            The file to be read
    
        Returns
        -------
        headers : list
            A list of headers read from the file
        """
        with open(self._filepath, 'rb') as ifile:
            xml_txt = find_xml(ifile, 0, stat(self._filepath).st_size)
        pydict = xml_to_dict(xml_txt)["meta"]["fits"]

        # Fix types
        for k, v in pydict.items():
            if v.isdigit():
                pydict[k] = int(v)
            elif self._is_float(v):
                pydict[k] = float(v)
    
        # Remove newlines from comment
        if 'comment' in pydict:
            pydict['comment'] = pydict['comment'].replace("\n", "")
    
        self._data = pydict
        return [FileHeader(pydict)]
    
    
    def _is_float(self, s):
        """Check to see if a string value is a valid float"""
        try:
            float(s)
            return True
        except ValueError:
            return False

    def _is_string(self, s):
        # if we use Python 3
        if (sys.version_info[0] >= 3):
            return isinstance(s, str)
        # we use Python 2
        return isinstance(s, basestring)

    def getDSun(self):
        """Returns the distance to the sun in meters
        For images where dsun is not specified it can be determined using:
            dsun = (rsun_1au / rsun_image) * dsun_1au
        """
        maxDSUN = 2.25e11 # A reasonable max for solar observatories ~1.5 AU

        try:
            dsun = self._data['DSUN_OBS'] # AIA, EUVI, COR, SWAP, SXT
        except Exception as e:
            try:
                rsun = self._data['SOLAR_R'] # EIT
            except Exception as e:
                try:
                    rsun = self._data['RADIUS'] # MDI
                except Exception as e:
                    #
                    value = 1
            
            try:
                rsun
            except NameError:
                #skipping
                value = 1
            else:
                scale = self._data['CDELT1']
                if scale == 0 :
                    print('JP2 WARNING! Invalid value for CDELT1 (' + self._filepath + '): ' + scale)
                if rsun == 0 :
                    print('JP2 WARNING! Invalid value for RSUN (' + self._filepath + '): ' + rsun)
                    
                dsun = (__HV_CONSTANT_RSUN__ / (rsun * scale)) * __HV_CONSTANT_AU__

        # HMI continuum images may have DSUN = 0.00
        # LASCO/MDI may have rsun=0.00
        try:
            dsun
        except NameError:
            dsun = __HV_CONSTANT_AU__
        
        if dsun <= 0:
            dsun = __HV_CONSTANT_AU__

        # Check to make sure header information is valid
        if self._is_float(dsun) == False or dsun <= 0 or dsun >= maxDSUN:
            print('JP2 WARNING! Invalid value for DSUN (' + self._filepath + '): ' + dsun)

        return dsun

    def getImageDimensions(self):
        """Returns the dimensions for a given image
        @return array JP2 width and height
        """
        ret = dict()
        
        try:
            ret['width']  = self._data['NAXIS1']
            ret['height'] = self._data['NAXIS2']
        except Exception as e:
            print('JP2 WARNING! Unable to locate image dimensions in header tags! (' + self._filepath + ')')

        return ret

    
    def getImagePlateScale(self):
        """Returns the plate scale for a given image
        @return string JP2 image scale
        """
        try:
            scale = self._data['CDELT1']
        except Exception as e:
            print( 'JP2 WARNING! Unable to locate image scale in header tags! (' + self._filepath + ')')

        # Check to make sure header information is valid
        if self._is_float(scale) == False or scale <= 0:
            print('JP2 WARNING! Invalid value for CDELT1 (' + self._filepath + '): ' + scale)

        return scale

    def getRefPixelCoords(self):
        """Returns the coordinates for the image's reference pixel.
           NOTE: The values for CRPIX1 and CRPIX2 reflect the x and y coordinates
                 with the origin at the bottom-left corner of the image, not the
                 top-left corner.
           @return array Pixel coordinates of the reference pixel
        """
        ret = dict()
        
        try:
            if self._data['INSTRUME'] == 'XRT':
                ret['x'] = -(self._data['CRVAL1'] / self._data['CDELT1'] - self._data['CRPIX1'])
                ret['y'] = -(self._data['CRVAL2'] / self._data['CDELT2'] - self._data['CRPIX2'])
            else:
                ret['x'] = self._data['CRPIX1']
                ret['y'] = self._data['CRPIX2']
        except Exception as e:
            print( 'JP2 WARNING! Unable to locate reference pixel coordinates in header tags! (' + self._filepath + ')')

        return ret

    
    def getSunCenterOffsetParams(self):
        """Returns the Header keywords containing any Sun-center location
           information
           @return array Header keyword/value pairs from JP2 file XML
        """
        sunCenterOffsetParams = dict()

        try:
            if self._data['INSTRUME'] == 'XRT':
                sunCenterOffsetParams['XCEN'] = self._data['XCEN']
                sunCenterOffsetParams['YCEN'] = self._data['YCEN']
                sunCenterOffsetParams['CDELT1'] = self._data['CDELT1']
                sunCenterOffsetParams['CDELT2'] = self._data['CDELT2']
        except Exception as e:
            print('JP2 WARNING! Unable to locate Sun center offset params in header! (' + self._filepath + ')')

        return sunCenterOffsetParams

    def getLayeringOrder(self):
        """Returns layering order based on data source
           NOTE: In the case of Hinode XRT, layering order is decided on an
                 image-by-image basis
           @return integer layering order
        """
        try:
            telescope = self._data['TELESCOP']
            if telescope == 'SOHO':
                layeringOrder = 2     # SOHO LASCO C2
                if self._data['INSTRUME'] == 'EIT':
                    layeringOrder = 1  # SOHO EIT
                elif self._data['INSTRUME'] == 'MDI':
                    layeringOrder = 1  # SOHO MDI
                elif self._data['DETECTOR'] == 'C3':
                    layeringOrder = 3 # SOHO LASCO C3
            elif telescope == 'STEREO':
                layeringOrder = 2     # STEREO_A/B SECCHI COR1
                if self._data['DETECTOR'] == 'COR2':
                    layeringOrder = 3 # STEREO_A/B SECCHI COR2
            elif telescope == 'HINODE':
                layeringOrder = 1     # Hinode XRT full disk
                if self._data['NAXIS1'] * self._data['CDELT1'] < 2048.0 and self._data['NAXIS2'] * self._data['CDELT2'] < 2048.0:
                    layeringOrder = 2 # Hinode XRT sub-field
            else:
                # All other data sources
                layeringOrder = 1
        except Exception as e:
            print('JP2 WARNING! Unable to determine layeringOrder from header tags! (' + self._filepath + ')')

        return layeringOrder

    def getImageRotationStatus(self):
        """Returns true if the image was rotated 180 degrees
          
           Note that while the image data may have been rotated to make it easier
           to line up different data sources, the meta-information regarding the
           sun center, etc. are not adjusted, and thus must be manually adjusted
           to account for any rotation.
          
           @return boolean True if the image has been rotated
        """
        try:
            rotation = self._data['CROTA1']
            if abs(rotation) > 170:
                return true
        except Exception as e:
            # AIA, EIT, and MDI do their own rotation
            return False
