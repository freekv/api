<?php 
/* vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4: */
/**
 * Image_ImageType_NSO-GONGImage class definition
 * There is one xxxImage for each type of detector Helioviewer supports.
 *
 * PHP version 5
 *
 * @category Image
 * @package  Helioviewer
 * @author   Keith Hughitt <keith.hughitt@nasa.gov>
 * @license  http://www.mozilla.org/MPL/MPL-1.1.html Mozilla Public License 1.1
 * @link     http://launchpad.net/helioviewer.org
 */
require_once HV_ROOT_DIR.'/../src/Image/HelioviewerImage.php';
/**
 * Image_ImageType_NSO-GONGImage class definition
 * There is one xxxImage for each type of detector Helioviewer supports.
 *
 * PHP version 5
 *
 * @category Image
 * @package  Helioviewer
 * @author   Keith Hughitt <keith.hughitt@nasa.gov>
 * @license  http://www.mozilla.org/MPL/MPL-1.1.html Mozilla Public License 1.1
 * @link     http://launchpad.net/helioviewer.org
 */
class Image_ImageType_NSO_GONG_MAGImage extends Image_HelioviewerImage
{
    /**
     * Creates a new NSO-GONGImage
     * 
     * @param string $jp2      Source JP2 image
     * @param string $filepath Location to output the file to
     * @param array  $roi      Top-left and bottom-right pixel coordinates on the image
     * @param string $inst     Instrument
     * @param string $det      Detector
     * @param string $meas     Measurement
     * @param int    $offsetX  Offset of the sun center from the image center
     * @param int    $offsetY  Offset of the sun center from the iamge center
     * @param array  $options  Optional parameters
     * 
     * @return void
     */     
    public function __construct($jp2, $filepath, $roi, $obs, $inst, $det, $meas, $offsetX, $offsetY, $options)
    {
        $colorTable = HV_ROOT_DIR . "/api/resources/images/color-tables/SDO_AIA_171.png";
        $this->setColorTable($colorTable);

        parent::__construct($jp2, $filepath, $roi, $obs, $inst, $det, $meas, $offsetX, $offsetY, $options);
    }
    
    /**
     * Gets a string that will be displayed in the image's watermark
     * 
     * @return string watermark name
     */
    public function getWaterMarkName() 
    {
        return "NSO-GONG $this->measurement\n";
    }
}
