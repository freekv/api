<?php
/* vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4: */
/**
 * Movie_HelioviewerMovie Class Definition
 *
 * PHP version 5
 *
 * @category Movie
 * @package  Helioviewer
 * @author   Keith Hughitt <keith.hughitt@nasa.gov>
 * @author   Jaclyn Beck <jaclyn.r.beck@gmail.com>
 * @license  http://www.mozilla.org/MPL/MPL-1.1.html Mozilla Public License 1.1
 * @link     http://launchpad.net/helioviewer.org
 */
require_once 'src/Image/Composite/HelioviewerMovieFrame.php';
require_once 'src/Helper/DateTimeConversions.php';
require_once 'src/Database/ImgIndex.php';
require_once 'src/Movie/FFMPEGEncoder.php';
/**
 * Represents a static (e.g. mp4/webm) movie generated by Helioviewer
 *
 * Note: For movies, it is easiest to work with Unix timestamps since that is what is returned
 *       from the database. To get from a javascript Date object to a Unix timestamp, simply
 *       use "date.getTime() * 1000." (getTime returns the number of miliseconds)
 *
 * @category Movie
 * @package  Helioviewer
 * @author   Keith Hughitt <keith.hughitt@nasa.gov>
 * @author   Jaclyn Beck <jaclyn.r.beck@gmail.com>
 * @license  http://www.mozilla.org/MPL/MPL-1.1.html Mozilla Public License 1.1
 * @link     http://launchpad.net/helioviewer.org
 */
class Movie_HelioviewerMovie
{
	protected $id;
    private $_db;
    private $_layers;
    private $_roi;
    private $_directory;
    private $_filename;
    private $_startTimestamp;
    private $_startDateString;
    private $_endTimestamp;
    private $_endDateString;
    private $_frames;
    private $_frameRate;
    private $_numFrames;
    
    /**
     * Prepares the parameters passed in from the api call and makes a movie from them.
     *
     * @return {String} a url to the movie, or the movie will display.
     */
    public function __construct($layers, $startDateString, $endDateString, $roi, $options)
    {
        $defaults = array(
            'format'      => "mp4",
            'frameRate'   => false,
            'maxFrames'   => HV_MAX_MOVIE_FRAMES,
            'watermarkOn' => true
        );
        $options = array_replace($defaults, $options);

        $this->_db        = new Database_ImgIndex();
        $this->_layers    = $layers;
        $this->_roi       = $roi;
        
        $this->_startDateString = $startDateString;
        $this->_endDateString   = $endDateString;

        $this->id = $this->_getMovieId($options['watermarkOn']);
        
        $this->_directory = $this->_buildDir();
        $this->_filename  = $this->_buildFilename($options['format']);
        
        // Also store as timestamps
        $this->_startTimestamp = toUnixTimestamp($startDateString);
        $this->_endTimestamp   = toUnixTimestamp($endDateString);
        
        var_dump($this);
        die();

        // Get timestamps for frames in the key movie layer
        $this->_timestamps = $this->_getTimeStamps($options['maxFrames']);

        $this->_numFrames  = sizeOf($this->_timestamps);

        if ($this->_numFrames == 0) {
        	$this->_abort("No images available for the requested time range");
        }

        //$this->_filename  = $this->_buildFilename($options['filename']);
        
        $this->_frameRate = $this->_determineOptimalFrameRate($options['frameRate']);

        $this->_setMovieDimensions();
        
        // Build movie frames
        $images = $this->_buildMovieFrames($options['watermarkOn']);

        // Compile movie
        $this->_build($images);
    }
    
    /**
     * Adds the movie to the database and returns its assigned identifier
     * 
     * @return int Movie id
     */
    private function _getMovieId($watermarkOn)
    {
        return $this->_db->insertMovie(
            $this->_startDateString,
            $this->_endDateString,
            $this->_roi->imageScale(),
            $this->_roi->getPolygonString(),
            $watermarkOn,
            $this->_layers->serialize(),
            $this->_layers->getBitMask()
        );
    }
    
    /**
     * Returns an array of the timestamps for the key movie layer
     * 
     * For single layer movies, the number of frames will be either HV_MAX_MOVIE_FRAMES, or the number of
     * images available for the requested time range. For multi-layer movies, the number of frames included
     * may be reduced to ensure that the total number of SubFieldImages needed does not exceed HV_MAX_MOVIE_FRAMES
     */
    private function _getTimeStamps($maxFrames)
    {
        $layerCounts = array();

        // Determine the number of images that are available for the request duration for each layer
        foreach ($this->_layers->toArray() as $layer) {
            $n = $this->_db->getImageCount($this->_startDateString, $this->_endDateString, $layer['sourceId']);
            $layerCounts[$layer['sourceId']] = $n;
        }

        // Choose the maximum number of frames that can be generated without exceeded the server limits defined
        // by HV_MAX_MOVIE_FRAMES
        $numFrames       = 0;
        $imagesRemaining = $maxFrames;
        $layersRemaining = $this->_layers->length();
        
        // Sort counts from smallest to largest
        asort($layerCounts);
        
        // Determine number of frames to create
        foreach($layerCounts as $dataSource => $count) {
            $numFrames = min($count, ($imagesRemaining / $layersRemaining));
            $imagesRemaining -= $numFrames;
            $layersRemaining -= 1;
        }
        
        // Number of frames to use
        $numFrames = floor($numFrames);

        // Get the entire range of available images between the movie start and end time 
        $entireRange = $this->_db->getImageRange($this->_startDateString, $this->_endDateString, $dataSource);
        
        // Sub-sample range so that only $numFrames timestamps are returned
        $timestamps = array();
        for ($i = 0; $i < $numFrames; $i++) {
        	$index = round($i * (sizeOf($entireRange) / $numFrames));
        	array_push($timestamps, $entireRange[$index]['date']);
        }
        return $timestamps;        
    }

    /**
     * Determines the directory to store the movie in.
     * 
     * @return string Directory
     */
    private function _buildDir ()
    {
        return sprintf("%s/movies/%s/%s/", HV_CACHE_DIR, date("Y/m/d"), $this->id);
    }

    /**
     * Determines filename to use for the movie
     * 
     * @param string $extension Extension of the movie format to be created 
     *
     * @return string Movie filename
     */
    private function _buildFilename($extension) {
        $start = str_replace(array(":", "-", "T", "Z"), "_", $this->_startDateString);
        $end   = str_replace(array(":", "-", "T", "Z"), "_", $this->_endDateString);

        return sprintf("%s_%s_%s.%s", $start, $end, $this->_layers->toString(), $extension);
    }

    /**
     * Takes in meta and layer information and creates movie frames from them.
     *
     * @param {String} $tmpDir     the directory where the frames will be stored
     *
     * @return $images an array of built movie frames
     */
    private function _buildMovieFrames($watermarkOn)
    {
        $movieFrames  = array();

        $frameNum = 0;

        // Movie frame parameters
        $options = array(
            'compress'   => false,
            'interlace'  => false,
            'format'     => 'bmp',
            'watermarkOn'=> $watermarkOn,
            'outputDir'  => $this->_directory . "/frames"
        );

        // Compile frames
        foreach ($this->_timestamps as $time) {
            $options = array_merge($options, array(
                'filename' => "frame" . $frameNum
            ));

            try {
	            $screenshot = new Image_Composite_HelioviewerMovieFrame($this->_layers, $time, $this->_roi, $options);
	            $filepath   = $screenshot->getFilepath();
	            $frameNum++;
	            array_push($movieFrames, $filepath);
            } catch (Exception $e) {
            	// Recover if failure occurs on a single frame
            	$this->_numFrames--;
            }
        }

        // TODO 2011/02/14: Verify that this is still necessary
        // Copy the last frame so that it actually shows up in the movie for the same amount of time
        // as the rest of the frames.
        $lastImage = dirname($filepath) . "/frame" . $frameNum . "." . $options['format'];

        copy($filepath, $lastImage);
        array_push($movieFrames, $lastImage);
        
        // Create preview image
        // TODO: Use middle frame instead last one...
        // TODO: Create standardized thumbnail sizes (e.g. thumbnail-med.png = 480x320, etc)
        $imagickImage = $screenshot->getIMagickImage();
        $imagickImage->setImageCompression(IMagick::COMPRESSION_LZW);
        $imagickImage->setImageCompressionQuality(PNG_LOW_COMPRESSION);
        $imagickImage->setInterlaceScheme(IMagick::INTERLACE_PLANE);
        $imagickImage->writeImage($this->_directory . "/" . $this->_filename . ".png");
        $imagickImage->destroy();

        return $movieFrames;
    }

    /**
     * Uses numFrames to calculate the frame rate that should
     * be used when encoding the movie.
     *
     * @return Int optimized frame rate
     */
    private function _determineOptimalFrameRate($requestedFrameRate)
    {
        // Subtract 1 because we added an extra frame to the end
        $frameRate = ($this->_numFrames - 1 ) / HV_DEFAULT_MOVIE_PLAYBACK_IN_SECONDS;

        // Take the smaller number in case the user specifies a larger frame rate than is practical.
        if ($requestedFrameRate) {
            $frameRate = min($frameRate, $requestedFrameRate);
        }

        return max(1, $frameRate);
    }

    /**
     * Builds the requested movie
     *
     * Makes a temporary directory to store frames in, calculates a timestamp for every frame, gets the closest
     * image to each timestamp for each layer. Then takes all layers belonging to one timestamp and makes a movie frame
     * out of it. When done with all movie frames, phpvideotoolkit is used to compile all the frames into a movie.
     *
     * @param array  $builtImages An array of built movie frames (in the form of HelioviewerCompositeImage objects)
     *
     * @return void
     */
    private function _build($builtImages)
    {
        $this->_frames = $builtImages;

        // Create and FFmpeg encoder instance
        $ffmpeg = new Movie_FFMPEGEncoder($this->_frameRate);

        // TODO 11/18/2010: add 'ipod' option to movie requests in place of the 'hqFormat' param
        $ipod = false;

        if ($ipod) {
            $ffmpeg->createIpodVideo($this->_directory, $this->_filename, "mp4", $this->_width, $this->_height);
        }
        
        // Create a high-quality H.264 video using an MPEG-4 (mp4) container format
        $ffmpeg->createVideo($this->_directory, $this->_filename . "-hq", "mp4", $this->_width, $this->_height, "ultrafast", 15);

        // Create a medium-quality H.264 video using an MPEG-4 (mp4) container format
        $ffmpeg->createVideo($this->_directory, $this->_filename, "mp4", $this->_width, $this->_height);

        //Create alternative container format options for medium-quality video (.flv)
        $ffmpeg->createAlternativeVideoFormat($this->_directory, $this->_filename, "mp4", "flv");
    }

    /**
     * Determines dimensions to use for movie and stores them
     * 
     * @return void
     */
    private function _setMovieDimensions() {
        $this->_width  = round($this->_roi->getPixelWidth());
        $this->_height = round($this->_roi->getPixelHeight());

        // Width and height must be divisible by 2 or ffmpeg will throw an error.
        if ($this->_width % 2 === 1) {
            $this->_width += 1;
        }
        
        if ($this->_height % 2 === 1) {
            $this->_height += 1;
        } 
    }

    /**
     * Cancels movie request
     * 
     * TODO 11/24/2010: Cleanup files
     * TODO 11/24/2010: Instead of using files to mark movie status, could instead use presence of 'frames' directory
     *                  and expected movies (frames directory is delete after movies are finished). In the longer run, 
     *                  movie status should be tracked in a database accessible to both Helioviewer and Helioqueuer.
     */
    private function _abort($msg) {
        touch($this->_directory . "/INVALID");
        throw new Exception("Unable to create movie: " . $msg, 1);
    }

    /**
     * Adds black border to movie frames if neccessary to guarantee a 16:9 aspect ratio
     *
     * Checks the ratio of width to height and adjusts each dimension so that the
     * ratio is 16:9. The movie will be padded with a black background in JP2Image.php
     * using the new width and height.
     *
     * @return array Width and Height of padded movie frames
     */
    private function _setAspectRatios()
    {
        $width  = $this->_roi->getPixelWidth();
        $height = $this->_roi->getPixelHeight();

        $ratio = $width / $height;

        // Commented out because padding the width looks funny.
        /*
        // If width needs to be adjusted but height is fine
        if ($ratio < 16/9) {
        $adjust = (16/9) * $height / $width;
        $width *= $adjust;
        }
        */
        // Adjust height if necessary
        if ($ratio > 16/9) {
            $adjust = (9/16) * $width / $height;
            $height *= $adjust;
        }

        $dimensions = array("width" => $width, "height" => $height);
        return $dimensions;
    }
    
    /**
     * Returns the base filepath for movie without any file extension
     */
    public function getFilepath()
    {
        return $this->_directory . "/" . $this->_filename;
    }
    
    /**
     * Returns the Base URL to the most recently created movie (without a file extension)
     */
    public function getURL()
    {
        return str_replace(HV_ROOT_DIR, HV_WEB_ROOT_URL, $this->getFilepath());
    }
    
    /**
     * Returns the movie frame rate
     */
    public function getFrameRate()
    {
        return $this->_frameRate;
    }
    
    /**
     * Returns the number of frames in the movie
     */
    public function getNumFrames()
    {
        return $this->_numFrames;
    }
    
    public function getDuration()
    {
        return $this->_numFrames / $this->_frameRate;
    }
    
    /**
     * Returns HTML for a video player with the requested movie loaded
     */
    public function getMoviePlayerHTML()
    {
        $filepath = str_replace(HV_ROOT_DIR, "../", $this->getFilepath());
        $css      = "width: {$this->_width}px; height: {$this->_height}px;";
        $duration = $this->_numFrames / $this->_frameRate;
        ?>
<!DOCTYPE html> 
<html> 
<head> 
    <title>Helioviewer.org - <?php echo $this->_filename;?></title>            
    <script type="text/javascript" src="http://html5.kaltura.org/js"></script> 
    <script src="http://ajax.googleapis.com/ajax/libs/jquery/1.4.2/jquery.js" type="text/javascript"></script>
</head> 
<body>
<div style="text-align: center;">
    <div style="margin-left: auto; margin-right: auto; <?php echo $css;?>";>
        <video style="margin-left: auto; margin-right: auto;" poster="<?php echo "$filepath.bmp"?>" durationHint="<?php echo $duration?>">
            <source src="<?php echo "$filepath.mp4"?>" /> 
            <source src="<?php echo "$filepath.mov"?>" />
            <source src="<?php echo "$filepath.flv"?>" /> 
        </video>
    </div>
</div>
</body> 
</html> 
        <?php        
    }
    
    /**
     * Destructor
     * 
     * @return void
     */
    public function __destruct()
    {
        // Clean up movie frame images that are no longer needed
        foreach ($this->_frames as $image) {
            if (file_exists($image)) {
                unlink($image);
            }
        }

        rmdir($this->_directory . "/frames");
        touch($this->_directory . "/READY");
    }
}
?>
