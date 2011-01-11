/**
 * @fileOverview Contains the class definition for an UserVideoGallery class.
 * @author <a href="mailto:keith.hughitt@nasa.gov">Keith Hughitt</a>
 * 
 * TODO 2011/01/10: Add a loading indicator
 * 
 */
/*jslint browser: true, white: true, onevar: true, undef: true, nomen: false, eqeqeq: true, plusplus: true, 
bitwise: true, regexp: true, strict: true, newcap: true, immed: true, maxlen: 120, sub: true */
/*global $, window, Class, addIconHoverEventListener */
"use strict";
var UserVideoGallery = Class.extend(
    /** @lends UserVideoGallery.prototype */
    {
    /**
     * @constructs
     * @description Creates a new UserVideoGallery component
     */
    init : function () {
        this._container   = $("#user-video-gallery-main");
        this._loader      = $("#user-video-gallery-spinner");
        this._nextPageBtn = $("#user-video-gallery-next");
        this._prevPageBtn = $("#user-video-gallery-prev");
        
        this._maxPageNum  = 50;
        
        this._working     = false;
        
        this._pageSize  = this._choosePageSize();
        this._pageNum   = 1;

        this._setupEventHandlers();
        this._fetchVideos();
      
        // TODO 2011/01/10: Add resize handler
    },
    
    /**
     * Retrieves a single page of video results and displays them to the user
     */
    _fetchVideos: function () {
        // Query parameters
        var params = {
            "action"   : "getUserVideos",
            "pageSize" : this._pageSize,
            "pageNum"  : this._pageNum
        };
        
        // Show loading indicator
        this._container.find("a").empty();
        this._loader.show();
        
        this._working = true;

        // Fetch videos
        $.getJSON("api/index.php", params, $.proxy(this._buildHTML, this));
    },
    
    /**
     * Builds video gallery HTML
     */
    _buildHTML: function (videos) {
        var html = "",
            now  = new Date().getTime();
        
        $.each(videos, function (i, vid) {
            var when = toFuzzyTime((now - getUTCTimestamp(vid.published)) / 1000) + " ago";
            
            html += "<a target='_blank' href='" + vid.watch + "' alt='video thumbnail'>" +
                    "<div class='user-video-thumbnail-container'>" +
                    "<div style='text-align: left; margin-left: 25px;'>" + when + "</div>" +
                    "<img src='" + vid.thumbnails[4].url + "' alt='user video thumbnail' /></div></a><br />";                        
        });
        
        // Drop tailing line break
        html = html.slice(0, -6);

        this._loader.hide();
        this._container.append(html);
        
        this._working = false;
    },
    
    /**
     * Go to the previous page
     */
    _prevPage: function () {
        if (this._working) {
            return false;
        }        
        
        if (this._pageNum < this._maxPageNum) {
            this._pageNum += 1;
        }
        
        this._fetchVideos();
    },
    
    /**
     * Go to the next page
     */
    _nextPage: function () {
        if (this._working) {
            return false;
        }

        if (this._pageNum > 1) {
            this._pageNum -= 1;
        }
        
        this._fetchVideos();
    },
    
    /**
     * Performs a quick check on the browser height to decide how many video thumbnails to show at once
     */
    _choosePageSize: function () {
        var height = $(window).height();

        if (height > 1035) {
            return 3;
        } else if (height > 780) {
            return 2;
        }        
        return 1;
    },
    
    /**
     * Setup event handlers
     */
    _setupEventHandlers: function () {
        // TODO 2011/01/10 Apply hover listen at .ui-icon level?
        addIconHoverEventListener(this._nextPageBtn);
        addIconHoverEventListener(this._prevPageBtn);
        
        this._nextPageBtn.click($.proxy(this._nextPage, this));
        this._prevPageBtn.click($.proxy(this._prevPage, this));
    }    
});