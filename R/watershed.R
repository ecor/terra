# Author: Emanuele Cordano
# Date : October 2023
# Version 1.0
# License GPL v3

setMethod("watershed", signature(x="SpatRaster"), 
    function(x, pourpoint, filename="", ...) { 
        opt <- spatOptions(filename, ...)		
		cell <- cellFromXY(x, pourpoint)
		if (is.na(cell)) error("watershed", "pourpoint not on raster")
        x@ptr <- x@ptr$watershed2(as.integer(cell-1), opt)
        messages(x, "watershed") ## EC 20210318
    }
)

setMethod("pitfinder", signature(x="SpatRaster"), 
    function(x, filename="", ...) { 
        opt <- spatOptions(filename, ...)
        x@ptr <- x@ptr$pitfinder2(opt)
        messages(x, "pitfinder") ## EC 20210318
    }
)

setMethod("NIDP", signature(x="SpatRaster"), 
    function(x, filename="", ...) { 
        opt <- spatOptions(filename, ...)
        x@ptr <- x@ptr$NIDP2(opt)
        messages(x, "NIDP") ## EC 20231031
    }
)

setMethod("flowAccumulation", signature(x="SpatRaster"), 
    function(x, weight=NULL, filename="", ...) { 
		opt <- spatOptions(filename, ...)
		if (is.null(weight)) {      
			x@ptr <- x@ptr$flowAccu2(opt)
	    } else {
			x@ptr <- x@ptr$flowAccu2_weight(weight@ptr, opt)
		} 
		messages(x, "flowAccumulation") 
    }      
)


