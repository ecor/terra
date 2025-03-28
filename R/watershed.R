# Author: Emanuele Cordano
# Date : October 2023
# Version 1.0
# License GPL v3

setMethod("watershed", signature(x="SpatRaster"), 
    function(x, pourpoint, filename="", ...) { 
        opt <- spatOptions(filename, ...)		
		cell <- cellFromXY(x, pourpoint)
		if (is.na(cell)) error("watershed", "pourpoint not on raster")
        x@pntr <- x@pntr$watershed2(as.integer(cell-1), opt)
        messages(x, "watershed") ## EC 20210318
    }
)

setMethod("pitfinder", signature(x="SpatRaster"), 
    function(x,pits_on_boundary=TRUE,filename="", ...) { 
        opt <- spatOptions(filename, ...)
        x@pntr <- x@pntr$pitfinder2(as.integer(pits_on_boundary),opt)

        messages(x, "pitfinder") ## EC 20210318
    }
)

setMethod("NIDP", signature(x="SpatRaster"), 
    function(x, filename="", ...) { 
        opt <- spatOptions(filename, ...)
        x@pntr <- x@pntr$NIDP2(opt)
        messages(x, "NIDP") ## EC 20231031
    }
)

setMethod("flowAccumulation", signature(x="SpatRaster"), 
    function(x, weight=NULL, filename="", ...) { 
		opt <- spatOptions(filename, ...)
		if (is.null(weight)) {      
			x@pntr <- x@pntr$flowAccu2(opt)
	    } else {
			x@pntr <- x@pntr$flowAccu2_weight(weight@pntr, opt)
		} 
		messages(x, "flowAccumulation") 
    }      
)

setMethod("flowdirD8ltd", signature(x="SpatRaster"), 
          function(x,lambda=0.5,deviation_type=c("ltd","lad"),filename="", ...) { 
            ## http://www.idrologia.unimore.it/orlandini/web-archive/seminars/nyc-2008-2.pdf
            ## ltd least transverse deviation
            ## lad least angular deviation
            deviation_type=deviation_type[1]
            use_lad=0
            if (deviation_type=="lad") use_lad=1
            opt <- spatOptions(filename, ...)
          ##  uselad=0
            x@pntr <- x@pntr$d8ltd(lambda,use_lad,opt)
            messages(x, "flowdirD8ltd") ## EC 20210318
          }
)

setMethod("flowdirD8lad", signature(x="SpatRaster"), 
          function(x,lambda=0.5,deviation_type="lad",filename="",...) { 
            ## http://www.idrologia.unimore.it/orlandini/web-archive/seminars/nyc-2008-2.pdf
            ## ltd least transverse deviation
            ## lad least angular deviation
           flowdirD8ltd(x=x,lambda=lambda,deviation_type=deviation_type,filename=filename,...)
          }
)

## EC 20241027


setMethod("pitfiller", signature(x="SpatRaster"), 
          function(x,pit=NULL,flowdir=NULL,niter=10,lambda=0,deviation_type="lad",U=1,D=300,beta=0.9,theta_exp=0.5,filename="",...) { 
           
            if (is.null(flowdir)) flowdir <- terrain(x,"flowdir") 
            if (is.null(pitfinder)) pit <- pitfinder(flowdir) 
            use_lad=1
            if (deviation_type=="ltd") use_lad=0
            flowdir[pit>0] <- 0 
            opt <- spatOptions(filename, ...)
            ##  uselad=0
            print(pit)
            x@pntr <- x@pntr$pitfillerm(pit@pntr,flowdir@pntr,niter,lambda,use_lad,U,D,beta,theta_exp,opt)
            messages(x, "pitfiller") ## EC 20210318
          }
)


# SpatRaster  SpatRaster::pitfillerm(SpatRaster pits,SpatRaster flowdirs,int niter, double lambda,int use_lad,
#                                    double U,double D,double beta,double theta_exp, // see // see reference doi:10.1016/j.advwatres.2006.11.016)    
# SpatOptions &opt) {