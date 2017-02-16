# ovl.R

# kdeOvl() was written to work specifically for the ABC data from Bi et al.
# and so, while general, does not provide any kind of data checking or exception handling

kdeOvl <- function(x1, x2, bandwidth)
{
	# estimates Weitzman's overlap coefficient (OVL)
	# 'x1' and 'x2' should both be numeric vectors representing D_ML_obs and D_ML_pseudo
	# 'bandwidth' is the bandwidth provided to density() for estimating the KDE

	# kde estimation
	cutsz <- 3;
	alldat <- c(x1,x2)
	minx <- min(alldat)
	maxx <- max(alldat)
	lb <- minx - cutsz*bandwidth
	ub <- maxx + cutsz*bandwidth

	kde1 <- density(x1, bw=bandwidth, kernel="gaussian", from=lb, to=ub, n=1000)
	kde2 <- density(x2, bw=bandwidth, kernel="gaussian", from=lb, to=ub, n=1000)

	# approximate functions for the KDEs
	f1 <- approxfun(kde1$x,kde1$y,yleft=0,yright=0) # could use rule=2 instead of yleft/yright
	f2 <- approxfun(kde2$x,kde2$y,yleft=0,yright=0)

	# find where kde1 intersects with kde2

	# determine values over which to calculate kde estimates
	xvals <- unique(sort(c(kde1$x,kde2$x)))

	# find where kde1 lies above kde2
	above <- f1(xvals)>f2(xvals)

	# find points of intersection (i.e. where kde1 switches from being > kde2 to < kde2)
	intersect <- which(diff(above)!=0)

	# set limits of integration
	limits <- unique(c(min(xvals), xvals[intersect], max(xvals)))

	# integrate to get area of overlap

	i <- 1
	area <- 0;
	while(i+1 <= length(limits))
	{
		interval_sample <- runif(1,limits[i],limits[i+1])
		ifelse(f1(interval_sample)<f2(interval_sample), area<-area+integrate(f1,limits[i],limits[i+1])$value, area<-area+integrate(f2,limits[i],limits[i+1])$value)
		i <- i+1
	}

	return(area)
}

