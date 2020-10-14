
barcodeRanks2 <- function(m, lower=100, fit.bounds=NULL, exclude.from=50, df=20, ...) {

m <- res_mat
    totals <- unname(colSums(m))
    o <- order(totals, decreasing=TRUE)

    stuff <- rle(totals[o])
    run.rank <- cumsum(stuff$lengths) - (stuff$lengths-1)/2 # Get mid-rank of each run.
    run.totals <- stuff$values

    keep <- run.totals > lower
    if (sum(keep)<3) { 
        stop("insufficient unique points for computing knee/inflection points")
    } 
    y <- log10(run.totals[keep])
    x <- log10(run.rank[keep])
    
    # Numerical differentiation to identify bounds for spline fitting.
    edge.out <- .find_curve_bounds(x=x, y=y, exclude.from=exclude.from) 
    left.edge <- edge.out["left"]
    right.edge <- edge.out["right"]

    # As an aside: taking the right edge to get the total for the inflection point.
    # We use the numerical derivative as the spline is optimized for the knee.
    inflection <- 10^(y[right.edge])

    # We restrict curve fitting to this region, thereby simplifying the shape of the curve.
    # This allows us to get a decent fit with low df for stable differentiation.
    if (is.null(fit.bounds)) {
        new.keep <- left.edge:right.edge
    } else {
        new.keep <- y > log10(fit.bounds[1]) & y < log10(fit.bounds[2])
    }

    # Smoothing to avoid error multiplication upon differentiation.
    # Minimizing the signed curvature and returning the total for the knee point.
    fitted.vals <- rep(NA_real_, length(keep))

    if (length(new.keep) >= 4) {
	#df=20
        fit <- smooth.spline(x[new.keep], y[new.keep], df=df)  #####################modify remove.............
        fitted.vals[keep][new.keep] <- 10^fitted(fit)

        d1 <- predict(fit, deriv=1)$y
        d2 <- predict(fit, deriv=2)$y
        curvature <- d2/(1 + d1^2)^1.5
        knee <- 10^(y[new.keep][which.min(curvature)])
    } else {
        # Sane fallback upon overly aggressive filtering by 'exclude.from', 'lower'.
        knee <- 10^(y[new.keep[1]]) 
    }

    # Returning a whole stack of useful stats.
    out <- DataFrame(
        rank=.reorder(run.rank, stuff$lengths, o), 
        total=.reorder(run.totals, stuff$lengths, o),
        fitted=.reorder(fitted.vals, stuff$lengths, o)  #######################
    )
    rownames(out) <- colnames(m)
    metadata(out) <- list(knee=knee, inflection=inflection)
    out
}

.reorder <- function(vals, lens, o) {
    out <- rep(vals, lens)
    out[o] <- out
    return(out)
}

.find_curve_bounds <- function(x, y, exclude.from) 
# The upper/lower bounds are defined at the plateau and inflection, respectively.
# Some exclusion of the LHS points avoids problems with discreteness.
{
    d1n <- diff(y)/diff(x)

    skip <- min(length(d1n) - 1, sum(x <= log10(exclude.from)))
    d1n <- tail(d1n, length(d1n) - skip)

    right.edge <- which.min(d1n)
    left.edge <- which.max(d1n[seq_len(right.edge)])

    c(left=left.edge, right=right.edge) + skip
}


