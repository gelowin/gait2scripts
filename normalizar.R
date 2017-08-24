
##' Given a numeric vector return its normalized non-parametric transformation
##'
##' NAs are allowed
##' 
##' @title Normal transformation of a numeric vector
##' @param X numeric vector to transform
##' @return Normalized vector
##' @author Angel Martinez-Perez, \email{angelmartinez@protonmail.com}
##' @examples
##' a <- c(runif(100,min=0,max=50))
##' a[sample(1:50,2)] <- NA
##' a.norm <- normalizar(a)
##' plot(density(a,na.rm=TRUE))
##' plot(density(a.norm,na.rm=TRUE))
##' 
normalizar <-function(X)
{
    if(class(X) != 'numeric')
    {
        warning('Need to be numeric')
        warning(paste('Returning the same',class(X)))
        return(X)
    }
    
    X <- rank(X,ties.method='average',na.last='keep')
    X <- qnorm(X /(length(na.omit(X)) + 1))
    return(X)
}



