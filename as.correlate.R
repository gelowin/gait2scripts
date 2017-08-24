##' transform a correlation data.frame into de class "correlate"
##'
##' .. content for \details{} ..
##' @title 
##' @param my_cor a correlation data.frame
##' @return object of class correlate
##' @author Angel Martinez-Perez, \email{angelmartinez@protonmail.com}
as.correlate <- function(my_cor)
{
    
    stopifnot(require(corrr,warn.conflicts=FALSE,quietly=TRUE))
    
    if(class(my_cor) !='data.frame')
        stop('my_cor need to be a data.frame')
    my_cor<- corrr::first_col(my_cor, rowname=colnames(my_cor))
    class(my_cor) <- c('cor_df',class(my_cor))
    return(my_cor)
}
