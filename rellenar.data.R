##' given a data.frame with NA, complete it with the transpose of this element
##'
##' .. content for \details{} ..
##' @title 
##' @param m squared numeric data.frame
##' @return 
##' @author Angel Martinez-Perez, \email{angelmartinez@protonmail.com}
rellenar.data <- function(m) {
    ## function complete a squared data frame in the following manner:
    ## if M_ij != NA then return M_ij
    ## if M_ij == NA then return M_ji   

    
    if(class(m) != 'data.frame'){
        warning('Argument need to be a data.frame')
        return(m)
    }
    if(ncol(m) != nrow(m))
        stop('m must be a squared data.frame')
       
            
  aux.1 <- is.na(m[lower.tri(m)])
  m[lower.tri(m)][aux.1] <- t(m)[lower.tri(m)][aux.1]

  aux.2 <- is.na(m[upper.tri(m)])
  m[upper.tri(m)][aux.2] <- t(m)[upper.tri(m)][aux.2]
  
  return(m)
}

