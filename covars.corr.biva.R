
##' build the covariate part of the formula in a bivariate model 
##'
##' 
##' @title 
##' @param tr1 first trait of the bivariate model 
##' @param tr2 second trait of the bivariate model
##' @param covarsList list of the covariates for each trait
##' @return character that contain the right part of the formula in the bivariate model
##' @author Angel Martinez-Perez, \email{angelmartinez@protonmail.com}
covars.corr.biva <- function(tr1, tr2, covarsList){
    ## este script nos devuelve las covariables tal y como las necesita solar para pasarlas
    ## para un modelo bivariado
    ## para los traits: tr1 y tr2
    ## dada la lista de covariables: covarsList (normalmente será "g2.ph.covariates"

    
    comun.cov <- na.omit(intersect(covarsList[[tr1]],covarsList[[tr2]]))
    tr1.cov <- na.omit(setdiff(covarsList[[tr1]],covarsList[[tr2]]))
    tr2.cov <- na.omit(setdiff(covarsList[[tr2]],covarsList[[tr1]]))

    covars.fin <- character(0)

    if(length(comun.cov)>0)
        covars.fin <- paste(comun.cov, collapse = ' + ')

  if(length(tr1.cov)>0)
        covars.fin <- paste(c(covars.fin, paste0(tr1.cov,'(',tr1,')')), collapse = ' + ')

  if(length(tr2.cov)>0)
        covars.fin <- paste(c(covars.fin, paste0(tr2.cov,'(',tr2,')')), collapse = ' + ')

    ## Si ambos traits no tienen covariables, devolvemos "SEX" para que no den error los modelos siguientes.
    if(length(covars.fin) == 0 )
        covars.fin <- 'SEX'
        
    return(covars.fin)
}
    
