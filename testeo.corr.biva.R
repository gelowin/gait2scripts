
##' Sanitize solarPolygenic objects and extract interesant data
##'
##' The "solarPolygenic" object must be created with  polygenic.options = '-testrhoe -testrhog -testrhop -rhopse'
##' @title 
##' @param kkk a solarPolygenic object
##' @return a data.frame with the sanitized results of the solarPolygenic object
##' @author Angel Martinez-Perez, \email{angelmartinez@protonmail.com}
testeo.corr.biva <- function(kkk)
### kkk es un objeto de la clase solarPolygenic que se ha creado con las opciones:
### polygenic.options = '-testrhoe -testrhog -testrhop -rhopse'
### Este objeto nos devuelve los coeficientes interesantes, o NA si el solar$ok es FALSE
{
    if(!kkk$solar$ok){
        warning(paste('modelo ha fallado: gen ',kkk$traits[2]))
        kk <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,tr1=kkk$traits[1],tr2=kkk$traits[2],NA,SOLARok=kkk$solar$ok)
        colnames(kk) <- c('rhog','Estimaterhog0','SErhog0','zscorerhog0','rhoe','Estimaterhoe','SErhoe','zscorerhoe','rhog0pval','rhoePval','tr1','tr2','N','SOLARok')
        return(kk)
    }

     
    kk <- kkk$vcf[kkk$vcf[,1]=='rhog',]
    kk <- cbind(kk, kkk$vcf[kkk$vcf[,1]=='rhoe',])
    
    kk$rhog0pval <- kkk$lf[kkk$lf$model=='rhog0',][,'pval']
                                     
    kk$rhoePval <- kkk$lf[kkk$lf$model=='rhoe0',][,'pval']
    kk$tr1 <- kkk$traits[1]
    kk$tr2 <- kkk$traits[2]
    kk$N <- as.numeric(strsplit(x=grep(pattern='Individuals:',x=kkk$solar$ret,value=TRUE),split='Individuals:' )[[1]][2]   )
    kk$SOLARok <- kkk$solar$ok
    colnames(kk) <- c('rhog','Estimaterhog0','SErhog0','zscorerhog0','rhoe','Estimaterhoe','SErhoe','zscorerhoe','rhog0pval','rhoePval','tr1','tr2','N','SOLARok')
    return(kk)
}


