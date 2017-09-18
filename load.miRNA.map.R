##' load into memory the map of the miRNA in the chromosome specified
##'
##' load miRNA map
##' @title load.miRNA.map
##' @param chr The chromosome of the map to return, if NULL, all chromosomes are returned, 1:22,X,Y
##' @return a data.frame with the required attributes of the map
##' @author Angel Martinez-Perez, \email{angelmartinez@protonmail.com}
##' @keywords GAIT2,GWAS,map
##' @examples
##' miRNA.map.chr5 <- load.miRNA.map(chr=5)
##'
##' @seealso \code{load.gwaGAIT2}
##' 
load.miRNA.map <- function( chr=NULL)
{
    stopifnot(require(dplyr))
    stopifnot(chr %in% c(NULL,1:22,'X','Y'))
    
    nodename <- Sys.info()[["nodename"]]
    switch(nodename,
           debian = {
               dir.getobj <- '/home/datasets/GAIT2/scripts'
               dir.map <- '/home/datasets/GAIT2/miRNA/mapas'
           },
           stop('*** Unknown machine ***')
           )
    
    mapfile <- file.path(dir.map, 'miRNA.map.RData')
    stopifnot(file.exists(mapfile))
    getobjfile <- file.path(dir.getobj,'getobj.R')
    stopifnot(file.exists(getobjfile))
    source(getobjfile)
    
    miRNA.map <- getobj(mapfile)

    if(!is.null(chr)){
        miRNA.map <- miRNA.map %>% filter(seqnames %in% paste0('chr',chr))
    }

    return(miRNA.map)       
    
}








    
