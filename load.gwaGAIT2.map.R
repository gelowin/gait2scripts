##' load into memory the map of the snps in the chromosome specified, with the attributes specified
##'
##' The rsIDref is the key that is unique
##' @title load.gwaGAIT2.map
##' @param chr The chromosome of the map to return
##' @param attributes One of the following: c("SNP", "rsIDref", "type", "source", "afr.aaf", "amr.aaf", "asn.aaf", "eur.aaf", "afr.maf", "amr.maf", "asn.maf", "eur.maf", "snpKey", "snpID", "snp", "rsID", "position", "alleleA", "alleleB", "chromosome", "exp_maf", "info", "LOCATION", "GENEID", "GENESYMBOL", "GENEdescript")
##' @return a data.frame with the required attributes of the map
##' @author Angel Martinez-Perez, \email{angelmartinez@protonmail.com}
##' @keywords GAIT2,GWAS,map
##' @examples
##' chr22 <- load.gwaGAIT2.map(chr=22,attributes=c('rsIDref','chromosome','position','LOCATION','info','GENESYMBOL'))
##'
##' @seealso \code{load.gwaGAIT2}
##' 
load.gwaGAIT2.map <- function( chr,attributes = c("SNP","rsIDref","chromosome","position","exp_maf","info","LOCATION","GENEID","GENESYMBOL","GENEdescript") )
{
    stopifnot(require(fst,warn.conflicts=FALSE,quietly=TRUE))
    stopifnot(chr %in% 1:22)
    
    nodename <- Sys.info()[["nodename"]]
    switch(nodename,
           debian = {
               dir.map <- '/home/datasets/GAIT2/04annotation/v4'
           },
           stop('*** Unknown machine ***')
           )
    
    mapfile <- file.path(dir.map, paste0('imputed.chr',chr,'.gwastools.snp.map.fst'))
    stopifnot(file.exists(mapfile))
        
    snpMap <- fstread(mapfile,columns = attributes )
    return(snpMap)
}








    
