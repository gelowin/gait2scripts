
##' load the GWAS data of GAIT2
##'
##' load the genotype data of the GAIT2
##' @title load.gwaGAIT2
##' @param chr numeric chromosome to retrieve, 1-22
##' @param snps NULL for all snps in the chromosome, or the "rsIDref" of the desired variants to retrieve
##' @return a data.frame with the values of the genetic variants for all the subject in the GAIT2 study
##' @author Angel Martinez-Perez, \email{angelmartinez@protonmail.com}
##' @keywords GAIT2,GWAS
##' @note Time to load all chr22, 151k snps = 33seg. Time to load chr2, 907k snps = 5.6min. Time to load 8 arbitrary snps in chr2 = 0.8seg. ALL times measured while system was copying files.
##' @examples
##' map <- load.gwaGAIT2.map(chr=2, attributes=c('rsIDref','chromosome','position'))
##' require(dplyr)
##' algunos <- map %>% filter(position>120990000, position<121000000) %>% select(rsIDref) %>% unlist
##' chr2 <- load.gwaGAIT2(chr=2, snps=c('id',algunos))
##' 
##' @seealso \code{load.gwaGAIT2.map}
##' 
load.gwaGAIT2 <- function( chr , snps = NULL)
{

    stopifnot(require(fst,warn.conflicts=FALSE,quietly=TRUE))

    nodename <- Sys.info()[["nodename"]]
    switch(nodename,
           debian = {
               dir.fst <- '/home/datasets/GAIT2/05postImputedSOLARdosages/fst'
               ## dir.map <- '/home/datasets/GAIT2/04annotation/v4'
           },
           stop('*** Unknown machine ***')
           )

    fileAUX <- file.path(dir.fst,paste0('imputed.chr',chr,'.fst'))
    stopifnot(file.exists(fileAUX))
    
    dat <- fstread(fileAUX,columns=snps)
    return(dat)
}



