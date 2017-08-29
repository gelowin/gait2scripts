##' Given the RNAseq gwas calculated with gwa.gait2.MatrixEQTL.RNAseq.R this script collect the results and put into a list
##'
##' This function applies to the output of the function "gwa.gait2.MatrixEQTL.RNAseq"
##' @title collect results of the GWAS calculate with MatrixEQTL for all RNAseq data
##' @param chr valid chromosomes are 1-22
##' @param inDir if missing "trait.resGWA.fst" is used
##' @param cis if TRUE results in "cis" GWA are collected, else results in "trans" are collected
##' @return list, one data.frame per gen, in each data.frame, the SNPs sugested and the name of the gen
##' @author Angel Martinez-Perez, \email{angelmartinez@protonmail.com}
##' @keywords GAIT2,GWAS
##' @note Update 20170829
##' @seealso \code{gwa.gait2.MatrixEQTL.RNAseq}
##' @examples
##' aux <- gwa.gait2.solariusRNA.getData(chr=2,inDir='gwaRNAseq')
##' 
gwa.gait2.solariusRNA.getData <- function(chr, inDir=NA,cis=FALSE){
    stopifnot(require(data.table,warn.conflicts=FALSE,quietly=TRUE))
    stopifnot(require(dplyr,warn.conflicts=FALSE,quietly=TRUE))
    
    if(cis)
    {
        inFile   <- file.path(inDir,paste0('gwaRNA.chr.',chr,'.MatrixEQTL.cis.txt'))
    }else{
        inFile <- file.path(inDir,paste0('gwaRNA.chr.',chr,'.MatrixEQTL.txt'))
    }
        
    stopifnot(file.exists(inFile))
    aux.gwa <- read.table(inFile,header=TRUE,stringsAsFactors = FALSE)
    aux.gwa <- aux.gwa %>% select(SNP,gene)   
    aux <- split( aux.gwa, f=aux.gwa$gene)
    return(aux)
}

