##' Load the RNAseq data of GAIT2
##'
##' Put into memory RNAseq data.frame of the GAIT2  .. Update on 20170821 ..
##'
##' @title load.RNAseqGAIT2: Load the map and data of RNAseq of GAIT2
##' @param rnaFile File to import, available files are: c("RNAseq.fst","RNAseq.normal.fst","RNAseq.normal.5pc.fst", "RNAseq.res.fst","RNAseq.res.noAGE.fst","RNAseq.res.noCounts.fst","RNAseq.res.noAGE.noCounts.fst")
##' @param gens If gens=c('id',gens) only "gens" are loaded. gens=NULL (default) to load all. 
##' @return The following object
##' \item{data.frame}{with the data of RNAseq}
##' @author Angel Martinez-Perez, \email{angelmartinez@protonmail.com}
##' @keywords GAIT2
##' @note c("RNAseq.fst","RNAseq.normal.fst","RNAseq.normal.5pc.fst") from Andrew use genNames=gene, from RNAseq.map, while c("RNAseq.res.fst","RNAseq.res.noAGE.fst","RNAseq.res.noCounts.fst","RNAseq.res.noAGE.noCounts.fst") ie, residuals from Sergui, use genNames=ENSEMBL from RNAseq.map. Only works in the allowed machines
##' 
##' @examples
##' RNAseq.dat <- load.RNAseqGAIT2( gens=c('id','ENSG00000238009','ENSG00000228463') )
##' RNAseq.dat.normal <- load.RNAseqGAIT2(rnaFile='RNAseq.normal.5pc.fst', gens=c('id','ENSG00000238009.2','ENSG00000228463.4') )
##'
##' @seealso \code{load.RNAseqGAIT2.map}
##' 
load.RNAseqGAIT2 <- function(rnaFile='RNAseq.res.noAGE.noCounts.fst',gens=NULL)
{

    stopifnot(require(fst,warn.conflicts=FALSE,quietly=TRUE))

    stopifnot(rnaFile %in% c("RNAseq.fst","RNAseq.normal.fst","RNAseq.normal.5pc.fst", "RNAseq.res.fst","RNAseq.res.noAGE.fst","RNAseq.res.noCounts.fst","RNAseq.res.noAGE.noCounts.fst"))
    
    
    nodename <- Sys.info()[["nodename"]]
    switch(nodename,
           debian = {
               dir.dat <- '/home/datasets/GAIT2/Expression'
           },
           stop('*** Unknown machine ***')
           )

    print(paste('loading',rnaFile,' of GAIT2, from',nodename))
    file.dat <- file.path(dir.dat,rnaFile)
    stopifnot(file.exists(file.dat))
    RNAdat <- fstread(file.dat,columns=gens)
    return(RNAdat)    
}



##' RNAseq.map with the gen map positions
##'
##' Only available in certain computers. It was enrich with "R package org.Hs.eg.db, Version 3.4.0". Last update in 20170821
##' @title load the RNAseq.map of GAIT2
##' @return RNAseq.map data.frame with the RNAseq map of GAIT2
##' \item {ENSEMBL}{gene id in ENSEMBL format}
##' \item{chr}{chromosome}
##' \item{gene}{gene id in original format}
##' \item{SYMBOL}{gene id in symbol format}
##' \item{descript}{Sort description of the name of the gene}
##' @author Angel Martinez-Perez, \email{angelmartinez@protonmail.com}
##' @keywords GAIT2,RNAseq,map
##' @examples
##' load.RNAseqGAIT2.map()
##'
##' @seealso \code{load.RNAseqGAIT2}
##' 
load.RNAseqGAIT2.map <- function()
{
    
    nodename <- Sys.info()[["nodename"]]
    switch(nodename,
           debian = {
               dir.dat <- '/home/datasets/GAIT2/Expression'
           },
           stop('*** Unknown machine ***')
           )
    
    file.map <- file.path(dir.dat,"RNAseq.map.RData")
    stopifnot(file.exists(file.map))
    
    load(file.map, envir = .GlobalEnv)
}


