##' Load the phenotypes of GAIT2
##'
##' Put into memory a R6 object with the phenotypes and the functions to extract the phenotypes.. Update on 20170818 ..
##'
##' @title load.fenoGAIT2: Load the phenotypes of GAIT2
##' @return The following object
##' \item{f1}{a R6 object, with the the phenotypes and functions to access the phenotypes of the GAIT2 project}
##' \item{as.correlate}{function}
##' \item{covars.corr.biva}{function}
##' \item{normalizar}{function}
##' \item{testeo.corr.biva}{function}
##' \item{definitions.gait2.500}{data.frame with the heritabilyties of the phenotypes}
##' @author Angel Martinez-Perez, \email{angelmartinez@protonmail.com}
##' @keywords GAIT2
##' @note only works in the allowed machines
##' 
##' @examples
##' load.fenoGAIT2()
##' f1$infoCategories()
##' traits <- c( f1$getTraits('solar'), f1$getTraits('fn.proteinALL'))
##' f1$setTraits( unique(traits) )
##' dat <- f1$getData()
##' f1$plot.pedigree('fam31')
##' 
load.fenoGAIT2 <- function()
{

    stopifnot(require(R6,warn.conflicts=FALSE,quietly=TRUE))
    
    nodename <- Sys.info()[["nodename"]]
    switch(nodename,
           debian = {
               dir.dat <- '/home/datasets/GAIT2/phen/classR6'
               dir.script <- '/home/datasets/GAIT2/scripts/'
           },
           salambo2 = {
               dir.dat <- '/home/amartinezp/Documents/GAIT2/fenoGAIT2'
               dir.script <-'/home/amartinezp/lib'
           },
           stop('*** Unknown machine ***')
           )

    print(paste('cargando clase ***fenoGAIT*** desde',nodename))
    file.dat <- file.path(dir.dat,'gait2.phe.all.RData')
    stopifnot(file.exists(file.dat))
    ## load(file.dat,envir = parent.frame())  
    load(file.dat,envir = .GlobalEnv) 
    
    print('*** Loading functions ***')
    source( file.path(dir.script,'normalizar.R') )
    source( file.path(dir.script,'testeo.corr.biva.R') )
    source( file.path(dir.script,'covars.corr.biva.R') )
    source( file.path(dir.script,'as.correlate.R') )
    
}

