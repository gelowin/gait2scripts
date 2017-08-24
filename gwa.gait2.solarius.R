
##' Calculate the association of SNPs that come out the first screening with MatrixEQTL in GAIT2. Only continuous traits.
##'
##' This function applies to the output of the function "gwa.gait2.MatrixEQTL"
##' @title Association study with solarius in GAIT2
##' @param chr valid chromosomes are 1-22
##' @param trait Name of the trait to calculate
##' @param normalize if TRUE then the trait is normalized
##' @param inDir if missing "trait.resGWA.fst" is used
##' @param outDir if missing same as inDir
##' @param cores number of cores to use
##' @return A file in the outDir directory
##' @author Angel Martinez-Perez, \email{angelmartinez@protonmail.com}
##' @keywords GAIT2,GWAS
##' @note only works if you have access to the genotypes. Update 20170823
##' @seealso \code{gwa.gait2.MatrixEQTL}
##' @examples
##' gwa.gait2.MatrixEQTL(chr=22,trait='sysVTuP')
##' gwa.gait2.solarius(chr=22,trait='sysVTuP')
##' 
gwa.gait2.solarius <- function(chr, trait, normalize = TRUE, inDir=NA, outDir=NA,cores=60){
    
    stopifnot(require(plyr,warn.conflicts=FALSE,quietly=TRUE))
    stopifnot(require(dplyr,warn.conflicts=FALSE,quietly=TRUE))
    stopifnot(require(solarius,warn.conflicts=FALSE,quietly=TRUE))

    stopifnot(is.logical(normalize))
    stopifnot(cores>0)
    stopifnot(chr %in% 1:22)
    
    ## load the phenotype and solar covariates:
    fuente <- '/home/datasets/GAIT2/scripts/load.fenoGAIT2.R'
    stopifnot(file.exists(fuente))
    source(fuente)
    ## source('/home/amartinezp/lib/load.fenoGAIT2.R')
    load.fenoGAIT2()
    source('/home/datasets/GAIT2/scripts/load.gwaGAIT2.R')
    source('/home/datasets/GAIT2/scripts/load.gwaGAIT2.map.R')
    
    f1$setTraits(unique( c(f1$getTraits('solar'),f1$getTraits('covariates'),trait)))
    df <- f1$getData()
    covlist <- f1$infoTrait(trait)$covariates
    conCovariables <- !is.na(covlist[1])
    
    ## first we collect all the results from "matrixEQTL"
    if(missing(inDir))
        inDir <- paste0(trait,'.resGWA.fst')
    inFile <- file.path(inDir,paste0(trait,'.chr.',chr,'.MatrixEQTL.txt'))
    stopifnot(file.exists(inFile))

    SNPs <- read.table(inFile,header=TRUE,stringsAsFactors = FALSE)$SNP

    ## load the genotype data from this SNPs
    df.snps <- load.gwaGAIT2(chr=chr, snps = c('id',SNPs))
    row.names(df.snps) <- df.snps$id
    df.snps <- df.snps %>% select(-id) %>% as.matrix
    
    df.map.all <- load.gwaGAIT2.map(chr=chr,attributes=c('rsIDref','chromosome','position','info','LOCATION','GENEID','GENESYMBOL','GENEdescript'))
    row.names(df.map.all) <- df.map.all$rsIDref
    df.map.all <- df.map.all[SNPs, ]
    df.map <- df.map.all %>% select(c('rsIDref','chromosome','position'))
    df.map <- rename(df.map, replace=c('rsIDref'='SNP'))



    ## Normalizo los traits para calcular los gwas
    if(normalize)
        df[,trait] <- normalizar(df[,trait])
    

    ## el solar da problemas con variables con nombre largo
    cambiar <- nchar(SNPs) > 19
    if(any(cambiar)){
        originales <- SNPs[cambiar]
        sinteticos <- paste0(paste0(sample(c(letters,LETTERS),6,replace=TRUE),collapse=''),1:sum(cambiar))
        colnames(df.snps) <- plyr::mapvalues(colnames(df.snps), from=originales, to=sinteticos)
        df.map$SNP <- plyr::mapvalues(df.map$SNP, from=originales, to=sinteticos)
    }
    
    ## build the formula
    if(conCovariables)
    {
        formula.sol <- as.formula(paste(trait, '~' ,paste(covlist,collapse=' + ')))
    }else{
        formula.sol <- as.formula(paste(trait, '~ 1'))
    }
    

    ## call solarius
    assoc.solarius <- solarAssoc(formula.sol, data= df ,snpcovdata=df.snps,snpmap=df.map, cores=cores)
    ## y finalmente deshacemos el cambio de nombre:
    resultados <- assoc.solarius$snpf
    if(any(cambiar))
    {
        resultados$SNP <- plyr::mapvalues(resultados$SNP, from=sinteticos, to=originales)
    }
    

    ##  "solarius"
    if(missing(outDir))
        outDir <- inDir
    
    if(!dir.exists(outDir))
    {
        aux <- dir.create(outDir,showWarnings = FALSE)
        if(!aux)
            stop('Not able to create the outDir')
    }
    file.out <- file.path(outDir,paste('solarius',trait,'chr',chr,'RData',sep='.'))

    resultados <- merge(resultados,df.map.all[,-which(colnames(df.map.all) %in% c('chromosome','position'))],by.x='SNP',by.y='rsIDref')
    ## rename a column.name
    resultados <- rename(resultados, replace=c('SNP' = 'rsIDref'))
    resultados <- arrange(resultados,pSNP)

    save(resultados, file=file.out)
}



    





