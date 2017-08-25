
##' Calculate aproximate GWAS in GAIT2 of all RNAseq data
##'
##' Fast screening of SNPs with MatrixEQTL
##' @title Association study with MatrixEQTL in GAIT2
##' @param chr valid chromosomes are 1-22
##' @param trait Name of the trait to calculate
##' @param normalize if TRUE then the trait is normalized
##' @param outDir folder to put the results
##' @param cores number of cores to use
##' @param plot available option "qqplot"
##' @param kin2 if default then load, else calculated with solarius
##' @param pvOutputThreshold p-value threshold to retain SNPs
##' @return A file in the outDir directory
##' @author Angel Martinez-Perez, \email{angelmartinez@protonmail.com}
##' @keywords GAIT2,GWAS
##' @note only works if you have access to the genotypes. Update 20170824
##' @seealso \code{gwa.gait2.solarius}
##' @examples
##' gwa.gait2.MatrixEQTL(chr=22,trait='sysVTuP'
##' gwa.gait2.solarius(chr=22,trait='sysVTuP')
##'
gwa.gait2.MatrixEQTL.RNAseq <- function(chr, normalize = TRUE, outDir=NA,cores=60,plot=FALSE,kin2='precalculatedGAIT2',pvOutputThreshold=1e-3){

   
    stopifnot(require(fst,warn.conflicts=FALSE,quietly=TRUE))
    stopifnot(require(dplyr,warn.conflicts=FALSE,quietly=TRUE))
    stopifnot(require(solarius,warn.conflicts=FALSE,quietly=TRUE))
    stopifnot(require(MatrixEQTL,warn.conflicts=FALSE,quietly=TRUE))
    
    useModel <- modelLINEAR # modelANOVA or modelLINEAR or modelLINEAR_CROSS


    ## load the data
    source('/home/datasets/GAIT2/scripts/load.gwaGAIT2.R')
    source('/home/datasets/GAIT2/scripts/load.gwaGAIT2.map.R')
    source('/home/datasets/GAIT2/scripts/load.exprGAIT2.R')
    
    ## source('/home/amartinezp/lib/load.fenoGAIT2.R')
  ##   ##     source('/home/datasets/GAIT2/scripts/load.fenoGAIT2.R')
  ##     load.fenoGAIT2()
    
  ##     f1$setTraits(unique( c(f1$getTraits('solar'),f1$getTraits('covariates'),'VT')))
  ##     df <- f1$getData()

########################################################################################
    ## cargamos los genes, ie. esto es el fenotipo cuantitativo que vamos a analizar
    RNAseq <- load.RNAseqGAIT2()
    id.1 <-  RNAseq$id
    row.names(RNAseq) <- id.1
    RNAseq <- RNAseq[,-1]
    ## > dim(RNAseq)
    ## [1]   915 16748
    load.RNAseqGAIT2.map()
    RNAseq.map <- RNAseq.map[,c('ENSEMBL','chr','start','stop')]
    RNAseq.map <- dplyr::rename(RNAseq.map ,geneid = ENSEMBL,left=start,right=stop)

    ## ordeno los datos de RNAseq para que correspondan con su mapa
    RNAseq <- RNAseq[,RNAseq.map$geneid]
    
########################################################################################
    ## preparamos la carga de los SNPs
    gdat <- load.gwaGAIT2(chr=chr)
    id.2 <- gdat$id
########################################################################################
    ## cargamos las covariables
    ## NO COVARIATES
 
#####################################
## hacemos la interseccion de todo ##
#####################################
    id.global <- intersect(id.1,id.2)



#####################################
##       ordenamos todo            ##
#####################################
    RNAseq <- RNAseq[id.global, , drop=FALSE]
    row.names(gdat) <- gdat$id
    gdat <- gdat[id.global, , drop=FALSE]
    gdat <- gdat[, -1, drop=FALSE]
    gdat <- t(gdat)

#####################################
##     cargamos el mapa            ##
#####################################

    snp.pos <- load.gwaGAIT2.map(chr=chr,attributes=c('rsIDref','chromosome','position'))
    snp.pos <- dplyr::rename(snp.pos , snpid = rsIDref,chr = chromosome,pos = position)
    row.names(snp.pos) <- snp.pos$snpid
    snp.pos <- snp.pos[row.names(gdat),] ## ordeno los snps para que esten en el mismo orden que gdat
    
########################################################################################
## cargamos los genes, ie. esto es el fenotipo cuantitativo que vamos a analizar
    genes <- SlicedData$new()
    genes$CreateFromMatrix( t(RNAseq) )

########################################################################################
## Cargamos los SNPs
    snps  <-  SlicedData$new();
    snps$CreateFromMatrix( gdat )
    

########################################################################################
    ## calculate the component of variance
    if(kin2=='precalculatedGAIT2')
    {
        kin2 <- fstread('/home/datasets/GAIT2/phen/kin2.GAIT2.fst',columns=c('id',id.global))
        row.names(kin2) <- kin2$id
        kin2 <- as.matrix(kin2[id.global,-1, drop=FALSE])
    }else{
        row.names(df) <- df$id
        kin2.global <- solarKinship2(df) ## ESTO SE PUEDE CALCULAR SOLO UNA VEZ PARA AHORRAR TIEMPO DE CALCULO
        ## ids.global <- df$id[!is.na(df[,traits[i]])]
        kin2 <- kin2.global[id.global, id.global]
    }

    h2r <- 0.5
    varcov <- kin2 * h2r + diag(nrow(kin2))

########################################################################################
## COMPROBAR QUE LOS DATOS ESTAN EN EL ORDEN CORRECTO:
        
    stopifnot(identical( snps$columnNames, genes$columnNames ))
     
########################################################################################
### create the directories and options of Matrix_eQTL_main
    if(missing(outDir))
        outDir <- file.path(paste0(trait,'.resGWA.fst'))
    dir.create(outDir,showWarnings=FALSE)
    output_file_name <- file.path(outDir,paste0('gwaRNA.chr.',chr,'.MatrixEQTL.txt'))
    output_file_name.cis <- file.path(outDir,paste0('gwaRNA.chr.',chr,'.MatrixEQTL.cis.txt'))
    pvOutputThreshold <- 2e-4
    pvOutputThreshold.cis <- 2e-4
    cisDist <- 1e6
    verbose <- TRUE
########################################################################################
    ## RUN
    
    assoc <- Matrix_eQTL_main(snps = snps,
                              gene = genes,
                              output_file_name = output_file_name,
                              output_file_name.cis = output_file_name.cis,                            
                              pvOutputThreshold = pvOutputThreshold,
                              pvOutputThreshold.cis = pvOutputThreshold.cis,                           
                              useModel = useModel,
                              errorCovariance = varcov,
                              snpspos = snp.pos,
                              genepos = RNAseq.map,
                              cisDist = cisDist,
                              verbose = verbose,
                              pvalue.hist = FALSE,
                              min.pv.by.genesnp = FALSE,
                              noFDRsaveMemory = FALSE)
}

















## una vez calculado el primero los datos del RNAseq no hace falta cargarlos de nuevo


for(chr in 1:21)
{
    
  
########################################################################################
    ## preparamos la carga de los SNPs
    gdat <- load.gwaGAIT2(chr=chr)
#####################################
##       ordenamos todo            ##
#####################################
    row.names(gdat) <- gdat$id
    gdat <- gdat[id.global, , drop=FALSE]
    gdat <- gdat[, -1, drop=FALSE]
    gdat <- t(gdat)

#####################################
##     cargamos el mapa            ##
#####################################

    snp.pos <- load.gwaGAIT2.map(chr=chr,attributes=c('rsIDref','chromosome','position'))
    snp.pos <- dplyr::rename(snp.pos , snpid = rsIDref,chr = chromosome,pos = position)
    row.names(snp.pos) <- snp.pos$snpid
    snp.pos <- snp.pos[row.names(gdat),] ## ordeno los snps para que esten en el mismo orden que gdat
    
########################################################################################
## cargamos los genes, ie. esto es el fenotipo cuantitativo que vamos a analizar
    genes <- SlicedData$new()
    genes$CreateFromMatrix( t(RNAseq) )

########################################################################################
## Cargamos los SNPs
    snps  <-  SlicedData$new();
    snps$CreateFromMatrix( gdat )
    




########################################################################################
## COMPROBAR QUE LOS DATOS ESTAN EN EL ORDEN CORRECTO:
        
    stopifnot(identical( snps$columnNames, genes$columnNames ))
     
########################################################################################
### create the directories and options of Matrix_eQTL_main
    if(missing(outDir))
        outDir <- file.path(paste0('chr',chr,'.resGWA.fst'))
    dir.create(outDir,showWarnings=FALSE)
    output_file_name <- file.path(outDir,paste0('gwaRNA.chr.',chr,'.MatrixEQTL.txt'))
    output_file_name.cis <- file.path(outDir,paste0('gwaRNA.chr.',chr,'.MatrixEQTL.cis.txt'))
    pvOutputThreshold <- 2e-4
    pvOutputThreshold.cis <- 2e-4
    cisDist <- 1e6
    verbose <- TRUE
########################################################################################
    ## RUN
    
    assoc <- Matrix_eQTL_main(snps = snps,
                              gene = genes,
                              output_file_name = output_file_name,
                              output_file_name.cis = output_file_name.cis,                            
                              pvOutputThreshold = pvOutputThreshold,
                              pvOutputThreshold.cis = pvOutputThreshold.cis,                           
                              useModel = useModel,
                              errorCovariance = varcov,
                              snpspos = snp.pos,
                              genepos = RNAseq.map,
                              cisDist = cisDist,
                              verbose = verbose,
                              pvalue.hist = FALSE,
                              min.pv.by.genesnp = FALSE,
                              noFDRsaveMemory = FALSE)


    rm(genes)
    rm(gdat)
    rm(assoc)

}











