
##' Calculate a proximate GWAS in GAIT2. Only continuous traits.
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
gwa.gait2.MatrixEQTL <- function(chr, trait, normalize = TRUE, outDir=NA,cores=60,plot=FALSE,kin2='precalculatedGAIT2',pvOutputThreshold=1e-3){

   
    stopifnot(require(fst,warn.conflicts=FALSE,quietly=TRUE))
    stopifnot(require(dplyr,warn.conflicts=FALSE,quietly=TRUE))
    stopifnot(require(solarius,warn.conflicts=FALSE,quietly=TRUE))
    stopifnot(require(MatrixEQTL,warn.conflicts=FALSE,quietly=TRUE))
    
    useModel <- modelLINEAR # modelANOVA or modelLINEAR or modelLINEAR_CROSS


    ## load the data
    source('/home/datasets/GAIT2/scripts/load.gwaGAIT2.R')
    ## source('/home/amartinezp/lib/load.fenoGAIT2.R')
    source('/home/datasets/GAIT2/scripts/load.fenoGAIT2.R')
    load.fenoGAIT2()
    
    f1$setTraits(unique( c(f1$getTraits('solar'),f1$getTraits('covariates'),trait)))
    df <- f1$getData()

########################################################################################
    ## cargamos los genes, ie. esto es el fenotipo cuantitativo que vamos a analizar
    gene <- df %>% select(trait) %>% na.omit %>% as.matrix
    id.1 <- row.names(gene)
########################################################################################
    ## preparamos la carga de los SNPs
    gdat <- load.gwaGAIT2(chr=chr)
    id.2 <- gdat$id
########################################################################################
    ## cargamos las covariables
    covlist <- f1$infoTrait(trait)$covariates
    conCovariables <- !is.na(covlist[1])

    if(conCovariables)
    {
        cov <- df %>% select(covlist) %>% na.omit %>% as.matrix
        id.3 <- row.names(cov)
    }else{
        id.3 <- id.1
    }

#####################################
## hacemos la interseccion de todo ##
#####################################
    id.global <- intersect(intersect(id.1,id.2),id.3)



#####################################
##       ordenamos todo            ##
#####################################
    gene <- gene[id.global, , drop=FALSE]
    row.names(gdat) <- gdat$id
    gdat <- gdat[id.global, , drop=FALSE]
    gdat <- gdat[, -1, drop=FALSE]
    gdat <- t(gdat)
    
########################################################################################
## cargamos los genes, ie. esto es el fenotipo cuantitativo que vamos a analizar
    genes <- SlicedData$new()
    genes$CreateFromMatrix( t(gene) )

########################################################################################
## cargamos las covariables en caso de existir
    if(conCovariables)
    {
        cov <- cov[id.global, ,drop=FALSE]
        covs <- SlicedData$new();
        covs$CreateFromMatrix( t(cov))
    }
########################################################################################
## Cargamos los SNPs
    snps  <-  SlicedData$new();
    snps$CreateFromMatrix( gdat )
    

########################################################################################
    ## Normalizo los traits para calcular los gwas
    if(normalize)
        df[,trait] <- normalizar(df[,trait])
    

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

    ## build the formula
    if(conCovariables)
    {
        formula.sol <- as.formula(paste(trait, '~' ,paste(covlist,collapse=' + ')))
    }else{
        formula.sol <- as.formula(paste(trait, '~ 1'))
    }
    
    poly <- solarPolygenic(formula = formula.sol, data=df,covtest=TRUE,household=TRUE)
    
    h2r <- poly$vcf %>% filter(varcomp=='h2r') %>% select(Var) %>% as.numeric
    varcov <- kin2 * h2r + diag(nrow(kin2))
    varcov <- varcov[as.character(id.global), as.character(id.global)]

########################################################################################
## COMPROBAR QUE LOS DATOS ESTAN EN EL ORDEN CORRECTO:
        
    stopifnot(identical( snps$columnNames, genes$columnNames ))
    if(conCovariables)
        stopifnot(identical( snps$columnNames , covs$columnNames ))
    
########################################################################################
### run
    if(missing(outDir))
        outDir <- file.path(paste0(trait,'.resGWA.fst'))
    dir.create(outDir,showWarnings=FALSE)
    output_file_name <- file.path(outDir,paste0(trait,'.chr.',chr,'.MatrixEQTL.txt'))
    


########################################################################################
        if(conCovariables)
        {
            assoc <- Matrix_eQTL_main(snps = snps, gene = genes, cvrt = covs,
                                      output_file_name = output_file_name,
                                      pvOutputThreshold = pvOutputThreshold,
                                      useModel = useModel,
                                      errorCovariance = varcov,
                                      verbose = TRUE,
                                      pvalue.hist = plot,
                                      min.pv.by.genesnp = FALSE,
                                      noFDRsaveMemory = FALSE)
        }else{
            assoc <- Matrix_eQTL_main(snps = snps, gene = genes,
                                      output_file_name = output_file_name,
                                      pvOutputThreshold = pvOutputThreshold,
                                      useModel = useModel,
                                      errorCovariance = varcov,
                                      verbose = TRUE,
                                      pvalue.hist = plot,
                                      min.pv.by.genesnp = FALSE,
                                      noFDRsaveMemory = FALSE)
        }
}









