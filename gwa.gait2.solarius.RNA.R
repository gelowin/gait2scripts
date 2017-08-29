##' Calculate the association of SNPs that come out the first screening with gwa.gait2.MatrixEQTL.RNA in GAIT2
##'
##' This function applies to the output of the function "gwa.gait2.solariusRNA.getData"
##' @title Refine the results with solarius of RNAseq GWAs
##' @param out the output of function "gwa.gait2.solariusRNA.getData"
##' @param chr the chromosome
##' @return data.frame with the results of the gwas, with all gens
##' @author Angel Martinez-Perez, \email{angelmartinez@protonmail.com}
##' @examples
##' out <- gwa.gait2.solariusRNA.getData(chr=22, inDir='gwaRNAseq')
##' gwa.gait2.solarius.RNA(out, chr=22)
##' 
gwa.gait2.solarius.RNA <- function(out , chr)
{
    trait <- out[1,2]
    stopifnot(chr %in% 1:22)
    coresMax <- 60
    stopifnot(class(out)=='data.frame')



    df.RNA <- load.RNAseqGAIT2(gens=c('id',trait))
    
    df <- plyr::join(df.Core,df.RNA, by = 'id', type = "left")
    row.names(df) <- df$id



    SNPs <- out$SNP

    ## load the genotype data from this SNPs
    df.snps <- load.gwaGAIT2(chr=chr, snps = c('id',SNPs))
    row.names(df.snps) <- df.snps$id
    df.snps <- df.snps %>% select(-id) %>% as.matrix
    
    df.map.all <- df.map.all2[SNPs, ]
    df.map <- df.map.all %>% select(c('rsIDref','chromosome','position'))
    df.map <- dplyr::rename(df.map, SNP=rsIDref)



    ## el solar da problemas con variables con nombre largo
    cambiar <- nchar(SNPs) > 19
    if(any(cambiar)){
        originales <- SNPs[cambiar]
        sinteticos <- paste0(paste0(sample(c(letters,LETTERS),6,replace=TRUE),collapse=''),1:sum(cambiar))
        colnames(df.snps) <- plyr::mapvalues(colnames(df.snps), from=originales, to=sinteticos)
        df.map$SNP <- plyr::mapvalues(df.map$SNP, from=originales, to=sinteticos)
    }
    
    ## build the formula
        formula.sol <- as.formula(paste(trait, '~ 1'))

    
    cores <- min(coresMax,ncol(df.snps))
    ## call solarius
    assoc.solarius <- solarAssoc(formula.sol, data= df ,snpcovdata=df.snps,snpmap=df.map, cores=cores,kinship=kin2)
    ## y finalmente deshacemos el cambio de nombre:
    resultados <- assoc.solarius$snpf
    if(any(cambiar))
    {
        resultados$SNP <- plyr::mapvalues(resultados$SNP, from=sinteticos, to=originales)
    }
    

    resultados <- merge(resultados,df.map.all[,-which(colnames(df.map.all) %in% c('chromosome','position'))],by.x='SNP',by.y='rsIDref')
    ## rename a column.name
    resultados <- dplyr::rename(resultados, rsIDref= SNP)
    return(resultados)
}

