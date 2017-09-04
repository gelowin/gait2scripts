#########################3
## genotipos
##########################
## Archivos con los VCF en el directorio:
## datasets@bioinf02/home/datasets/GAIT2/03VCFs/chr1.vcf.gz
## datasets@bioinf02/home/datasets/GAIT2/03VCFs/chr1.vcf.gz.tbi

## YA HECHO, CON UNA SOLA VEZ VALE.
## Comprimimos y construimos el índice de los archivos correctamente con
## tabix
## en R
## for ( i in 1:22)
## {
##
## orden <- paste0( 'bgzip chr',i,'.vcf && tabix -p vcf chr',i,'.vcf.gz')
## system(orden)
## }



#########################3
## Fenotipos
##########################
## Documentation
## http://fastqtl.sourceforge.net/


## Exportamos los valores crudos del RNAseq

## a formato UCSC BED


## ##Chr start end ID id1 id2 id3
## chr1 685395 685396 mi-hsa-45-3p 45 54 63  
## chr2 4685395 4685396 mi-hsa-345-3p 5 4 3

## en la fila de comentario ##chr  los espacios son espacios
## en el resto de filas los espacios son tabulaciones
## mejor poner todo tabulaciones

sourceDir <- '/home/datasets/GAIT2/scripts'

source(file.path(sourceDir, 'load.exprGAIT2.R'))

library(dplyr)

load.RNAseqGAIT2.map()



map <- dplyr::select(RNAseq.map,c(chr,start,stop,ENSEMBL))
miRNA.data$Chr <- sub('chr','',miRNA.data$Chr)

## kkmap %>% count(chr) %>% as.data.frame %>% .[c(gtools::mixedorder(.[,'chr'])),]


RNAdat <- load.RNAseqGAIT2()
row.names(RNAdat) <- RNAdat$id
RNAdat <- select(RNAdat, -id)
RNAdat <- t(RNAdat)


RNAdat2 <- merge(map,RNAdat,by.x='ENSEMBL',by.y='row.names', all.x=FALSE,all.y=TRUE)
RNAdat2 <- RNAdat2[,c(2:4,1,5:ncol(RNAdat2))]
RNAdat2 <- dplyr::arrange(RNAdat2,chr,start)
RNAdat2 <- dplyr::rename(RNAdat2, `##Chr` = chr, end = stop, ID=ENSEMBL)


library(Rsamtools)

salida <- 'RNAdat.bed'
write.table(RNAdat2,file=salida,sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE,na='')
salida.gz <- paste0(salida,'.gz')
bgzip(salida,salida.gz)
indexTabix(salida.gz,format='bed')


dir.create('results')
rutaFastQTL <- '/home/amartinezp/analisis/protocolo.miRNA.GWAS.fastQTL/FastQTL/bin/fastQTL.static'
file.exists(rutaFastQTL)



for(chr in 2:22)
{
    regionMin <- min(RNAdat2$start)
    regionMax <- max(RNAdat2$end)
    
    orden <- paste0(rutaFastQTL,' --vcf /home/datasets/GAIT2/03VCFs/chr',chr,'.vcf.gz --bed ',salida.gz,' --threshold 0.001 --out results/RNAcis.chr',chr,'.fastQTL.gz --log results/RNAcis.chr',chr,'.log --include-samples include.samples.txt --region ',chr,':',regionMin,'-',regionMax)
### FALTA POR LANZARLO PERO TODO ESTA LISTO.
    system(orden)
}

rm(list=ls())
q()
n

FALTA por hacer el X, MT y el Y




