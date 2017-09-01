require(plyr)
require(dplyr)
require(solarius)



chr <- 7
inDir <- "gwaRNAseqSolarius"
sourceDir <- '/home/datasets/GAIT2/scripts'
inFile <- file.path(inDir,paste0('LIST.gwaRNAseq.chr',chr,'.RData'))

stopifnot(file.exists(inFile))
load(inFile)


source(file.path(sourceDir,'load.gwaGAIT2.map.R'))
source(file.path(sourceDir,'load.exprGAIT2.R'))
source(file.path(sourceDir,'load.gwaGAIT2.R'))

source(file.path(sourceDir,'load.fenoGAIT2.R'))

load.fenoGAIT2()
df.Core <- f1$getData()



df.map.all2 <- load.gwaGAIT2.map(chr=chr,attributes=c('rsIDref','chromosome','position','exp_maf','info','LOCATION','GENEID','GENESYMBOL','GENEdescript'))
row.names(mapALL) <- mapALL$rsIDref
source(file.path(sourceDir,'gwa.gait2.solarius.RNA.R'))

kin2 <- fstread('/home/datasets/GAIT2/phen/kin2.GAIT2.fst')
row.names(kin2) <- kin2$id
kin2 <- kin2 %>% select(-id) %>% as.matrix




resultados <- ldply(aux, gwa.gait2.solarius.RNA,chr=chr)

df.map.add <- df.map.all2 %>% select(c('rsIDref','exp_maf','info','LOCATION','GENEID','GENESYMBOL','GENEdescript'))


resultados <- merge(resultados,df.map.add,by='rsIDref',all.x=TRUE,all.y=FALSE)
resultados <- arrange(resultados,pSNP)


outFile <- file.path(inDir,paste0('gwaRNAseq.chr',chr,'.solarius.RData'))

save(resultados, file=outFile)





