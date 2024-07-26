library(tidyverse)
library(data.table)
library(TwoSampleMR)
library(AnnotationHub)
library(locuscomparer)
source("./coloc_diy.R")
load("./05gene.anno.RData")

#### 1.下载exposure
mr.res <- fread("./01mr.res.positive.csv") %>% 
  dplyr::mutate(adj.pvalue=p.adjust(pval),
                file=exposure) %>% 
  dplyr::filter(adj.pvalue<0.05) %>% 
  dplyr::left_join(fileDf,by="file")
write.csv(mr.res,file = "02exposure.csv",row.names = F)

#### 2.decode注释文件
anno <- fread("./assocvariants.annotated.txt.gz")
anno <- anno %>% dplyr::select(Name,rsids,efa.anno=effectAllele,ota.anno=otherAllele,effectAlleleFreq)
exclude <- fread("./assocvariants.excluded.txt.gz")

#### 3.加载暴露 35559
exposure.list <- list()
for (x in 1:nrow(mr.res)) {
  print(x)
  df <- fread(str_c("./exposure/",mr.res$file[x])) %>% 
    dplyr::mutate(minp=ifelse(Pval==0,NA,Pval),
                  minp=min(na.omit(minp)),
                  Pval=ifelse(Pval==0,minp,Pval)) %>% 
    dplyr::filter(!rsids%in%exclude$rsids,
                  Chrom==mr.res$chromosome_name[x],
                  Pos>=mr.res$start_position[x],
                  Pos<=mr.res$end_position[x],
                  ImpMAF>=0.01)%>% 
    dplyr::filter(!(Chrom=="chr6" & Pos>=26000000 & Pos<=34000000)) %>% 
    dplyr::select(Chrom,Pos,Name,effectAllele,otherAllele,Beta,Pval,SE) %>% 
    dplyr::left_join(anno[,c("Name","rsids","effectAlleleFreq")],by="Name") %>% 
    dplyr::filter(rsids!="",!is.na(rsids),!str_detect(rsids,","),rsids!=".",
                  str_length(effectAllele)==1,str_length(otherAllele)==1) %>% 
    dplyr::mutate(maf=ifelse(effectAlleleFreq<0.5,effectAlleleFreq,(1-effectAlleleFreq))) %>%
    dplyr::select(rsids,chrom=Chrom,position=Pos, pValue=Pval, maf, beta=Beta, se=SE)
  exposure.list <- rlist::list.append(exposure.list,df)
}

save(exposure.list,file = "02exposure.list.RData")

#### 4.加载结局 342569
outcome <- fread("./finngen_R10_I9_NONISCHCARDMYOP_STRICT.gz") %>% 
  dplyr::select(rsids,chrom="#chrom",position=pos, pValue=pval, maf=af_alt, beta=beta, se=sebeta) %>% 
  dplyr::mutate(maf=ifelse(maf<0.5,maf,(1-maf)))

#### 5.coloc分析
coloc.res <- lapply(1:length(exposure.list),function(ep){
  qtl <- exposure.list[[ep]] %>% 
    dplyr::distinct(rsids,.keep_all = T)
  gwas <- outcome %>% 
    dplyr::filter(rsids%in%qtl$rsids) %>% 
    dplyr::select(-position) %>% 
    dplyr::left_join(qtl[,c("rsids","position")],by="rsids") %>% 
    dplyr::select(rsids,chrom,position, pValue, maf, beta, se)
  res <- coloc_diy(gwas,qtl,
                   gwasSampleNum = 342569,qtlSampleNum = 35559,
                   gwasType = "cc",qtlType = "quant",
                   gwasCase = 1754)
  return(res[[1]])
})
coloc.df <- do.call(rbind,coloc.res) %>% 
  dplyr::mutate(decode=mr.res$id.exposure)
write.csv(coloc.df,file = "02coloc.res.csv",row.names = F)

#### 6.绘图
plot.list <- lapply(1:length(exposure.list),function(ep){
  qtl <- exposure.list[[ep]] %>% 
    dplyr::distinct(rsids,.keep_all = T) 
  gwas <- outcome %>% 
    dplyr::filter(rsids%in%qtl$rsids) %>% 
    dplyr::select(-position) %>% 
    dplyr::left_join(qtl[,c("rsids","position")],by="rsids") %>% 
    dplyr::select(rsid=rsids,pval=pValue)
  qtl <- qtl %>% dplyr::select(rsid=rsids,pval=pValue)
  return(list(out=gwas,gene=qtl))
})

for (x in 1:length(plot.list)) {
  plt <- locuscompare(in_fn1 = plot.list[[x]]$out, in_fn2 = plot.list[[x]]$gene, 
                      title1 = "NISCM", title2 = mr.res$id.exposure[x])
  pdf(str_c("./coloc.figure/",mr.res$id.exposure[x],".pdf"))
  print(plt)
  dev.off()
}
