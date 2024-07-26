library(tidyverse)
library(data.table)
library(TwoSampleMR)

#### 1.加载暴露
load("./05r2.0.1.1mb.biomaRt.for.res.cis.RData")

#### 2.读取结局
outcome <- read_outcome_data(
  filename = "./finngen_R10_I9_NONISCHCARDMYOP_STRICT.gz",
  snps = decode.pqtl$SNP,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  eaf_col = "af_alt",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval",
  chr_col = "#chrom",
  pos_col = "pos"
)
outcome$outcome <- "NISCM"
outcome$id.outcome <- "NISCM"
outcome$samplesize.outcome <- 342569
outcome <- outcome %>% dplyr::filter(pval.outcome>1e-05)

#### 3.mr
harm.res <- harmonise_data(decode.pqtl,outcome)
mr.res <- mr(harm.res)
mr.res.filter <- mr.res %>% 
  dplyr::filter((nsnp==1&method=="Wald ratio")|(nsnp>1&method=="Inverse variance weighted")) %>% 
  dplyr::filter(pval<0.05)

write.csv(harm.res,file = "01harm.res.all.csv",row.names = F)
write.csv(mr.res,file = "01mr.res.all.csv",row.names = F)
write.csv(mr.res.filter,file = "01mr.res.positive.csv",row.names = F)

#### 4.敏感性分析
mr.res.filter <- mr.res.filter %>% 
  dplyr::mutate(index=str_c(id.exposure,"_",id.outcome))
sensitive.mr <- mr.res %>% 
  dplyr::mutate(index=str_c(id.exposure,"_",id.outcome)) %>% 
  dplyr::filter(index%in%mr.res.filter$index)
sensitive.harm <- harm.res %>% 
  dplyr::mutate(index=str_c(id.exposure,"_",id.outcome)) %>% 
  dplyr::filter(index%in%mr.res.filter$index)
## scatter
mr.scatter.plot <- mr_scatter_plot(sensitive.mr,sensitive.harm)
for (x in 1:length(mr.scatter.plot)) {
  p <- mr.scatter.plot[[x]]+
    theme_bw()+
    theme(legend.position = "top",
          aspect.ratio = 1)
  pdf(str_c("./figure/01scatter.",names(mr.scatter.plot)[x],".pdf"),height = 8,width = 8)
  print(p)
  dev.off()
}
## 多效性
heterogeneity.all <- mr_heterogeneity(sensitive.harm)
write.csv(heterogeneity.all,file = "01heterogeneity.res.csv",row.names = F)
mr.hetero.plot <- mr_funnel_plot(mr_singlesnp(sensitive.harm))
for (x in 1:length(mr.hetero.plot)) {
  p <- mr.hetero.plot[[x]]+
    theme_bw()+
    theme(legend.position = "top",
          aspect.ratio = 1)
  pdf(str_c("./figure/01funnel.",names(mr.hetero.plot)[x],".pdf"),height = 8,width = 8)
  print(p)
  dev.off()
}
## 多效性检测(多效性有问题不建议再用，pvalue>0.05没有多效性)
mr.pleiotropy <- mr_pleiotropy_test(sensitive.harm)
write.csv(mr.pleiotropy,file = "01pleiotropy.res.csv",row.names = F)
## 留一法分析
leaveoneout.res <- mr_leaveoneout(sensitive.harm)
leaveoneout.res.plot <- mr_leaveoneout_plot(leaveoneout.res)
for (x in 1:length(leaveoneout.res.plot)) {
  p <- leaveoneout.res.plot[[x]]+
    theme_bw()+
    theme(legend.position = "top",
          aspect.ratio = 1)
  pdf(str_c("./figure/01leaveOne.",names(leaveoneout.res.plot)[x],".pdf"),height = 8,width = 8)
  print(p)
  dev.off()
}
## 森林图
mr.forest.plot <- mr_forest_plot(mr_singlesnp(sensitive.harm))
for (x in 1:length(mr.forest.plot)) {
  p <- mr.forest.plot[[x]]+
    theme_bw()+
    theme(legend.position = "top",
          aspect.ratio = 1)
  pdf(str_c("./figure/01forest.",names(mr.forest.plot)[x],".pdf"),height = 8,width = 8)
  print(p)
  dev.off()
}
