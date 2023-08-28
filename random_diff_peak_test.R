Args=commandArgs(trailingOnly = T)
options(stringsAsFactors = F)
row.wilcox.test=function(common_peak,tumor_sample_in_peak,normal_sample_in_peak,paired.set=F){
  i=1
  n=nrow(common_peak)
  test=c()
  while(i<=n){
    pvalue.paired=wilcox.test(as.numeric(common_peak[i,tumor_sample_in_peak]),
                              as.numeric(common_peak[i,normal_sample_in_peak]),paired = paired.set)
    log2FC.paired=log2(mean(as.numeric(common_peak[i,tumor_sample_in_peak]))/
                         mean(as.numeric(common_peak[i,normal_sample_in_peak])))
    test=rbind(test,c(rownames(common_peak)[i],pvalue.paired$p.value,log2FC.paired))
    i=i+1
  }
  test=data.frame(test)
  colnames(test)=c("Peak.id","P.value.paired","log2FC.paired")
  test$P.value.paired=as.numeric(test$P.value.paired)
  test=cbind(test,FDR.paired=p.adjust(test$P.value.paired,method = "fdr"))
  test$log2FC.paired=as.numeric(test$log2FC.paired)
  return(test)
  
}

i=as.numeric(Args[1])
#n=1000
#diff.common.peak.random_group=c()
#while(i<=n){
rpkm_ip_norm_m.test=read.table("All_sample_m6A_data.txt",header = T,sep="\t")
  set.seed(i)
  random_group_1=sample(colnames(rpkm_ip_norm_m.test), 33, replace = FALSE, prob = NULL)
  set.seed(i+1001)
  random_group_2=sample(setdiff(colnames(rpkm_ip_norm_m.test),random_group_1), 33, replace = FALSE, prob = NULL)
  diff.common.peak.random_group.temp=row.wilcox.test(rpkm_ip_norm_m.test,random_group_1,random_group_2,paired.set = F)
  write.table(diff.common.peak.random_group.temp,paste("Random_diff_peak_test_",i,".txt",sep=""),quote = F,row.names = F,sep="\t")
  # diff.common.peak.random_group=rbind(diff.common.peak.random_group,
  #                                     cbind(diff.common.peak.random_group.temp,type=i))
  print(i)
  #i=i+1
#}
