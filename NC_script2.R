#basic loading
options(stringsAsFactors = F)
options(scipen = 200)
library(CaSpER)
library(limma)
library(VennDiagram)
library(ggplot2)
library(pheatmap)
library(DESeq2)
library(edgeR)
library(viridis)
library(ConsensusClusterPlus)
library(survival)
library(survminer)
library(reshape2)
library(cluster)
library(ggsci)
library(Hmisc)
library(ROCR)
library(mice)
library(pROC)
library(glmnet)
library(caret)
library(papeR)
library(ggthemes)
library(risksetROC)
library(rms)
library(clusterProfiler)
library(data.table)
library(pathview)
source("ggcdf.R")
source("function.R")
source("pie_plot.R")
library(plyr)
library(ggrepel)
library(RColorBrewer)
library(GenomicFeatures)
library(Biostrings)
library(sampling)

#color
tumor.color="#F7413C"
normal.color="#0005FC"
random.color="#808080"
tumor.color2="firebrick"
normal.color2="navy"
overlap.color="#A83690"

annotation_color1="#00B5BC"
annotation_color2="#F25B51"
annotation_color3="#F1EA0D"
annotation_color4="#058142"
annotation_color5="#F57F25"
annotation_color6="#D359A1"
annotation_color7="#3E6AB2"
annotation_color8="#9CBB5C"
annotation_color9="#8063A0"
annotation_color10="IndianRed1"
annotation_color11="Wheat"
annotation_color12="VioletRed"
annotation_color13="LightGreen"


#Peaks overlap
system("cat analysis/macs2/*.narrowPeak > analysis/macs2/peak_sample_peaks.bed")
system("intersectBed -a analysis/macs2/peak_sample_peaks.bed -b analysis/macs2/merged_Peak.bed -wa -wb > analysis/macs2/all_sample_in_merged_Peak.bed")
system("intersectBed -a analysis/macs2/all_sample_in_merged_Peak.bed -b analysis/metpeak/metpeak_merged_Peak.bed -wa -u > analysis/macs2/all_sample_in_macs2_and_metpeak.bed")

peak_in_macs2=system("awk -F \"\t\" '{print $14}' analysis/macs2/all_sample_in_merged_Peak.bed|sort|uniq|wc -l",intern = T)
peak_in_metpeak=system("awk -F \"\t\" '{print $3}' analysis/metpeak/metpeak_merged_Peak.bed|sort|uniq|wc -l",intern = T)

peak_in_macs2=as.numeric(peak_in_macs2)
peak_in_metpeak=as.numeric(peak_in_metpeak)

allpeak=read.table("analysis/macs2/all_sample_in_macs2_and_metpeak.bed")
allpeak=cbind(allpeak,all=strsplit2(strsplit2(allpeak$V4,split="/")[,2],split="[.]")[,1])
allpeak=allpeak[!duplicated(paste(allpeak$V14,allpeak$all,sep="_")),]
overlap_peak_bed=unique(allpeak$V14)
overlap_peak_bed=cbind(strsplit2(overlap_peak_bed,split=":|-"),overlap_peak_bed)
write.table(overlap_peak_bed,"/data/xingyang/m6A_zhengjian/analysis/annotation/overlap_peak.bed",sep="\t",quote = F,row.names = F,col.names = F)


#Overlap peak annotation
setwd("/data/xingyang/m6A_zhengjian/analysis/annotation/")
system("perl m6A_annotate_forGTF_xingyang_v2.pl /data/database/hg38/GENCODE/gencode.v25.annotation.gtf /data/xingyang/m6A_zhengjian/analysis/annotation/overlap_peak.bed /data/xingyang/m6A_zhengjian/analysis/annotation/new_annotation/overlap")
setwd("/data/xingyang/m6A_zhengjian/")
overlap.anno=read.table("analysis/annotation/new_annotation/overlap.anno.txt",header = F,sep="\t")
colnames(overlap.anno)=c("Chr","Start","End","Peak.id","Chr.gene","Start.gene","End.gene","Transcript.id",
                 "Nouse","Strand","Gene","Gene.type","Gene.site","Peak.position","Ensembl.gene.id","Level.1.gene.type",
                 "Level.2.gene.type")
rownames(overlap.anno)=overlap.anno$Peak.id

overlap.unanno=read.table("analysis/annotation/new_annotation/overlap.unanno.txt",header = F,sep="\t")
colnames(overlap.unanno)[4]="Peak.id"
overlap.unanno=merge(overlap.unanno,overlap.anno,by="Peak.id",all.x=T)
overlap.unanno=overlap.unanno[,colnames(overlap.anno)]

overlap.anno=rbind(overlap.anno,overlap.unanno)

overlap.anno[which(is.na(overlap.anno$Level.2.gene.type)),"Level.2.gene.type"]="Unknown"
rownames(overlap.anno)=overlap.anno$Peak.id

overlap.anno[which(overlap.anno$Gene.type=="coding" & overlap.anno$Level.2.gene.type!="mRNA"),"Level.2.gene.type"]="mRNA"
overlap.anno[which(overlap.anno$Gene.type=="coding" & overlap.anno$Level.2.gene.type!="mRNA"),"Level.1.gene.type"]="mRNA"

rm(overlap.unanno)


#Find m6Am 5'UTR peaks
anno_5UTR=overlap.anno[which(overlap.anno$Gene.site=="5UTR"),]

gtf.temp=fread("/data/database/hg38/GENCODE/gencode.v25.annotation.gtf",sep="\t",skip = 5,data.table = F)
gtf.temp=cbind(gtf.temp,Transcript.id=strsplit2(strsplit2(gtf.temp$V9,split = "transcript_id ")[,2],split = ";")[,1])
gtf.temp$Transcript.id=gsub("\"","",gtf.temp$Transcript.id)
gtf.temp$Gene.site=gtf.temp$V3
gtf.temp$UTR.temp=strsplit2(gtf.temp[,9],split = ";")[,9]
gtf.temp[which(gtf.temp$Gene.site=="UTR" & gtf.temp$UTR.temp==" exon_number 1"),"Gene.site"]="5UTR"
gtf.temp[which(gtf.temp$Gene.site=="UTR" & gtf.temp$UTR.temp!=" exon_number 1"),"Gene.site"]="3UTR"
gtf.temp$UTR.temp=NULL

write.table(strsplit2(anno_5UTR$Peak.id,split=":|-"),"UTR5.peak.bed",row.names = F,col.names = F,quote = F,sep="\t")
system("fastaFromBed -fi /data/database/hg38/genome.fa -bed UTR5.peak.bed -fo UTR5.peak.fa")
system("/data/software/homer/bin/homer2 find -i UTR5.peak.fa -m /data/xingyang/m6A_zhengjian/BCA.motif -p 50 > /data/xingyang/m6A_zhengjian/analysis/BCA_peak_offset.txt")
BCA_in_5UTR_offset=read.table("analysis/BCA_peak_offset.txt",header=F)

anno_5UTR=overlap.anno[unique(BCA_in_5UTR_offset$V1),]
anno_5UTR=merge(gtf.temp,anno_5UTR,by="Transcript.id",all.y=T)
anno_5UTR=anno_5UTR[which(anno_5UTR$V3=="UTR"),]
anno_5UTR$temp.start=anno_5UTR$V4-anno_5UTR$Start
anno_5UTR$temp.end=anno_5UTR$V5-anno_5UTR$Start
anno_5UTR[which(anno_5UTR$temp.start>0),"temp.start"]=1
anno_5UTR[which(anno_5UTR$temp.start<0),"temp.start"]=(-1)
anno_5UTR[which(anno_5UTR$temp.end>0),"temp.end"]=1
anno_5UTR[which(anno_5UTR$temp.end<0),"temp.end"]=(-1)
anno_5UTR=anno_5UTR[which((anno_5UTR$temp.start*anno_5UTR$temp.end)<=0),]
anno_5UTR.bed=cbind(anno_5UTR$V1,anno_5UTR$V4,anno_5UTR$V5,anno_5UTR$Peak.id,".",anno_5UTR$V7)
write.table(anno_5UTR.bed,"UTR5.peak.bed",row.names = F,col.names = F,quote = F,sep="\t")
system("fastaFromBed -fi /data/database/hg38/genome.fa -bed UTR5.peak.bed -s -name -fo UTR5.peak.fa")

temp.utr5=read.table("UTR5.peak.fa",sep="\n")
temp.utr5=cbind(temp.utr5,substr(temp.utr5[,1],1,1))

i=2
n=nrow(temp.utr5)
temp.utr5=cbind(temp.utr5,type=NA)
while(i<=n){
  if(temp.utr5[i,2]=="A"){
    temp.utr5[c(i-1,i),"type"]="m6Am"
  }
  i=i+2
}
m6Am=na.omit(temp.utr5)
m6Am=m6Am[grep(">",m6Am$V1),]
m6Am=gsub(">","",m6Am$V1)
m6Am=strsplit2(m6Am,split="[(]")[,1]

overlap.anno=overlap.anno[setdiff(rownames(overlap.anno),m6Am),]

#Peaks filtered by sample
fuck_peak=table(allpeak$V14)
fuck_peak=fuck_peak[which(fuck_peak>4.9)]
fuck_peak=intersect(names(fuck_peak),overlap.anno$Peak.id)
fuck_peak=cbind(strsplit2(fuck_peak,split=":|-"),fuck_peak)

write.table(fuck_peak,"/data/xingyang/m6A_zhengjian/analysis/macs2/fix_merged_Peak.bed",quote=F,row.names = F,col.names = F,sep="\t")
system("cat /data/xingyang/m6A_zhengjian/analysis/macs2/*T.macs2_peaks.narrowPeak > /data/xingyang/m6A_zhengjian/analysis/macs2/all.tumor.sample.macs2.peak.bed")
system("cat /data/xingyang/m6A_zhengjian/analysis/macs2/*N.macs2_peaks.narrowPeak > /data/xingyang/m6A_zhengjian/analysis/macs2/all.normal.sample.macs2.peak.bed")
system("intersectBed -a /data/xingyang/m6A_zhengjian/analysis/macs2/all.tumor.sample.macs2.peak.bed -b /data/xingyang/m6A_zhengjian/analysis/macs2/fix_merged_Peak.bed -wa -wb > /data/xingyang/m6A_zhengjian/analysis/annotation/Tumor_vs_Normal/fix_Tumor_all_overlap.bed")
system("intersectBed -a /data/xingyang/m6A_zhengjian/analysis/macs2/all.normal.sample.macs2.peak.bed -b /data/xingyang/m6A_zhengjian/analysis/macs2/fix_merged_Peak.bed -wa -wb > /data/xingyang/m6A_zhengjian/analysis/annotation/Tumor_vs_Normal/fix_Normal_all_overlap.bed")


#peak annotation
setwd("/data/xingyang/m6A_zhengjian/analysis/annotation/")
system("perl m6A_annotate_forGTF_xingyang_v2.pl /data/database/hg38/GENCODE/gencode.v25.annotation.gtf /data/xingyang/m6A_zhengjian/analysis/macs2/fix_merged_Peak.bed /data/xingyang/m6A_zhengjian/analysis/annotation/new_annotation/fix")
setwd("/data/xingyang/m6A_zhengjian/")
anno=read.table("analysis/annotation/new_annotation/fix.anno.txt",header = F,sep="\t")
colnames(anno)=c("Chr","Start","End","Peak.id","Chr.gene","Start.gene","End.gene","Transcript.id",
                 "Nouse","Strand","Gene","Gene.type","Gene.site","Peak.position","Ensembl.gene.id","Level.1.gene.type",
                 "Level.2.gene.type")
rownames(anno)=anno$Peak.id

unanno=read.table("analysis/annotation/new_annotation/fix.unanno.txt",header = F,sep="\t")
colnames(unanno)[4]="Peak.id"
unanno=merge(unanno,anno,by="Peak.id",all.x=T)
unanno=unanno[,colnames(anno)]

anno=rbind(anno,unanno)

anno[which(is.na(anno$Level.2.gene.type)),"Level.2.gene.type"]="Unknown"
rownames(anno)=anno$Peak.id

anno[which(anno$Gene.type=="coding" & anno$Level.2.gene.type!="mRNA"),"Level.2.gene.type"]="mRNA"
anno[which(anno$Gene.type=="coding" & anno$Level.2.gene.type!="mRNA"),"Level.1.gene.type"]="mRNA"

rm(unanno)
write.table(anno,"anno.txt",quote = F,
            sep="\t")

#Find RRACH motif in peaks
system("fastaFromBed -fi /data/database/hg38/genome.fa -bed /data/xingyang/m6A_zhengjian/analysis/macs2/fix_merged_Peak.bed -fo analysis/macs2/two_method_overlap_peak.fa")
system("/data/software/homer/bin/homer2 find -i /data/xingyang/m6A_zhengjian/analysis/macs2/two_method_overlap_peak.fa -m /data/xingyang/m6A_zhengjian/GGACH.motif -p 50 > /data/xingyang/m6A_zhengjian/analysis/GGACH.in.peak.txt")
GGACH_in_peak=read.table("analysis/GGACH.in.peak.txt",header=F)
GGACH_in_peak[,1]=gsub("\\(NA\\)","",GGACH_in_peak[,1])


overlap.peak.frep=data.frame(table(allpeak$V14))
colnames(overlap.peak.frep)[1]="Peak.id"
overlap.peak.frep=merge(data.frame(Peak.id=overlap.anno$Peak.id),overlap.peak.frep,all.x=T,by="Peak.id")
overlap.peak.frep=cbind(strsplit2(overlap.peak.frep$Peak.id,split=":|-"),overlap.peak.frep)
overlap.peak.frep$`1`=gsub("chr","hs",overlap.peak.frep$`1`)
overlap.peak.frep$Peak.id=NULL
#overlap.peak.frep$Freq=paste("value=",overlap.peak.frep$Freq,sep="")
write.table(overlap.peak.frep,"analysis/new_figure/peak_freq_circos/overlap_peak_freq",sep=" ",row.names = F,col.names = F,quote = F)


#Stop codon annotation
ggtopo(all_topology)+scale_color_manual(values = c(normal.color))+guides(color=F)
ggsave("analysis/new_figure/figure1_allpeak_topology.pdf",width=5,heigh=4,dpi=500)

stop_codon_regoins_anno=coding_anno(annotation_table = anno,gtf.temp = gtf.temp)
pie_plot(stop_codon_regoins_anno[,"Gene.site"])+
  scale_fill_manual(values=c(annotation_color2,annotation_color4,annotation_color1,annotation_color9,annotation_color3))
ggsave("analysis/new_figure/All_coding_annotation_with_stop_codon.pdf",width = 5,height = 5)

coding_type_counting=gtf.temp[grep("protein_coding",gtf.temp[,9]),]
exon.temp=table(paste(coding_type_counting$Transcript.id,"-",coding_type_counting$Gene.site,sep=""))
exon.temp=intersect(grep("exon",names(exon.temp)),which(exon.temp==1))
coding_type_counting=table(coding_type_counting$Gene.site)
coding_type_counting["intron"]=coding_type_counting["exon"]
coding_type_counting["intron"]=coding_type_counting["intron"]-length(exon.temp)
coding_type_counting["intron"]=coding_type_counting["intron"]-coding_type_counting["transcript"]

tissue_peak_type_counting=table(stop_codon_regoins_anno$Gene.site)
coding_type_counting=coding_type_counting[names(tissue_peak_type_counting)]

i=1
n=length(tissue_peak_type_counting)
peak.enrichment=c()
while(i<=n){
  enrichment.temp=fisher.test(matrix(c(tissue_peak_type_counting[i],sum(tissue_peak_type_counting[-i]),
                                       coding_type_counting[i],sum(coding_type_counting[-i])),2,2))
  peak.enrichment=rbind(peak.enrichment,c(type=names(tissue_peak_type_counting)[i],score=as.numeric(enrichment.temp$estimate)))
  i=i+1
}
peak.enrichment=data.frame(peak.enrichment)
peak.enrichment$score=as.numeric(peak.enrichment$score)
peak.enrichment$type=ordered(peak.enrichment$type,levels=c("5UTR","CDS","intron","stop_codon","3UTR"))
ggplot(peak.enrichment)+
  geom_bar(aes(x=peak.enrichment$type,y=peak.enrichment$score,fill=peak.enrichment$type),stat="identity")+
  scale_fill_manual(values=c(annotation_color3,annotation_color2,annotation_color4,annotation_color9,annotation_color1),name="")+
  theme_classic()+scale_y_continuous(expand = c(0,0))+xlab("m6A peaks")+ylab("Peaks enrichment score")
ggsave("analysis/new_figure/All_coding_annotation_enrichment_bar.pdf",width = 5,height = 5)


#identified sample in RNA-seq
tumor_sample_in_gene=colnames(gene_count_data_t)
normal_sample_in_gene=colnames(gene_count_data_n)

#load peak enrichment matrix
rpkm_data=read.table("analysis/merged_Peak_abundance_normbytmn.txt",header=T,row.names = 1)
rpkm_ip_data=rpkm_data[,grep("ip",colnames(rpkm_data))]
rpkm_input_data=rpkm_data[,grep("input",colnames(rpkm_data))]
rpkm_ip_norm=(rpkm_ip_data)/(rpkm_input_data+1)
colnames(rpkm_ip_norm)=gsub(".ip","",colnames(rpkm_ip_norm))
rpkm_ip_norm_n=rpkm_ip_norm[,grep("N",colnames(rpkm_ip_norm))]
rpkm_ip_norm_t=rpkm_ip_norm[,grep("T",colnames(rpkm_ip_norm))]
rpkm_ip_norm_m=cbind(rpkm_ip_norm_n,rpkm_ip_norm_t)

#identified sample in m6A-seq
tumor_sample_in_peak=colnames(rpkm_ip_norm_t)
normal_sample_in_peak=colnames(rpkm_ip_norm_n)

#Paired sample peaks
paired_sample=read.table("analysis/enyzme/samplename_ID.txt")
colnames(paired_sample)=c("normal","tumor")

paired_sample_normal=paste("X",paired_sample$normal,sep="")
paired_sample_tumor=paste("X",paired_sample$tumor,sep="")

#diff peak analysis
diff.common.peak=row.wilcox.test(rpkm_ip_norm_m,paired_sample_tumor,paired_sample_normal,paired.set = T)
peak_for_write=cbind(diff.common.peak,anno[diff.common.peak$Peak.id,])
write.csv(peak_for_write,"analysis/new_figure/peak_for_write.csv",row.names = F)

rownames(diff.common.peak)=diff.common.peak$Peak.id
diff.common.peak=diff.common.peak[anno[which(!is.na(anno$Gene)),"Peak.id"],]

significan.diff.common.peak=diff.common.peak[which(diff.common.peak$FDR.paired<0.1),]

significan.diff.common.peak=significan.diff.common.peak[order(significan.diff.common.peak$log2FC.paired,decreasing = T),]
hyper_peaks=significan.diff.common.peak[which(significan.diff.common.peak$log2FC.paired>0),]
hypo_peaks=significan.diff.common.peak[which(significan.diff.common.peak$log2FC.paired<0),]

diff_peak_annotation=significan.diff.common.peak
diff_peak_annotation=cbind(diff_peak_annotation,Peak_type="Hyper")
diff_peak_annotation[intersect(diff_peak_annotation$Peak.id,hypo_peaks$Peak.id),"Peak_type"]="Hypo"
diff_peak_annotation=diff_peak_annotation[order(diff_peak_annotation$Peak_type,
                                                diff_peak_annotation$log2FC.paired,decreasing = F),]
rownames(diff_peak_annotation)=diff_peak_annotation$Peak.id
diff_peak_annotation$Peak.id=NULL
diff_peak_annotation$log2FC.paired=NULL
diff_peak_annotation$P.value.paired=NULL
diff_peak_annotation$FDR.paired=NULL

pheatmap(rpkm_ip_norm_m[rownames(diff_peak_annotation),
                        c(normal_sample_in_peak,paired_sample_tumor,single_sample_tumor)],scale = "row",
         color=colorRampPalette(c(rep("navy"),"white",rep("firebrick3")))(50),
         annotation_col = tumor_normal_annotation,cluster_rows = F,
         cluster_cols = F,
         annotation_row = diff_peak_annotation,
         annotation_colors = list(Sample_type=c(Tumor_paired=tumor.color,Normal_paired=normal.color,Tumor_single=overlap.color),
                                  Peak_type=c(Hyper=tumor.color2,Hypo=normal.color2)),
         show_rownames = F,
         show_colnames = F,
         border_color = NULL,
         filename = "analysis/new_figure/heatmap_m6A_scale.pdf",
         width = 10,height = 8)
dev.off()

#PDAC tissue differential methylated m6A validation
qPCR_data=read.table("peak_PCR",header = T)
validation_point=cbind(qPCR_data,m6A_FC=diff.common.peak[qPCR_data$m6A_peak,"log2FC.paired"])
validation_point$m6A_Fold_change=log2(validation_point$m6A_Fold_change)

cor.temp=cor.test(validation_point$m6A_Fold_change,validation_point$m6A_FC,method="spearman")

ggplot(validation_point)+
  geom_point(aes(x=validation_point$m6A_Fold_change,y=validation_point$m6A_FC),
             color=annotation_color7,size=4)+
  geom_smooth(aes(x=validation_point$m6A_Fold_change,y=validation_point$m6A_FC),
              method = "lm",se = F,linetype="dashed",color="black",size=1.2)+
  geom_label(aes(x=2,y=-0.5,label=paste("Correlation: ",round(cor.temp$estimate,3),"\n",
                                        "P value: ",round(cor.temp$p.value,6))))+
  theme_classic()+xlab("log2(fold change) in PDAC tissue identified by m6A specific PCR")+
  ylab("log2(fold change) in PDAC tissue identified by m6A-seq")
ggsave("analysis/new_figure/PDAC_m6A_validation_point.pdf",width = 5.2,height = 5)

#random sampling test
i=1
n=1000
run.script=c()
while(i<=n){
  run.script=append(run.script,
                   paste("Rscript random_diff_peak_test.R ",i,sep=""))
  i=i+1
}
write.table(run.script,"run_random_diff_peak_test.sh",quote = F,row.names = F,col.names = F)

i=1
n=1000
random_diff_peak=c()
while(i<=n){
  data.temp=fread(paste("Random_diff_peak_test_",i,".txt",sep=""),sep="\t",header=T)
  data.temp=data.frame(data.temp)
  random_diff_peak=rbind(random_diff_peak,cbind(data.temp,type=i))
  print(i)
  i=i+1
}

sig.random_diff_peak=random_diff_peak[which(random_diff_peak$FDR.paired<0.1),]
sig.random_diff_peak.temp=sig.random_diff_peak
sig.random_diff_peak.temp$type=factor(sig.random_diff_peak.temp$type,levels=c(4500:5500))
random_diff_count=table(table(sig.random_diff_peak.temp$type))
P.value.temp=sum(random_diff_count[which(as.numeric(names(random_diff_count))>nrow(significan.diff.common.peak))])/sum(random_diff_count)
random_diff_count=data.frame(random_diff_count)
random_diff_count$Var1=as.character(random_diff_count$Var1)
random_diff_count$type="Random"
random_diff_count=rbind(random_diff_count,c(nrow(significan.diff.common.peak),1,"Sig"))
random_diff_count$Freq=as.numeric(random_diff_count$Freq)
random_diff_count=random_diff_count[order(as.numeric(random_diff_count$Var1)),]
random_diff_count$Var1=ordered(random_diff_count$Var1,levels=random_diff_count$Var1)

i=1
n=nrow(random_diff_count)
mean.temp=c()
while(i<=n){
  mean.temp=append(mean.temp,rep(as.numeric(as.character(random_diff_count[i,1])),
                                 as.numeric(as.character(random_diff_count[i,2]))))
  i=i+1
}
mean.temp=data.frame(mean.temp)

ggplot(mean.temp)+geom_density(aes(x=mean.temp),color="skyblue3",size=1.5)+
  geom_vline(xintercept = mean(mean.temp$mean.temp),linetype="dashed",color="Darkgreen",size=1.5)+
  geom_vline(xintercept = 288,linetype="dashed",color="firebrick",size=1.5)+
  theme_classic()+
  xlab("Differentially methylated m6A peaks number")+
  ylab("Density")+scale_y_continuous(expand = c(0,0))
ggsave("analysis/new_figure/Diff_peak_permutation_test.pdf",width = 4,height = 3)


#annotation for peaks
pie_plot(anno[,"Level.2.gene.type"])+
  scale_fill_manual(values=c(annotation_color7,annotation_color2,annotation_color1,annotation_color3,
                             random.color))
ggsave("analysis/new_figure/Annotation_all_gene.pdf",width = 5,heigh=4,dpi=500)

#Precentage of peak contain GGACH motif
peak.with.GGACH.motif=unique(GGACH_in_peak$V1)
peak.wothout.GGACH.motif=setdiff(anno$Peak.id,strsplit2(peak.with.GGACH.motif,split = "[(]")[,1])
peak.with.GGACH.motif=strsplit2(peak.with.GGACH.motif,split = "[(]")[,1]

all_regoins=na.omit(anno)
all_regoins[which(all_regoins$Gene.site!="intron"),"Gene.site"]="Exon"
pie_plot(all_regoins[,"Gene.site"])+
  scale_fill_manual(values=c(annotation_color4,annotation_color5))
ggsave("analysis/new_figure/Annotation_all_regoins.pdf",width = 5,heigh=4,dpi=500)

all_regoins[intersect(peak.with.GGACH.motif,rownames(all_regoins)),"type"]="Contain_GGACH"
all_regoins[which(is.na(all_regoins$type)),"type"]="Withou_GGACH"
exon_intron_contain_motif=dcast(all_regoins,Gene.site~type)

ggplot()+
  geom_histogram(data=all_regoins,aes(all_regoins$Gene.site,fill=all_regoins$type),position="fill",stat = "count")+
  #geom_text(data=diff.peak.type.label,aes(x=x,y=y,label=label),color="white",size=5)+
  theme_classic(base_size = 15)+scale_fill_manual(values = c(annotation_color7,annotation_color2))+
  scale_y_continuous(expand = c(0,0))+coord_flip()+xlab("")+ylab("Precentage")+theme(legend.title = element_blank())
ggsave("analysis/new_figure/exon_intron_contain_GGACH.pdf",width = 8,height = 4.5)


#distanc between exon or intron and m6A sites
intron_distribution=site_distance(anno[which(anno$Gene.site=="intron"),]
                                             ,exon_intron_splicing.anno,gtf.temp)
exon_distribution=site_distance(anno[which(anno$Gene.site=="exon"),]
                                           ,exon_intron_splicing.anno,gtf.temp)

plot_distribution=rbind(cbind(exon_distribution,type="Exon"),
                                cbind(intron_distribution,type="Intron"))

plot_distribution[which(plot_distribution$type=="Exon"&plot_distribution$Distance>= 0),"type"]="Exon_down"
plot_distribution[which(plot_distribution$type=="Exon"&plot_distribution$Distance<= 0),"type"]="Exon_up"
plot_distribution[which(plot_distribution$type=="Intron"&plot_distribution$Distance>= 0),"type"]="Intron_up"
plot_distribution[which(plot_distribution$type=="Intron"&plot_distribution$Distance<= 0),"type"]="Intron_down"

plot_distribution[which(plot_distribution$type=="Exon_down"),"Distance"]=4000+
  plot_distribution[which(plot_distribution$type=="Exon_down"),"Distance"]
plot_distribution[which(plot_distribution$type=="Intron_down"),"Distance"]=4000+
  plot_distribution[which(plot_distribution$type=="Intron_down"),"Distance"]

plot_distribution[which(plot_distribution$type=="Exon_up"&
                                  (plot_distribution$Distance< (-2000)|plot_distribution$Distance>0)),"type"]=NA
plot_distribution[which(plot_distribution$type=="Intron_up"&
                                  (plot_distribution$Distance< (0)|plot_distribution$Distance>2000)),"type"]=NA
plot_distribution[which(plot_distribution$type=="Intron_down"&
                                  (plot_distribution$Distance< (2000)|plot_distribution$Distance>4000)),"type"]=NA
plot_distribution[which(plot_distribution$type=="Exon_down"&
                                  (plot_distribution$Distance< (4000)|plot_distribution$Distance>6000)),"type"]=NA
#plot_distribution[is.na(plot_distribution$Distance),"Distance"]=2000
plot_distribution=plot_distribution[!is.na(plot_distribution$type),]

ggplot(plot_distribution)+
  geom_line(aes(x=plot_distribution$Distance,color=plot_distribution$type),
            stat = "density",size=1.2,trim=T)+
  theme_classic()+
  #xlim(-2000,6000)+
  scale_color_manual(values=c(normal.color,normal.color,tumor.color,tumor.color),name="")+
  xlab("Distance of m6A site from exon splicing")+ylab("Density")+
  scale_x_continuous(breaks = c(-2000,0,2000,4000,6000),labels = c("2000","0","2000","0","2000"))
ggsave("analysis/new_figure/m6A_sites_vs_exon_splicing.pdf",width = 7,height = 4.5)

#hyper
pie_plot(anno[hyper_peaks$Peak.id, "Level.2.gene.type"])+
  scale_fill_manual(values=c(annotation_color7,annotation_color2,annotation_color1,annotation_color3,random.color))
ggsave("analysis/new_figure/Annotation_hyper_gene.pdf",width = 5,heigh=4,dpi=500)

pie_plot(anno[which(anno$Gene.type=="coding"),][hyper_peaks$Peak.id,"Gene.site"])+
  scale_fill_manual(values=c(annotation_color2,annotation_color4,annotation_color1,annotation_color3))
ggsave("analysis/new_figure/Annotation_hyper_coding.pdf",width = 5,heigh=4,dpi=500)

#Stop codon annotation-hyper
hyper_anno=anno[hyper_peaks$Peak.id,]
stop_codon_regoins_anno=coding_anno(annotation_table = hyper_anno,gtf.temp = gtf.temp)
pie_plot(stop_codon_regoins_anno[,"Gene.site"])+
  scale_fill_manual(values=c(annotation_color2,annotation_color4,annotation_color1,annotation_color9,annotation_color3))
ggsave("analysis/new_figure/Hyper_coding_annotation_with_stop_codon.pdf",width = 5,height = 5)

tissue_peak_type_counting=table(stop_codon_regoins_anno$Gene.site)
coding_type_counting=coding_type_counting[names(tissue_peak_type_counting)]

i=1
n=length(tissue_peak_type_counting)
peak.enrichment=c()
while(i<=n){
  enrichment.temp=fisher.test(matrix(c(tissue_peak_type_counting[i],sum(tissue_peak_type_counting[-i]),
                                       coding_type_counting[i],sum(coding_type_counting[-i])),2,2))
  peak.enrichment=rbind(peak.enrichment,c(type=names(tissue_peak_type_counting)[i],score=as.numeric(enrichment.temp$estimate)))
  i=i+1
}
peak.enrichment=data.frame(peak.enrichment)
peak.enrichment$score=as.numeric(peak.enrichment$score)
peak.enrichment$type=ordered(peak.enrichment$type,levels=c("5UTR","CDS","intron","stop_codon","3UTR"))
ggplot(peak.enrichment)+
  geom_bar(aes(x=peak.enrichment$type,y=peak.enrichment$score,fill=peak.enrichment$type),stat="identity")+
  scale_fill_manual(values=c(annotation_color3,annotation_color2,annotation_color4,annotation_color9,annotation_color1),name="")+
  theme_classic()+scale_y_continuous(expand = c(0,0))+xlab("m6A peaks")+ylab("Peaks enrichment score")
ggsave("analysis/new_figure/Hyper_coding_annotation_enrichment_bar.pdf",width = 6,height = 5)



#hypo
pie_plot(anno[hypo_peaks$Peak.id,"Level.2.gene.type"])+
  scale_fill_manual(values=c(annotation_color7,annotation_color2,annotation_color1,annotation_color3,random.color))
ggsave("analysis/new_figure/Annotation_hypo_gene.pdf",width = 5,heigh=4,dpi=500)

pie_plot(anno[which(anno$Gene.type=="coding"),][hypo_peaks$Peak.id,"Gene.site"])+
  scale_fill_manual(values=c(annotation_color2,annotation_color1,annotation_color4,annotation_color3))
ggsave("analysis/new_figure/Annotation_hypo_coding.pdf",width = 5,heigh=4,dpi=500)

#Stop codon annotation-hypo
hypo_anno=anno[hypo_peaks$Peak.id,]
stop_codon_regoins_anno=coding_anno(annotation_table = hypo_anno,gtf.temp = gtf.temp)
pie_plot(stop_codon_regoins_anno[,"Gene.site"])+
  scale_fill_manual(values=c(annotation_color2,annotation_color4,annotation_color1,annotation_color9,annotation_color3))
ggsave("analysis/new_figure/hypo_coding_annotation_with_stop_codon.pdf",width = 5,height = 5)


tissue_peak_type_counting=table(stop_codon_regoins_anno$Gene.site)
coding_type_counting=coding_type_counting[names(tissue_peak_type_counting)]

i=1
n=length(tissue_peak_type_counting)
peak.enrichment=c()
while(i<=n){
  enrichment.temp=fisher.test(matrix(c(tissue_peak_type_counting[i],sum(tissue_peak_type_counting[-i]),
                                       coding_type_counting[i],sum(coding_type_counting[-i])),2,2))
  peak.enrichment=rbind(peak.enrichment,c(type=names(tissue_peak_type_counting)[i],score=as.numeric(enrichment.temp$estimate)))
  i=i+1
}
peak.enrichment=data.frame(peak.enrichment)
peak.enrichment$score=as.numeric(peak.enrichment$score)
peak.enrichment$type=ordered(peak.enrichment$type,levels=c("5UTR","CDS","intron","stop_codon","3UTR"))
ggplot(peak.enrichment)+
  geom_bar(aes(x=peak.enrichment$type,y=peak.enrichment$score,fill=peak.enrichment$type),stat="identity")+
  scale_fill_manual(values=c(annotation_color3,annotation_color2,annotation_color4,annotation_color9,annotation_color1),name="")+
  theme_classic()+scale_y_continuous(expand = c(0,0))+xlab("m6A peaks")+ylab("Peaks enrichment score")
ggsave("analysis/new_figure/hypo_coding_annotation_enrichment_bar.pdf",width = 6,height = 5)

#hyper pathway and hypo pathway (figure 1F)
hyper.eg <- bitr(unique(na.omit(anno[hyper_peaks$Peak.id,"Gene"])),
                 fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
#go
hyper.go.bp=enrichGO(hyper.eg$ENTREZID, 
                     OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 1, 
                     qvalueCutoff = 1,keytype = 'ENTREZID')

#kegg
hyper.kegg=enrichKEGG(hyper.eg$ENTREZID,organism = "hsa", 
                      pAdjustMethod = 'BH',pvalueCutoff = 1, 
                      qvalueCutoff = 1,keyType = 'kegg')

hypo.eg <- bitr(unique(na.omit(anno[hypo_peaks$Peak.id,"Gene"])),
                fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
#go
hypo.go.bp=enrichGO(hypo.eg$ENTREZID, 
                    OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 1, 
                    qvalueCutoff = 1,keytype = 'ENTREZID')

#kegg
hypo.kegg=enrichKEGG(hypo.eg$ENTREZID,organism = "hsa", 
                     pAdjustMethod = 'BH',pvalueCutoff = 1, 
                     qvalueCutoff = 1,keyType = 'kegg')

#load gene expression matrix
id2name=read.table("ensemble_id2gene_name.txt")
rownames(id2name)=id2name$V1
gene_count_data=read.table("analysis/gencode.rsem.count.txt",header=T,row.names = 1)
gene_count_data=cbind(gene=id2name[rownames(gene_count_data),2],gene_count_data)
rownames(gene_count_data)=gene_count_data$gene
gene_count_data$gene=NULL
gene_count_data_n=gene_count_data[,grep("N",colnames(gene_count_data))]
gene_count_data_t=gene_count_data[,grep("T",colnames(gene_count_data))]
gene_count_data_m=cbind(gene_count_data_n,gene_count_data_t)

setwd("/data/xingyang/m6A_zhengjian/")
gene_exp_data=read.table("analysis/gencode.rsem.fpkm.txt",header=T,row.names = 1)
gene_exp_data=cbind(gene=id2name[rownames(gene_exp_data),2],gene_exp_data)
rownames(gene_exp_data)=gene_exp_data$gene
gene_exp_data$gene=NULL
gene_exp_data_n=gene_exp_data[,grep("N",colnames(gene_exp_data))]
gene_exp_data_t=gene_exp_data[,grep("T",colnames(gene_exp_data))]
gene_exp_data_m=cbind(gene_exp_data_n,gene_exp_data_t)

#overlap with diffexpression
diff.gene=diff.exp.gene(gene_count_data_m[,c(paired_sample_normal,paired_sample_tumor)],
                        group_1 = paired_sample_normal,
                        group_2 = paired_sample_tumor)
rownames(diff.gene)=diff.gene$Gene_ID
colnames(diff.gene)[1]="gene.x"

write.table(diff.gene,"analysis/new_figure/diff.gene.txt",sep="\t",quote = F,row.names = F)

#subtype
TCGA_subtype=read.table("TCGA_subtype",sep="\t")
clinic_data=read.table("clinic_data.new3.txt",sep="\t",header = T)
rownames(clinic_data)=paste("X",clinic_data$t_id,sep="")
clinic_data=clinic_data[tumor_sample_in_peak,]
clinic_data[which(clinic_data$PFS_Status=="Yes"),"PFS_Status"]=1
clinic_data[which(clinic_data$PFS_Status=="No"),"PFS_Status"]=0 
Clinic.factor.temp=c("gender","age","smoke","drink","fenhua_class","nerve","xueguan","node","TNM_3")
clinic_data_for_analysis=clinic_data
peak_for_analysis=anno[diff_peak_annotation$Peak.id,]
peak_for_analysis=peak_for_analysis[which(peak_for_analysis$Level.2.gene.type=="mRNA"),]
peak_for_analysis=peak_for_analysis$Peak.id
dc=as.matrix(peak_matrix[peak_for_analysis,])
dc=t(scale(t(dc)))
rownames(clinic_data_for_analysis)=paste(sep="","X",clinic_data_for_analysis$t_id)
clinic_data_for_analysis=clinic_data_for_analysis[tumor_sample_in_peak,]

ConsensusClusterPlus(dc,maxK=3,reps=1000,pItem=0.8,pFeature=1,
                                   clusterAlg="pam",distance="spearman",plot='pdf',
                                   title="pam_spearman/")
    
subtype_sample_id=results[[2]][['consensusClass']]
subtype_sample_id=data.frame(subtype=subtype_sample_id)
subtype_sample_id=cbind(item=rownames(subtype_sample_id),subtype_sample_id)
subtype_sample_id$subtype=paste("s",subtype_sample_id$subtype,sep="")

subtype_annotation=subtype_sample_id
colnames(subtype_annotation)[1]="Sample"
subtype_annotation=merge(subtype_annotation,TCGA_subtype,by="Sample",all.x=T)
clinic_data_for_subtype_annotation=clinic_data_for_analysis[,c("gender","age","fenhua_class",
                                                               "smoke","drink","nerve","xueguan",
                                                               "node","TNM_3")]
clinic_data_for_subtype_annotation=data.frame(cbind(Sample=rownames(clinic_data_for_subtype_annotation),
                                                    clinic_data_for_subtype_annotation))
subtype_annotation=merge(subtype_annotation,clinic_data_for_subtype_annotation,by="Sample",all.x=T)
rownames(subtype_annotation)=subtype_annotation$Sample
subtype_annotation$Sample=NULL
surv_data=data.frame(cbind(month=clinic_data_for_analysis[subtype_sample_id$item,"month"],
                           status=clinic_data_for_analysis[subtype_sample_id$item,"status"],
                           subtype=subtype_sample_id$subtype))
surv_data$month=as.numeric(surv_data$month)
surv=Surv(surv_data$month,surv_data$status==1)
pvalue.survival=coxph(surv~surv_data$subtype,data=surv_data)
pvalue.survival=summary(pvalue.survival)
pvalue.survival=pvalue.survival$sctest[3]

ggsurvplot(survfit(surv~subtype,data=surv_data),data=surv_data,pval = T,
           palette =  c("#ce090d", "#142585","#188514","#d2b346","#d246ce"),
           title=paste("coxph p value:",pvalue.survival,sep=""))
write.table(subtype_annotation,"pam_spearman/subtype_sample_id.txt",
            quote=F,row.name=T,col.names = T,sep='\t')

i=1
n=nrow(peak_matrix)
pvalue.subtype.tumor=c()
while(i<=n){
  s1=as.numeric(peak_matrix[i,subtype_sample_id[which(subtype_sample_id$subtype=="s1"),"item"]])
  s2=as.numeric(peak_matrix[i,subtype_sample_id[which(subtype_sample_id$subtype=="s2"),"item"]])
  pvalue.temp=wilcox.test(s1,s2)
  pvalue.temp=pvalue.temp$p.value
  log2FC.temp=log2(mean(s2)/mean(s1))
  pvalue.subtype.tumor=rbind(pvalue.subtype.tumor,c(rownames(peak_matrix)[i],
                                                    log2FC.temp,pvalue.temp))
  i=i+1
}
pvalue.subtype.tumor=data.frame(pvalue.subtype.tumor)
colnames(pvalue.subtype.tumor)=c("Peak.id","log2FC","P.value")
pvalue.subtype.tumor$log2FC=as.numeric(pvalue.subtype.tumor$log2FC)
pvalue.subtype.tumor$P.value=as.numeric(pvalue.subtype.tumor$P.value)
pvalue.subtype.tumor$FDR=p.adjust(pvalue.subtype.tumor$P.value,method = "fdr")
pvalue.subtype.tumor$subtype="s1"
pvalue.subtype.tumor[which(pvalue.subtype.tumor$log2FC>0),"subtype"]="s2"
pvalue.subtype.tumor=pvalue.subtype.tumor[order(pvalue.subtype.tumor$log2FC),]

ann_colors = list(
  TNM_3 = c(`0`='#66C2A5', `1`='#BD0026'),
  node =  c(`0`='#66C2A5', `1`='#BD0026'),
  xueguan = c(`0`='#66C2A5', `1`='#BD0026'),
  nerve = c(`0`='#66C2A5', `1`='#BD0026'),
  drink = c(`0`='#66C2A5', `1`='#BD0026'),
  smoke = c(`0`='#66C2A5', `1`='#BD0026'),
  fenhua_class = c(`0`='#00A2E8', `1`='#13238B',`3`='#FB8072'),
  age = c(`0`='#66C2A5', `1`='#BD0026'),
  gender = c( `1`='#66C2A5',`2`='#BD0026'),
  moffitt = c(classic='#FB8072', basal_like='#3690C0'),
  collisson = c(Exocrine_like_PDA='#253494', Classical_PDA='#BD0026',QM_PDA='#006837'),
  bailey = c(ADEX='#FF7F27', Immunogenic='#ADD2E5',Pancreatic_Progenitor='#B5E61D',Squamous='#4A9471'),
  subtype = c(s2='#006AB1', s1='#BB3A27',s3='#1c7911'),
  Stage= c(I='#BB3A27', II='#006AB1',III='#E18625',IV="#006837"),
  Signature_in_subtype=c(s2='#BB3A27', s1='#006AB1',s3='#1c7911',s4="#006837",s5="#B5E61D")
)            

subtype_2_plot=subtype_sample_id
subtype_2_plot=subtype_2_plot[order(subtype_2_plot$subtype),]
subtype_annotation.clanc=subtype_2_plot
subtype_annotation.clanc$item=NULL

write.table(cbind(pvalue.subtype.tumor,anno[pvalue.subtype.tumor$Peak.id,]),
            "pam_spearman/subtype.diff.peak.txt",
            quote = F,row.names = F,sep="\t")

pvalue.subtype.tumor.plot=pvalue.subtype.tumor[which(pvalue.subtype.tumor$FDR<0.1),]

signature.type_plot.annotation=pvalue.subtype.tumor$subtype
signature.type_plot.annotation=data.frame(signature.type_plot.annotation)
rownames(signature.type_plot.annotation)=pvalue.subtype.tumor$Peak.id
colnames(signature.type_plot.annotation)="subtype"

pheatmap(peak_matrix[pvalue.subtype.tumor.plot$Peak.id,subtype_2_plot$item],scale='row',border_color = FALSE,
         color = colorRampPalette(c(rep('#0a268a',9),
                                    'white',
                                    rep('#830502',9)))(50),
         annotation_col = subtype_annotation.clanc,annotation_row=signature.type_plot.annotation,
         annotation_colors = ann_colors,
         show_rownames = F,cluster_cols = F,cluster_rows = F,fontsize = 9,
         filename =  "pam_spearman/subtype.signature.heatmapt.png",
         width = 10,height = 10)

#subtype 1
S1.eg <- bitr(unique(na.omit(anno[pvalue.subtype.tumor.plot[which(pvalue.subtype.tumor.plot$log2FC<0),"Peak.id"],
                                  "Gene"])),
              fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
#go
S1.go.bp=enrichGO(S1.eg$ENTREZID, 
                  OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 1, 
                  qvalueCutoff = 1,keytype = 'ENTREZID')

#kegg
S1.kegg=enrichKEGG(S1.eg$ENTREZID,organism = "hsa", 
                   pAdjustMethod = 'BH',pvalueCutoff = 0.17, 
                   qvalueCutoff = 1,keyType = 'kegg')

#subtype 2
S2.eg <- bitr(unique(na.omit(anno[pvalue.subtype.tumor.plot[which(pvalue.subtype.tumor.plot$log2FC>0),"Peak.id"],
                                  "Gene"])),
              fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
#go
S2.go.bp=enrichGO(S2.eg$ENTREZID, 
                  OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 1, 
                  qvalueCutoff = 1,keytype = 'ENTREZID')

#kegg
S2.kegg=enrichKEGG(S2.eg$ENTREZID,organism = "hsa", 
                   pAdjustMethod = 'BH',pvalueCutoff = 1, 
                   qvalueCutoff = 1,keyType = 'kegg')


#barplot
ClusterProfile_pathway=rbind(cbind(as.data.frame(S1.kegg)[1:15,],type="s1"),
                             cbind(as.data.frame(S2.kegg)[1:15,],type="s2"))


ClusterProfile_pathway=ClusterProfile_pathway[which(ClusterProfile_pathway$pvalue<0.05 & ClusterProfile_pathway$Count>1),]
ClusterProfile_pathway=ClusterProfile_pathway[order(ClusterProfile_pathway$type,ClusterProfile_pathway$pvalue,decreasing = T),]
ClusterProfile_pathway$Description=paste(ClusterProfile_pathway$Description,": ",ClusterProfile_pathway$type,sep="")
ClusterProfile_pathway$Description=ordered(ClusterProfile_pathway$Description,level=ClusterProfile_pathway$Description)

ggplot(ClusterProfile_pathway)+geom_bar(aes(x=ClusterProfile_pathway$Description,y=-log10(ClusterProfile_pathway$pvalue),
                                            fill=ClusterProfile_pathway$type),
                                        stat = "identity",position = "stack")+
  coord_flip()+theme_classic()+xlab("")+ylab("-log10(P value)")+scale_fill_manual(values=c("#CE0A0E","#142585"),
                                                                                  name="")+
  theme(text = element_text(size=20))
ggsave("pam_spearman/subtype_pathway_kegg.png",width = 10,height = 7)

ClusterProfile_pathway=rbind(cbind(as.data.frame(S1.go.bp)[1:15,],type="s1"),
                             cbind(as.data.frame(S2.go.bp)[1:15,],type="s2"))


ClusterProfile_pathway=ClusterProfile_pathway[which(ClusterProfile_pathway$pvalue<0.05 & ClusterProfile_pathway$Count>1),]
ClusterProfile_pathway=ClusterProfile_pathway[order(ClusterProfile_pathway$type,ClusterProfile_pathway$pvalue,decreasing = T),]
ClusterProfile_pathway$Description=paste(ClusterProfile_pathway$Description,": ",ClusterProfile_pathway$type,sep="")
ClusterProfile_pathway$Description=ordered(ClusterProfile_pathway$Description,level=ClusterProfile_pathway$Description)

ggplot(ClusterProfile_pathway)+geom_bar(aes(x=ClusterProfile_pathway$Description,y=-log10(ClusterProfile_pathway$pvalue),
                                            fill=ClusterProfile_pathway$type),
                                        stat = "identity",position = "stack")+
  coord_flip()+theme_classic()+xlab("")+ylab("-log10(P value)")+scale_fill_manual(values=c("#CE0A0E","#142585"),
                                                                                  name="")+
  theme(text = element_text(size=20))
ggsave("pam_spearman/subtype_pathway_go.png",width = 10,height = 7)

j=1
clinic_factor=colnames(subtype_annotation)
subtype_clinic=c()
while(j<=length(clinic_factor)){
  if(clinic_factor[j]=="subtype"){
    j=j+1
    next
  }
  plot_data=dcast(data.frame(cbind(factor=subtype_annotation[,clinic_factor[j]],
                                   subtype=subtype_annotation[,"subtype"])),factor~subtype)
  rownames(plot_data)=plot_data$factor
  plot_data$factor=NULL
  pvalue=fisher.test(plot_data)
  pvalue.title=pvalue$p.value
  subtype_clinic=rbind(subtype_clinic,c(clinic_factor[j],pvalue.title))
  if(pvalue.title<0.6){
    jj=1
    kk=dim(plot_data)[1]
    pvalue.matrix=matrix(rep(1,dim(plot_data)[2]*dim(plot_data)[1]),dim(plot_data)[1],dim(plot_data)[2])
    while(jj<=kk){
      jjj=1
      kkk=dim(plot_data)[2]
      while(jjj<=kkk){
        test.matrix=matrix(c(plot_data[jj,jjj],sum(plot_data[-jj,jjj]),sum(plot_data[jj,-jjj]),sum(plot_data[-jj,-jjj])),2,2)
        pvalue=fisher.test(test.matrix,alternative = "greater")
        pvalue=pvalue$p.value
        pvalue.matrix[jj,jjj]=pvalue
        jjj=jjj+1
      }
      jj=jj+1
    }
    rownames(pvalue.matrix)=rownames(plot_data)
    colnames(pvalue.matrix)=colnames(plot_data)
    pheatmap(pvalue.matrix,scale='none',border_color = FALSE,
             display_numbers = matrix(ifelse(pvalue.matrix < 0.05, '*', ''),nrow(pvalue.matrix)),
             color = colorRampPalette(c('#cb3120','white','white','white'))(50),
             show_rownames = T,cluster_cols = F,fontsize_number=50,main=paste(sep="","P value: ",round(pvalue.title,4)),
             cluster_rows = F,fontsize = 9,
             filename = paste("pam_spearman/",clinic_factor[j],"_overlap_heatmap.png",sep=""))
  }
  
  j=j+1
}
subtype_clinic=data.frame(subtype_clinic)
colnames(subtype_clinic)=c("Clinic.factor","P.value")
write.table(subtype_clinic,"pam_spearman/subtype_clinic.fisher.txt",row.name=F,quote=F,sep='\t')

subtype=read.table("pam_spearman/subtype_sample_id.txt",header = T)

#Subtype diff peak
diff.subtype.peak=row.wilcox.test(rpkm_ip_norm_t[c(hyper_peaks$Peak.id,hypo_peaks$Peak.id),],
                                  rownames(subtype[which(subtype$subtype=="s1"),]),
                                  rownames(subtype[which(subtype$subtype=="s2"),]),paired.set = F)


rownames(diff.subtype.peak)=diff.subtype.peak$Peak.id
diff.subtype.peak=diff.subtype.peak[anno[which(!is.na(anno$Gene)),"Peak.id"],]

significan.diff.subtype.peak=diff.subtype.peak[which(diff.subtype.peak$FDR.paired<0.1),]

significan.diff.subtype.peak=significan.diff.subtype.peak[order(significan.diff.subtype.peak$log2FC.paired,decreasing = T),]

diff_peak_annotation=significan.diff.subtype.peak
diff_peak_annotation[intersect(diff_peak_annotation$Peak.id,hyper_peaks$Peak.id),"Peak_type"]="Hyper"
diff_peak_annotation[intersect(diff_peak_annotation$Peak.id,hypo_peaks$Peak.id),"Peak_type"]="Hypo"
diff_peak_annotation=diff_peak_annotation[order(diff_peak_annotation$Peak_type,diff_peak_annotation$log2FC.paired,decreasing = T),]
rownames(diff_peak_annotation)=diff_peak_annotation$Peak.id
diff_peak_annotation$Peak.id=NULL
diff_peak_annotation$log2FC.paired=NULL
diff_peak_annotation$P.value.paired=NULL
diff_peak_annotation$FDR.paired=NULL

subtype_annotation=data.frame(subtype$subtype)
rownames(subtype_annotation)=rownames(subtype)

pheatmap(rpkm_ip_norm_t[rownames(diff_peak_annotation),
                        c(rownames(subtype[order(subtype$subtype),]))],scale = "row",
         color=colorRampPalette(c(rep("navy",10),"white",rep("firebrick3",10)))(50),
         annotation_col = subtype_annotation,cluster_rows = F,cluster_cols = F,
         annotation_row = diff_peak_annotation,
         annotation_colors = list(subtype.subtype=c(s2=tumor.color2,s1=normal.color2),
                                  Peak_type=c(Hyper=tumor.color2,Hypo=normal.color2)),
         show_rownames = F,show_colnames = F,border_color = F,filename = "analysis/new_figure/subtype_heatmap.pdf",
         width = 10,height = 8)
dev.off()

diff_peak_annotation_n=subtype_annotation
diff_peak_annotation_n=rbind(diff_peak_annotation_n,c("s2"))
diff_peak_annotation_n$subtype.subtype=paste0(diff_peak_annotation_n$subtype.subtype,"_normal")
rownames(diff_peak_annotation_n)=gsub("T","N",rownames(diff_peak_annotation_n))
diff_peak_annotation_n=rbind(diff_peak_annotation_n,subtype_annotation)
diff_peak_annotation_n=data.frame(diff_peak_annotation_n)

plot_tmp=scale(rpkm_ip_norm_m[rownames(diff_peak_annotation),
                              c(intersect(c(gsub("T","N",rownames(subtype[order(subtype$subtype),]))),
                                          colnames(rpkm_ip_norm_n)),
                                rownames(subtype[order(subtype$subtype),]))][,-1])
plot_tmp=as.matrix(plot_tmp)
plot_tmp[which(plot_tmp>2)]=2
plot_tmp[which(plot_tmp< -2)]= -2
pheatmap(plot_tmp,scale = "row",
         color=colorRampPalette(c(rep("navy",3),"white",rep("firebrick3",3)))(50),
         annotation_col = diff_peak_annotation_n,cluster_rows = F,cluster_cols = F,
         annotation_row = diff_peak_annotation,
         annotation_colors = list(subtype.subtype=c(s2=tumor.color2,s1=normal.color2,
                                                    s1_normal="skyblue",s2_normal="yellowgreen"),
                                  Peak_type=c(Hyper=tumor.color2,Hypo=normal.color2)),
         show_rownames = F,show_colnames = F,border_color = F,filename = "analysis/new_figure/subtype_heatmap_m.pdf",
         width = 10,height = 8)

subtyp_normal_test=row.wilcox.test(rpkm_ip_norm_n,
                tumor_sample_in_peak = intersect(rownames(diff_peak_annotation_n)[which(diff_peak_annotation_n$subtype.subtype=="s2")],colnames(rpkm_ip_norm_n)),
                normal_sample_in_peak =intersect(rownames(diff_peak_annotation_n)[which(diff_peak_annotation_n$subtype.subtype=="s1")],colnames(rpkm_ip_norm_n)))
write.table(subtyp_normal_test,"analysis/new_figure/subtyp_normal_test.txt",quote = F,row.names = F,sep="\t")
  
#subtype pathway
#subtype 1
S1.eg <- bitr(unique(anno[significan.diff.subtype.peak[which(significan.diff.subtype.peak$log2FC.paired>0),"Peak.id"],
                          "Gene"]),
              fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")

#kegg
S1.kegg=enrichKEGG(S1.eg$ENTREZID,organism = "hsa", 
                   pAdjustMethod = 'BH',pvalueCutoff = 1, 
                   qvalueCutoff = 1,keyType = 'kegg')

#subtype 2
S2.eg <- bitr(unique(anno[significan.diff.subtype.peak[which(significan.diff.subtype.peak$log2FC.paired<0),"Peak.id"],
                          "Gene"]),
              fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")

#kegg
S2.kegg=enrichKEGG(S2.eg$ENTREZID,organism = "hsa", 
                   pAdjustMethod = 'BH',pvalueCutoff = 1, 
                   qvalueCutoff = 1,keyType = 'kegg')
ClusterProfile_pathway=rbind(cbind(as.data.frame(S1.kegg),type="s1"),
                             cbind(as.data.frame(S2.kegg),type="s2"))


ClusterProfile_pathway=ClusterProfile_pathway[which(ClusterProfile_pathway$pvalue<0.05 & ClusterProfile_pathway$Count>1),]
ClusterProfile_pathway=ClusterProfile_pathway[order(ClusterProfile_pathway$type,ClusterProfile_pathway$pvalue,decreasing = T),]
ClusterProfile_pathway$Description=paste(ClusterProfile_pathway$Description,": ",ClusterProfile_pathway$type,sep="")
ClusterProfile_pathway$Description=ordered(ClusterProfile_pathway$Description,level=ClusterProfile_pathway$Description)

ggplot(ClusterProfile_pathway)+geom_bar(aes(x=ClusterProfile_pathway$Description,y=-log10(ClusterProfile_pathway$pvalue),
                                            fill=ClusterProfile_pathway$type),
                                        stat = "identity",position = "stack")+
  coord_flip()+theme_classic()+xlab("")+ylab("-log10(P value)")+scale_fill_manual(values=c(s2=tumor.color2,s1=normal.color2),
                                                                                  name="")+
  scale_y_continuous(expand = c(0,0))+
  theme(text = element_text(size=20))
ggsave("analysis/new_figure/subtype_pathway.pdf",width = 10.2,height = 8)

#m6A seq mutation
maf.file=list.files("annova/result/")
maf.file=maf.file[grep(".hg38_multianno.txt",maf.file)]
i=1
n=length(maf.file)
maf.res=c()
while(i<=n){
  maf.temp=fread(paste("annova/result/",maf.file[i],sep=""),sep="\t",fill=T)
  maf.temp=data.frame(maf.temp)
  maf.temp$Sample=maf.file[i]
  maf.res=rbind(maf.res,maf.temp)
  i=i+1
}
rm(maf.temp)

maf.res=maf.res[which(maf.res$X1000g2015aug_all<0.001|maf.res$X1000g2015aug_all=="."),]
maf.res=maf.res[which(maf.res$V89=="PASS"),]
maf.res=maf.res[which(maf.res$Func.refGene=="exonic"|maf.res$Func.refGene=="splicing"|maf.res$Func.refGene=="exonic;splicing"),]
AD.temp=strsplit2(strsplit2(maf.res$V92,split = "[:]")[,2],split=",")
AD.temp[which(AD.temp=="")]=0
AD.temp=data.frame(AD.temp)
AD.temp=apply(AD.temp,2,as.numeric)
AD.temp=rowSums(AD.temp)
maf.res$AD=AD.temp
maf.res=maf.res[which(maf.res$AD>5),]
maf.res$Sample=strsplit2(maf.res$Sample,split="[.]")[,1]
maf.res$Sample=paste("X",maf.res$Sample,sep="")
maf.res=maf.res[which(maf.res$ExonicFunc.refGene!="synonymous SNV"&
                        maf.res$ExonicFunc.refGene!="unknown"&
                        maf.res$ExonicFunc.refGene!="."),]

unpaired.maf.file=list.files("annova/result_unpaired//")
unpaired.maf.file=unpaired.maf.file[grep(".hg38_multianno.txt",unpaired.maf.file)]
i=1
n=length(unpaired.maf.file)
unpaired.maf.res=c()
while(i<=n){
  unpaired.maf.temp=fread(paste("annova/result_unpaired/",unpaired.maf.file[i],sep=""),sep="\t",fill=T)
  unpaired.maf.temp=data.frame(unpaired.maf.temp)
  unpaired.maf.temp$Sample=unpaired.maf.file[i]
  unpaired.maf.res=rbind(unpaired.maf.res,unpaired.maf.temp)
  i=i+1
}
rm(unpaired.maf.temp)

unpaired.maf.res$ID=paste(unpaired.maf.res$Chr,unpaired.maf.res$Start,sep=":")

normal.panel.site=read.table("annova/result_unpaired/normal.unique.output.site")
normal.panel.site$ID=paste(normal.panel.site[,1],normal.panel.site[,2],sep=":")

unpaired.maf.res=unpaired.maf.res[!(unpaired.maf.res$ID %in% normal.panel.site$ID),]

unpaired.maf.res=unpaired.maf.res[which(unpaired.maf.res$X1000g2015aug_all<0.001|unpaired.maf.res$X1000g2015aug_all=="."),]
unpaired.maf.res=unpaired.maf.res[which(unpaired.maf.res$V89=="PASS"),]
unpaired.maf.res=unpaired.maf.res[which(unpaired.maf.res$Func.refGene=="exonic"|unpaired.maf.res$Func.refGene=="splicing"|unpaired.maf.res$Func.refGene=="exonic;splicing"),]
AD.temp=strsplit2(strsplit2(unpaired.maf.res$V92,split = "[:]")[,2],split=",")
AD.temp[which(AD.temp=="")]=0
AD.temp=data.frame(AD.temp)
AD.temp=apply(AD.temp,2,as.numeric)
AD.temp=rowSums(AD.temp)
unpaired.maf.res$AD=AD.temp
unpaired.maf.res=unpaired.maf.res[which(unpaired.maf.res$AD>5),]
unpaired.maf.res$Sample=strsplit2(unpaired.maf.res$Sample,split="[.]")[,1]
unpaired.maf.res$Sample=paste("X",unpaired.maf.res$Sample,sep="")
unpaired.maf.res=unpaired.maf.res[which(unpaired.maf.res$ExonicFunc.refGene!="synonymous SNV"&
                        unpaired.maf.res$ExonicFunc.refGene!="unknown"&
                        unpaired.maf.res$ExonicFunc.refGene!="."),]
unpaired.maf.res$ID=NULL

maf.res=rbind(unpaired.maf.res,maf.res)

write.table(maf.res,"maf.res.tsv",row.names = F,col.names = T,sep="\t")

#Mut state
KRAS.mut.stat=unique(maf.res[which(maf.res$Gene.refGene=="KRAS"),"Sample"])
CDKN2A.mut.stat=unique(maf.res[which(maf.res$Gene.refGene=="CDKN2A"),"Sample"])
TP53.mut.stat=unique(maf.res[which(maf.res$Gene.refGene=="TP53"),"Sample"])
SMAD4.mut.stat=unique(maf.res[which(maf.res$Gene.refGene=="SMAD4"),"Sample"])

#mutation survival analysis
high_freq_mut_gene=maf.res
high_freq_mut_gene$ID=paste(high_freq_mut_gene$Sample,high_freq_mut_gene$Gene.refGene,sep="_")
high_freq_mut_gene=high_freq_mut_gene[!duplicated(high_freq_mut_gene$ID),]
high_freq_mut_gene_count=table(high_freq_mut_gene$Gene.refGene)
high_freq_mut_gene_count=high_freq_mut_gene_count[which(high_freq_mut_gene_count>6)]

mut_gene_count=high_freq_mut_gene[high_freq_mut_gene$Gene.refGene %in% c("KRAS","TP53","CDKN2A","SMAD4"),]
mut_gene_count_survival=matrix(nrow = length(tumor_sample_in_peak),ncol = 4+1)
mut_gene_count_survival=data.frame(mut_gene_count_survival)
mut_gene_count_survival[,1]=tumor_sample_in_peak
colnames(mut_gene_count_survival)=c("Sample",c("KRAS","TP53","CDKN2A","SMAD4"))
rownames(mut_gene_count_survival)=mut_gene_count_survival$Sample
mut_gene_count_survival$Sample=NULL

i=1
n=4
while(i<=n){
  mut_gene_count_survival[,c("KRAS","TP53","CDKN2A","SMAD4")[i]]=0
  mut_gene_count_survival[mut_gene_count[which(mut_gene_count$Gene.refGene==c("KRAS","TP53","CDKN2A","SMAD4")[i]),"Sample"],
                          c("KRAS","TP53","CDKN2A","SMAD4")[i]]=1
  i=i+1
}
mut_gene_count_survival=t(mut_gene_count_survival)

#Survival table
survival_table=data.frame(Sample=rownames(subtype),Subtype=subtype$subtype)
survival_table$new_ID=clinic_data[survival_table$Sample,"New_ID"]
survival_table$new_PFS=clinic_data[survival_table$Sample,"PFS_month"]
survival_table$new_PFS_Status=clinic_data[survival_table$Sample,"PFS_Status"]
survival_table$new_zhuanyi.fufa=clinic_data[survival_table$Sample,"zhuanyi.fufa"]
survival_table$new_OS=clinic_data[survival_table$Sample,"month"]
survival_table$new_vital.status=clinic_data[survival_table$Sample,"status"]
survival_table$new_vital.status=clinic_data[survival_table$Sample,"status"]
survival_table=merge(survival_table,TCGA_subtype,by="Sample",all.x=T)

mut_gene_count_survival.temp=data.frame(t(mut_gene_count_survival))
mut_gene_count_survival.temp$Sample=rownames(mut_gene_count_survival.temp)

survival_table=merge(survival_table,mut_gene_count_survival.temp[,c("KRAS","TP53","SMAD4","CDKN2A","Sample")],
                     by="Sample",all.x=T)
write.table(survival_table,"analysis/new_figure/survival_table.txt",sep="\t",row.names = F,quote = F)

#subtype survival
ggsurvplot(survfit(Surv(survival_table$new_PFS,survival_table$new_PFS_Status==1)~survival_table$Subtype), 
           data=survival_table,conf.int = F,pval = T,palette = c("#0852aa","#b40816","yellow","green"),fontsize=20)
ggsave("analysis/new_figure/subtype_survival_new_PFS.pdf",width = 5,height = 5)

ggsurvplot(survfit(Surv(survival_table$new_OS,survival_table$new_vital.status==1)~survival_table$Subtype), 
           data=survival_table,conf.int = F,pval = T,palette = c("#0852aa","#b40816","yellow","green"),fontsize=20)
ggsave("analysis/new_figure/subtype_survival_new_OS.pdf",width = 5,height = 5)

pvalue=summary(coxph(Surv(survival_table$new_PFS,survival_table$new_PFS_Status==1)~
                       survival_table$Subtype,
                     data = survival_table))

pvalue.multi=summary(coxph(Surv(survival_table$new_PFS,survival_table$new_PFS_Status==1)~
                             survival_table$Subtype+gender+age+fenhua_class+smoke+drink+nerve+xueguan+node+TNM_3+
                             as.numeric(factor(survival_table$Sample %in% KRAS.mut.stat,levels = c(FALSE,TRUE)))+
                             as.numeric(factor(survival_table$Sample %in% TP53.mut.stat,levels = c(FALSE,TRUE)))+
                             factor(subtype[combine_survival_plot$ID,"bailey"])+
                             factor(subtype[combine_survival_plot$ID,"Collisson"])+
                             factor(subtype[combine_survival_plot$ID,"Moffitt"]),
                           data = clinic_data_for_combind_survival[survival_table$Sample,]))

ggsurvplot(survfit(Surv(survival_table$new_PFS,survival_table$new_PFS_Status==1)~survival_table$Subtype),
           data=survival_table,
           conf.int = F,pval = T,palette = c("#0852aa","#b40816"),fontsize=20,
           title=paste(sep="","\nHR: ",round(pvalue$coefficients[2],5),
                       "\nHR pvalue: ",round(pvalue$coefficients[5],5),
                       "\nHR CI: ",round(pvalue$conf.int[3],5),"-",round(pvalue$conf.int[4],5),
                       "\nClinic scale HR: ",round(pvalue.multi$coefficients[1,2],5),
                       "\nClinic scale HR pvalue: ",round(pvalue.multi$coefficients[1,5],5),
                       "\nClinic scale HR CI: ",round(pvalue.multi$conf.int[1,3],5),"-",round(pvalue.multi$conf.int[1,4],5)))
ggsave("/data/xingyang/m6A_zhengjian/analysis/new_figure/subtype.survival.new.PFS.with.4mut.subtype.pdf",width = 5,height = 6.5)

median(survival_table[which(survival_table$Subtype=="s1"),"new_PFS"])
median(survival_table[which(survival_table$Subtype=="s2"),"new_PFS"])

pvalue=summary(coxph(Surv(survival_table$new_OS,survival_table$new_vital.status==1)~
                       survival_table$Subtype,
                     data = survival_table))

pvalue.multi=summary(coxph(Surv(survival_table$new_OS,survival_table$new_vital.status==1)~
                             survival_table$Subtype+gender+age+fenhua_class+smoke+drink+nerve+xueguan+node+TNM_3+
                             as.numeric(factor(survival_table$Sample %in% KRAS.mut.stat,levels = c(FALSE,TRUE)))+
                             as.numeric(factor(survival_table$Sample %in% TP53.mut.stat,levels = c(FALSE,TRUE)))+
                             factor(subtype[combine_survival_plot$ID,"bailey"])+
                             factor(subtype[combine_survival_plot$ID,"Collisson"])+
                             factor(subtype[combine_survival_plot$ID,"Moffitt"]),
                           data = clinic_data_for_combind_survival[survival_table$Sample,]))

ggsurvplot(survfit(Surv(survival_table$new_OS,survival_table$new_vital.status==1)~survival_table$Subtype),
           data=survival_table,
           conf.int = F,pval = T,palette = c("#0852aa","#b40816"),fontsize=20,
           title=paste(sep="","\nHR: ",round(pvalue$coefficients[2],5),
                       "\nHR pvalue: ",round(pvalue$coefficients[5],5),
                       "\nHR CI: ",round(pvalue$conf.int[3],5),"-",round(pvalue$conf.int[4],5),
                       "\nClinic scale HR: ",round(pvalue.multi$coefficients[1,2],5),
                       "\nClinic scale HR pvalue: ",round(pvalue.multi$coefficients[1,5],5),
                       "\nClinic scale HR CI: ",round(pvalue.multi$conf.int[1,3],5),"-",round(pvalue.multi$conf.int[1,4],5)))
ggsave("/data/xingyang/m6A_zhengjian/analysis/new_figure/subtype.survival.new.OS.with.4mut.subtype.pdf",width = 5,height = 6.5)

median(survival_table[which(survival_table$Subtype=="s1"),"new_OS"])
median(survival_table[which(survival_table$Subtype=="s2"),"new_OS"])

#TCGA subtype survival
combine_survival_plot_bailey=survival_table
combine_survival_plot_bailey$bailey=subtype[combine_survival_plot_bailey$ID,"bailey"]
ggsurvplot(survfit(Surv(combine_survival_plot_bailey$Time,combine_survival_plot_bailey$Status==1)~combine_survival_plot_bailey$bailey), 
           data=combine_survival_plot_bailey,conf.int = F,pval = T,palette = c("#b40816","#0852aa","yellow","green"),fontsize=20)
ggsave("analysis/new_figure/bailey_survival.new.PFS.pdf",width = 5,height = 5)

combine_survival_plot_Collisson=survival_table
combine_survival_plot_Collisson$Collisson=subtype[combine_survival_plot_Collisson$ID,"Collisson"]
ggsurvplot(survfit(Surv(combine_survival_plot_Collisson$Time,combine_survival_plot_Collisson$Status==1)~combine_survival_plot_Collisson$Collisson), 
           data=combine_survival_plot_Collisson,conf.int = F,pval = T,palette = c("#b40816","#0852aa","yellow","green"),fontsize=20)
ggsave("analysis/new_figure/Collisson_survival.new.PFS.pdf",width = 5,height = 5)


combine_survival_plot_Moffitt=survival_table
combine_survival_plot_Moffittdd$Moffitt=subtype[combine_survival_plot_Moffitt$ID,"Moffitt"]
ggsurvplot(survfit(Surv(combine_survival_plot_Moffitt$Time,combine_survival_plot_Moffitt$Status==1)~combine_survival_plot_Moffitt$Moffitt), 
           data=combine_survival_plot_Moffitt,conf.int = F,pval = T,palette = c("#b40816","#0852aa","yellow","green"),fontsize=20)
ggsave("analysis/new_figure/Moffitt_survival.new.PFS.pdf",width = 5,height = 5)

combine_survival_plot_bailey=survival_table
combine_survival_plot_bailey$bailey=subtype[combine_survival_plot_bailey$ID,"bailey"]
ggsurvplot(survfit(Surv(combine_survival_plot_bailey$new_OS,combine_survival_plot_bailey$Status==1)~combine_survival_plot_bailey$bailey), 
           data=combine_survival_plot_bailey,conf.int = F,pval = T,palette = c("#b40816","#0852aa","yellow","green"),fontsize=20)
ggsave("analysis/new_figure/bailey_survival.new.OS.pdf",width = 5,height = 5)

combine_survival_plot_Collisson=survival_table
combine_survival_plot_Collisson$Collisson=subtype[combine_survival_plot_Collisson$ID,"Collisson"]
ggsurvplot(survfit(Surv(combine_survival_plot_Collisson$new_OS,combine_survival_plot_Collisson$Status==1)~combine_survival_plot_Collisson$Collisson), 
           data=combine_survival_plot_Collisson,conf.int = F,pval = T,palette = c("#b40816","#0852aa","yellow","green"),fontsize=20)
ggsave("analysis/new_figure/Collisson_survival.new.OS.pdf",width = 5,height = 5)


combine_survival_plot_Moffitt=survival_table
combine_survival_plot_Moffittdd$Moffitt=subtype[combine_survival_plot_Moffitt$ID,"Moffitt"]
ggsurvplot(survfit(Surv(combine_survival_plot_Moffitt$new_OS,combine_survival_plot_Moffitt$Status==1)~combine_survival_plot_Moffitt$Moffitt), 
           data=combine_survival_plot_Moffitt,conf.int = F,pval = T,palette = c("#b40816","#0852aa","yellow","green"),fontsize=20)
ggsave("analysis/new_figure/Moffitt_survival.new.OS.pdf",width = 5,height = 5)



subtype_bar=rbind(cbind(subtype=subtype$subtype,type="Subtype",value=subtype$subtype),
                  cbind(subtype=subtype$subtype,type="Age",value=subtype$age),
                  cbind(subtype=subtype$subtype,type="Gender",value=subtype$gender),
                  cbind(subtype=subtype$subtype,type="Smoke",value=subtype$smoke),
                  cbind(subtype=subtype$subtype,type="Drink",value=subtype$drink),
                  cbind(subtype=subtype$subtype,type="Neural invasion",value=subtype$nerve),
                  cbind(subtype=subtype$subtype,type="Vascular invasion",value=subtype$xueguan),
                  cbind(subtype=subtype$subtype,type="Lymph node metastasis",value=subtype$node),
                  cbind(subtype=subtype$subtype,type="Stage",value=subtype$TNM_3),
                  cbind(subtype=subtype$subtype,type="Differentiation",value=subtype$fenhua_class),
                  cbind(subtype=subtype$subtype,type="Bailey",value=subtype$bailey),
                  cbind(subtype=subtype$subtype,type="Collisson",value=subtype$Collisson),
                  cbind(subtype=subtype$subtype,type="Moffitt",value=subtype$Moffitt),
                  cbind(subtype=subtype.temp$subtype,type="KRAS.mut",value=subtype.temp$KRAS.mut),
                  cbind(subtype=subtype.temp$subtype,type="TP53.mut",value=subtype.temp$TP53.mut),
                  cbind(subtype=subtype.temp$subtype,type="CDKN2A.mut",value=subtype.temp$CDKN2A.mut),
                  cbind(subtype=subtype.temp$subtype,type="SMAD4.mut",value=subtype.temp$SMAD4.mut),
                  cbind(subtype=stroma_score_subtype$Subtype,type="Stroma contain",value=round(stroma_score_subtype$V1/100,1)))


subtype_bar=data.frame(subtype_bar)
subtype_bar$value=paste(subtype_bar$subtype,subtype_bar$type,subtype_bar$value,sep="_")
subtype_bar=subtype_bar[order(subtype_bar$value,decreasing = T),]
subtype_bar$value=ordered(subtype_bar$value,levels=unique(subtype_bar$value))
subtype_bar$type=ordered(subtype_bar$type,levels=rev(c("Subtype","Age","Gender","Smoke","Drink",
                                                       "Neural invasion","Vascular invasion","Lymph node metastasis",
                                                       "Stage","Differentiation","Bailey","Collisson","Moffitt",
                                                       "KRAS.mut","TP53.mut","CDKN2A.mut","SMAD4.mut","Stroma contain")))

dcast(subtype_bar[which(subtype_bar$type=="KRAS.mut"),],subtype~value)
fisher.test(matrix(c(19,15,5,26),2,2),alternative = "greater")

dcast(subtype_bar[which(subtype_bar$type=="TP53.mut"),],subtype~value)
fisher.test(matrix(c(11,23,1,30),2,2))

dcast(subtype_bar[which(subtype_bar$type=="CDKN2A.mut"),],subtype~value)
fisher.test(matrix(c(7,27,3,28),2,2))

dcast(subtype_bar[which(subtype_bar$type=="SMAD4.mut"),],subtype~value)
fisher.test(matrix(c(3,31,0,31),2,2))

ggplot()+
  geom_histogram(data=subtype_bar,aes(subtype_bar$type,fill=subtype_bar$value),position="fill",stat = "count")+
  #geom_text(data=diff.peak.type.label,aes(x=x,y=y,label=label),color="white",size=5)+
  theme_classic(base_size = 15)+scale_fill_manual(values = c(s1_Subtype_s1=normal.color2,s2_Subtype_s2=tumor.color2,
                                                             s1_Age_1="#B4E9BF",s1_Age_0="#009C54",
                                                             s2_Age_1="#B4E9BF",s2_Age_0="#009C54",
                                                             s1_Gender_1="#A2E1EC",s1_Gender_2="#FFCAD4",
                                                             s2_Gender_1="#A2E1EC",s2_Gender_2="#FFCAD4",
                                                             s1_Smoke_0="#00AED8",s1_Smoke_1="#FF4746",
                                                             s2_Smoke_0="#00AED8",s2_Smoke_1="#FF4746",
                                                             s1_Drink_1="#EED840",s1_Drink_0="#2090FE",
                                                             s2_Drink_1="#EED840",s2_Drink_0="#2090FE",
                                                             `s1_Neural invasion_0`="#0065A0",`s1_Neural invasion_1`="#00BC57",
                                                             `s2_Neural invasion_0`="#0065A0",`s2_Neural invasion_1`="#00BC57",
                                                             `s1_Vascular invasion_0`="#86CDF9",`s1_Vascular invasion_1`="#3D8754",
                                                             `s2_Vascular invasion_0`="#86CDF9",`s2_Vascular invasion_1`="#3D8754",
                                                             `s1_Lymph node metastasis_1`="#DAC388",`s1_Lymph node metastasis_0`="#79BDA6",
                                                             `s2_Lymph node metastasis_1`="#DAC388",`s2_Lymph node metastasis_0`="#79BDA6",
                                                             s1_Stage_1="#36A3AB",s1_Stage_2="#1B2088",
                                                             s2_Stage_1="#36A3AB",s2_Stage_2="#1B2088",
                                                             s1_Differentiation_1="#C93FF8",s1_Differentiation_2="#8699FF",s1_Differentiation_3="#0057A0",
                                                             s2_Differentiation_1="#C93FF8",s2_Differentiation_2="#8699FF",s2_Differentiation_3="#0057A0",
                                                             s1_Moffitt_Classical="#FD35A1",s1_Moffitt_BasalLike="#002EDC",
                                                             s2_Moffitt_Classical="#FD35A1",s2_Moffitt_BasalLike="#002EDC",
                                                             s1_Collisson_Exocrine_like_PDA="#00C500",s1_Collisson_Classical_PDA="#0097DA",s1_Collisson_QM_PDA="#FF0FD7",
                                                             s2_Collisson_Exocrine_like_PDA="#00C500",s2_Collisson_Classical_PDA="#0097DA",s2_Collisson_QM_PDA="#FF0FD7",
                                                             s1_Bailey_ADEX="#DE30FF",s1_Bailey_Immunogenic="#E80000",s1_Bailey_Squamous="#FFC64F",s1_Bailey_Pancreatic_Progenitor="#8699FF",
                                                             s2_Bailey_ADEX="#DE30FF",s2_Bailey_Immunogenic="#E80000",s2_Bailey_Squamous="#FFC64F",s2_Bailey_Pancreatic_Progenitor="#8699FF",
                                                             s1_KRAS.mut_WT="#0097DA",s1_KRAS.mut_Mut="#FF8C00",
                                                             s2_KRAS.mut_WT="#0097DA",s2_KRAS.mut_Mut="#FF8C00",
                                                             s1_TP53.mut_WT="#0097DA",s1_TP53.mut_Mut="#FF8C00",
                                                             s2_TP53.mut_WT="#0097DA",s2_TP53.mut_Mut="#FF8C00",
                                                             s1_CDKN2A.mut_WT="#0097DA",s1_CDKN2A.mut_Mut="#FF8C00",
                                                             s2_CDKN2A.mut_WT="#0097DA",s2_CDKN2A.mut_Mut="#FF8C00",
                                                             s1_SMAD4.mut_WT="#0097DA",s1_SMAD4.mut_Mut="#FF8C00",
                                                             s2_SMAD4.mut_WT="#0097DA",s2_SMAD4.mut_Mut="#FF8C00",
                                                             `S2_Stroma contain_0.4`="#323232",`S2_Stroma contain_0.3`="#5E5E5E",`S2_Stroma contain_0.2`="#858585",`S2_Stroma contain_0.1`="#CACACA",
                                                             `S1_Stroma contain_0.4`="#323232",`S1_Stroma contain_0.3`="#5E5E5E",`S1_Stroma contain_0.2`="#858585",`S1_Stroma contain_0.1`="#CACACA"))+
  scale_y_continuous(expand = c(0,0))+coord_flip()+xlab("")+ylab("Precentage")+theme(legend.title = element_blank())
ggsave("analysis/new_figure/subtype_bar2.pdf",width = 24,height = 7)


#diff peak host gene corr
hyper_peaks_host_gene_corr.tumor=host_gene_corr_peak(
  gene_for_analysis = anno[hyper_peaks$Peak.id,'Gene'],
  peak_for_analysis =hyper_peaks$Peak.id,
  peak_matrix = rpkm_ip_norm_t[,paired_sample_tumor],gene_matrix = gene_exp_data_t[,paired_sample_tumor])

hyper_peaks_host_gene_corr.normal=host_gene_corr_peak(
  gene_for_analysis = anno[hyper_peaks$Peak.id,'Gene'],
  peak_for_analysis =hyper_peaks$Peak.id,
  peak_matrix = rpkm_ip_norm_n[,paired_sample_normal],gene_matrix = gene_exp_data_n[,paired_sample_normal])


hyper_peaks_host_gene_corr=merge(hyper_peaks_host_gene_corr.tumor,hyper_peaks_host_gene_corr.normal,by="Peak.id",all.x=T)

hypo_peaks_host_gene_corr.tumor=host_gene_corr_peak(
  gene_for_analysis = anno[hypo_peaks$Peak.id,'Gene'],
  peak_for_analysis =hypo_peaks$Peak.id,
  peak_matrix = rpkm_ip_norm_t[,paired_sample_tumor],gene_matrix = gene_exp_data_t[,paired_sample_tumor])

hypo_peaks_host_gene_corr.normal=host_gene_corr_peak(
  gene_for_analysis = anno[hypo_peaks$Peak.id,'Gene'],
  peak_for_analysis =hypo_peaks$Peak.id,
  peak_matrix = rpkm_ip_norm_n[,paired_sample_normal],gene_matrix = gene_exp_data_n[,paired_sample_normal])

hypo_peaks_host_gene_corr=merge(hypo_peaks_host_gene_corr.tumor,hypo_peaks_host_gene_corr.normal,by="Peak.id",all.x=T)

point_plot=rbind(cbind(merge(hyper_peaks_host_gene_corr,diff.gene,by="gene.x",all.x=T),type="Non"),
                 cbind(merge(hypo_peaks_host_gene_corr,diff.gene,by="gene.x",all.x=T),type="Non"))
point_plot=merge(point_plot,significan.diff.common.peak,by="Peak.id",all.x=T)

point_plot[which(point_plot$padj<0.1 & point_plot$FDR<0.1 & point_plot$log2FoldChange>0 & 
                   point_plot$log2FC.paired>0&
                   abs(point_plot$Corr.x)>0.25 & point_plot$Pvalue.x<0.05),"type"]="Hyper_up"
point_plot[which(point_plot$padj<0.1 & point_plot$FDR<0.1 & point_plot$log2FoldChange<0 & 
                   point_plot$log2FC.paired<0&
                   abs(point_plot$Corr.x)>0.25& point_plot$Pvalue.x<0.05),"type"]="Hypo_down"
point_plot[which(point_plot$padj<0.1 & point_plot$FDR<0.1 & point_plot$log2FoldChange<0 & 
                   point_plot$log2FC.paired>0&
                   abs(point_plot$Corr.x)>0.25& point_plot$Pvalue.x<0.05),"type"]="Hyper_down"
point_plot[which(point_plot$padj<0.1 & point_plot$FDR<0.1 & point_plot$log2FoldChange>0 & 
                   point_plot$log2FC.paired<0&
                   abs(point_plot$Corr.x)>0.25& point_plot$Pvalue.x<0.05),"type"]="Hypo_up"

point_plot.sig=point_plot[which(point_plot$type!="Non"),]
point_plot.non=point_plot[which(point_plot$type=="Non"),]
point_plot.non$Corr.x=0

point_label=data.frame(x.lo=c(4,4,-4,-4),y.lo=c(4,-4,4,-4),
                       label=c(dim(point_plot.sig[which(point_plot.sig$log2FC.paired>0 &
                                                          point_plot.sig$log2FoldChange>0),])[1],
                               dim(point_plot.sig[which(point_plot.sig$log2FC.paired<0 &
                                                          point_plot.sig$log2FoldChange>0),])[1],
                               dim(point_plot.sig[which(point_plot.sig$logFC>0 &
                                                          point_plot.sig$log2FoldChange<0),])[1],
                               dim(point_plot.sig[which(point_plot.sig$logFC<0 &
                                                          point_plot.sig$log2FoldChange<0),])[1]))

ggplot()+geom_point(aes(y=point_plot.non$log2FC.paired,
                        x=point_plot.non$log2FoldChange,
                        color=point_plot.non$type),size=5,alpha=.9)+
  geom_point(aes(y=point_plot.sig$log2FC.paired,
                 x=point_plot.sig$log2FoldChange,
                 color=point_plot.sig$type),size=5,alpha=.9)+
  
  geom_text(data=point_label,
            aes(x=point_label$x.lo,y=point_label$y.lo,label=point_label$label),size=5)+
  theme_classic()+
  geom_hline(yintercept = 0,size=1.5,linetype="dashed")+
  geom_vline(xintercept = 0,size=1.5,linetype="dashed")+
  xlab(paste("Gene log2Foldchange",sep=""))+
  ylab("m6A log2Foldchange")+theme(text = element_text(size=20))+
  scale_color_manual(values = c(tumor.color,normal.color,annotation_color1,random.color),name="")
# scale_color_gradient2(midpoint=0,mid=random.color,limits=c(-1,1),high = "firebrick",low = "navy",name="Correlation")
ggsave("analysis/new_figure/Diff.point.pdf",width = 6.5,height = 5)


#POSTAR hyper peak analysis

significan.diff.common.peak.for.POSTAR=point_plot[which(point_plot$log2FC.paired>0 ),"Peak.id"]
significan.diff.common.peak.for.POSTAR=anno[significan.diff.common.peak.for.POSTAR,]
significan.diff.common.peak.for.POSTAR.offset.2000=significan.diff.common.peak.for.POSTAR
significan.diff.common.peak.for.POSTAR.offset.2000$Start=significan.diff.common.peak.for.POSTAR.offset.2000$Start-1000
significan.diff.common.peak.for.POSTAR.offset.2000$End=significan.diff.common.peak.for.POSTAR.offset.2000$End+1000
write.table(significan.diff.common.peak.for.POSTAR.offset.2000[,1:4],"temp.bed",quote = F,row.names = F,col.names = F,sep="\t")
system("intersectBed -a temp.bed -b human_RBP_binding_sites.txt -wa -wb > diff_peak_POSTAR_overlap.bed")
diff_peak_POSTAR_overlap=read.table("diff_peak_POSTAR_overlap.bed",sep="\t")
diff_peak_POSTAR_overlap=diff_peak_POSTAR_overlap[!duplicated(paste(diff_peak_POSTAR_overlap$V4,diff_peak_POSTAR_overlap$V11,sep="_")),]
diff_peak_POSTAR_overlap=data.frame(table(diff_peak_POSTAR_overlap$V11))
diff_peak_POSTAR_overlap$Var1=as.character(diff_peak_POSTAR_overlap$Var1)
diff_peak_POSTAR_overlap=diff_peak_POSTAR_overlap[order(diff_peak_POSTAR_overlap$Freq,decreasing = T),]
rownames(diff_peak_POSTAR_overlap)=diff_peak_POSTAR_overlap$Var1
diff_peak_POSTAR_overlap=diff_peak_POSTAR_overlap[setdiff(rownames(diff_peak_POSTAR_overlap),reader),]
diff_peak_POSTAR_overlap_plot=diff_peak_POSTAR_overlap[unique(c(diff_peak_POSTAR_overlap[1:30,"Var1"],WE)),]
diff_peak_POSTAR_overlap_plot$Var1=unique(c(diff_peak_POSTAR_overlap[1:30,"Var1"],WE))
diff_peak_POSTAR_overlap_plot[which(is.na(diff_peak_POSTAR_overlap_plot$Freq)),"Freq"]=0
diff_peak_POSTAR_overlap_plot=diff_peak_POSTAR_overlap_plot[order(diff_peak_POSTAR_overlap_plot$Freq,decreasing = T),]
diff_peak_POSTAR_overlap_plot$Var1=ordered(diff_peak_POSTAR_overlap_plot$Var1,levels=diff_peak_POSTAR_overlap_plot$Var1)
ggplot(diff_peak_POSTAR_overlap_plot)+
  geom_bar(aes(x=diff_peak_POSTAR_overlap_plot$Var1,y=diff_peak_POSTAR_overlap_plot$Freq),
           fill=tumor.color,stat="identity")+
  scale_y_continuous(expand = c(0,0))+theme_classic()+theme(axis.text.x = element_text(angle=90))+
  ylab("POSTAR hyper peak overlap frequancy")+xlab("RBP")
ggsave("analysis/new_figure/POSTAR_clipdb_binding_overlap.pdf",width = 6,height = 5)


#Hyper peaks RBP correlation
RBP_hyper_peak_corr=gene_corr_peak(peak_for_analysis = point_plot[which(point_plot$log2FC.paired>0 ),"Peak.id"],
                                   gene_for_analysis = as.character(diff_peak_POSTAR_overlap_plot$Var1),
                                   peak_matrix = rpkm_ip_norm_t[,paired_sample_tumor],
                                   gene_matrix = gene_exp_data_t[,paired_sample_tumor])
RBP_hyper_peak_corr_plot=cbind(RBP_hyper_peak_corr$gene,RBP_hyper_peak_corr$Corr)

#Non diff peaks RBP correlation
RBP_non_peak_corr=gene_corr_peak(peak_for_analysis = setdiff(rownames(rpkm_ip_norm_t),point_plot$Peak.id),
                                 gene_for_analysis =as.character(diff_peak_POSTAR_overlap_plot$Var1),
                                 peak_matrix = rpkm_ip_norm_t[,paired_sample_tumor],
                                 gene_matrix = gene_exp_data_t[,paired_sample_tumor])
RBP_non_peak_corr_plot=cbind(RBP_non_peak_corr$gene,RBP_non_peak_corr$Corr)

#fisher test
sig_score_hyper_peak_corr=RBP_hyper_peak_corr
sig_score_non_peak_corr=RBP_non_peak_corr
temp=table(RBP_hyper_peak_corr[which(abs(RBP_hyper_peak_corr$Corr)>0.25 & 
                                       RBP_hyper_peak_corr$Pvalue<0.05),"gene"])
temp=temp[order(temp,decreasing = T)]
sig_RBP=names(temp)
i=1
n=length(sig_RBP)
enrich_test_RBP=c()
while(i<=n){
  p.temp=fisher.test(
    matrix(c(nrow(sig_score_hyper_peak_corr[which(sig_score_hyper_peak_corr$gene==sig_RBP[i]&
                                                    abs(sig_score_hyper_peak_corr$Corr)>0.25&
                                                    sig_score_hyper_peak_corr$Pvalue<0.05),])[1],
             nrow(sig_score_non_peak_corr[which(sig_score_non_peak_corr$gene==sig_RBP[i]&
                                                  abs(sig_score_non_peak_corr$Corr)>0.25&
                                                  sig_score_hyper_peak_corr$Pvalue<0.05),])[1],
             nrow(sig_score_hyper_peak_corr[which(sig_score_hyper_peak_corr$gene==sig_RBP[i]&
                                                    (abs(sig_score_hyper_peak_corr$Corr)<0.25|
                                                       sig_score_hyper_peak_corr$Pvalue>0.05)),])[1],
             nrow(sig_score_non_peak_corr[which(sig_score_non_peak_corr$gene==sig_RBP[i]&
                                                  (abs(sig_score_non_peak_corr$Corr)<0.25|
                                                     sig_score_hyper_peak_corr$Pvalue>0.05)),])[1]),2,2),
    alternative = "greater")
  enrich_test_RBP=rbind(enrich_test_RBP,c(P=p.temp$p.value,OR=p.temp$estimate,RBP=sig_RBP[i]))
  i=i+1
}
enrich_test_RBP=data.frame(enrich_test_RBP)
enrich_test_RBP$P=as.numeric(enrich_test_RBP$P)
enrich_test_RBP$OR.odds.ratio=as.numeric(enrich_test_RBP$OR.odds.ratio)

enrich_test_RBP=enrich_test_RBP[order(enrich_test_RBP$P),]
enrich_test_RBP$RBP=ordered(enrich_test_RBP$RBP,levels=enrich_test_RBP$RBP)
enrich_test_RBP$logP=-log10(enrich_test_RBP$P)

enrich_test_RBP_plot_hyper=enrich_test_RBP
enrich_test_RBP_plot=enrich_test_RBP

#Hyper peaks RBP RF
   hyper_peak_RF=run_RF(clinic_data.for.analysis = clinic_data,
                        target.peak.matrix = rpkm_ip_norm_t[hyper_peaks$Peak.id,
                                                            paired_sample_tumor],
                        gene_matrix = gene_exp_data_t[,paired_sample_tumor],
                        target.gene = as.character(diff_peak_POSTAR_overlap_plot$Var1),
                        Clinic.factor.temp = Clinic.factor.temp)
   hyper_peak_RF=merge(hyper_peak_RF,data.frame(X1=as.character(diff_peak_POSTAR_overlap$Var1[1:30])),by="X1",all.y=T)
   contribution.WER_plot=hyper_peak_RF
   write.table(contribution.WER_plot,"contribution.WER_plot_hyper.txt")

   #Contribution summary
   score_summary=ddply(hyper_peak_RF,"X1",summarise,Mean=mean(X2),Median=median(X2))
   score_summary=score_summary[order(score_summary$Median,decreasing = T),]
   score_summary$X1=ordered(score_summary$X1,levels=score_summary$X1)
   score_summary$type=NA
   score_summary[1:5,"type"]="Top5"
   score_summary[-c(1:5),"type"]="Others"
   colnames(score_summary)[1]="RBP"
    
#Plot
combine_plot=merge(enrich_test_RBP_plot_hyper,score_summary,by="RBP")
combine_plot$type="Others"
combine_plot[which(combine_plot$logP>20),"type"]=as.character(combine_plot[which(combine_plot$logP>20),"RBP"])
combine_plot[which(combine_plot$Median>0.03),"type"]=as.character(combine_plot[which(combine_plot$Median>0.03),"RBP"])
combine_plot[which(combine_plot$Median>0.03),"label"]=as.character(combine_plot[which(combine_plot$Median>0.03),"RBP"])
combine_plot$label=NA
combine_plot[which(combine_plot$type!="Others"),"label"]=combine_plot[which(combine_plot$type!="Others"),"type"]
ggplot(combine_plot)+geom_point(aes(x=combine_plot$logP,y=combine_plot$Median,color=combine_plot$type),size=3)+
  geom_text(aes(x=combine_plot$logP,y=combine_plot$Median+0.005,label=combine_plot$label),size=3)+
  xlab("-log(P value)\nRegulation enrichment")+ylab("Contribution score\nRandom forest")+
  guides(color=F)+
  theme_classic()
ggsave(paste("analysis/new_figure/combine_point.png",sep=""),width = 6,height = 5)

#CSTF2 boxplot
xyz <- data.frame(Tumor = as.numeric(gene_exp_data_t["CSTF2",paired_sample_tumor]),
                  Normal = as.numeric(gene_exp_data_n["CSTF2",paired_sample_normal])) 
p.test=wilcox.test(as.numeric(gene_exp_data_t["CSTF2",paired_sample_tumor]),
                   as.numeric(gene_exp_data_n["CSTF2",paired_sample_normal]),paired = T) 
xyz=melt(xyz)
xyz$variable=ordered(xyz$variable,levels=c("Normal","Tumor"))
ggplot(xyz)+geom_boxplot(aes(x=xyz$variable,y=xyz$value,color=xyz$variable))+
  theme_classic()+scale_color_manual(values = c(normal.color,tumor.color),name="")+
  xlab(paste("P value: ",round(p.test$p.value,4),sep=""))+
  ylab("CSTF2 Gene expression (FPKM)")+guides(fill=F)
ggsave("analysis/new_figure/CSTF2.exp.boxplot.pdf",width = 4,height = 3.5,dpi=500)

ggplot(xyz)+geom_boxplot(aes(x=xyz$variable,y=log2(xyz$value+1),color=xyz$variable))+
  theme_classic()+scale_color_manual(values = c(normal.color,tumor.color),name="")+
  xlab(paste("P value: ",round(p.test$p.value,4),sep=""))+
  ylab("CSTF2 Gene expression (logFPKM+1)")+guides(fill=F)
ggsave("analysis/new_figure/CSTF2.exp.boxplot_log.pdf",width = 4,height = 3.5,dpi=500)

#CSTF2 in subtype
CSTF2_exp=data.frame(t(gene_exp_data_t["CSTF2",]),type=NA)
CSTF2_exp[rownames(subtype[which(subtype$subtype=="s1"),]),"type"]="S1"
CSTF2_exp[rownames(subtype[which(subtype$subtype=="s2"),]),"type"]="S2"
p.temp=wilcox.test(CSTF2_exp[which(CSTF2_exp$type=="S1"),"CSTF2"],
                   CSTF2_exp[which(CSTF2_exp$type=="S2"),"CSTF2"])

ggplot(CSTF2_exp)+geom_boxplot(aes(x=CSTF2_exp$type,y=log2(CSTF2_exp$CSTF2+1),color=CSTF2_exp$type))+
  theme_classic()+scale_color_manual(values = c(normal.color2,tumor.color2),name="")+
  xlab(paste("P value: ",round(p.temp$p.value,5)))+
  ylab("CSTF2 gene expression (logFPKM+1)")+guides(fill=F)
ggsave("analysis/new_figure/CSTF2_subtype_log.pdf",width = 4.5,height = 3.5)

#CSTF2 peak annotation
setwd("/data/xingyang/m6A_zhengjian/analysis/annotation/")
system("perl m6A_annotate_forGTF_xingyang_v2.pl /data/database/hg38/GENCODE/gencode.v25.annotation.gtf /data/xingyang/m6A_zhengjian/analysis/final.cellline.bed /data/xingyang/m6A_zhengjian/analysis/annotation/new_annotation/CSTF2_cellline")
setwd("/data/xingyang/m6A_zhengjian/")
CSTF2.peak.anno=read.table("analysis/annotation/new_annotation/CSTF2_cellline.anno.txt",header = F,sep="\t")
colnames(CSTF2.peak.anno)=c("Chr","Start","End","Peak.id","Chr.gene","Start.gene","End.gene","Transcript.id",
                            "Nouse","Strand","Gene","Gene.type","Gene.site","Peak.position","Ensembl.gene.id","Level.1.gene.type",
                            "Level.2.gene.type")
rownames(CSTF2.peak.anno)=CSTF2.peak.anno$Peak.id

CSTF2.peak.unanno=read.table("analysis/annotation/new_annotation/CSTF2_cellline.unanno.txt",header = F,sep="\t")
colnames(CSTF2.peak.unanno)[4]="Peak.id"
CSTF2.peak.unanno=merge(CSTF2.peak.unanno,CSTF2.peak.anno,by="Peak.id",all.x=T)
CSTF2.peak.unanno=CSTF2.peak.unanno[,colnames(CSTF2.peak.anno)]

CSTF2.peak.anno=rbind(CSTF2.peak.anno,CSTF2.peak.unanno)

CSTF2.peak.anno[which(is.na(CSTF2.peak.anno$Level.2.gene.type)),"Level.2.gene.type"]="Unknown"
rownames(CSTF2.peak.anno)=CSTF2.peak.anno$Peak.id

CSTF2.peak.anno[which(CSTF2.peak.anno$Gene.type=="coding" & CSTF2.peak.anno$Level.2.gene.type!="mRNA"),"Level.2.gene.type"]="mRNA"
CSTF2.peak.anno[which(CSTF2.peak.anno$Gene.type=="coding" & CSTF2.peak.anno$Level.2.gene.type!="mRNA"),"Level.1.gene.type"]="mRNA"

rm(CSTF2.peak.unanno)

#remove m6Am 5'UTR peaks in CSTF2 si
CSTF2_anno_5UTR=CSTF2.peak.anno[which(CSTF2.peak.anno$Gene.site=="5UTR"),]

write.table(strsplit2(CSTF2_anno_5UTR$Peak.id,split=":|-"),"UTR5.peak.bed",row.names = F,col.names = F,quote = F,sep="\t")
system("fastaFromBed -fi /data/database/hg38/genome.fa -bed UTR5.peak.bed -fo UTR5.peak.fa")
system("/data/software/homer/bin/homer2 find -i UTR5.peak.fa -m /data/xingyang/m6A_zhengjian/BCA.motif -p 50 > /data/xingyang/m6A_zhengjian/analysis/BCA_peak_offset.txt")
BCA_in_5UTR_offset=read.table("analysis/BCA_peak_offset.txt",header=F)

CSTF2_anno_5UTR=CSTF2.peak.anno[unique(BCA_in_5UTR_offset$V1),]
CSTF2_anno_5UTR=merge(gtf.temp,CSTF2_anno_5UTR,by="Transcript.id",all.y=T)
CSTF2_anno_5UTR=CSTF2_anno_5UTR[which(CSTF2_anno_5UTR$V3=="UTR"),]
CSTF2_anno_5UTR$temp.start=CSTF2_anno_5UTR$V4-CSTF2_anno_5UTR$Start
CSTF2_anno_5UTR$temp.end=CSTF2_anno_5UTR$V5-CSTF2_anno_5UTR$Start
CSTF2_anno_5UTR[which(CSTF2_anno_5UTR$temp.start>0),"temp.start"]=1
CSTF2_anno_5UTR[which(CSTF2_anno_5UTR$temp.start<0),"temp.start"]=(-1)
CSTF2_anno_5UTR[which(CSTF2_anno_5UTR$temp.end>0),"temp.end"]=1
CSTF2_anno_5UTR[which(CSTF2_anno_5UTR$temp.end<0),"temp.end"]=(-1)
CSTF2_anno_5UTR=CSTF2_anno_5UTR[which((CSTF2_anno_5UTR$temp.start*CSTF2_anno_5UTR$temp.end)<=0),]
CSTF2_anno_5UTR.bed=cbind(CSTF2_anno_5UTR$V1,CSTF2_anno_5UTR$V4,CSTF2_anno_5UTR$V5,CSTF2_anno_5UTR$Peak.id,".",CSTF2_anno_5UTR$V7)
write.table(CSTF2_anno_5UTR.bed,"UTR5.peak.bed",row.names = F,col.names = F,quote = F,sep="\t")
system("fastaFromBed -fi /data/database/hg38/genome.fa -bed UTR5.peak.bed -s -name -fo UTR5.peak.fa")

temp.utr5=read.table("UTR5.peak.fa",sep="\n")
temp.utr5=cbind(temp.utr5,substr(temp.utr5[,1],1,1))

i=2
n=nrow(temp.utr5)
temp.utr5=cbind(temp.utr5,type=NA)
while(i<=n){
  if(temp.utr5[i,2]=="A"){
    temp.utr5[c(i-1,i),"type"]="m6Am"
  }
  i=i+2
}
m6Am_CSTF2=na.omit(temp.utr5)
m6Am_CSTF2=m6Am_CSTF2[grep(">",m6Am_CSTF2$V1),]
m6Am_CSTF2=gsub(">","",m6Am_CSTF2$V1)
m6Am_CSTF2=strsplit2(m6Am_CSTF2,split="[(]")[,1]

CSTF2.peak.anno=CSTF2.peak.anno[setdiff(rownames(CSTF2.peak.anno),m6Am_CSTF2),]

CSTF2_m6A_matrix=read.table("CSTF2_m6A.txt",header = T,row.names = 1)
CSTF2_m6A_matrix=CSTF2_m6A_matrix[,grep("ip",colnames(CSTF2_m6A_matrix))]/
  (CSTF2_m6A_matrix[,grep("input",colnames(CSTF2_m6A_matrix))]+1)
CSTF2_m6A_matrix=CSTF2_m6A_matrix[CSTF2.peak.anno$Peak.id,]


#si-CSTF2 m6A PANC
si_CSTF2_FC_matrix=CSTF2_m6A_matrix[,grep("siCS_",colnames(CSTF2_m6A_matrix))]/
  CSTF2_m6A_matrix[,grep("NC_",colnames(CSTF2_m6A_matrix))]

si_CSTF2_FC_matrix$siCS_1.ip=log2(si_CSTF2_FC_matrix$siCS_1.ip)
si_CSTF2_FC_matrix$siCS_2.ip=log2(si_CSTF2_FC_matrix$siCS_2.ip)

summary(si_CSTF2_FC_matrix)

si_CSTF2_FC_matrix$Color=NA

si_CSTF2_FC_matrix[which(si_CSTF2_FC_matrix$siCS_1.ip<(-0.26)&si_CSTF2_FC_matrix$siCS_2.ip<(-0.26)),"Color"]="Down"
si_CSTF2_FC_matrix[which(si_CSTF2_FC_matrix$siCS_1.ip>(0.26)&si_CSTF2_FC_matrix$siCS_2.ip>(0.26)),"Color"]="Up"
si_CSTF2_FC_matrix[which(is.na(si_CSTF2_FC_matrix$Color)),"Color"]="Non"
si_CSTF2_FC_matrix[which(si_CSTF2_FC_matrix$siCS_2.ip==-Inf),"siCS_2.ip"]=
  min(si_CSTF2_FC_matrix[which(si_CSTF2_FC_matrix$siCS_2.ip!=-Inf),"siCS_2.ip"])
si_CSTF2_FC_matrix[which(si_CSTF2_FC_matrix$siCS_1.ip==-Inf),"siCS_1.ip"]=
  min(si_CSTF2_FC_matrix[which(si_CSTF2_FC_matrix$siCS_1.ip!=-Inf),"siCS_1.ip"])

text_plot=data.frame(table(si_CSTF2_FC_matrix$Color))
text_plot=text_plot[which(text_plot$Var1!="Non"),]
text_plot$X.temp=NA
text_plot$Y.temp=NA
text_plot[which(text_plot$Var1=="Up"),c("X.temp","Y.temp")]=c(0,5.5)
text_plot[which(text_plot$Var1=="Down"),c("X.temp","Y.temp")]=c(6,0)

si_CSTF2_FC_matrix$mean_NC=rowMeans(CSTF2_m6A_matrix[,grep("NC_",colnames(CSTF2_m6A_matrix))])
si_CSTF2_FC_matrix$mean_siCS=rowMeans(CSTF2_m6A_matrix[,grep("siCS_",colnames(CSTF2_m6A_matrix))])
si_CSTF2_FC_matrix$Color=ordered(si_CSTF2_FC_matrix$Color,level=c("Down","Up","Non"))
si_CSTF2_FC_matrix=si_CSTF2_FC_matrix[order(si_CSTF2_FC_matrix$Color,decreasing = T),]

ggplot(si_CSTF2_FC_matrix)+
  geom_point(aes(x=log2(si_CSTF2_FC_matrix$mean_NC+1),y=log2(si_CSTF2_FC_matrix$mean_siCS+1),
                 color=si_CSTF2_FC_matrix$Color),size=3)+
  scale_color_manual(values = c(Down="#17489C",Up="#E7423B",Non=random.color),name="")+
  theme_classic()+
  geom_abline(intercept = 0)+
  xlab("Log2(Mean m6A level in NC +1)")+
  ylab("Log2(Mean m6A level si-CSTF2 +1)")+
  geom_text(data = text_plot,aes(x=text_plot$X.temp,y=text_plot$Y.temp,label=text_plot$Freq))
ggsave("analysis/new_figure/PA_si_CSTF2_m6A_diff_point_1.pdf",width = 5.5,height = 4.5)

#SW and siCSTF2 m6A
siCSTF2_matrix2=CSTF2_m6A_matrix

siCSTF2_FC_SW_matrix2=(siCSTF2_matrix2[,c("SWsiCSTF21.ip","SWsiCSTF22.ip")])/
  (siCSTF2_matrix2[,c("SWcon1.ip","SWcon2.ip")])
siCSTF2_FC_SW_matrix2[which(is.na(siCSTF2_FC_SW_matrix2$SWsiCSTF21.ip)),1]=1
siCSTF2_FC_SW_matrix2[which(is.na(siCSTF2_FC_SW_matrix2$SWsiCSTF22.ip)),2]=1

siCSTF2_FC_SW_matrix2[which(siCSTF2_FC_SW_matrix2$SWsiCSTF21.ip==Inf),1]=max(siCSTF2_FC_SW_matrix2[which(siCSTF2_FC_SW_matrix2$SWsiCSTF21.ip!=Inf),1])
siCSTF2_FC_SW_matrix2[which(siCSTF2_FC_SW_matrix2$SWsiCSTF22.ip==Inf),2]=max(siCSTF2_FC_SW_matrix2[which(siCSTF2_FC_SW_matrix2$SWsiCSTF22.ip!=Inf),2])

siCSTF2_FC_SW_matrix2$SWsiCSTF21.ip=log2(siCSTF2_FC_SW_matrix2$SWsiCSTF21.ip)
siCSTF2_FC_SW_matrix2$SWsiCSTF22.ip=log2(siCSTF2_FC_SW_matrix2$SWsiCSTF22.ip)

siCSTF2_FC_SW_matrix2$type_SW="Non"
siCSTF2_FC_SW_matrix2[which(siCSTF2_FC_SW_matrix2$SWsiCSTF21.ip>0.26& siCSTF2_FC_SW_matrix2$SWsiCSTF22.ip >0.26),"type_SW"]="Up"
siCSTF2_FC_SW_matrix2[which(siCSTF2_FC_SW_matrix2$SWsiCSTF21.ip< -0.26 & siCSTF2_FC_SW_matrix2$SWsiCSTF22.ip < -0.26),"type_SW"]="Down"


siCSTF2_FC_SW_matrix2$Mean_SWcon=
  rowMeans(siCSTF2_matrix2[rownames(siCSTF2_FC_SW_matrix2),c("SWcon1.ip","SWcon2.ip")])

siCSTF2_FC_SW_matrix2$Mean_SWsi=
  rowMeans(siCSTF2_matrix2[rownames(siCSTF2_FC_SW_matrix2),c("SWsiCSTF21.ip","SWsiCSTF22.ip")])


text_plot=data.frame(table(siCSTF2_FC_SW_matrix2$type_SW))
text_plot=text_plot[which(text_plot$Var1!="Non"),]
text_plot$X.temp=NA
text_plot$Y.temp=NA
text_plot[which(text_plot$Var1=="Up"),c("X.temp","Y.temp")]=c(0.5,4)
text_plot[which(text_plot$Var1=="Down"),c("X.temp","Y.temp")]=c(4,0.5)

siCSTF2_FC_SW_matrix2$type_SW=ordered(siCSTF2_FC_SW_matrix2$type_SW,levels=c("Non","Down","Up"))
siCSTF2_FC_SW_matrix2=siCSTF2_FC_SW_matrix2[order(siCSTF2_FC_SW_matrix2$type_SW),]

ggplot(siCSTF2_FC_SW_matrix2)+geom_point(aes(x=log2(siCSTF2_FC_SW_matrix2$Mean_SWcon+1),
                                             y=log2(siCSTF2_FC_SW_matrix2$Mean_SWsi+1),
                                             color=siCSTF2_FC_SW_matrix2$type_SW),size=2.2)+
  geom_abline(intercept = 0,linetype="dashed")+theme_classic()+
  scale_color_manual(values=c(Down="#17489C",Non=random.color,Up="#E7423B"),name="type")+
  geom_text(data=text_plot,aes(x=text_plot$X.temp,y=text_plot$Y.temp,label=text_plot$Freq))
ggsave("analysis/new_figure/SW_siCSTF2_diff_m6A_point_1.pdf",width = 5,height = 4)

#PA CSTF2 OE m6A
CSTF2_OE_FC_SW_matrix=CSTF2_m6A_matrix
CSTF2_OE_FC_SW_matrix2=(CSTF2_OE_FC_SW_matrix[,c("SW_CS_OE_1.ip","SW_CS_OE_2.ip","PA_CS_OE_1.ip","PA_CS_OE_2.ip")])/
  (CSTF2_OE_FC_SW_matrix[,c("SW_CTR_1.ip","SW_CTR_2.ip","PA_CTR_1.ip","PA_CTR_2.ip")])

CSTF2_OE_FC_SW_matrix2[which(is.na(CSTF2_OE_FC_SW_matrix2[,1])),1]=1
CSTF2_OE_FC_SW_matrix2[which(is.na(CSTF2_OE_FC_SW_matrix2[,2])),2]=1
CSTF2_OE_FC_SW_matrix2[which(is.na(CSTF2_OE_FC_SW_matrix2[,3])),3]=1
CSTF2_OE_FC_SW_matrix2[which(is.na(CSTF2_OE_FC_SW_matrix2[,4])),4]=1

CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2[,1]==Inf),1]=max(CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2[,1]!=Inf),1])
CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2[,2]==Inf),2]=max(CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2[,2]!=Inf),2])
CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2[,3]==Inf),3]=max(CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2[,3]!=Inf),3])
CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2[,4]==Inf),4]=max(CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2[,4]!=Inf),4])

CSTF2_OE_FC_SW_matrix2[,1]=log2(CSTF2_OE_FC_SW_matrix2[,1])
CSTF2_OE_FC_SW_matrix2[,2]=log2(CSTF2_OE_FC_SW_matrix2[,2])
CSTF2_OE_FC_SW_matrix2[,3]=log2(CSTF2_OE_FC_SW_matrix2[,3])
CSTF2_OE_FC_SW_matrix2[,4]=log2(CSTF2_OE_FC_SW_matrix2[,4])

CSTF2_OE_FC_SW_matrix2$type_SW="Non"
CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2[,1]>0.26 & CSTF2_OE_FC_SW_matrix2[,2] >0.26),"type_SW"]="Up"
CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2[,1]< -0.26 & CSTF2_OE_FC_SW_matrix2[,2] < -0.26),"type_SW"]="Down"


CSTF2_OE_FC_SW_matrix2$type_PA="Non"
CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2[,3]>0.26 & CSTF2_OE_FC_SW_matrix2[,4] >0.26),"type_PA"]="Up"
CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2[,3]< -0.26 & CSTF2_OE_FC_SW_matrix2[,4] < -0.26),"type_PA"]="Down"

temp.OE1=table(CSTF2_OE_FC_SW_matrix2$type_SW)
pre_SW=temp.OE1["Up"]/(temp.OE1["Down"]+temp.OE1["Up"])
temp.OE2=table(CSTF2_OE_FC_SW_matrix2$type_PA)
pre_PA=temp.OE2["Up"]/(temp.OE2["Down"]+temp.OE2["Up"])

CSTF2_OE_FC_SW_matrix2$Mean_SWOE=
  rowMeans(CSTF2_OE_FC_SW_matrix[rownames(CSTF2_OE_FC_SW_matrix2),
                                 c("SW_CS_OE_1.ip","SW_CS_OE_2.ip")])

CSTF2_OE_FC_SW_matrix2$Mean_PAOE=
  rowMeans(CSTF2_OE_FC_SW_matrix[rownames(CSTF2_OE_FC_SW_matrix2),
                                 c("PA_CS_OE_1.ip","PA_CS_OE_2.ip")])

CSTF2_OE_FC_SW_matrix2$Mean_SWcon=
  rowMeans(CSTF2_OE_FC_SW_matrix[rownames(CSTF2_OE_FC_SW_matrix2),
                                 c("SW_CTR_1.ip","SW_CTR_2.ip")])

CSTF2_OE_FC_SW_matrix2$Mean_PAcon=
  rowMeans(CSTF2_OE_FC_SW_matrix[rownames(CSTF2_OE_FC_SW_matrix2),
                                 c("PA_CTR_1.ip","PA_CTR_2.ip")])

text_plot=data.frame(table(CSTF2_OE_FC_SW_matrix2$type_SW))
text_plot=text_plot[which(text_plot$Var1!="Non"),]
text_plot$X.temp=NA
text_plot$Y.temp=NA
text_plot[which(text_plot$Var1=="Up"),c("X.temp","Y.temp")]=c(0.5,6)
text_plot[which(text_plot$Var1=="Down"),c("X.temp","Y.temp")]=c(6,0.7)

CSTF2_OE_FC_SW_matrix2$type_SW=ordered(CSTF2_OE_FC_SW_matrix2$type_SW,levels=c("Non","Down","Up"))
CSTF2_OE_FC_SW_matrix2=CSTF2_OE_FC_SW_matrix2[order(CSTF2_OE_FC_SW_matrix2$type_SW),]

ggplot(CSTF2_OE_FC_SW_matrix2)+geom_point(aes(x=log2(CSTF2_OE_FC_SW_matrix2$Mean_SWcon+1),
                                              y=log2(CSTF2_OE_FC_SW_matrix2$Mean_SWOE+1),
                                              color=CSTF2_OE_FC_SW_matrix2$type_SW),size=2.2)+
  geom_abline(intercept = 0,linetype="dashed")+theme_classic()+
  scale_color_manual(values=c(Down="#17489C",Non=random.color,Up="#E7423B"),name="type")+
  geom_text(data=text_plot,aes(x=text_plot$X.temp,y=text_plot$Y.temp,label=text_plot$Freq))
ggsave("analysis/new_figure/SW_CSTF2_OE_diff_m6A_point_2v2.pdf",width = 5,height = 4)

text_plot=data.frame(table(CSTF2_OE_FC_SW_matrix2$type_PA))
text_plot=text_plot[which(text_plot$Var1!="Non"),]
text_plot$X.temp=NA
text_plot$Y.temp=NA
text_plot[which(text_plot$Var1=="Up"),c("X.temp","Y.temp")]=c(0.5,7)
text_plot[which(text_plot$Var1=="Down"),c("X.temp","Y.temp")]=c(5,0.5)

CSTF2_OE_FC_SW_matrix2$type_PA=ordered(CSTF2_OE_FC_SW_matrix2$type_PA,levels=c("Non","Down","Up"))
CSTF2_OE_FC_SW_matrix2=CSTF2_OE_FC_SW_matrix2[order(CSTF2_OE_FC_SW_matrix2$type_PA),]

ggplot(CSTF2_OE_FC_SW_matrix2)+geom_point(aes(x=log2(CSTF2_OE_FC_SW_matrix2$Mean_PAcon+1),
                                              y=log2(CSTF2_OE_FC_SW_matrix2$Mean_PAOE+1),
                                              color=CSTF2_OE_FC_SW_matrix2$type_PA),size=2.2)+
  geom_abline(intercept = 0,linetype="dashed")+theme_classic()+
  scale_color_manual(values=c(Down="#17489C",Non=random.color,Up="#E7423B"),name="type")+
  geom_text(data=text_plot,aes(x=text_plot$X.temp,y=text_plot$Y.temp,label=text_plot$Freq))
ggsave("analysis/new_figure/PA_CSTF2_OE_diff_m6A_point_2v2.pdf",width = 5,height = 4)

length(which(CSTF2.peak.anno[unique(c(rownames(siCSTF2_FC_SW_matrix2[which(siCSTF2_FC_SW_matrix2$type_SW=="Down"),]))),"Gene"] %in%
               CSTF2.peak.anno[rownames(CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2$type_SW=="Up"),]),"Gene"]))

length(which(!(CSTF2.peak.anno[unique(c(rownames(siCSTF2_FC_SW_matrix2[which(siCSTF2_FC_SW_matrix2$type_SW=="Down"),]))),"Gene"] %in%
                 CSTF2.peak.anno[rownames(CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2$type_SW=="Up"),]),"Gene"])))

length(which(!(CSTF2.peak.anno[rownames(CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2$type_SW=="Up"),]),"Gene"] %in% 
                 CSTF2.peak.anno[unique(c(rownames(siCSTF2_FC_SW_matrix2[which(siCSTF2_FC_SW_matrix2$type_SW=="Down"),]))),"Gene"])))

length(which(CSTF2.peak.anno[unique(c(rownames(si_CSTF2_FC_matrix[which(si_CSTF2_FC_matrix$Color=="Down"),]))),"Gene"] %in%
               CSTF2.peak.anno[rownames(CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2$type_PA=="Up"),]),"Gene"]))

length(which(!(CSTF2.peak.anno[unique(c(rownames(si_CSTF2_FC_matrix[which(si_CSTF2_FC_matrix$Color=="Down"),]))),"Gene"] %in%
                 CSTF2.peak.anno[rownames(CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2$type_PA=="Up"),]),"Gene"])))

length(which(!(CSTF2.peak.anno[rownames(CSTF2_OE_FC_SW_matrix2[which(CSTF2_OE_FC_SW_matrix2$type_PA=="Up"),]),"Gene"] %in% 
                 CSTF2.peak.anno[unique(c(rownames(si_CSTF2_FC_matrix[which(si_CSTF2_FC_matrix$Color=="Down"),]))),"Gene"])))

S2_hyper=significan.diff.subtype.peak[which(significan.diff.subtype.peak$log2FC.paired>0),]
rownames(S2_hyper)=S2_hyper$Peak.id
S2_hyper$Peak_type=NA
S2_hyper[intersect(S2_hyper$Peak.id,hyper_peaks$Peak.id),"Peak_type"]="Hyper"
S2_hyper=S2_hyper[which(!is.na(S2_hyper$Peak_type)),]

length(which(anno[S2_hyper$Peak.id,"Gene"] %in%
               c(CSTF2.peak.anno[rownames(si_CSTF2_FC_matrix[which(si_CSTF2_FC_matrix$Color=="Down"),]),"Gene"],
                 CSTF2.peak.anno[rownames(siCSTF2_FC_SW_matrix2[which(siCSTF2_FC_SW_matrix2$type_SW=="Down"),]),"Gene"])))

length(which(!(anno[S2_hyper$Peak.id,"Gene"] %in%
               c(CSTF2.peak.anno[rownames(si_CSTF2_FC_matrix[which(si_CSTF2_FC_matrix$Color=="Down"),]),"Gene"],
                 CSTF2.peak.anno[rownames(siCSTF2_FC_SW_matrix2[which(siCSTF2_FC_SW_matrix2$type_SW=="Down"),]),"Gene"]))))

length(which(!(c(CSTF2.peak.anno[rownames(si_CSTF2_FC_matrix[which(si_CSTF2_FC_matrix$Color=="Down"),]),"Gene"],
                 CSTF2.peak.anno[rownames(siCSTF2_FC_SW_matrix2[which(siCSTF2_FC_SW_matrix2$type_SW=="Down"),]),"Gene"])%in%
               anno[S2_hyper$Peak.id,"Gene"])))



#m6A iCLIP
setwd("/data/xingyang/m6A_zhengjian/analysis/annotation/")
system("perl m6A_annotate_forGTF_xingyang_v2.pl /data/database/hg38/GENCODE/gencode.v25.annotation.gtf /data/xingyang/m6A_zhengjian/analysis/iCLIP_zhengjian/m6A_iCLIP_3sample.bed /data/xingyang/m6A_zhengjian/analysis/annotation/new_annotation/m6A_iCLIP")
setwd("/data/xingyang/m6A_zhengjian/")
m6A.iCLIP.anno=read.table("analysis/annotation/new_annotation/m6A_iCLIP.anno.txt",sep="\t")
rownames(m6A.iCLIP.anno)=m6A.iCLIP.anno$V4
colnames(m6A.iCLIP.anno)=colnames(anno)
m6A.iCLIP.anno[which(m6A.iCLIP.anno$Gene.type=="coding"),"Level.2.gene.type"]="mRNA"

m6A.iCLIP.unanno=read.table("analysis/annotation/new_annotation/m6A_iCLIP.unanno.txt",header = F,sep="\t")
colnames(m6A.iCLIP.unanno)[4]="Peak.id"
m6A.iCLIP.unanno=merge(m6A.iCLIP.unanno,m6A.iCLIP.anno,by="Peak.id",all.x=T)
m6A.iCLIP.unanno=m6A.iCLIP.unanno[,colnames(m6A.iCLIP.anno)]

m6A.iCLIP.anno=rbind(m6A.iCLIP.anno,m6A.iCLIP.unanno)


#CSTF2 HITS CLIP
setwd("/data/xingyang/m6A_zhengjian/analysis/annotation/")
system("perl m6A_annotate_forGTF_xingyang_v2.pl /data/database/hg38/GENCODE/gencode.v25.annotation.gtf /data/xingyang/m6A_zhengjian/analysis/iCLIP_zhengjian/PA_C_CIMS.FDR001.merge.bed /data/xingyang/m6A_zhengjian/analysis/annotation/new_annotation/CSTF2_HITS")
setwd("/data/xingyang/m6A_zhengjian/")
CSTF2.CLIP.anno=read.table("analysis/annotation/new_annotation/CSTF2_HITS.anno.txt",sep="\t")
rownames(CSTF2.CLIP.anno)=CSTF2.CLIP.anno$V4
colnames(CSTF2.CLIP.anno)=colnames(anno)
CSTF2.CLIP.anno[which(CSTF2.CLIP.anno$Gene.type=="coding"),"Level.2.gene.type"]="mRNA"

CSTF2.CLIP.unanno=read.table("analysis/annotation/new_annotation/CSTF2_HITS.unanno.txt",header = F,sep="\t")
colnames(CSTF2.CLIP.unanno)[4]="Peak.id"
CSTF2.CLIP.unanno=merge(CSTF2.CLIP.unanno,CSTF2.CLIP.anno,by="Peak.id",all.x=T)
CSTF2.CLIP.unanno=CSTF2.CLIP.unanno[,colnames(CSTF2.CLIP.anno)]

CSTF2.CLIP.anno=rbind(CSTF2.CLIP.anno,CSTF2.CLIP.unanno)

#CSTF2 CLIP
setwd("/data/xingyang/m6A_zhengjian/analysis/annotation/")
system("perl m6A_annotate_forGTF_xingyang_v2.pl /data/database/hg38/GENCODE/gencode.v25.annotation.gtf /data/xingyang/m6A_zhengjian/analysis/iCLIP_zhengjian/CSTF2_293T_POSTAR.bed /data/xingyang/m6A_zhengjian/analysis/annotation/new_annotation/CSTF2_293T_POSTAR")
setwd("/data/xingyang/m6A_zhengjian/")
CSTF2_eCLIP_293T.CLIP.anno=read.table("analysis/annotation/new_annotation/CSTF2_293T_POSTAR.anno.txt",sep="\t")
rownames(CSTF2_eCLIP_293T.CLIP.anno)=CSTF2_eCLIP_293T.CLIP.anno$V4
colnames(CSTF2_eCLIP_293T.CLIP.anno)=colnames(anno)

#HEK293T
setwd("/data/xingyang/m6A_zhengjian/analysis/annotation/")
system("perl m6A_annotate_forGTF_xingyang_v2.pl /data/database/hg38/GENCODE/gencode.v25.annotation.gtf /data/xingyang/m6A_zhengjian/analysis/iCLIP_zhengjian/hg38_GSE63753_miCLIP.bed /data/xingyang/m6A_zhengjian/analysis/annotation/new_annotation/miCLIP_293T")
setwd("/data/xingyang/m6A_zhengjian/")
miCLIP_293T.CLIP.anno=read.table("analysis/annotation/new_annotation/miCLIP_293T.anno.txt",sep="\t")
rownames(miCLIP_293T.CLIP.anno)=miCLIP_293T.CLIP.anno$V4
colnames(miCLIP_293T.CLIP.anno)=colnames(anno)


#HEK293T
HEK293T_miCLIP_vs_eCLIP_distribution=site_distance(miCLIP_293T.CLIP.anno,CSTF2_eCLIP_293T.CLIP.anno,gtf.temp)
HEK293T_miCLIP_vs_eCLIP_distribution[which(is.na(HEK293T_miCLIP_vs_eCLIP_distribution$Distance)),"Distance"]=
  as.integer(runif(length(HEK293T_miCLIP_vs_eCLIP_distribution[which(is.na(HEK293T_miCLIP_vs_eCLIP_distribution$Distance)),
                                                               "Distance"]),-100,100)*2000)

CSTF2_eCLIP_293T.CLIP.anno_2kb=cbind(strsplit2(HEK293T_miCLIP_vs_eCLIP_distribution$Peak.id,split=":|-"),
                                     HEK293T_miCLIP_vs_eCLIP_distribution$Peak.id)
CSTF2_eCLIP_293T.CLIP.anno_2kb[,2]=as.numeric(CSTF2_eCLIP_293T.CLIP.anno_2kb[,2])-2000
CSTF2_eCLIP_293T.CLIP.anno_2kb[,3]=as.numeric(CSTF2_eCLIP_293T.CLIP.anno_2kb[,3])+2000
CSTF2_eCLIP_293T.CLIP.anno_2kb[which(as.numeric(CSTF2_eCLIP_293T.CLIP.anno_2kb[,2])<0),2]=0
write.table(CSTF2_eCLIP_293T.CLIP.anno_2kb,"CSTF2_eCLIP_293T.CLIP.anno_2kb.bed",quote = F,
            row.names = F,col.names = F,sep="\t")

ggplot(HEK293T_miCLIP_vs_eCLIP_distribution)+
  geom_line(aes(x=HEK293T_miCLIP_vs_eCLIP_distribution$Distance),stat = "density",size=1.2,color=normal.color)+
  theme_classic()+
  xlim(-200,200)+
  xlab("Distance of 293T miCLIP from eCLIP site")+ylab("Density")
ggsave("analysis/new_figure/293T_miCLIP_vs_eCLIP.pdf",width = 5.5,height = 4.5)

fisher.test(matrix(c(6287,7584,1508,56269-6287-7584-1508),2,2),alternative = "greater")

venn.diagram(list(miCLIP=miCLIP_293T.CLIP.anno$Gene,
                  PAR_CLIP=CSTF2_eCLIP_293T.CLIP.anno$Gene),
             fill = c("red","blue"),cat.dist=c(0.05,0.05),cat.prompts=2,
             margin=0.1,ext.text=1,ext.line.lwd=0,ext.dist=c(-0.05),
             filename = 'analysis/new_figure/293T_miCLIP_vs_eCLIP_overlap.png',imagetype = 'png')


#m6A iCLIP vs HIST CLIP
CSTF2_m6A_m6A_iCLIP_HISTCLIP_distribution=site_distance(CSTF2.CLIP.anno,m6A.iCLIP.anno,gtf.temp)
CSTF2_m6A_m6A_iCLIP_HISTCLIP_distribution[which(is.na(CSTF2_m6A_m6A_iCLIP_HISTCLIP_distribution$Distance)),"Distance"]=
  as.integer(runif(length(CSTF2_m6A_m6A_iCLIP_HISTCLIP_distribution[which(is.na(CSTF2_m6A_m6A_iCLIP_HISTCLIP_distribution$Distance)),
                                                          "Distance"]),-100,100)*2000)



CSTF2.CLIP.anno_2kb=cbind(strsplit2(HEK293T_miCLIP_vs_eCLIP_distribution$Peak.id,split=":|-"),
                                     HEK293T_miCLIP_vs_eCLIP_distribution$Peak.id)
CSTF2.CLIP.anno_2kb[,2]=as.numeric(CSTF2.CLIP.anno_2kb[,2])-2000
CSTF2.CLIP.anno_2kb[,3]=as.numeric(CSTF2.CLIP.anno_2kb[,3])+2000
CSTF2.CLIP.anno_2kb[which(as.numeric(CSTF2.CLIP.anno_2kb[,2])<0),2]=0
write.table(CSTF2.CLIP.anno_2kb,"CSTF2.CLIP.anno_2kb.bed",quote = F,
            row.names = F,col.names = F,sep="\t")


ggplot(CSTF2_m6A_m6A_iCLIP_HISTCLIP_distribution)+geom_line(aes(x=CSTF2_m6A_m6A_iCLIP_HISTCLIP_distribution$Distance),stat = "density",size=1.2,color="navy")+
  theme_classic()+
  xlim(-2000,2000)+
  xlab("Distance of m6A iCLIP from CSTF2 HITS-CLIP site")+ylab("Density")
ggsave("analysis/new_figure/CSTF2_m6A_iCLIP_vs_HISTCLIP.pdf",width = 5.5,height = 4.5)

#Peaks overlap near m6A clip and HITS clip
fisher.test(matrix(c(2386,5535,1675,56269-5535-2386-1675),2,2),alternative = "greater")
venn.diagram(list(HITS_CLIP=na.omit(CSTF2.CLIP.anno$Gene),
                  iCLIP=na.omit(m6A.iCLIP.anno$Gene)
                  #m6A_seq=na.omit(CSTF2_down_reg_peaks_anno$Gene)
),
fill = c("red","blue"),cat.dist=c(0.05,0.05),cat.prompts=2,
margin=0.1,ext.text=1,ext.line.lwd=0,ext.dist=c(-0.05),total.population = 56269,hyper.test = T,
filename = 'analysis/new_figure/Peak_CLIP_overlap_all_2.png',imagetype = 'png')



#CSTF2 HITS CLIP topology
temp=CSTF2.CLIP.anno
colnames(temp)=colnames(anno)
temp=cbind(temp,type="CSTF2")
ggtopo(temp)+scale_color_manual(values = c("firebrick"))+guides(color=F)
ggsave("analysis/new_figure/CSTF2_HITS_CLIP_topology.pdf",width=5,heigh=4,dpi=500)

#m6A iCLIP topology
temp=m6A.iCLIP.anno
colnames(temp)=colnames(anno)
temp=cbind(temp,type="m6A")
ggtopo(temp)+scale_color_manual(values = c(normal.color))+guides(color=F)
ggsave("analysis/new_figure/m6A_iCLIP_topology.pdf",width=5,heigh=4,dpi=500)


#CSTF2 HITS CLIP annotation
pie_plot(CSTF2.CLIP.anno[,"Level.2.gene.type"])+
  scale_fill_manual(values=c(annotation_color7,annotation_color2,annotation_color3,annotation_color1,random.color))
ggsave("analysis/new_figure/Annotation_CSTF2_HITSCLIP_1_gene.pdf",width = 5,heigh=4,dpi=500)

pie_plot(CSTF2.CLIP.anno[which(CSTF2.CLIP.anno$Gene.type=="coding"),
                         ][,"Gene.site"])+
  scale_fill_manual(values=c(annotation_color4,annotation_color2,annotation_color1,annotation_color3))
ggsave("analysis/new_figure/Annotation_CSTF2_HITSCLIP_1_coding.pdf",width = 5,heigh=4,dpi=500)

#m6A iCLIP annotation
pie_plot(m6A.iCLIP.anno[,"Level.2.gene.type"])+
  scale_fill_manual(values=c(annotation_color7,annotation_color2,annotation_color1,annotation_color3,random.color))
ggsave("analysis/new_figure/Annotation_m6A_iCLIP_1_gene.pdf",width = 5,heigh=4,dpi=500)

pie_plot(m6A.iCLIP.anno[which(m6A.iCLIP.anno$Gene.type=="coding"),
                        ][,"Gene.site"])+
  scale_fill_manual(values=c(annotation_color2,annotation_color1,annotation_color4,annotation_color3))
ggsave("analysis/new_figure/Annotation_m6A_iCLIP_1_coding.pdf",width = 5,heigh=4,dpi=500)

#Stop codon annotation-CSTF2 CLIP
CSTF2_CLIP_anno=CSTF2.CLIP.anno
stop_codon_regoins_anno=coding_anno(annotation_table = CSTF2_CLIP_anno,gtf.temp = gtf.temp)
pie_plot(stop_codon_regoins_anno[,"Gene.site"])+
  scale_fill_manual(values=c(annotation_color4,annotation_color9,annotation_color1,annotation_color2,annotation_color3))
ggsave("analysis/new_figure/CSTF2_CLIP_coding_annotation_with_stop_codon.pdf",width = 5,height = 5)

tissue_peak_type_counting=table(stop_codon_regoins_anno$Gene.site)
coding_type_counting=coding_type_counting[names(tissue_peak_type_counting)]

i=1
n=length(tissue_peak_type_counting)
peak.enrichment=c()
while(i<=n){
  enrichment.temp=fisher.test(matrix(c(tissue_peak_type_counting[i],sum(tissue_peak_type_counting[-i]),
                                       coding_type_counting[i],sum(coding_type_counting[-i])),2,2))
  peak.enrichment=rbind(peak.enrichment,c(type=names(tissue_peak_type_counting)[i],score=as.numeric(enrichment.temp$estimate)))
  i=i+1
}
peak.enrichment=data.frame(peak.enrichment)
peak.enrichment$score=as.numeric(peak.enrichment$score)
peak.enrichment$type=ordered(peak.enrichment$type,levels=c("5UTR","CDS","intron","stop_codon","3UTR"))
ggplot(peak.enrichment)+
  geom_bar(aes(x=peak.enrichment$type,y=peak.enrichment$score,fill=peak.enrichment$type),stat="identity")+
  scale_fill_manual(values=c(annotation_color3,annotation_color2,annotation_color4,annotation_color9,annotation_color1),name="")+
  theme_classic()+scale_y_continuous(expand = c(0,0))+xlab("m6A peaks")+ylab("Peaks enrichment score")
ggsave("analysis/new_figure/CSTF2_CLIP_coding_annotation_enrichment_bar.pdf",width = 6,height = 5)


#Stop codon annotation
m6A_CLIP_anno=m6A.iCLIP.anno
stop_codon_regoins_anno=coding_anno(annotation_table = m6A_CLIP_anno,gtf.temp = gtf.temp)
pie_plot(stop_codon_regoins_anno[,"Gene.site"])+
  scale_fill_manual(values=c(annotation_color4,annotation_color9,annotation_color2,annotation_color1,annotation_color3))
ggsave("analysis/new_figure/m6A_CLIP_coding_annotation_with_stop_codon.pdf",width = 5,height = 5)


tissue_peak_type_counting=table(stop_codon_regoins_anno$Gene.site)
coding_type_counting=coding_type_counting[names(tissue_peak_type_counting)]

i=1
n=length(tissue_peak_type_counting)
peak.enrichment=c()
while(i<=n){
  enrichment.temp=fisher.test(matrix(c(tissue_peak_type_counting[i],sum(tissue_peak_type_counting[-i]),
                                       coding_type_counting[i],sum(coding_type_counting[-i])),2,2))
  peak.enrichment=rbind(peak.enrichment,c(type=names(tissue_peak_type_counting)[i],score=as.numeric(enrichment.temp$estimate)))
  i=i+1
}
peak.enrichment=data.frame(peak.enrichment)
peak.enrichment$score=as.numeric(peak.enrichment$score)
peak.enrichment$type=ordered(peak.enrichment$type,levels=c("5UTR","CDS","intron","stop_codon","3UTR"))
ggplot(peak.enrichment)+
  geom_bar(aes(x=peak.enrichment$type,y=peak.enrichment$score,fill=peak.enrichment$type),stat="identity")+
  scale_fill_manual(values=c(annotation_color3,annotation_color2,annotation_color4,annotation_color9,annotation_color1),name="")+
  theme_classic()+scale_y_continuous(expand = c(0,0))+xlab("m6A peaks")+ylab("Peaks enrichment score")
ggsave("analysis/new_figure/m6A_CLIP_coding_annotation_enrichment_bar.pdf",width = 6,height = 5)


#Pol II tag density
temp.bed=cbind(CSTF2_eCLIP_293T.CLIP.anno$Chr,CSTF2_eCLIP_293T.CLIP.anno$Start-2000,CSTF2_eCLIP_293T.CLIP.anno$End+2000,
               CSTF2_eCLIP_293T.CLIP.anno$Peak.id)
temp.bed[,2]=as.numeric(temp.bed[,2])
temp.bed[,3]=as.numeric(temp.bed[,3])
write.table(temp.bed,"temp.bed",row.names = F,col.names = F,sep="\t",quote = F)
system("/data/software/liftover/liftOver  temp.bed /data/software/liftover/hg38ToHg19.over.chain temp.hg19.bed unmapping")
#system("/data/software/bwtool-master/bwtool extract jsp temp.bed GSM3272322_Pol_II_hg38.bw Pol_II.txt")
system("/data/software/bwtool-master/bwtool dist GSM3272322_Pol_II.bw Pol_II_dist.txt -regions=temp.hg19.bed")
Pol_II_tag=fread("Pol_II_dist.txt",sep="\t",header = F)
Pol_II_tag=data.frame(Pol_II_tag)

ggplot(Pol_II_tag)+
  geom_line(aes(x=Pol_II_tag[,1],y=Pol_II_tag[,2]),stat = "identity",color="navy",size=1)+
  ylab("Reads count")+xlab("Distance of CSTF2 target sites from Pol II sites in HEK293 cell line")+
  theme_classic()
ggsave("analysis/new_figure/CSTF2_Pol_II.pdf",width = 5.5,height = 4.5)

#Spearman corr
S2_hyper_IGF2BP2_peak_corr=gene_corr_peak(peak_for_analysis = S2_hyper$Peak.id,
                                  gene_for_analysis = "IGF2BP2",
                                  peak_matrix = rpkm_ip_norm_t[,paired_sample_tumor],
                                  gene_matrix = gene_exp_data_t[,paired_sample_tumor])

IGF2BP2_reg_hyper_up.uniq=unique(S2_hyper_IGF2BP2_peak_corr[which(abs(S2_hyper_IGF2BP2_peak_corr$Corr)>0.25 & 
                                                              S2_hyper_IGF2BP2_peak_corr$Pvalue<0.05),"Peak.id"])
