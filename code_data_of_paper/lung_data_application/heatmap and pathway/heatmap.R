# ------------------------- K=3 ---------------------------------------#
load("D:/research/Disease subtyping/application/lung_disease_our_top2000_elastic_fast_K3.Rdata")
#----------------------------------------------------------------------#
# choose lambda=0.06897959
# grp_assign<-apply(pai_est,1,which.max)
# grp_assign<-ifelse(grp_assign==2,3,ifelse(grp_assign==3,2,1))
# heatmap ----------------------------------------------------------
table(grp_diag, grp_assign)
#          grp_assign
#grp_diag   1   2   3
#       1   7  68  59
#       2 113  61   7
genenames<-names(g.sd[g.sd>=cut.sd])
colnames(gdata2)<-genenames
selected.genes<-genenames[apply(gamma_est_matrix[-1,],1,function(x) sum(abs(x)))>0.0024]
selected.genes.expr<-gdata2[,selected.genes]
sample.index<-c(which(grp_assign==1),which(grp_assign==2),which(grp_assign==3))
gdata.hm<-t(selected.genes.expr[sample.index,])
annotation_col<-data.frame(Groups=factor(grp_assign[sample.index]))
rownames(annotation_col)<-colnames(gdata.hm)
#gaps=cumsum(c(
#  sum(annotation_col$Groups==1),
#  sum(annotation_col$Groups==2),
#  sum(annotation_col$Groups==3)
#))
#colors
Labels<-annotation_col$Groups
colors<-rainbow(length(unique(Labels)))
names(colors)<-levels(Labels)

#heatmap
library(pheatmap)
library(gplots)

gdata.hm[gdata.hm>3]<-3
gdata.hm[gdata.hm< -3]<- -3


pdf("D:/research/Disease subtyping/application/genes selected/heatmap_our_K3.pdf", width = 10, height= 6)
out<-pheatmap(gdata.hm, cluster_rows = T,cluster_cols = F, show_rownames = F, show_colnames = F,border_color = NA, 
              color=greenred(75),treeheight_row = 0, treeheight_col = 0, labels_col="Samples",
              annotation =annotation_col, annotation_colors=list(Groups=colors))
dev.off()

# par gamma --------------------------------------------------------------------
genes_gamma_matrix<-gamma_est_matrix[-1,]
rownames(genes_gamma_matrix)<-genenames
#selected.genes.ordered<-rev(c("DIO2","TNN","PRMT8","SEC14L3","SLITRK6","HP","CA3","DAPL1","ZBTB16","DEFA3","PF4","SERPINE1","JUNB","DUSP2","ARC","S100A12","CMTM2","IL8RA","CLC","TMEM130","FGG","AZGP1"))
selected.genes.ordered<-rownames(genes_gamma_matrix)[apply(gamma_est_matrix,1,sum)!=0]
genes_selected_gamma_matrix<-genes_gamma_matrix[selected.genes.ordered,]
gamma1=data.frame(genes=selected.genes.ordered, gamma=genes_selected_gamma_matrix[,1])
gamma2=data.frame(genes=selected.genes.ordered, gamma=genes_selected_gamma_matrix[,2])
gamma.data<-rbind(gamma1,gamma2)
gamma.data<-data.frame(gamma.data,gamma_n=c(rep(1,length(selected.genes.ordered)),rep(2,length(selected.genes.ordered))))
gamma.data$gamma_n<-as.factor(gamma.data$gamma_n)
gamma.data$genes<-factor(gamma.data$genes, levels=selected.genes.ordered)
gamma.data$gamma<-round(gamma.data$gamma,digits = 2)
library(ggplot2)
theme_set(theme_bw()) 

# Plot
pdf("D:/research/Disease subtyping/application/genes selected/heatmap_our_K3_gamma_dotplot.pdf", width =4, height= 10)
ggplot(gamma.data, aes(x=genes, y=gamma, label=gamma)) + 
  geom_point(stat='identity', aes(col=gamma_n), size=6)  +
  scale_color_manual(name="gamma_k", 
                     labels = c("k=1", "k=2"), 
                     values = c("1"="#00ba38", "2"="#f8766d")) + 
  geom_text(color="white", size=2) +
  #labs(title="Diverging Dot Plot", 
  #     subtitle="Normalized mileage from 'mtcars': Dotplot") + 
  ylim(-0.5, 0.5) +
  coord_flip()
dev.off()
pdf("D:/research/Disease subtyping/application/genes selected/heatmap_our_K3_gamma_barplot.pdf", width =4, height= 10)
ggplot(gamma.data, aes(x=genes, y=gamma, label=gamma)) + 
  geom_bar(stat='identity', aes(fill=gamma_n), width=.5)  +
  scale_fill_manual(name="gamma_k", 
                    labels = c("k=1", "k=2"), 
                    values = c("1"="#00ba38", "2"="#f8766d")) + 
  #  labs(subtitle="Normalised mileage from 'mtcars'", 
  #       title= "Diverging Bars") + 
  coord_flip()
dev.off()
# pie chart --------------------------------------------------------------------
count.table<-table(grp_diag, grp_assign)
prop.table<-sapply(1:3,function(k) count.table[,k]/sum(count.table[,k]))

par(mfrow=c(1,1))
# Pie Chart with Percentages
for(k in 1:3){
  slices <- prop.table[,k]
  lbls <- c("COPD", "ILD")
  pct <- round(slices/sum(slices)*100)
  lbls <- paste(lbls, pct) # add percents to labels
  lbls <- paste(lbls,"%",sep="") # ad % to labels
  pie(slices, col=c("#00ba38", "#f8766d")) 
}

library(dplyr)
pie<-list()
for(k in 1:3){
  count <- count.table[,k]
  slices <- prop.table[,k]
  pct <- round(slices/sum(slices),digits=2)
  diagnosis <- c("COPD", "ILD")
  data=data.frame(diagnosis,count,pct)
  data <- data %>%
    arrange(desc(diagnosis)) %>%
    mutate(lab.ypos = cumsum(pct) - 0.5*pct)
  mycols=c("#0073C2FF","#CD534CFF")
  pie[[k]]<-ggplot(data, aes(x=1, y=pct,fill=diagnosis))+
    geom_bar(width = 1,stat = "identity",color="white")+
    coord_polar("y")+
    geom_text(aes(y = lab.ypos, label = pct), color = "white")+
    scale_fill_manual(name = "Diagnosis",values = mycols) +
    theme_void()+
    theme(legend.position = "none")
}
gridExtra::marrangeGrob(pie, nrow=1, ncol=3)
grp_diag2<-LungData.all$diagnosis$diagnosis[index.Y]
count.table<-table(grp_diag2, grp_assign)
prop.table<-sapply(1:3,function(k) count.table[,k]/sum(count.table[,k]))

par(mfrow=c(1,3))
# Pie Chart with Percentages
for(k in 1:3){
  slices <- prop.table[,k]
  lbls <- c("COPD", "Fibrosis-Uncharacterized","Hypersensitivity Pneumonitis","Idiopathic UIP","NSIP","Other ILD","Respiratory bronchiolitis")
  pct <- round(slices/sum(slices)*100)
  lbls <- paste(lbls, pct) # add percents to labels
  lbls <- paste(lbls,"%",sep="") # ad % to labels
  pie(slices,labels = lbls, col=rainbow(length(slices))) 
}

grp_diag3<-LungData.all$diagnosis$Major_Diagnosis_Final_Pathological[index.Y]
count.table<-table(grp_diag3, grp_assign)
prop.table<-sapply(1:3,function(k) count.table[,k]/sum(count.table[,k]))

par(mfrow=c(1,3))
# Pie Chart with Percentages
for(k in 1:3){
  slices <- prop.table[,k]
  lbls <- c("ILD", "COPD/Emphysema","Control","Other Lung")
  pct <- round(slices/sum(slices)*100)
  lbls <- paste(lbls, pct) # add percents to labels
  lbls <- paste(lbls,"%",sep="") # ad % to labels
  pie(slices,labels = lbls, col=rainbow(length(slices))) 
}

grp_diag4<-LungData.all$diagnosis$Severity_classification_of_DLCO_abnormality_Categorical[index.Y]
count.table<-table(grp_diag4, grp_assign)
prop.table<-sapply(1:3,function(k) count.table[,k]/sum(count.table[,k]))

par(mfrow=c(1,3))
# Pie Chart with Percentages
for(k in 1:3){
  slices <- prop.table[,k]
  lbls <- c("missing","Advanced disease", "Limited disease")
  pct <- round(slices/sum(slices)*100)
  lbls <- paste(lbls, pct) # add percents to labels
  lbls <- paste(lbls,"%",sep="") # ad % to labels
  pie(slices,labels = lbls, col=rainbow(length(slices))) 
}

#----------------------------------------------------------------------#
load("D:/research/Disease subtyping/application/lung_disease_skm_top2000_K3.Rdata")
#----------------------------------------------------------------------#
# choose lambda=0.06897959
# grp_assign<-apply(pai_est,1,which.max)
grp_assign<-ifelse(grp_assign==3,1,ifelse(grp_assign==1,2,3))
# heatmap ----------------------------------------------------------
table(grp_diag, grp_assign)
#           grp_assign
#grp_diag   1   2   3
#       1   6  81  47
#       2 123  58   0
genenames<-names(g.sd[g.sd>=cut.sd])
colnames(gdata2)<-genenames
selected.genes<-weight!=0
selected.genes.expr<-gdata2[,selected.genes]
sample.index<-c(which(grp_assign==1),which(grp_assign==2),which(grp_assign==3))
gdata.hm<-t(selected.genes.expr[sample.index,])
annotation_col<-data.frame(Groups=factor(grp_assign[sample.index]))
rownames(annotation_col)<-colnames(gdata.hm)
#gaps=cumsum(c(
#  sum(annotation_col$Groups==1),
#  sum(annotation_col$Groups==2),
#  sum(annotation_col$Groups==3)
#))
#colors
Labels<-annotation_col$Groups
colors<-rainbow(length(unique(Labels)))
names(colors)<-levels(Labels)

#heatmap
library(pheatmap)
library(gplots)

gdata.hm[gdata.hm>3]<-3
gdata.hm[gdata.hm< -3]<- -3


pdf("D:/research/Disease subtyping/application/genes selected/heatmap_skm_K3.pdf", width = 10, height= 6)
out<-pheatmap(gdata.hm, cluster_rows = T,cluster_cols = F, show_rownames = F, show_colnames = F,border_color = NA, 
              color=greenred(75),treeheight_row = 0, treeheight_col = 0, labels_col="Samples",
              annotation =annotation_col, annotation_colors=list(Groups=colors))
dev.off()

#----------------------------------------------------------------------#
load("D:/research/Disease subtyping/application/lung_disease_scluster_top2000_K3.Rdata")
#----------------------------------------------------------------------#
# choose lambda=0.06897959
# grp_assign<-apply(pai_est,1,which.max)
grp_assign<-ifelse(grp_assign==2,3,ifelse(grp_assign==3,2,1))
# heatmap ----------------------------------------------------------
table(grp_diag, grp_assign)
#        grp_assign
#grp_diag  1  2  3
#       1 14 47 73
#       2 76 72 33
genenames<-names(g.sd[g.sd>=cut.sd])
colnames(gdata2)<-genenames
selected.genes<-rep(FALSE,length(genenames))
selected.genes[order(ps)[1:final.ngs]]<-TRUE
selected.genes.expr<-gdata2[,selected.genes]
sample.index<-c(which(grp_assign==1),which(grp_assign==2),which(grp_assign==3))
gdata.hm<-t(selected.genes.expr[sample.index,])
annotation_col<-data.frame(Groups=factor(grp_assign[sample.index]))
rownames(annotation_col)<-colnames(gdata.hm)
#gaps=cumsum(c(
#  sum(annotation_col$Groups==1),
#  sum(annotation_col$Groups==2),
#  sum(annotation_col$Groups==3)
#))
#colors
Labels<-annotation_col$Groups
colors<-rainbow(length(unique(Labels)))
names(colors)<-levels(Labels)

#heatmap
library(pheatmap)
library(gplots)

gdata.hm[gdata.hm>3]<-3
gdata.hm[gdata.hm< -3]<- -3


pdf("D:/research/Disease subtyping/application/genes selected/heatmap_our_K3.pdf", width = 10, height= 6)
out<-pheatmap(gdata.hm, cluster_rows = T,cluster_cols = F, show_rownames = F, show_colnames = F,border_color = NA, 
              color=greenred(75),treeheight_row = 0, treeheight_col = 0, labels_col="Samples",
              annotation =annotation_col, annotation_colors=list(Groups=colors))
dev.off()
