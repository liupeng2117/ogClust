
#--------------------------------------------------------------------#
#----------------------- bootstrap result ---------------------------#
load("D:/research/Disease subtyping/application/bootstrap/lung_disease_our_top2000_elastic_fast_bootstrap_alpha0.3_B500.Rdata")
#--------------------------------------------------------------------#
genenames<-names(g.sd[g.sd>=cut.sd])
genes.id.matrix<-matrix(,ncol=2000,nrow=500) #column = number of genes
for(i in 1:500){
  gamma_est=boot.res[[i]]$theta_est[(np+1):((K-1)*(NG+1)+np)]
  gamma_est_matrix=matrix(gamma_est,ncol=K-1,byrow=TRUE)
  genes.id.matrix[i,]<-apply(gamma_est_matrix[-1,,drop=F],1, function(x) sum(x!=0)>0)
}
apply(genes.id.matrix,1,sum) # by repitation
apply(genes.id.matrix,2,sum) # by gene
counts<-apply(genes.id.matrix,2,sum)
bootCounts<-counts[order(counts,decreasing = T)]
sum(counts!=0)
names(bootCounts)<-genenames[order(counts,decreasing = T)]
bootCountsdf<-data.frame(names(bootCounts),bootCounts)
colnames(bootCountsdf)<-c("gene","count")
write.csv(bootCountsdf,file="D:/research/Disease subtyping/application/bootstrap/genes selected/bootCounts2000_B500_alpha0.3.csv",row.names=F)

#--------------------------------------------------------------------#
load("D:/research/Disease subtyping/application/bootstrap/lung_disease_skm_top2000_bootstrap_B500.Rdata")
genenames<-names(g.sd[g.sd>=cut.sd])
genes.id.matrix<-matrix(,ncol=2000,nrow=100) #column = number of genes
for(i in 1:100){
  weight=boot.res[[i]]$weight
  genes.id.matrix[i,]<- (weight!=0)
}
apply(genes.id.matrix,1,sum) # by repitation
apply(genes.id.matrix,2,sum) # by gene
counts<-apply(genes.id.matrix,2,sum)
bootCounts<-counts[order(counts,decreasing = T)]
sum(counts!=0)
names(bootCounts)<-genenames[order(counts,decreasing = T)]
bootCountsdf<-data.frame(names(bootCounts),bootCounts)
colnames(bootCountsdf)<-c("gene","count")
write.csv(bootCountsdf,file="D:/research/Disease subtyping/application/bootstrap/genes selected/bootCounts2000_B500_skm.csv",row.names=F)

#--------------------------------------------------------------------#
load("D:/research/Disease subtyping/application/bootstrap/lung_disease_scluster_top2000_bootstrap_B500.Rdata")
genenames<-names(g.sd[g.sd>=cut.sd])
genes.id.matrix<-matrix(,ncol=2000,nrow=500) #column = number of genes
for(i in 1:500){
  genes_selected=boot.res[[i]]$genes_selected
  genes.id.matrix[i,]<- genenames %in% genes_selected
}
apply(genes.id.matrix,1,sum) # by repitation
apply(genes.id.matrix,2,sum) # by gene
counts<-apply(genes.id.matrix,2,sum)
bootCounts<-counts[order(counts,decreasing = T)]
sum(counts!=0)
names(bootCounts)<-genenames[order(counts,decreasing = T)]
bootCountsdf<-data.frame(names(bootCounts),bootCounts)
colnames(bootCountsdf)<-c("gene","count")
write.csv(bootCountsdf,file="D:/research/Disease subtyping/application/bootstrap/genes selected/bootCounts2000_B500_scluster.csv",row.names=F)

#----------------------------------------------------------------------#