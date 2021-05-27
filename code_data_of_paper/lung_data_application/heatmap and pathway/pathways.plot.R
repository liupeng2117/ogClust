#-----------
#Plot of top10 pathway results
setwd("D:/research/Disease subtyping/application/bootstrap/pathways")
top=10
library(readxl)
pathways.skm<-read_excel("pathway_skm_top2000_B500.xls")
pathways.skm.top<-pathways.skm[top:1,c(1,2)]
pathways.scluster<-read_excel("pathway_scluster_top2000_B500.xls")
pathways.scluster.top<-pathways.scluster[top:1,c(1,2)]
pathways.our<-read_excel("pathway_our_top2000_B500.xls")
pathways.our.top<-pathways.our[top:1,c(1,2)]

break.data<-matrix(,2,3)
colnames(break.data)<-c("pathways","log10p","method")
break.data[,1]=c("","")

pathways.skm.top<-cbind(pathways.skm.top,method=rep("sparse K-means",top))
pathways.scluster.top<-cbind(pathways.scluster.top,method=rep("supervised clustering",top))
pathways.our.top<-cbind(pathways.our.top,method=rep("hiearchical mixture regression",top))
colnames(pathways.skm.top)<-c("pathways","log10p","method")
colnames(pathways.scluster.top)<-c("pathways","log10p","method")
colnames(pathways.our.top)<-c("pathways","log10p","method")
pathways.data<-rbind(pathways.skm.top, break.data, pathways.scluster.top, break.data, pathways.our.top)
pathways.data$log10p<-as.numeric(pathways.data$log10p)

ggplot(data=pathways.data, aes(x=1:34, y=log10p,fill=method)) + 
  geom_bar(stat = "identity",
           position = position_dodge(0.9)) +
  scale_fill_manual(values = alpha(c("blue","orange","purple"),0.7)) +
  scale_x_continuous(breaks=c(1:34),labels=pathways.data$pathways) +
  coord_flip() + 
  xlab("")+
  ylab(expression(paste("-log"[10],"(p-value)")))+
  theme(legend.position = "none",
        #axis.text.x = element_text(face = "bold", size=14),
        #axis.text.y = element_text(face = "bold", size=14),
        #axis.title.x = element_text(colour = "#6C8EBF", face = "bold", size =16),
        #axis.title.y = element_text(colour = "#6C8EBF", face = "bold", size =16),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(size = 1, colour="#DAE8FC"),
        panel.grid.minor.x = element_line(colour="#DAE8FC"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.text = element_text(size=12, face="bold")
        ) 