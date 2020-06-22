load(file="df.delta1.gamma1.rdata")
load(file="df.delta1.gamma3.rdata")
load(file="df.delta2.gamma1.rdata")
load(file="df.delta2.gamma3.rdata")
df<-rbind(data.frame(case="A",df.delta1.gamma1),
          data.frame(case="B",df.delta1.gamma3),
          data.frame(case="C",df.delta2.gamma1),
          data.frame(case="D",df.delta2.gamma3))

library(dplyr)
library(ggplot2)
colnames(df)[5]<-"Methods"
df<-df[df$measures!="FPs",]
levels(df$case)<-c(expression(paste("A: ",gamma,"=1, ",delta,"=1")),
                expression(paste("B: ",gamma,"=3, ",delta,"=1")),
                expression(paste("C: ",gamma,"=1, ",delta,"=2")),
                expression(paste("D: ",gamma,"=3, ",delta,"=2")))
pdf(file="plot.4cases.pdf", width = 10, height=10)
ggplot(df,aes(x,y,group)) +
  geom_line(aes(color=Methods,linetype=Methods), size=1)+
  #geom_point(aes(color=group, shape=group))+
  xlab("Number of selected genes") + ylab("Measures") +
  facet_grid(measures ~ case, scales="free", labeller = label_parsed) +
  theme_bw()
dev.off()
