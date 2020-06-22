load(file="df.normal3.rdata")
load(file="df.lognormal3.rdata")
load(file="df.outlier3.rdata")
df<-rbind(data.frame(case="A: normal",df.normal),
          data.frame(case="B: outliers",df.outlier),
          data.frame(case="C: lognormal",df.lognormal))
names(df)[5]<-"Methods"
df$Methods<-factor(df$Methods, levels=c("none","median","huber","adhuber"))

library(plyr)
library(ggplot2)
df$Methods<-mapvalues(df$Methods,from=c("none","median","huber","adhuber"), to=c("None-robust","Median-truncated","Huber","ad-Huber"))

df<-df[df$measures!="FPs",]

pdf(file="plot.3cases.pdf", width = 10, height=10)
ggplot(df,aes(x,y,group)) +
  geom_line(aes(color=Methods,linetype=Methods), size=1)+
  #geom_point(aes(color=group, shape=group))+
  xlab("Number of selected genes") + ylab("Measures") +
  facet_grid(measures ~ case, scales="free") +
  theme_bw()
dev.off()
