require(ggplot2)

HvsK = as.data.frame(read.table(paste(getwd(),"/HvsK.txt",sep =""),header=FALSE))
colnames(HvsK) = c("Date","Kilian Index", "IP Index")
HvsK$Date = as.character(HvsK$Date)

K = ggplot(data = HvsK, aes(x = Date, y= `Kilian Index`,group = 1))+
  geom_line()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs('Killian Index',y = NULL, x = NULL) + theme(plot.title = element_text(hjust=0.5))+
  scale_x_continuous(expand = c(0, 0))+
  geom_hline(yintercept=0,size = 0.5, alpha = 0.5, linetype = "dashed")
  