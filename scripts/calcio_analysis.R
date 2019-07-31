##############################
# Jessica Calcium analysis #
##############################

# Obs: Quantities

#Necessary libraries

library(readODS)
library(ggplot2)

###############
# Import data #
setwd("~/Documentos/Doutorado/NeuroCell/Jéssica/raw_data/")
output <- "~/Documentos/Doutorado/NeuroCell/Jéssica/output/"
chol <- read.table(url("http://assets.datacamp.com/blog_assets/chol.txt"),header = T)

###
# Raw Data table
table <- data.frame(read.ods("momento.ods",sheet = 1))
table <- data.frame(momento=as.numeric(gsub(",", ".", gsub("\\.", "", table[c(-1,-19,-21,-28),]))))

table_2 <- data.frame(read.ods("momento.ods",sheet = 1))
table_2 <- data.frame(momento=as.numeric(gsub(",", ".", gsub("\\.", "", table_2[-1,]))))
table_2$momento[table_2$momento > 45] <- 45.1
table_2$momento[table_2$momento < -35] <- -35

serif <- element_text(face = "bold", colour = "black", size = 12)
legenda <- element_text(face= "bold", colour = "black", size = 14)

ggplot(table_2, aes(momento)) + labs(y="% neurons",x=NULL) +
  theme(panel.background = element_rect(fill="white"),
        axis.line.x = element_line(color = "black", size = 1),
        axis.line.y = element_line(color = "black", size = 1),
        axis.text = serif, title= legenda) +
  geom_histogram(aes(y= (..density..)),breaks=c(seq(-45,55,by=10)),
                 col="black",
                 fill="white", size = 1.5) +
  geom_density() +
  stat_function(fun=dnorm, color="red",
                args=list(mean=mean(x), 
                          sd=sd(x))) +
  scale_x_continuous(limits = c(-45,55),breaks = c(seq(-45,55,by=10)),
                     labels = c(" ","> -35",seq(-25,35,by=10),"< 45"," ")) +
  scale_y_continuous(expand = c(0,0))

ggplot(table_2, aes(momento)) +labs(y="Cumulative distribution",x=NULL) +
  theme(panel.background = element_rect(fill="white"),
        axis.line.x = element_line(color = "black", size = 1),
        axis.line.y = element_line(color = "black", size = 1),
        axis.text = serif, title= legenda) +
  stat_bin(aes(y=cumsum(..count../sum(..count..))),geom = "step",size=1) +
  scale_x_continuous(limits = c(-45,55),breaks = c(seq(-45,55,by=10)),
                     labels = c(" ","> -35",seq(-25,35,by=10),"< 45"," ")) +
  scale_y_continuous(expand = c(0,0))
