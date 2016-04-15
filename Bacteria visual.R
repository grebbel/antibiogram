setwd("~/Documents/workspace/R-project/Antibiogram")


library(xlsx)
library(reshape2)
library(dplyr)
library(scales)
library(ggplot2)


## read data 
newDF3 <- read.csv("newDF3.csv", sep = ",")


Plot7 <- newDF3[newDF3$n > 6,]



### new dataset without SOURCE, aggregated again ####
Plot7 <- aggregate(Plot7[,c("S","R","n")], by=list(Plot7$Bacteria, Plot7$Antibiotics), "sum")

names(Plot7)[names(Plot7)=="Group.1"] <- "Bacteria"
names(Plot7)[names(Plot7)=="Group.2"] <- "Antibiotics"

Plot7$Resis <- ((Plot7$R/Plot7$n)*100)

Plot7$Resis <- formatC(Plot7$Resis, digits=3)
Plot7$Resis <- as.numeric(Plot7$Resis)
  

sum(newDF3$n)
sum(Plot7$n)

### MAKE PLOT

Plot7g <- ggplot(Plot7) + geom_point(aes(x = Bacteria, y = Antibiotics, colour = Resis, size = n)) +
  scale_colour_gradient(name = '% resistance', high = "red", low = "green") +
  scale_size(name = 'sample size', range = c(2,9)) +
  xlab('Species') + ylab('Antibiotic') +
  theme(text = element_text(size = 12), axis.text.x = element_text(angle=60, hjust=1, size=12) )

Plot7g



