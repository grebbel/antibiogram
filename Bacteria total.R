setwd("~/Documents/workspace/R-project/Antibiogram")


library(xlsx)
library(reshape2)
library(dplyr)
library(scales)
library(ggplot2)

## IMPORT DATA###
newDF2 <- read.csv("newDF2.csv", sep = ",")

newDF2 <- subset(newDF2, select=c(Bacteria, Source, Antibiotics, S, I, R, N))

sum(newDF2$N)
######## Do an analysis on the total data set, 
######## before bacteria species are grouped and filtered in file bacteria resis2.R


####### bacteria number per source ######

total_n <- aggregate(newDF2[,c("S","R","N")], by=list(newDF2$Bacteria, newDF2$Source), "sum")

names(total_n)[names(total_n)=="Group.1"] <- "Bacteria"
names(total_n)[names(total_n)=="Group.2"] <- "Source"
total_n = dcast(total_n, Bacteria ~ Source, value.var = "N")
View(total_n)
write.xlsx(total_n, "total_n.xlsx")


### new dataset without SOURCE, aggregated again ####
newDF4 <- aggregate(newDF2[,c("S","R","N")], by=list(newDF2$Bacteria, newDF2$Antibiotics), "sum")

names(newDF4)[names(newDF4)=="Group.1"] <- "Bacteria"
names(newDF4)[names(newDF4)=="Group.2"] <- "Antibiotics"

newDF4$n <- (newDF4$R + newDF4$S)
newDF4$Resis <- ((newDF4$R/newDF4$n)*100)
newDF4$Resis <- formatC(newDF4$Resis, digits=3)
newDF4$Resis <- as.numeric(newDF4$Resis)


# make pivot table Bacterie ~ Antibiotics for % Resistance 
basic_summ_t = dcast(newDF4, Bacteria ~ Antibiotics, value.var = "Resis")
View(basic_summ_t)

basic_summ_n2 = dcast(newDF4, Bacteria ~ Antibiotics, value.var = "n")
View(basic_summ_n2)
write.xlsx(basic_summ_n2, "basic_sum_n.xlsx")


#### To make table containing 'resistance' TOGETHER WITH 'n'
# first melt the data frame to put all the metrics in a single column
basic_summ_t2 = melt(newDF4, id.vars = c("Bacteria","Antibiotics"), 
                     measure.vars = c("Resis", "n"))
View(basic_summ_t2)

# then transpose the quality and variable co
basic_summ_t2 = dcast(basic_summ_t2, Bacteria ~ Antibiotics + variable, value.var = "value")
View(basic_summ_t2)

##write.xlsx(basic_summ_t2, "basic_summ_t2.xlsx")



### remove incomplete cases ###

Bactresis <- newDF4[complete.cases(newDF4), ]

# make pivot table Bacterie ~ Antibiotics
Bactresis_pivot = dcast(Bactresis, Bacteria ~ Antibiotics, value.var = "n")
View(Bactresis_pivot)

summary(Bactresis_pivot)

##write.xlsx(Bactresis_pivot, "Bactresis_pivot.xlsx")

### select cases n >6 ####

Bactresis6 <- Bactresis[Bactresis$n > 6,]


# make pivot table >6 Bacterie ~ Antibiotics
Bactresis_pivot6 = dcast(Bactresis6, Bacteria ~ Antibiotics, value.var = "Resis")
View(Bactresis_pivot6)

write.xlsx(Bactresis_pivot6, "Bactresis_pivot6.xlsx")


# first melt the data frame to put all the metrics in a single column n>6
basic_summ_n6 = melt(Bactresis6, id.vars = c("Bacteria","Antibiotics"), 
                     measure.vars = c("Resis", "n"))
View(basic_summ_n6)

# then transpose the quality and variable co
basic_summ_n6 = dcast(basic_summ_n6, Bacteria ~ Antibiotics + variable, value.var = "value")
View(basic_summ_n6)

write.xlsx(basic_summ_n6, "basic_summ_n6.xlsx")


#################### VISUALIZATIONS ###################

###### PLOT #####
ggplot(data = melt(Bactresis_pivot6), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable))

#### TABLE ######
## Select all obs.##
Bactoplot <- ggplot(Bactresis) + geom_point(aes(x = Bacteria, y = Antibiotics, colour = Resis, size = n)) +
  scale_colour_gradient(name = '% resistance', high = "red", low = "green") +
  scale_size(name = 'sample size', range = c(2,6)) +
  xlab('Species') + ylab('Antibiotic') +
  theme(text = element_text(size = 10), axis.text.x = element_text(angle=60, hjust=1, size=10) )

Bactoplot

## Select n>6 ##

Bactoplot6 <- ggplot(Bactresis6) + geom_point(aes(x = Bacteria, y = Antibiotics, colour = Resis, size = n)) +
  scale_colour_gradient(name = '% resistance', high = "red", low = "green") +
  scale_size(name = 'sample size', range = c(2,10)) +
  xlab('Species') + ylab('Antibiotic') +
  theme(text = element_text(size = 12), axis.text.x = element_text(angle=60, hjust=1, size=12) )

Bactoplot6



################## Bacterial resistance per source

## read data 
newDF3 <- read.csv("newDF3.csv", sep = ",")



###### Ear culture resistance

by_ear <- select(newDF3, Bacteria, Source, Antibiotics, S, I, R)
by_ear <- group_by(by_ear, Source)
by_ear <- by_ear[grep("Ear culture", by_ear$Source), ]

by_ear$n <- by_ear$S + by_ear$R
by_ear$total <- by_ear$n + by_ear$I

by_ear$resis <- (by_ear$R/by_ear$n)*100
by_ear$resis <- format(by_ear$resis, digits=2, nsmall=2)
by_ear$resis <- as.numeric(by_ear$resis)
sum(by_ear$total)

# make pivot table EAR Bacterie ~ Antibiotics

by_ear_tabl = dcast(by_ear, Bacteria ~ Antibiotics, value.var = "resis")
View(by_ear_tabl)

summary(by_ear_tabl)

min(by_ear$resis)
max(by_ear$resis)


####### STOOL culture resistance

by_stool <- select(newDF3, Bacteria, Source, Antibiotics, S, I, R)
by_stool <- group_by(by_stool, Source)
by_stool <- by_stool[grep("Stool", by_stool$Source), ]

View(by_stool)

######## EYE culture resistance

by_eye <- select(newDF3, Bacteria, Source, Antibiotics, S, I, R)
by_eye <- group_by(by_eye, Source)
by_eye <- by_eye[grep("Eye culture", by_eye$Source), ]

by_eye$n <- by_eye$S + by_eye$R
by_eye$total <- by_eye$n + by_eye$I

by_eye$resis <- (by_eye$R/by_eye$n)*100
by_eye$resis <- format(by_eye$resis, digits=2, nsmall=2)
by_eye$resis <- as.numeric(by_eye$resis)
sum(by_eye$total)
View(by_eye)

# make pivot table EYE Bacterie ~ Antibiotics

by_eye_tabl = dcast(by_eye, Bacteria ~ Antibiotics, value.var = "resis")
View(by_eye_tabl)

summary(by_eye_tabl)

min(by_eye$resis)
max(by_eye$resis)

# total nr eye infective bacteria
by_eye_aggr <- aggregate(by_eye[,c("total")], by=list(by_eye$Bacteria), "sum")
View(by_eye_aggr)

####### WOUND culture resistance

by_wnd <- select(newDF3, Bacteria, Source, Antibiotics, S, I, R)
by_wnd <- group_by(by_wnd, Source)
by_wnd <- by_wnd[grep("Wound", by_wnd$Source), ]

by_wnd$n <- by_wnd$S + by_wnd$R
by_wnd$total <- by_wnd$n + by_wnd$I

by_wnd$resis <- (by_wnd$R/by_wnd$n)*100
by_wnd$resis <- format(by_wnd$resis, digits=2, nsmall=2)
by_wnd$resis <- as.numeric(by_wnd$resis)
sum(by_wnd$total)

# make pivot table WOUND Bacterie ~ Antibiotics

by_wnd_tabl = dcast(by_wnd, Bacteria ~ Antibiotics, value.var = "resis")
View(by_wnd_tabl)

summary(by_wnd_tabl)

min(by_wnd$resis)
max(by_wnd$resis)

ggplot(data = melt(by_wnd_tabl), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable))



# total nr WOUND infective bacteria
by_wnd_aggr <- aggregate(by_wnd[,c("total")], by=list(by_wnd$Bacteria), "sum")
View(by_wnd_aggr)



###### URINE culture resistance

##Select >6 obs ###

by_urine <- newDF3[complete.cases(newDF3), ]
by_urine$n <- by_urine$S + by_urine$R
by_urine <- by_urine[by_urine$n > 6,]

by_urine <- select(by_urine, Bacteria, Source, Antibiotics, S, I, R, n, Resis)
by_urine <- group_by(by_urine, Source)
by_urine <- by_urine[grep("Urine", by_urine$Source), ]

sum(by_urine$n)

# make pivot table URINE Bacterie ~ Antibiotics

by_urine_tabl = dcast(by_urine, Antibiotics ~ Bacteria, value.var = "Resis")
View(by_urine_tabl)


by_urine_tabl2 = dcast(by_urine, Antibiotics ~ Bacteria, value.var = "n")
View(by_urine_tabl2)


summary(by_urine_tabl)
summary(by_urine_tabl$`Enterococcus spp.`)


write.xlsx(by_urine_tabl2, "by_urine_tabl2.xlsx")

# Export to excel
write.xlsx(by_urine_tabl, "by_urine_tabl.xlsx")

