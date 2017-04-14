setwd("~/Documents/workspace/R-project/Antibiogram")


library(xlsx)
library(reshape2)
library(dplyr)
library(scales)
library(ggplot2)





########################### Bacterial resistance per source ##############################

## read data 
newDF3 <- read.csv("newDF3.csv", sep = ",")

write.xlsx(newDF3,"newDF3.xlsx")
########################## SUMMARY per Source total ################

### Aggregate on Source.

DFsource <- aggregate(cbind(S,I,R,N,n)~Source, data=newDF3, sum, na.rm=TRUE)

sum(DFsource$S)
sum(DFsource$I)
sum(DFsource$R)
sum(DFsource$N)
sum(DFsource$n)

View(DFsource)
write.xlsx(DFsource, "DFsource.xlsx")


########################### SELECT cases n >6 #########################

Bactresis7 <- newDF3[newDF3$n > 6,]
sum(Bactresis7$S)
sum(Bactresis7$I)
sum(Bactresis7$R)
sum(Bactresis7$N)
sum(Bactresis7$n)

sum(newDF3$N) - sum(Bactresis7$N)
write.csv(Bactresis7, "Bactresis7.csv")

########################### SELECT Glycopeptides #########################

vanco <- Bactresis7[grep("VA30", Bactresis7$Antibiotics), ]




########################### SELECT beta-lactams AND Staph aureus #########################


toMatch1 <- Bactresis7[grep("Staphylococcus aureus",Bactresis7$Bacteria), ]


toMatch2 <-subset(toMatch1, Antibiotics == "AMC30" | Antibiotics == "AMP10" | Antibiotics == "CTX30" 
                  | Antibiotics == "FOX30" | Antibiotics == "CRO30" | Antibiotics == "CXM30" 
                  | Antibiotics == "CF30"  | Antibiotics == "OX1"   | Antibiotics == "P10" )

toMatch2 <- toMatch1[grep("AMC30", "AMP10","OX1","P10",Bactresis7$Antibiotics), ]


toMatch <- c("A1", "A9", "A6")
matches <- unique (grep(paste(toMatch,collapse="|"), 
                        myfile$Letter, value=TRUE))


###################### URINE culture resistance #############################


by_urine <- Bactresis7[complete.cases(Bactresis7), ]
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



###################### WOUND culture resistance #############################

by_wnd <- select(Bactresis7, Bacteria, Source, Antibiotics, N, n, Resis)
by_wnd <- group_by(by_wnd, Source)
by_wnd <- by_wnd[grep("Wound", by_wnd$Source), ]

sum(by_wnd$N)
sum(by_wnd$n)

# make pivot table WOUND Bacterie ~ Antibiotics

by_wnd_tabl = dcast(by_wnd, Antibiotics ~ Bacteria , value.var = "Resis")
View(by_wnd_tabl)


by_wnd_tabl2 = dcast(by_wnd, Antibiotics ~ Bacteria , value.var = "n")
View(by_wnd_tabl2)



summary(by_wnd_tabl)

min(by_wnd$Resis)
max(by_wnd$Resis)

sum(by_wnd_tabl2$CNS)
summary(by_wnd_tabl2$CNS)

CNS_r <- by_wnd_tabl$CNS[complete.cases(by_wnd_tabl$CNS)]
summary(CNS_r)

# ggplot(data = melt(by_wnd_tabl), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable))



# total nr WOUND infective bacteria
by_wnd_aggr <- aggregate(by_wnd[,c("N")], by=list(by_wnd$Bacteria), "sum")
View(by_wnd_aggr)



###################### EAR culture resistance #############################

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


###################### STOOL culture resistance #############################

by_stool <- select(newDF3, Bacteria, Source, Antibiotics, S, I, R)
by_stool <- group_by(by_stool, Source)
by_stool <- by_stool[grep("Stool", by_stool$Source), ]

View(by_stool)

###################### EYE culture resistance #############################

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







