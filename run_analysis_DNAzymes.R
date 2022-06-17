## LOAD FILE
DZ <- read.table("HPV_FullLength.tab", header = TRUE, row.name = 1, sep = '\t')

# Subset Cairns sequences
DZ.Cairns <- DZ[grepl("DT", row.names(DZ)),]

DZ.Cairns$Energy <- DZ.Cairns$Energy * (-1)
DZ.Cairns$Dimer <- as.numeric(as.character(DZ.Cairns$Dimer)) * (-1)
DZ.Cairns$Internal <- as.numeric(as.character(DZ.Cairns$Internal)) * (-1)

## FIGURE 2
#
png("figure_2.png", width=16,height=16,units="cm",res=1200)
palette(rainbow(nrow(DZ.Cairns)))
plot(DZ.Cairns$Energy, DZ.Cairns$time60, col=1:nrow(DZ.Cairns), pch=19, xlab = "DNAzyme/RNA deltaG (kcal/mol)", ylab = "Efficiency (%)", xaxt = 'n', yaxt = 'n')
axis(1, at = c(10,15,20,25,30), labels = c("-10", "-15","-20", "-25", "-30"), tick = TRUE)
axis(2, las = 2, at = c(0,20,40,60,80), labels = c("0", "20","40", "60", "80"), tick = TRUE)
legend(25.5, 71, legend=rownames(DZ.Cairns), pch = 19, col=1:nrow(DZ.Cairns), cex=0.8)
mtext('A)', las = 2, side=2, line=1, at=93, cex = 1.3)
dev.off()


## FIGURE 3
#
png("figure_3.png", width=16,height=16,units="cm",res=1200)
par(mfrow = c(2,2))
DZ.Cairns$good20<-sapply(DZ.Cairns$time60, function(x) ifelse(x >= 20, 1, 0))
fit_logit_DZ.Cairns20 <- glm(formula = good20 ~ Energy, family = binomial(link = "logit"), data = DZ.Cairns)
# MODEL STATS, UNCOMMENT TO REPORT RESULTS IN CONSOLE
# summary(fit_logit_DZ.Cairns20)
# anova(fit_logit_DZ.Cairns20,test="Chisq")
newdat20 <- data.frame(Energy=seq(min(DZ.Cairns$Energy), max(DZ.Cairns$Energy),len=100))
plot(good20 ~ Energy, data=DZ.Cairns, col=1:nrow(DZ.Cairns), pch = 19, xaxt = 'n', xlab = "DNAzyme/RNA deltaG (kcal/mol)", yaxt = 'n', ylab = "Efficiency > 20% (0/1)")
axis(2, las = 2, at = c(0,1), labels = c("0", "1"), tick = TRUE)
axis(1, at = c(10,15,20,25,30), labels = c("-10", "-15","-20", "-25", "-30"), tick = TRUE)
#legend(25, 0.5, legend=rownames(DZ.Cairns), pch = 19, col=1:nrow(DZ.Cairns), cex=0.8)
newdat20$pred = predict(fit_logit_DZ.Cairns20, newdata=newdat20, type="response")
lines(pred ~ Energy, newdat20, lwd=2)
points(good20 ~ Energy, data=DZ.Cairns, col=1:nrow(DZ.Cairns), pch = 19, xaxt = 'n')
mtext('A)', las = 2, side=2, line=1, at=1.2, cex = 1.3)

DZ.Cairns$good40<-sapply(DZ.Cairns$time60, function(x) ifelse(x >= 40, 1, 0))
fit_logit_DZ.Cairns40 <- glm(formula = good40 ~ Energy, family = binomial(link = "logit"), data = DZ.Cairns)
# summary(fit_logit_DZ.Cairns40)
# anova(fit_logit_DZ.Cairns40,test="Chisq")
newdat40 <- data.frame(Energy=seq(min(DZ.Cairns$Energy), max(DZ.Cairns$Energy),len=100))
plot(good40 ~ Energy, data=DZ.Cairns, col=1:nrow(DZ.Cairns), pch = 19, xaxt = 'n', xlab = "DNAzyme/RNA deltaG (kcal/mol)", yaxt = 'n', ylab = "Efficiency > 40% (0/1)")
axis(1, at = c(10,15,20,25,30), labels = c("-10", "-15","-20", "-25", "-30"), tick = TRUE)
axis(2, las = 2, at = c(0,1), labels = c("0", "1"), tick = TRUE)
#legend(25, 0.5, legend=rownames(DZ.Cairns), pch = 19, col=1:nrow(DZ.Cairns), cex=0.8)
newdat40$pred = predict(fit_logit_DZ.Cairns40, newdata=newdat40, type="response")
lines(pred ~ Energy, newdat40, lwd=2)
mtext('B)', las = 2, side=2, line=1, at=1.2, cex = 1.3)

fit_logit_DZ.Cairns40.internal <- glm(formula = good40 ~ Internal, family = binomial(link = "logit"), data = DZ.Cairns)
# summary(fit_logit_DZ.Cairns40.internal)
# anova(fit_logit_DZ.Cairns40.internal,test="Chisq")
newdat40.internal <- data.frame(Internal=seq(min(DZ.Cairns$Internal), max(DZ.Cairns$Internal),len=100))
plot(good40 ~ Internal, data=DZ.Cairns, col=1:nrow(DZ.Cairns), pch = 19, xaxt = 'n', xlab = "DNAzyme-internal deltaG (kcal/mol)", yaxt = 'n', ylab = "Efficiency > 40% (0/1)")
axis(1, at = c(0,1,2,3), labels = c("0", "-1","-2", "-3"), tick = TRUE)
axis(2, las = 2, at = c(0,1), labels = c("0", "1"), tick = TRUE)
#legend(25, 0.5, legend=rownames(DZ.Cairns), pch = 19, col=1:nrow(DZ.Cairns), cex=0.8)
newdat40.internal$pred = predict(fit_logit_DZ.Cairns40.internal, newdata=newdat40.internal, type="response")
lines(pred ~ Internal, newdat40.internal, lwd=2)
mtext('C)', las = 2, side=2, line=1, at=1.2, cex = 1.3)

fit_logit_DZ.Cairns40.dimer <- glm(formula = good40 ~ Dimer, family = binomial(link = "logit"), data = DZ.Cairns)
# summary(fit_logit_DZ.Cairns40.dimer)
# anova(fit_logit_DZ.Cairns40.dimer,test="Chisq")
newdat40.dimer <- data.frame(Dimer=seq(min(DZ.Cairns$Dimer), max(DZ.Cairns$Dimer),len=100))
plot(good40 ~ Dimer, data=DZ.Cairns, col=1:nrow(DZ.Cairns), pch = 19, xaxt = 'n', xlab = "DNAzyme-homodimer deltaG (kcal/mol)", yaxt = 'n', ylab = "Efficiency > 40% (0/1)")
axis(1, at = c(5,6,7,8), labels = c("-5", "-6","-7", "-8"), tick = TRUE)
axis(2, las = 2, at = c(0,1), labels = c("0", "1"), tick = TRUE)
#legend(25, 0.5, legend=rownames(DZ.Cairns), pch = 19, col=1:nrow(DZ.Cairns), cex=0.8)
newdat40.dimer$pred = predict(fit_logit_DZ.Cairns40.dimer, newdata=newdat40.dimer, type="response")
lines(pred ~ Dimer, newdat40.dimer, lwd=2)
mtext('D)', las = 2, side=2, line=1, at=1.2, cex = 1.3)

dev.off()


# internal energies computed with RNA Vienna Package
# RNAup -b -d2 --noLP -P dna_mathews2004.par --noconv -c 'S' < sequences.fa > RNAup.out
# RNAfold -p -d2 --noLP -P dna_mathews2004.par --noconv < sequence1.fa > sequence1.out


## PREDICT EFFICIENCIES FOR NOVEL DNAzymes
data.to.fit <- DZ[16:30,]
data.to.fit$Energy <- data.to.fit$Energy * (-1)
data.to.fit$Dimer <- as.numeric(as.character(data.to.fit$Dimer)) * (-1)
data.to.fit$Internal <- as.numeric(as.character(data.to.fit$Internal)) * (-1)

# PREDICT Efficiency >= 20%
single.prob <- predict(fit_logit_DZ.Cairns20, newdata=data.to.fit, type="response")
single.pred <- ifelse(single.prob > 0.5, 1, 0)
single.pred

#ASSESS
data.to.fit$good<-sapply(data.to.fit$time60, function(x) ifelse(x >= 20, 1, 0))
data.to.fit$fittedS<-single.pred
# Accuracy (TP+TN)/ALL
misClasificErrorS <- mean(data.to.fit$fittedS != data.to.fit$good)
TP.S <- sum(data.to.fit$fittedS==1 & data.to.fit$good==1)
FP.S <- sum(data.to.fit$fittedS==1 & data.to.fit$good==0)
FN.S <- sum(data.to.fit$fittedS==0 & data.to.fit$good==1)
# precision TP/(TP+FP) AND recall	TP/(TP+FN)
print(paste('Single Accuracy',1-misClasificErrorS))
print(paste('Single Precision',TP.S/(TP.S+FP.S)))
print(paste('Single Recall',TP.S/(TP.S+FN.S)))


## Multiple logistic regression
# PREDICT Efficiency >= 40%

fit_logit_DZ.Cairns_multiple <- glm(formula = good40 ~ Energy + Dimer + Internal, family = binomial(link = "logit"), data = DZ.Cairns)
summary(fit_logit_DZ.Cairns_multiple)
anova(fit_logit_DZ.Cairns_multiple,test="Chisq")
# Simple model
single.prob <- predict(fit_logit_DZ.Cairns40, newdata=data.to.fit, type="response")
single.pred <- ifelse(single.prob > 0.5, 1, 0)
single.pred
multiple.prob <- predict(fit_logit_DZ.Cairns_multiple, newdata=data.to.fit, type="response")
multiple.pred <- ifelse(multiple.prob > 0.5, 1, 0)
multiple.pred
#ASSESS
data.to.fit$good<-sapply(data.to.fit$time60, function(x) ifelse(x >= 40, 1, 0))
data.to.fit$fittedM<-multiple.pred
data.to.fit$fittedS<-single.pred
# Accuracy (TP+TN)/ALL
misClasificErrorS <- mean(data.to.fit$fittedS != data.to.fit$good)
# 0.625
misClasificErrorM <- mean(data.to.fit$fittedM != data.to.fit$good)
# 0.625
TP.M <- sum(data.to.fit$fittedM==1 & data.to.fit$good==1)
FP.M <- sum(data.to.fit$fittedM==1 & data.to.fit$good==0)
FN.M <- sum(data.to.fit$fittedM==0 & data.to.fit$good==1)
TP.S <- sum(data.to.fit$fittedS==1 & data.to.fit$good==1)
FP.S <- sum(data.to.fit$fittedS==1 & data.to.fit$good==0)
FN.S <- sum(data.to.fit$fittedS==0 & data.to.fit$good==1)
# precision TP/(TP+FP) AND recall	TP/(TP+FN)
print(paste('Single Accuracy',1-misClasificErrorS))
print(paste('Single Precision',TP.S/(TP.S+FP.S)))
print(paste('Single Recall',TP.S/(TP.S+FN.S)))
print(paste('Multiple Accuracy',1-misClasificErrorM))
print(paste('Multiple Precision',TP.M/(TP.M+FP.M)))
print(paste('Multiple Recall',TP.M/(TP.M+FN.M)))



## PREDICT over the second set of novel DNAzymes

DZ.expanded <- DZ[1:30,]
DZ.expanded$good<-sapply(DZ.expanded$time60, function(x) ifelse(x >= 40, 1, 0))

fit_logit_DZ.expanded_multiple <- glm(formula = good ~ Energy + Dimer + Internal, family = binomial(link = "logit"), data = DZ.expanded)
summary(fit_logit_DZ.expanded_multiple)
anova(fit_logit_DZ.expanded_multiple,test="Chisq")

data.to.fit.2 <- DZ[31:38,]
data.to.fit.2$Energy <- data.to.fit.2$Energy * (-1)
data.to.fit.2$Dimer <- as.numeric(as.character(data.to.fit.2$Dimer)) * (-1)
data.to.fit.2$Internal <- as.numeric(as.character(data.to.fit.2$Internal)) * (-1)

multiple.prob <- predict(fit_logit_DZ.Cairns_multiple, newdata=data.to.fit.2, type="response")
multiple.pred <- ifelse(multiple.prob > 0.5, 1, 0)
multiple.pred
#ASSESS
data.to.fit.2$good<-sapply(data.to.fit.2$time60, function(x) ifelse(x >= 40, 1, 0))
data.to.fit.2$fittedM<-multiple.pred
data.to.fit.2$fittedS<-single.pred
# Accuracy (TP+TN)/ALL
misClasificErrorS <- mean(data.to.fit.2$fittedS != data.to.fit.2$good)
misClasificErrorM <- mean(data.to.fit.2$fittedM != data.to.fit.2$good)
TP.M <- sum(data.to.fit.2$fittedM==1 & data.to.fit.2$good==1)
FP.M <- sum(data.to.fit.2$fittedM==1 & data.to.fit.2$good==0)
FN.M <- sum(data.to.fit.2$fittedM==0 & data.to.fit.2$good==1)
TP.S <- sum(data.to.fit.2$fittedS==1 & data.to.fit.2$good==1)
FP.S <- sum(data.to.fit.2$fittedS==1 & data.to.fit.2$good==0)
FN.S <- sum(data.to.fit.2$fittedS==0 & data.to.fit.2$good==1)
print(paste('Multiple Accuracy',1-misClasificErrorM))
print(paste('Multiple Precision',TP.M/(TP.M+FP.M)))
print(paste('Multiple Recall',TP.M/(TP.M+FN.M)))


## MODEL WITH ALL DNAzymes IN THE PAPER

DZ.ALL <- DZ[1:38,]
DZ.ALL$good<-sapply(DZ.ALL$time60, function(x) ifelse(x >= 40, 1, 0))

fit_logit_DZ.all_multiple <- glm(formula = good ~ Energy + Dimer + Internal, family = binomial(link = "logit"), data = DZ.ALL)
summary(fit_logit_DZ.all_multiple)
anova(fit_logit_DZ.all_multiple,test="Chisq")


quit()
