require(Hmisc)
N <- 112
OR_CD <- c(4, 5)
OR_UC <- c(1, 0.6)

dosage <- read.table("dat/example.raw", header=T)
dosage <- dosage[,-c(1:6)]

set.seed(1)
sum(temp<0.2&temp>0.15)
signal <- sample((1:(dim(dosage)[2]))[temp<0.2&temp>0.15], 2)
risk_CD <-  as.matrix(dosage[, signal])  %*% log(OR_CD) 
risk_UC <-  as.matrix(dosage[, signal])  %*% log(OR_UC) 
risk_total <- 1+exp(risk_CD) + exp(risk_UC)
prob <- cbind(1/risk_total, exp(risk_UC)/risk_total, exp(risk_CD)/risk_total) 
pheno <- rMultinom(prob, 1)-1

write.table(pheno,"dat/example.pheno", col.names=F, row.names=F, quote=F)

PCA <- matrix(rnorm(2*N, sd=0.01), nrow=N)
write.table(PCA, "dat/pcs.txt", col.names=F, row.names=F, quote=F)

names(dosage)[signal]
apply(dosage[, signal], 2, mean)/2
cor(dosage[,signal])
summary(glm(I(pheno[pheno!=1]>0)~as.matrix(dosage[pheno!=1, signal]), family=binomial))
apply(dosage, 2, mean)/2 -> temp


