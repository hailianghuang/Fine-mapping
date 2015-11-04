source("./code/math.r")

MIN_LD_REPOS <- 0.1

#start of main program
dosage <- read.table("dat/example.raw", header=T, stringsAsFactors=F)
dosage <- dosage[,-c(1:6)]
pheno <- scan("dat/example.pheno")
PCA <- read.table("dat/pcs.txt", header=F)
Neff <- getT(cov(dosage))
result <- modelSearch(dosage, pheno, PCA, Neff)

#start repositioning
traits <- c()
hits1 <- c()
nParam <- c()
hit_pval <- c()
prev_max <- result[1,]

for (i in c(1:max(result$model_size))){
  result_sel <- result[result$model_size==i, ]
  imax <- which.max(result_sel$rel_logProb)
  max <- result_sel[imax,]
  if( (i>1 && max$logProb-prev_max$logProb <= 0 )  ){
    break
  }
  traits <- c(traits, as.character(max$add_trait))
  hits1 <- c(hits1,  as.character(max$add_SNP))
  hit_pval <- c(hit_pval, max$pval)
  nParam <- c(nParam, max$n_param_add)
  prev_max <- max
}

repos <- c()
if(length(hits1)==0){
  
  print("model size of 0, quit")
  quit(save="no")
  
}else if(length(hits1)==1){
  
  print("model size of 1, printing result")
  repos <- result[result$model_size==1, -c(1,2,3, 20, 24, 26, 27)]

}else{
  # set model
  constrain <- rep(0, 1 + length(hits1) + dim(PCA)[2])
  constrain[1 + which(traits == "UC")] <- 1
  constrain[1 + which(traits == "CD")] <- 2

  iterhits <- matrix(NA,nrow=30,ncol=length(hits1))
  iterhits[1,] <- hits1 
  
  iter <- 2
  while(TRUE){
    
    # extract genotypes for the current best model
    xhits <- dosage[, iterhits[iter - 1,]]
    xs <- cbind(1,xhits,PCA)
      
    # find the likelihood for the current best model
    bestL <- -multinomConstrain(pheno,xs,constrain)$value
    pval_list <- c()
      
    # for each of the signals
    for (hiti in 1:length(hits1)){
      if(constrain[hiti+1]==3){
        pval_list <- c(pval_list, 1);
        iterhits[iter,hiti] <- NA
        next;
      }else{
        constrain_null <- constrain
        constrain_null[hiti + 1] <- 3
        NULL_L <- -multinomConstrain(pheno,xs, constrain_null)$value
          
          # get data for the current best SNP
        xhit <- xs[,hiti + 1]
          
        x_other_to_keep <- rep(T,length(hits1))
        x_other_to_keep[hiti] <- F
        x_other <- xhits[, x_other_to_keep ]
          
        print(paste("working on ", iter-1, "-",  hiti))
          
        r2  <- sapply( names(dosage), function(x){ r2 <- cor(dosage[, x],xhit)^2 } ) 
        rr <- lapply(names(subset(dosage,  select=r2>MIN_LD_REPOS) ), reposTest, dosage, x_other, xs, hiti, pheno, constrain)
        ret <- do.call("rbind", rr)
          
        imax <- which.max(ret$newL)
        bestsnp <- as.character(ret$snp[imax])
        pval_list <- c(pval_list, pchisq((ret$newL[imax]-NULL_L)*2, nParam[hiti], low=F))
          
        # record the best SNP
        iterhits[iter,hiti] <- bestsnp
        xs[,hiti + 1] <- dosage[, bestsnp]
      }
    }
      
     # if nothing changed this iterations, break the loop
    if ( all(iterhits[iter,] == iterhits[iter - 1,] ) | iter > 9) {
      break
    }else{
      iter <- iter + 1
    }
    
  }
  print("done re-positioning")
  
  nhit <- length(hits1)
  null_logProb <- result[1,]$logProb
  cnt <- 1
  
  for (hiti in c(1:nhit)){
    
    print(paste("hit #", hiti))
    x_in_current <- rep(T, length(hits1))
    x_in_current[hiti] <- F
    current_model <- modelProb(pheno, PCA, buildModel ( iterhits[iter,x_in_current], traits[x_in_current]), dosage, Neff)		
    rr <- lapply(names(dosage), testSNP, pheno, PCA, Neff, current_model, null_logProb, dosage)
    ret <- do.call("rbind", rr)
    repos <- rbind(repos,  data.frame(ret[, c(4, 5)], model_size=cnt, ret[,c(7:10, 12:14, 17)]))
    cnt <- cnt +1
  }
  
}

write.table(repos, file="result/repos.txt", col.names=TRUE, row.names=F, quote=F, sep="\t")

traits1 <- c()
hits1 <- c()
for (i in c(1:max(repos$model_size))){
  result_sel <- repos[repos$model_size==i, ]
  max <- result_sel[which.max(result_sel$rel_logProb),]
  traits1 <- c(traits1, as.character(max$add_trait))
  hits1 <- c(hits1,  as.character(max$add_SNP))
}

output <- data.frame()
for (model_size in 1:length(hits1)){
  
  best_trait <- traits1[model_size]
  best_snp <- hits1[model_size]
  result_sel <- repos[repos$model_size==model_size, ]
  
  log10_pval_IBD <- -log10(result_sel$pval_IBD[result_sel$add_trait=="IBD"])
  rel_logProb_IBD <- result_sel$rel_logProb[result_sel$add_trait=="IBD"] - apply(cbind(result_sel$rel_logProb[result_sel$add_trait=="CD"] , result_sel$rel_logProb[result_sel$add_trait=="UC"]), 1, max)
  result_sel <- result_sel[result_sel$add_trait==best_trait, ]
  result_best <- result_sel[match(best_snp,result_sel$add_SNP ), ]
  
  if(dim(result_sel)[1] == 0){
    break;
  }
  
  ret <- getCredibleSNP(best_snp, result_sel$rel_logProb, result_sel$add_trait, as.character(result_sel$add_SNP) )
  
  temp <- data.frame(SNP=ret$credible_set_ichip, 
                     signal=rep(model_size, ret$nSNP), 
                     trait=as.character(ret$credible_trait), 
                     logOR_CD=sprintf("%.3f",  result_sel$logOR_CD[match(ret$credible_set_ichip,  result_sel$add_SNP)]), 
                     logOR_UC=sprintf("%.3f", result_sel$logOR_UC[match(ret$credible_set_ichip,  result_sel$add_SNP)]), 
                     deviance = sprintf("%.1f", result_sel$deviance[match(ret$credible_set_ichip,  result_sel$add_SNP)]), 
                     nParam=result_sel$n_param[match(ret$credible_set_ichip,  result_sel$add_SNP)], 
                     log10Pval=sprintf("%.2f", (pchisq(result_sel$deviance, result_sel$n_param_add, low=F,log.p=T)/log(10))[match(ret$credible_set_ichip,  result_sel$add_SNP)] ), 
                     logProb=sprintf("%.2f",result_sel$rel_logProb[match(ret$credible_set_ichip,  result_sel$add_SNP)] ), 
                     probNormed=sprintf("%.3f",ret$prob_normed[match(ret$credible_set_ichip,  result_sel$add_SNP)] ), 
                     R=sprintf("%.3f",ret$R[match(ret$credible_set_ichip,  result_sel$add_SNP)]),
                     log10_pval_IBD=sprintf("%.3f", log10_pval_IBD[match(ret$credible_set_ichip,  result_sel$add_SNP)] )
  )
  
  output <- rbind(output, temp)
}
write.table(output, file="result/credible.txt", col.names=TRUE, row.names=F, quote=F, sep="\t")

