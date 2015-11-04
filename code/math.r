require(nnet)

getT <- function(cov){
  w <- rep(1, dim(cov)[1])
  T <-0
  while(max(w)>0){
    T <- T+max(w)
    v1 <- cov[which.max(w), which.max(w)]
    v2 <- diag(cov)
    cor2 <- 0
    v2[v2==0]=0.001
    if(v1>0){
      cor2 <- cov[,which.max(w)]^2/v2/v1
    }
    cor2[v2==0]=0
    w <- w - cor2*max(w)
  }
  T
}

multinomConstrain <- function(y,xs,constrain){
 # fit a constrained multinomial model
 Npheno <- length(unique(y)) - 1
 Npred <- length(xs[1,])
 X <- xs
 Y <- t(sapply(y,function(l) l == (0:Npheno)))

 mask <- rep(FALSE,(Npheno + 1)*(Npred + 1))
 for (i in 1:Npheno){
	mask[(Npred + 1)*i + 2:(2+Npred - 1)] <- TRUE
	}

	for (i in 1:Npred){
	if (constrain[i] == 1) mask[(Npred + 1)*2 + 1 + i] <- FALSE
	if (constrain[i] == 2) mask[(Npred + 1)*1 + 1 + i] <- FALSE
		if (constrain[i] == 3) {
			mask[(Npred + 1)*1 + 1 + i] <- FALSE
			mask[(Npred + 1)*2 + 1 + i] <- FALSE
		}
	}

	nnet.default(X, Y, mask = mask, size = 0,skip = TRUE, softmax = TRUE, rang=0, trace=FALSE)
}

getlogLik <- function(y, covar, model_spec, dat ){

	ret <- NA
	ncovar <- dim(covar)[2]
	
	# set model
	var_total <- length(model_spec$IBD_SNP) + length(model_spec$CD_SNP) + length(model_spec$UC_SNP)

	constrain <- rep(0,var_total + 1 + ncovar)
	if(length(model_spec$CD_SNP) > 0){
		constrain[1 + length(model_spec$IBD_SNP) + 1:length(model_spec$CD_SNP)] <- 2
	}
	if(length(model_spec$UC_SNP) > 0){
		constrain[1 + length(model_spec$IBD_SNP) + length(model_spec$CD_SNP)  + 1:length(model_spec$UC_SNP)] <- 1
	}

	model_element <- dat[,  c(model_spec$IBD_SNP, model_spec$CD_SNP, model_spec$UC_SNP) ]

	xs <- cbind(1,model_element,covar)

	temp <- multinomConstrain(y,xs,constrain)
	bestL <- -temp$value
	param <- temp$wts

	N_uc <- sum(y==1 | y==0, na.rm=T) 
	N_cd <- sum(y==2 | y==0, na.rm=T)
	N_ibd <- sum(!is.na(y))

	P_cd <- length(model_spec$CD_SNP)
	P_uc <- length(model_spec$UC_SNP)
	P_ibd <- length(model_spec$IBD_SNP)
	P = P_cd + P_uc +P_ibd


	if(P_ibd+P_uc+P_cd==0){
		param_CD <- NULL
		param_UC <- NULL
	}else{
		param_CD <- param[((1+1+P+ncovar)*2+1+1+1): ((1+1+P+ncovar)*2+1+1+1 + P-1) ]
		param_UC <- param[(1+1+P+ncovar+1+1+1): (1+1+P+ncovar+1+1+1 + P-1) ]
	}	

	ret <- list(N_cd = N_cd, 
							N_uc = N_uc,
							N_ibd = N_ibd,
							P_cd = P_cd, 
							P_uc = P_uc,
							P_ibd = P_ibd, 
							CD_SNP = model_spec$CD_SNP, 
							UC_SNP = model_spec$UC_SNP,
							IBD_SNP = model_spec$IBD_SNP, 
							param_CD = param_CD,
							param_UC = param_UC,
							logLik = bestL,
							nnet=temp
							)
	ret	
}

growModel <- function(model_spec, i, trait, dat){

	existing_snp <- c(model_spec$CD_SNP, model_spec$UC_SNP, model_spec$IBD_SNP)

	existing_snp_geno <- dat[, existing_snp ] 
	  
	new_pos <- 1;
	if(i %in% existing_snp){
		model_spec$good <- 0;
	}else if(length(existing_snp) >0 && summary(lm(dat[,i]~ as.matrix(existing_snp_geno) ))$r.squared >0.8){
		model_spec$good <- 0;
	}else{ 
		if(trait=="CD"){
			model_spec$CD_SNP <- c(model_spec$CD_SNP, i)
			new_pos <- length(model_spec$IBD_SNP) + length(model_spec$CD_SNP)
		}else if(trait=="UC"){
			model_spec$UC_SNP <- c(model_spec$UC_SNP, i)
			new_pos <- length(model_spec$IBD_SNP) + length(model_spec$CD_SNP)+ length(model_spec$UC_SNP)
		}else if(trait=="IBD"){
			model_spec$IBD_SNP <- c(model_spec$IBD_SNP, i)
			new_pos <- length(model_spec$IBD_SNP)
		}
		model_spec$good <- 1;
		model_spec$new_pos <- new_pos;
	}
	model_spec
}

buildModel <- function( snp, trait){

	model_spec <- c()
	model_spec$CD_SNP <- as.character(snp[trait=="CD"])
	model_spec$UC_SNP <- as.character(snp[trait=="UC"])
	model_spec$IBD_SNP <- as.character(snp[trait=="IBD"])
	model_spec$good  <- 1
	model_spec$new_pos  <- -1

	model_spec

}

getSpec <- function(test_ret){

	model_spec <- list(CD_SNP=NULL, UC_SNP=NULL, IBD_SNP=NULL, good=1)

	snps <- unlist(strsplit(as.character(test_ret$CD_SNP), '@'))
	if(length(snps)>0){
		model_spec$CD_SNP <- snps
	}

	snps <- unlist(strsplit(as.character(test_ret$UC_SNP), '@'))
	if(length(snps)>0){
		model_spec$UC_SNP <- snps
	}

	snps <- unlist(strsplit(as.character(test_ret$IBD_SNP), '@'))
	if(length(snps)>0){
		model_spec$IBD_SNP <- snps
	}

	model_spec
}

modelProb<-function(y, covar, model_spec, dat, Neff){
  
  ret <- getlogLik(y, covar, model_spec, dat)
  
  if(sum(!is.na(ret)) > 0 ){
    BIC_penalty <- - log(ret$N_cd)/2*(ret$P_ibd+ret$P_cd) - log(ret$N_uc)/2*(ret$P_ibd+ret$P_uc)
    model_penalty <- -lchoose(Neff, length(ret$IBD_SNP)) - lchoose(Neff-length(ret$IBD_SNP), length(ret$CD_SNP)+length(ret$UC_SNP))
    
    logProb <- ret$logLik + BIC_penalty + model_penalty
    ret <- c(ret, BIC_penalty=BIC_penalty , model_penalty=model_penalty, logProb =logProb )
  }
  
  ret
}

testSNP <- function(i, y, covar, Neff, prev_best_model, null_logProb, dat){
  
  ret <- c();
  bestDeviance <- -1
  bestlogProb <- -Inf
  prev_best_model_spec <- getSpec(prev_best_model)
  
  for (d in c("CD", "UC", "IBD")){
    
    current_model_spec <- growModel(prev_best_model_spec, i, d, dat)
    if(current_model_spec$good==0){
      next;
    }
    current_model <- modelProb(y, covar, current_model_spec, dat, Neff)
    
    deviance <- 2*(current_model$logLik-prev_best_model$logLik)
    n_param_add <- 1;
    pval_IBD <- 1;
    logBF_IBD <- 0;
    
    if(d=="IBD"){
      n_param_add <- 2;
      pval_IBD <- pchisq(deviance-bestDeviance, 1, low=F)
      logBF_IBD <- current_model$logProb-bestlogProb
    }else{
      if(bestDeviance < deviance){
        bestDeviance <- deviance
      }
      if(bestlogProb < current_model$logProb){
        bestlogProb <- current_model$logProb
      }
    }
    pval <- pchisq(deviance, n_param_add, low=F)
    
    temp <- data.frame(CD_SNP=paste(current_model$CD_SNP, collapse="@"), 
                       UC_SNP=paste(current_model$UC_SNP, collapse="@"), 
                       IBD_SNP=paste(current_model$IBD_SNP, collapse="@"),
                       add_SNP=as.character(i), 
                       add_trait=d,
                       model_size=length(current_model$CD_SNP) + length(current_model$UC_SNP) + length(current_model$IBD_SNP),
                       n_param=length(current_model$CD_SNP) + length(current_model$UC_SNP) + 2*length(current_model$IBD_SNP),
                       n_param_add=n_param_add,
                       logOR_CD=current_model$param_CD[current_model_spec$new_pos],
                       logOR_UC=current_model$param_UC[current_model_spec$new_pos],
                       logLik=current_model$logLik,
                       deviance=deviance,
                       pval=pval,
                       pval_IBD=pval_IBD,
                       logBF_IBD=logBF_IBD,
                       logProb=current_model$logProb, 
                       rel_logProb=current_model$logProb- null_logProb,
                       BIC_penalty=current_model$BIC_penalty,
                       model_penalty=current_model$model_penalty
    )			
    
    ret <- rbind(ret, temp)
  }
  ret
}

modelSearch <- function(dosage, pheno, PCA, Neff){	
  
  iter <- 0
  
  null_model <- modelProb(pheno, PCA, list(CD_SNP=NULL, UC_SNP=NULL, IBD_SNP=NULL, good=1), dosage,  Neff)
  
  result <- data.frame(CD_SNP="", 
                       UC_SNP="", 
                       IBD_SNP="",
                       add_SNP="", 
                       add_trait="-",
                       model_size=0,
                       n_param=0,
                       n_param_add=0,
                       logOR_CD=0,
                       logOR_UC=0,
                       logLik=null_model$logLik,
                       deviance=0,
                       pval=1,
                       logBF_IBD=0,
                       pval_IBD=1,
                       logProb=null_model$logProb, 
                       rel_logProb=0,
                       BIC_penalty=0,
                       model_penalty=0
  )
  
  while(TRUE){
    
    if(iter==0){
      prev_best_model <- result
      print(paste(sep="", "Start searching model of size ", iter+1, ". logProb of the null model is ", sprintf("%.2f", null_model$logProb)  ));
    }else{
      
      print(paste("Completed searching for model of size", iter))
      print(paste("Previous relative logProb is", sprintf("%.2f", prev_best_model$logProb - null_model$logProb )))
      print(paste("Prev best model is CD:", prev_best_model$CD_SNP, "UC:", prev_best_model$UC_SNP,  "IBD:", prev_best_model$IBD_SNP) )		
      print(paste("Current best model is CD:", best_model$CD_SNP, "UC:", best_model$UC_SNP,  "IBD:", best_model$IBD_SNP) )		
      print(paste("Current relative logProb is", sprintf("%.2f", (best_model$logProb- null_model$logProb)) ))
      print(paste("Start searching model of size", iter+1))
      prev_best_model <- best_model
    }
    
    rr <- lapply(names(dosage), testSNP, pheno, PCA, Neff, prev_best_model, null_model$logProb, dosage)
    ret <- do.call("rbind", rr)
    result <- rbind(result, ret)
    
    imax <- which.max(ret$logProb)
    best_model <- ret[imax,]
    
    if( best_model$logProb - prev_best_model$logProb <= 0  ){
      print("Model search stopped")
      print(paste("Current relative logProb is", sprintf("%.2f",  best_model$logProb - null_model$logProb) ))
      print(paste("Current best model is CD:", best_model$CD_SNP, "UC:", best_model$UC_SNP,  "IBD:", best_model$IBD_SNP) )		
      break; 
    }
    
    iter <- iter+1
  }
  
  result
}

reposTest <- function(snp, dosage,  x_other, xs, hiti, y, constrain){
  
  xnew <- dosage[, snp ]
  
  xs_new <- xs
  xs_new[,(hiti + 1)] <- xnew
  
  if( length(x_other)> 0 && summary(lm(xnew~ x_other))$r.squared > 0.8){
    newL <- NA
  }else{
    newL <- -multinomConstrain(y, xs_new, constrain)$value
  }
  
  data.frame(snp, newL)
}

getCredibleSNP <- function(best_snp, logProb,model, snp_name, threshold=0.99){
  
  R <- sapply(snp_name, function(x){y <- dosage[, x]; cor(y, dosage[,best_snp])})
  select <- R^2 <0.3
  prob <- exp(logProb-max(logProb))
  prob[select] <- 0
  prob_normed <-  prob/sum(prob)
  
  prob_cumsum <- cumsum(sort(prob_normed,  decreasing=TRUE))
  nSNP <- which.max(prob_cumsum>threshold )
  
  if(nSNP<=0){
    nSNP=1
  }
  
  credible_set <- snp_name[order(prob_normed, decreasing=TRUE)[1:nSNP]]
  credible_trait <- model[order(prob_normed, decreasing=TRUE)[1:nSNP]]
  credible_set_ichip <- as.character(credible_set)
  credible_set_ii <- snp_name %in% credible_set_ichip
  
  ret <- list(nSNP = nSNP,
              prob_normed = prob_normed,
              credible_set_ichip=credible_set_ichip,
              credible_set_ii=credible_set_ii,
              credible_trait=credible_trait,
              R = R
  )
  
  ret
}
