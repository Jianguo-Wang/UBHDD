
# 02/13/22
# R script in:
# The trait coding rule in phenotype space
# 
# Jianguo Wang (wangjg36@mail.sysu.edu.cn)
# 


######################### library ###########################################

setwd("/home/wjg/UBHDD/yeastSegregants")

library(glmnet)
library(parallel)


######################### functions ###########################################



## notice: input x should be without missing values; if exist, function will make an error
## notice: input y is a vector; if y is a matrix or data.frame, function will make an error; linear_learn_matrix can deal with matrix-wide learning

############### linear.learn ###############
linear_learn<-function(x_ll,y_ll,alpha_ll=1,nlambda_ll=200,testSize_ll=100){

    performance_ll<-matrix(NA,nrow=1,ncol=3)
    colnames(performance_ll)<-c("Rsquare","alpha","lambda")

    testIndex_ll<-c()
	
    if(class(x_ll)=="numeric"|class(x_ll)=="integer"){
	    
		naIndex_ll<-which(is.na(y_ll))
        if(length(naIndex_ll)!=0){
            x_ll<-x_ll[-naIndex_ll]
            y_ll<-y_ll[-naIndex_ll]
        }

	    testIndex_ll_temp<-sample(1:length(y_ll),testSize_ll)
		trainIndex_ll_temp<-setdiff(1:length(y_ll),testIndex_ll_temp)
		x_train_ll<-x_ll[trainIndex_ll_temp]
		y_train_ll<-y_ll[trainIndex_ll_temp]
		
		x_test_ll<-x_ll[testIndex_ll_temp]
		y_test_ll<-y_ll[testIndex_ll_temp]

		
        lm_ll<-lm(y_train_ll~x_train_ll)	
		
		performance_ll[1,1]<-cor(y_test_ll,x_test_ll*lm_ll$coefficients[2]+lm_ll$coefficients[1])^2
		performance_ll[1,2]<-NA
		performance_ll[1,3]<-NA
		
		coefficients_ll<-matrix(lm_ll$coefficients,nrow=2,ncol=1)
		
		if(length(naIndex_ll)!=0){
            testIndex_ll<-c(1:(length(y_ll)+length(naIndex_ll)))[-naIndex_ll][testIndex_ll_temp]
        }else{
            testIndex_ll<-testIndex_ll_temp	    
	    }		
	}
	else{
	
        naIndex_ll<-which(is.na(y_ll))
		## ****
        if(length(naIndex_ll)!=0){
            x_ll<-x_ll[-naIndex_ll,]
            y_ll<-y_ll[-naIndex_ll]
        }

        testIndex_ll_temp<-sample(1:length(y_ll),testSize_ll)
        trainIndex_ll_temp<-setdiff(1:length(y_ll),testIndex_ll_temp)
		
        x_train_ll<-x_ll[trainIndex_ll_temp,]
        y_train_ll<-y_ll[trainIndex_ll_temp]
		
        x_test_ll<-x_ll[testIndex_ll_temp,]
        y_test_ll<-y_ll[testIndex_ll_temp]

        cvglm_ll<-cv.glmnet(x_train_ll, y_train_ll, family="gaussian", alpha = alpha_ll, nlambda =nlambda_ll,
                            nfolds=10,thresh = 1e-07,  type.measure="deviance")
							
        performance_ll[1,1]<-cor(x_test_ll%*%coef(cvglm_ll,s="lambda.min")[-1,1],y_test_ll)^2
        performance_ll[1,2]<-alpha_ll
		performance_ll[1,3]<-cvglm_ll$lambda.min
		
		coefficients_ll<-coef(cvglm_ll,s="lambda.min")
		
        if(length(naIndex_ll)!=0){
            testIndex_ll<-c(1:(length(y_ll)+length(naIndex_ll)))[-naIndex_ll][testIndex_ll_temp]
        }else{
            testIndex_ll<-testIndex_ll_temp	    
	    }
    }
	

	
    return(list(performance_ll,coefficients_ll,testIndex_ll))
}


#now this program is designed for corresponding relationship between the traits of X and Y 

linear_learn_down_sampling<-function(X_llds,Y_llds,coincide_llds=F,unCorCutoff_llds=1.5,paraComp_llds=T,cores_llds=2,cycle_llds=1,alpha_llds=1,nlambda_llds=200,testSize_llds=100,downRatio_llds=c(1:10)*0.1){

    perfList_llds<-list()
    coefList_llds<-list()
    testIndexList_llds<-list()    
	
	for(i_llds in 1:length(downRatio_llds)){
	    set.seed(1)
	    focalSample_llds<-sample(1:dim(X_llds)[2],round(downRatio_llds[i_llds]*dim(X_llds)[2]))
        cvglm_llds<-linear_learn_matrix(X_llm=X_llds[,focalSample_llds],Y_llm=Y_llds[,focalSample_llds],coincide_llm=coincide_llds,unCorCutoff_llm=unCorCutoff_llds,paraComp_llm=paraComp_llds,cores_llm=cores_llds,cycle_llm=cycle_llds,alpha_llm=alpha_llds,nlambda_llm=nlambda_llds,testSize_llm=testSize_llds)
        perfList_llds[[i_llds]]<-cvglm_llds[[1]]
        coefList_llds[[i_llds]]<-cvglm_llds[[2]]
        testIndexList_llds[[i_llds]]<-cvglm_llds[[3]]

        rownames(perfList_llds[[i_llds]])<-colnames(Y_llds)[focalSample_llds]
		names(coefList_llds[[i_llds]])<-colnames(Y_llds)[focalSample_llds]
		names(testIndexList_llds[[i_llds]])<-colnames(Y_llds)[focalSample_llds]
		
        for(tt_llds in 1:length(focalSample_llds)){
		    rownames(coefList_llds[[i_llds]][[tt_llds]])<-c("intercept",colnames(X_llds)[focalSample_llds])
		}	

		
	}
	
	names(perfList_llds)<-paste0("X",downRatio_llds)
	names(coefList_llds)<-paste0("X",downRatio_llds)
	names(testIndexList_llds)<-paste0("X",downRatio_llds)
	

    return(list(perfList_llds,coefList_llds,testIndexList_llds))

}



### X_llm & Y_llm are both matrix ;

linear_learn_matrix<-function(X_llm,Y_llm,coincide_llm=F,unCorCutoff_llm=1.5,paraComp_llm=T,cores_llm=2,cycle_llm=1,alpha_llm=1,nlambda_llm=200,testSize_llm=100){
      
	#1: linear learning based on all variables ；without parallel computation
	if(coincide_llm==F&unCorCutoff_llm>1&paraComp_llm==F){
	
	    perfMat_llm<-matrix(nrow=dim(Y_llm)[2],ncol=cycle_llm)
		coefList_llm<-list()
		testIndexList_llm<-list()
		
		for(tt_llm in 1:dim(Y_llm)[2]){
	        coefMat_llm<-matrix(0,nrow=dim(X_llm)[2]+1,ncol=cycle_llm)        
			testIndexMat_llm<-matrix(nrow=testSize_llm,ncol=cycle_llm)
		    for(i_llm in 1:cycle_llm){
				cvglm_llm<-linear_learn(X_llm,Y_llm[,tt_llm],alpha_ll=alpha_llm,nlambda_ll=nlambda_llm,testSize_ll=testSize_llm)
	            perfMat_llm[tt_llm,i_llm]<-cvglm_llm[[1]][1]
				coefMat_llm[,i_llm]<-cvglm_llm[[2]][,1]
				testIndexMat_llm[,i_llm]<-cvglm_llm[[3]]
			}
		    coefList_llm[[tt_llm]]<-coefMat_llm
			testIndexList_llm[[tt_llm]]<-testIndexMat_llm
		}
		
	}
	
	#2 : linear learning based on all variables ；each column of X corresponds to each column of Y; without parallel computation
	if(coincide_llm==T&unCorCutoff_llm>1&paraComp_llm==F){
	
	    perfMat_llm<-matrix(nrow=dim(Y_llm)[2],ncol=cycle_llm)
		coefList_llm<-list()
		testIndexList_llm<-list()
		
		for(tt_llm in 1:dim(Y_llm)[2]){
		
		    coefMat_llm<-matrix(0,nrow=dim(X_llm)[2]+1,ncol=cycle_llm) 
			testIndexMat_llm<-matrix(nrow=testSize_llm,ncol=cycle_llm)
			
			unCorInd_llm<-setdiff(1:dim(X_llm)[2],tt_llm)		
			
		    for(i_llm in 1:cycle_llm){
				cvglm_llm<-linear_learn(X_llm[,unCorInd_llm],Y_llm[,tt_llm],alpha_ll=alpha_llm,nlambda_ll=nlambda_llm,testSize_ll=testSize_llm)
	            perfMat_llm[tt_llm,i_llm]<-cvglm_llm[[1]][1]
				coefMat_llm[c(1,unCorInd_llm+1),i_llm]<-cvglm_llm[[2]][,1]
				testIndexMat_llm[,i_llm]<-cvglm_llm[[3]]
			}
		    coefList_llm[[tt_llm]]<-coefMat_llm
			testIndexList_llm[[tt_llm]]<-testIndexMat_llm
		}
		
	}
	
	
	
	#3: linear learning based on correlation-constrained variables ；without parallel computation
	if(coincide_llm==F&unCorCutoff_llm<=1&paraComp_llm==F){
	    corMat_llm<-abs(cor(X_llm,Y_llm))
		
		
	    perfMat_llm<-matrix(nrow=dim(Y_llm)[2],ncol=cycle_llm)
		coefList_llm<-list()
		testIndexList_llm<list()
		
		for(tt_llm in 1:dim(Y_llm)[2]){
	        coefMat_llm<-matrix(0,nrow=dim(X_llm)[2]+1,ncol=cycle_llm)        
			testIndexMat_llm<-matrix(nrow=testSize_llm,ncol=cycle_llm)
			
			unCorInd_llm<-which(corMat_llm[,tt_llm]<unCorCutoff_llm)
			if(length(unCorInd_llm)==0){
			    perfMat_llm[tt_llm,]<-NA
				coefMat_llm[,i_llm]<-NA
				testIndexMat_llm[,i_llm]<-NA
			}
			else{
		        for(i_llm in 1:cycle_llm){
				    cvglm_llm<-linear_learn(X_llm[,unCorInd_llm],Y_llm[,tt_llm],alpha_ll=alpha_llm,nlambda_ll=nlambda_llm,testSize_ll=testSize_llm)
	                perfMat_llm[tt_llm,i_llm]<-cvglm_llm[[1]][1]
				    coefMat_llm[c(1,unCorInd_llm+1),i_llm]<-cvglm_llm[[2]][,1]
				    testIndexMat_llm[,i_llm]<-cvglm_llm[[3]]
			    }			    
			}

		    coefList_llm[[tt_llm]]<-coefMat_llm
			testIndexList_llm[[tt_llm]]<-testIndexMat_llm
		}
		
	}
		
	#4: linear learning based on correlation-constrained variables；each column of X corresponds to each column of Y ；without parallel computation
	if(coincide_llm==T&unCorCutoff_llm<=1&paraComp_llm==F){
	    corMat_llm<-abs(cor(X_llm,Y_llm))
		
		
	    perfMat_llm<-matrix(nrow=dim(Y_llm)[2],ncol=cycle_llm)
		coefList_llm<-list()
		testIndexList_llm<-list()
		
		for(tt_llm in 1:dim(Y_llm)[2]){
		
            coefMat_llm<-matrix(0,nrow=dim(X_llm)[2]+1,ncol=cycle_llm)        
			testIndexMat_llm<-matrix(nrow=testSize_llm,ncol=cycle_llm)
			
			unCorInd_llm<-setdiff(which(corMat_llm[,tt_llm]<unCorCutoff_llm),tt_llm)
			if(length(unCorInd_llm)==0){
			    perfMat_llm[tt_llm,]<-NA
				coefMat_llm[,i_llm]<-NA
				testIndexMat_llm[,i_llm]<-NA
			}
			else{
		        for(i_llm in 1:cycle_llm){
				    cvglm_llm<-linear_learn(X_llm[,unCorInd_llm],Y_llm[,tt_llm],alpha_ll=alpha_llm,nlambda_ll=nlambda_llm,testSize_ll=testSize_llm)
	                perfMat_llm[tt_llm,i_llm]<-cvglm_llm[[1]][1]
				    coefMat_llm[c(1,unCorInd_llm+1),i_llm]<-cvglm_llm[[2]][,1]
				    testIndexMat_llm[,i_llm]<-cvglm_llm[[3]]
			    }
			}
			
		    coefList_llm[[tt_llm]]<-coefMat_llm
			testIndexList_llm[[tt_llm]]<-testIndexMat_llm
			
		}
		
	}
		


    
	#5: linear learning based on all variables ；with parallel computation
	if(coincide_llm==F&unCorCutoff_llm>1&paraComp_llm==T){
	
	    perfMat_llm<-matrix(nrow=dim(Y_llm)[2],ncol=cycle_llm)
		coefList_llm<-list()
		testIndexList_llm<-list()
		
        cl_lmm <- makeCluster(cores_llm,type="FORK")
		
	    para_llm<-parLapply(cl_lmm,1:dim(Y_llm)[2],function(tt_llm){
            
			perf_llm<-c()
	        coefMat_llm<-matrix(0,nrow=dim(X_llm)[2]+1,ncol=cycle_llm)        
			testIndexMat_llm<-matrix(nrow=testSize_llm,ncol=cycle_llm)
			
		    for(i_llm in 1:cycle_llm){
				cvglm_llm<-linear_learn(X_llm,Y_llm[,tt_llm],alpha_ll=alpha_llm,nlambda_ll=nlambda_llm,testSize_ll=testSize_llm)
	            perf_llm[i_llm]<-cvglm_llm[[1]][1]
				coefMat_llm[,i_llm]<-cvglm_llm[[2]][,1]
				testIndexMat_llm[,i_llm]<-cvglm_llm[[3]]
			}
			
			return(list(perf_llm,coefMat_llm,testIndexMat_llm))
		
		})
		
		stopCluster(cl_lmm)		
		
		for(tt_llm in 1:dim(Y_llm)[2]){
		    perfMat_llm[tt_llm,]<-para_llm[[tt_llm]][[1]]
		    coefList_llm[[tt_llm]]<-para_llm[[tt_llm]][[2]]
			testIndexList_llm[[tt_llm]]<-para_llm[[tt_llm]][[3]]		    
			
		}
				
	}
	
	
	#6: linear learning based on all variables ；each column of X corresponds to each column of Y; with parallel computation
	if(coincide_llm==T&unCorCutoff_llm>1&paraComp_llm==T){
	
	    perfMat_llm<-matrix(nrow=dim(Y_llm)[2],ncol=cycle_llm)
		coefList_llm<-list()
		testIndexList_llm<-list()

        cl_lmm <- makeCluster(cores_llm,type="FORK")
		
	    para_llm<-parLapply(cl_lmm,1:dim(Y_llm)[2],function(tt_llm){
		
			perf_llm<-c()		
	        coefMat_llm<-matrix(0,nrow=dim(X_llm)[2]+1,ncol=cycle_llm)        
			testIndexMat_llm<-matrix(nrow=testSize_llm,ncol=cycle_llm)
			
			unCorInd_llm<-setdiff(1:dim(X_llm)[2],tt_llm)		
			
		    for(i_llm in 1:cycle_llm){
				cvglm_llm<-linear_learn(X_llm[,unCorInd_llm],Y_llm[,tt_llm],alpha_ll=alpha_llm,nlambda_ll=nlambda_llm,testSize_ll=testSize_llm)
	            perf_llm[i_llm]<-cvglm_llm[[1]][1]
				coefMat_llm[c(1,unCorInd_llm+1),i_llm]<-cvglm_llm[[2]][,1]
				testIndexMat_llm[,i_llm]<-cvglm_llm[[3]]
			}
			
			return(list(perf_llm,coefMat_llm,testIndexMat_llm))
			
		})
		
		stopCluster(cl_lmm)		
		
		for(tt_llm in 1:dim(Y_llm)[2]){
		    perfMat_llm[tt_llm,]<-para_llm[[tt_llm]][[1]]
		    coefList_llm[[tt_llm]]<-para_llm[[tt_llm]][[2]]
			testIndexList_llm[[tt_llm]]<-para_llm[[tt_llm]][[3]]		    
			
		}
		
	}
	
	
	
	#7: linear learning based on correlation-constrained variables ；with parallel computation
	if(coincide_llm==F&unCorCutoff_llm<=1&paraComp_llm==T){
	    corMat_llm<-abs(cor(X_llm,Y_llm))
				
	    perfMat_llm<-matrix(nrow=dim(Y_llm)[2],ncol=cycle_llm)
		coefList_llm<-list()
		testIndexList_llm<-list()

        cl_lmm <- makeCluster(cores_llm,type="FORK")
		
	    para_llm<-parLapply(cl_lmm,1:dim(Y_llm)[2],function(tt_llm){
            
			perf_llm<-c()
	        coefMat_llm<-matrix(0,nrow=dim(X_llm)[2]+1,ncol=cycle_llm)        
			testIndexMat_llm<-matrix(nrow=testSize_llm,ncol=cycle_llm)
			
			unCorInd_llm<-which(corMat_llm[,tt_llm]<unCorCutoff_llm)
			if(length(unCorInd_llm)==0){
			    perf_llm<-rep(NA,cycle_llm)
				coefMat_llm[,i_llm]<-NA
				testIndexMat_llm[,i_llm]<-NA
			}
			else{
		        for(i_llm in 1:cycle_llm){
				    cvglm_llm<-linear_learn(X_llm[,unCorInd_llm],Y_llm[,tt_llm],alpha_ll=alpha_llm,nlambda_ll=nlambda_llm,testSize_ll=testSize_llm)
	                perf_llm[i_llm]<-cvglm_llm[[1]][1]
				    coefMat_llm[c(1,unCorInd_llm+1),i_llm]<-cvglm_llm[[2]][,1]
				    testIndexMat_llm[,i_llm]<-cvglm_llm[[3]]
			    }
			}
			return(list(perf_llm,coefMat_llm,testIndexMat_llm))
			
		})
		
		stopCluster(cl_lmm)		
		
		for(tt_llm in 1:dim(Y_llm)[2]){
		    perfMat_llm[tt_llm,]<-para_llm[[tt_llm]][[1]]
		    coefList_llm[[tt_llm]]<-para_llm[[tt_llm]][[2]]
			testIndexList_llm[[tt_llm]]<-para_llm[[tt_llm]][[3]]		    
			
		}
		
	}
		
	#8: linear learning based on correlation-constrained variables；each column of X corresponds to each column of Y ；with parallel computation
	if(coincide_llm==T&unCorCutoff_llm<=1&paraComp_llm==T){
	    corMat_llm<-abs(cor(X_llm,Y_llm))
				
	    perfMat_llm<-matrix(nrow=dim(Y_llm)[2],ncol=cycle_llm)
		coefList_llm<-list()
		testIndexList_llm<-list()

        cl_lmm <- makeCluster(cores_llm,type="FORK")
		
	    para_llm<-parLapply(cl_lmm,1:dim(Y_llm)[2],function(tt_llm){
            
			perf_llm<-c()
	        coefMat_llm<-matrix(0,nrow=dim(X_llm)[2]+1,ncol=cycle_llm)        
			testIndexMat_llm<-matrix(nrow=testSize_llm,ncol=cycle_llm)
			
			unCorInd_llm<-setdiff(which(corMat_llm[,tt_llm]<unCorCutoff_llm),tt_llm)
			if(length(unCorInd_llm)==0){
			    perf_llm<-rep(NA,cycle_llm)
				coefMat_llm[,i_llm]<-NA
				testIndexMat_llm[,i_llm]<-NA
			}
			else{
		        for(i_llm in 1:cycle_llm){
				    cvglm_llm<-linear_learn(X_llm[,unCorInd_llm],Y_llm[,tt_llm],alpha_ll=alpha_llm,nlambda_ll=nlambda_llm,testSize_ll=testSize_llm)
	                perf_llm[i_llm]<-cvglm_llm[[1]][1]
				    coefMat_llm[c(1,unCorInd_llm+1),i_llm]<-cvglm_llm[[2]][,1]
				    testIndexMat_llm[,i_llm]<-cvglm_llm[[3]]
			    }
			}
			
			return(list(perf_llm,coefMat_llm,testIndexMat_llm))
			
		})
		
		stopCluster(cl_lmm)		
		
		for(tt_llm in 1:dim(Y_llm)[2]){
		    perfMat_llm[tt_llm,]<-para_llm[[tt_llm]][[1]]
		    coefList_llm[[tt_llm]]<-para_llm[[tt_llm]][[2]]
			testIndexList_llm[[tt_llm]]<-para_llm[[tt_llm]][[3]]		    
			
		}
		
	}


    rownames(perfMat_llm)<-colnames(Y_llm)
	names(coefList_llm)<-colnames(Y_llm)
	names(testIndexList_llm)<-colnames(Y_llm)
	for(tt_llm in 1:length(coefList_llm)){
	    rownames(coefList_llm[[tt_llm]])<-c("intercept",colnames(X_llm))
	}
		
	
    return(list(perfMat_llm,coefList_llm,testIndexList_llm))

}



### X_llm & Y_llm are both matrix ;

linear_learn_matrix_shuffle<-function(X_llms,Y_llms,coincide_llms=F,unCorCutoff_llms=1.5,paraComp_llms=T,cores_llms=2,cycle_llms=1,alpha_llms=1,nlambda_llms=200,testSize_llms=100){
		
	#8: linear learning based on correlation-constrained variables；each column of X corresponds to each column of Y ；with parallel computation
	if(coincide_llms==T&unCorCutoff_llms<=1&paraComp_llms==T){
	    corMat_llms<-abs(cor(X_llms,Y_llms))
				
	    perfMat_llms<-matrix(nrow=dim(Y_llms)[2],ncol=cycle_llms)
		coefList_llms<-list()
		testIndexList_llms<-list()

        cl_lmms <- makeCluster(cores_llms,type="FORK")
		
	    para_llms<-parLapply(cl_lmms,1:dim(Y_llms)[2],function(tt_llms){
            
			perf_llms<-c()
	        coefMat_llms<-matrix(0,nrow=dim(X_llms)[2]+1,ncol=cycle_llms)        
			testIndexMat_llms<-matrix(nrow=testSize_llms,ncol=cycle_llms)
			
			unCorInd_llms<-setdiff(which(corMat_llms[,tt_llms]<unCorCutoff_llms),tt_llms)
			if(length(unCorInd_llms)==0){
			    perf_llms<-rep(NA,cycle_llms)
				coefMat_llms[,i_llms]<-NA
				testIndexMat_llms[,i_llms]<-NA
			}
			else{
		        for(i_llms in 1:cycle_llms){
				    set.seed(i_llms)
				    focalShuffleInd<-sample(1:dim(Y_llms)[1],dim(Y_llms)[1])
				    cvglm_llms<-linear_learn(X_llms[,unCorInd_llms],Y_llms[focalShuffleInd,tt_llms],alpha_ll=alpha_llms,nlambda_ll=nlambda_llms,testSize_ll=testSize_llms)
	                perf_llms[i_llms]<-cvglm_llms[[1]][1]
				    coefMat_llms[c(1,unCorInd_llms+1),i_llms]<-cvglm_llms[[2]][,1]
				    testIndexMat_llms[,i_llms]<-cvglm_llms[[3]]
			    }
			}
			
			return(list(perf_llms,coefMat_llms,testIndexMat_llms))
			
		})
		
		stopCluster(cl_lmms)		
		
		for(tt_llms in 1:dim(Y_llms)[2]){
		    perfMat_llms[tt_llms,]<-para_llms[[tt_llms]][[1]]
		    coefList_llms[[tt_llms]]<-para_llms[[tt_llms]][[2]]
			testIndexList_llms[[tt_llms]]<-para_llms[[tt_llms]][[3]]		    
			
		}
		
	}


    rownames(perfMat_llms)<-colnames(Y_llms)
	names(coefList_llms)<-colnames(Y_llms)
	names(testIndexList_llms)<-colnames(Y_llms)
	for(tt_llms in 1:length(coefList_llms)){
	    rownames(coefList_llms[[tt_llms]])<-c("intercept",colnames(X_llms))
	}
		
	
    return(list(perfMat_llms,coefList_llms,testIndexList_llms))

}


############################### modelling ##################################

########################################################################################


################# UBHDD ###########################



load("data/basic/pdata.02.RData")


pcMatCutoff<-sqrt(qt(1-0.005/404,815-2)^2/(815-2+qt(1-0.005/404,815-2)^2))


cvglm_unCorr_one_rep<-linear_learn_matrix(pdata.02,pdata.02,coincide_llm=T,unCorCutoff_llm=pcMatCutoff,paraComp_llm=T,cores_llm=100,cycle_llm=20,alpha_llm=1,nlambda_llm=200,testSize_llm=100)



perfMat_unCorr_one_rep<-cvglm_unCorr_one_rep[[1]]
coefList_unCorr_one_rep<-cvglm_unCorr_one_rep[[2]]
testIndexList_unCorr_one_rep<-cvglm_unCorr_one_rep[[3]]



save.image("data/oneRep/perfMat_unCorr_one_rep_20cycles.RData")


#################### Shuffling #################################



pcMatCutoff<-sqrt(qt(1-0.005/404,815-2)^2/(815-2+qt(1-0.005/404,815-2)^2))


cvglm_unCorr_one_rep_shuffle<-linear_learn_matrix_shuffle(pdata.02,pdata.02,coincide_llms=T,unCorCutoff_llms=pcMatCutoff,paraComp_llms=T,cores_llms=100,cycle_llms=20,alpha_llms=1,nlambda_llms=200,testSize_llms=100)



perfMat_unCorr_one_rep_shuffle<-cvglm_unCorr_one_rep_shuffle[[1]]
coefList_unCorr_one_rep_shuffle<-cvglm_unCorr_one_rep_shuffle[[2]]
testIndexList_unCorr_one_rep_shuffle<-cvglm_unCorr_one_rep_shuffle[[3]]



save.image("data/oneRep/perfMat_unCorr_one_rep_shuffle_20cycles.RData")






#################### TOTAL #####################################


###

cvglm_total_one_rep<-linear_learn_matrix(pdata.02,pdata.02,coincide_llm=T,unCorCutoff_llm=1,paraComp_llm=T,cores_llm=50,cycle_llm=20,alpha_llm=1,nlambda_llm=200,testSize_llm=100)


perfMat_total_one_rep<-cvglm_total_one_rep[[1]]
coefList_total_one_rep<-cvglm_total_one_rep[[2]]
testIndexList_total_one_rep<-cvglm_total_one_rep[[3]]


save.image("data/oneRep/perfMat_total_one_rep_20cycles.RData")


####################### robust of Ru ###################################



### 0.1

cvglm_unCorr01_one_rep<-linear_learn_matrix(pdata.02,pdata.02,coincide_llm=T,unCorCutoff_llm=0.1,paraComp_llm=T,cores_llm=100,cycle_llm=20,alpha_llm=1,nlambda_llm=200,testSize_llm=100)


perfMat_unCorr01_one_rep<-cvglm_unCorr01_one_rep[[1]]
coefList_unCorr01_one_rep<-cvglm_unCorr01_one_rep[[2]]
testIndexList_unCorr01_one_rep<-cvglm_unCorr01_one_rep[[3]]


save.image("data/oneRep/perfMat_unCorr01_one_rep_20cycles.RData")

### 0.11

cvglm_unCorr011_one_rep<-linear_learn_matrix(pdata.02,pdata.02,coincide_llm=T,unCorCutoff_llm=0.11,paraComp_llm=T,cores_llm=50,cycle_llm=20,alpha_llm=1,nlambda_llm=200,testSize_llm=100)


perfMat_unCorr011_one_rep<-cvglm_unCorr011_one_rep[[1]]
coefList_unCorr011_one_rep<-cvglm_unCorr011_one_rep[[2]]
testIndexList_unCorr011_one_rep<-cvglm_unCorr011_one_rep[[3]]


save.image("data/oneRep/perfMat_unCorr011_one_rep_20cycles.RData")


### 0.12
rm(cvglm_unCorr011_one_rep)
rm(perfMat_unCorr011_one_rep)
rm(coefList_unCorr011_one_rep)
rm(testIndexList_unCorr011_one_rep)


cvglm_unCorr012_one_rep<-linear_learn_matrix(pdata.02,pdata.02,coincide_llm=T,unCorCutoff_llm=0.12,paraComp_llm=T,cores_llm=50,cycle_llm=20,alpha_llm=1,nlambda_llm=200,testSize_llm=100)


perfMat_unCorr012_one_rep<-cvglm_unCorr012_one_rep[[1]]
coefList_unCorr012_one_rep<-cvglm_unCorr012_one_rep[[2]]
testIndexList_unCorr012_one_rep<-cvglm_unCorr012_one_rep[[3]]


save.image("data/oneRep/perfMat_unCorr012_one_rep_20cycles.RData")



### 0.13
rm(cvglm_unCorr012_one_rep)
rm(perfMat_unCorr012_one_rep)
rm(coefList_unCorr012_one_rep)
rm(testIndexList_unCorr012_one_rep)


cvglm_unCorr013_one_rep<-linear_learn_matrix(pdata.02,pdata.02,coincide_llm=T,unCorCutoff_llm=0.13,paraComp_llm=T,cores_llm=50,cycle_llm=20,alpha_llm=1,nlambda_llm=200,testSize_llm=100)


perfMat_unCorr013_one_rep<-cvglm_unCorr013_one_rep[[1]]
coefList_unCorr013_one_rep<-cvglm_unCorr013_one_rep[[2]]
testIndexList_unCorr013_one_rep<-cvglm_unCorr013_one_rep[[3]]


save.image("data/oneRep/perfMat_unCorr013_one_rep_20cycles.RData")


### 0.14
rm(cvglm_unCorr013_one_rep)
rm(perfMat_unCorr013_one_rep)
rm(coefList_unCorr013_one_rep)
rm(testIndexList_unCorr013_one_rep)


cvglm_unCorr014_one_rep<-linear_learn_matrix(pdata.02,pdata.02,coincide_llm=T,unCorCutoff_llm=0.14,paraComp_llm=T,cores_llm=50,cycle_llm=20,alpha_llm=1,nlambda_llm=200,testSize_llm=100)


perfMat_unCorr014_one_rep<-cvglm_unCorr014_one_rep[[1]]
coefList_unCorr014_one_rep<-cvglm_unCorr014_one_rep[[2]]
testIndexList_unCorr014_one_rep<-cvglm_unCorr014_one_rep[[3]]


save.image("data/oneRep/perfMat_unCorr014_one_rep_20cycles.RData")



### 0.15

rm(cvglm_unCorr01_one_rep)
rm(perfMat_unCorr01_one_rep)
rm(coefList_unCorr01_one_rep)
rm(testIndexList_unCorr01_one_rep)

cvglm_unCorr015_one_rep<-linear_learn_matrix(pdata.02,pdata.02,coincide_llm=T,unCorCutoff_llm=0.15,paraComp_llm=T,cores_llm=50,cycle_llm=20,alpha_llm=1,nlambda_llm=200,testSize_llm=100)


perfMat_unCorr015_one_rep<-cvglm_unCorr015_one_rep[[1]]
coefList_unCorr015_one_rep<-cvglm_unCorr015_one_rep[[2]]
testIndexList_unCorr015_one_rep<-cvglm_unCorr015_one_rep[[3]]


save.image("data/oneRep/perfMat_unCorr015_one_rep_20cycles.RData")



### 0.16
rm(cvglm_unCorr014_one_rep)
rm(perfMat_unCorr014_one_rep)
rm(coefList_unCorr014_one_rep)
rm(testIndexList_unCorr014_one_rep)


cvglm_unCorr016_one_rep<-linear_learn_matrix(pdata.02,pdata.02,coincide_llm=T,unCorCutoff_llm=0.16,paraComp_llm=T,cores_llm=50,cycle_llm=20,alpha_llm=1,nlambda_llm=200,testSize_llm=100)


perfMat_unCorr016_one_rep<-cvglm_unCorr016_one_rep[[1]]
coefList_unCorr016_one_rep<-cvglm_unCorr016_one_rep[[2]]
testIndexList_unCorr016_one_rep<-cvglm_unCorr016_one_rep[[3]]


save.image("data/oneRep/perfMat_unCorr016_one_rep_20cycles.RData")


### 0.17
rm(cvglm_unCorr016_one_rep)
rm(perfMat_unCorr016_one_rep)
rm(coefList_unCorr016_one_rep)
rm(testIndexList_unCorr016_one_rep)


cvglm_unCorr017_one_rep<-linear_learn_matrix(pdata.02,pdata.02,coincide_llm=T,unCorCutoff_llm=0.17,paraComp_llm=T,cores_llm=50,cycle_llm=20,alpha_llm=1,nlambda_llm=200,testSize_llm=100)


perfMat_unCorr017_one_rep<-cvglm_unCorr017_one_rep[[1]]
coefList_unCorr017_one_rep<-cvglm_unCorr017_one_rep[[2]]
testIndexList_unCorr017_one_rep<-cvglm_unCorr017_one_rep[[3]]


save.image("data/oneRep/perfMat_unCorr017_one_rep_20cycles.RData")


### 0.18
rm(cvglm_unCorr017_one_rep)
rm(perfMat_unCorr017_one_rep)
rm(coefList_unCorr017_one_rep)
rm(testIndexList_unCorr017_one_rep)


cvglm_unCorr018_one_rep<-linear_learn_matrix(pdata.02,pdata.02,coincide_llm=T,unCorCutoff_llm=0.18,paraComp_llm=T,cores_llm=50,cycle_llm=20,alpha_llm=1,nlambda_llm=200,testSize_llm=100)


perfMat_unCorr018_one_rep<-cvglm_unCorr018_one_rep[[1]]
coefList_unCorr018_one_rep<-cvglm_unCorr018_one_rep[[2]]
testIndexList_unCorr018_one_rep<-cvglm_unCorr018_one_rep[[3]]


save.image("data/oneRep/perfMat_unCorr018_one_rep_20cycles.RData")


### 0.19
rm(cvglm_unCorr018_one_rep)
rm(perfMat_unCorr018_one_rep)
rm(coefList_unCorr018_one_rep)
rm(testIndexList_unCorr018_one_rep)


cvglm_unCorr019_one_rep<-linear_learn_matrix(pdata.02,pdata.02,coincide_llm=T,unCorCutoff_llm=0.19,paraComp_llm=T,cores_llm=50,cycle_llm=20,alpha_llm=1,nlambda_llm=200,testSize_llm=100)


perfMat_unCorr019_one_rep<-cvglm_unCorr019_one_rep[[1]]
coefList_unCorr019_one_rep<-cvglm_unCorr019_one_rep[[2]]
testIndexList_unCorr019_one_rep<-cvglm_unCorr019_one_rep[[3]]


save.image("data/oneRep/perfMat_unCorr019_one_rep_20cycles.RData")



### 0.2

rm(cvglm_unCorr015_one_rep)
rm(perfMat_unCorr015_one_rep)
rm(coefList_unCorr015_one_rep)
rm(testIndexList_unCorr015_one_rep)

cvglm_unCorr02_one_rep<-linear_learn_matrix(pdata.02,pdata.02,coincide_llm=T,unCorCutoff_llm=0.2,paraComp_llm=T,cores_llm=50,cycle_llm=20,alpha_llm=1,nlambda_llm=200,testSize_llm=100)


perfMat_unCorr02_one_rep<-cvglm_unCorr02_one_rep[[1]]
coefList_unCorr02_one_rep<-cvglm_unCorr02_one_rep[[2]]
testIndexList_unCorr02_one_rep<-cvglm_unCorr02_one_rep[[3]]


save.image("data/oneRep/perfMat_unCorr02_one_rep_20cycles.RData")


######################### down sampling ###############################



pcMatCutoff<-sqrt(qt(1-0.005/404,815-2)^2/(815-2+qt(1-0.005/404,815-2)^2))


cvglm_unCorr_one_rep_down_sampling<-linear_learn_down_sampling(pdata.02,pdata.02,coincide_llds=T,unCorCutoff_llds=pcMatCutoff,paraComp_llds=T,cores_llds=100,cycle_llds=1,alpha_llds=1,nlambda_llds=200,testSize_llds=100)



perfList_unCorr_one_rep_down_sampling<-cvglm_unCorr_one_rep_down_sampling[[1]]
coefList_unCorr_one_rep_down_sampling<-cvglm_unCorr_one_rep_down_sampling[[2]]
testIndexList_unCorr_one_rep_down_sampling<-cvglm_unCorr_one_rep_down_sampling[[3]]



save.image("data/downSampling/perfList_unCorr_one_rep_down_sampling_1cycle.RData")


########################################### robustness between conditions ###################################

setwd("/home/wjg/UBHDD/yeastDeletome")

load("data/mt4718.RData")
load("data/pdata.02_ori.RData")


###
commonTrait<-colnames(pdata.02_ori)

mt4718_common<-mt4718[,match(commonTrait,colnames(mt4718))]

mt4718_common_scale<-foreach(i=1:405,.combine=cbind)%do%{
    (mt4718_common[,i]-mean(pdata.02_ori[,i],na.rm=T))/sd(pdata.02_ori[,i],na.rm=T)
}

colnames(mt4718_common_scale)<-colnames(pdata.02_ori)


pcMatCutoff<-sqrt(qt(1-0.005/404,815-2)^2/(815-2+qt(1-0.005/404,815-2)^2))


cvglm_mt4718_common_scale<-linear_learn_matrix(mt4718_common_scale,mt4718_common_scale,coincide_llm=T,unCorCutoff_llm=pcMatCutoff,paraComp_llm=T,cores_llm=100,cycle_llm=20,alpha_llm=1,nlambda_llm=200,testSize_llm=100)



perfMat_unCorr_mt4718_common_scale<-cvglm_mt4718_common_scale[[1]]
coefList_unCorr_mt4718_common_scale<-cvglm_mt4718_common_scale[[2]]
testIndexList_unCorr_mt4718_common_scale<-cvglm_mt4718_common_scale[[3]]



save.image("data/perfMat_unCorr_mt4718_common_scale_20cycles.RData")








