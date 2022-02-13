
# 02/13/22
# R script in:
# The trait coding rule in phenotype space
# 
# Jianguo Wang (wangjg36@mail.sysu.edu.cn)
# 



################### library ###########################################

setwd("/home/wjg/UBHDD/theory")



library(glmnet)
library(MASS)
library(parallel)
library(foreach)
library(doParallel)


################### functions ##########################################


### calculate uncorrelation cutoff
calUnCorCutoff<-function(samSize,traitNum,alpha){
    return(sqrt(qt(1-alpha*0.5/(traitNum-1),samSize-2)^2/(samSize-2+qt(1-alpha*0.5/(traitNum-1),samSize-2)^2)))
}



### get simulated traits
get_sim_trait<-function(sim_PG_dim,sim_PNG_dim,sim_PG_dim_per_trait,sim_PNG_dim_per_trait=NULL,sim_cluster,sim_sam_size=1000,sim_H2_min=0,sim_H2_max=1){
	
    sim_trait_num<-sum(sim_cluster)
	
    sim_cluster_len<-length(sim_cluster)

    if(is.null(sim_PNG_dim_per_trait)){
	    sim_PNG_dim_per_trait=sim_PG_dim_per_trait
	}

    if(length(table(sim_cluster))==1){
	    sim_struct<-FALSE
	}else{
	    sim_struct<-TRUE
	}


    
    #sim basic traits satisfying standard multiple normal distribution


    ## mu and sigma of basic traits

    #sim_bas_mu<-rep(0,sim_PNG_dim+sim_PG_dim)
    #sim_bas_sigma<-matrix(0,nrow=sim_PNG_dim+sim_PG_dim,ncol=sim_PNG_dim+sim_PG_dim)
    #diag(sim_bas_sigma)<-1

    ## basic traits

    #set.seed(10)
    #sim_bas<-mvrnorm(sim_sam_size,mu=sim_bas_mu,Sigma=sim_bas_sigma,tol=1)

    #sim_bas_PG<-sim_bas[,1:sim_PG_dim]
    #sim_bas_PNG<-sim_bas[,(sim_PG_dim+1):(sim_PG_dim+sim_PNG_dim)]

    #save time by generating in advance
    sim_bas_PG<-mvnormPre[,1:sim_PG_dim]
    sim_bas_PNG<-mvnormPre[,(sim_PG_dim+1):(sim_PG_dim+sim_PNG_dim)]

    ### random sampling  

    ## PG component (random sampling)

    set.seed(10)
    sim_PG_coef<-matrix(rnorm(sim_PG_dim*sim_cluster_len),nrow=sim_PG_dim,ncol=sim_cluster_len)
    sim_PG_coef<-foreach(i=1:sim_cluster_len,.combine=cbind)%do%{
        set.seed(i)
	    temp<-sim_PG_coef[,i]
        temp[sample(1:sim_PG_dim,sim_PG_dim-sim_PG_dim_per_trait)]<-0
	    return(temp)
    }


    sim_PG_random<-sim_bas_PG%*%sim_PG_coef
    sim_PG_random<-scale(sim_PG_random)

    ## PNG component (random sampling)

    set.seed(10)
    sim_PNG_coef<-matrix(rnorm(sim_PNG_dim*sim_cluster_len),nrow=sim_PNG_dim,ncol=sim_cluster_len)
    sim_PNG_coef<-foreach(i=1:sim_cluster_len,.combine=cbind)%do%{
        set.seed(i)
	    temp<-sim_PNG_coef[,i]
        temp[sample(1:sim_PNG_dim,sim_PNG_dim-sim_PNG_dim_per_trait)]<-0
	    return(temp)
    }


    sim_PNG_random<-sim_bas_PNG%*%sim_PNG_coef
    sim_PNG_random<-scale(sim_PNG_random)

    ## setting H2 [random] obeying uniform distribution
	if(sim_struct){
	    set.seed(1)
		#here can be as a parameter
		sim_struct_len<-length(which(sim_cluster>min(sim_cluster)))
		sim_H2_struct<-(c(1:sim_struct_len)-0.5)/sim_struct_len*(sim_H2_max-sim_H2_min)+sim_H2_min
		sim_H2_unstruct<-runif(sim_cluster_len-sim_struct_len,sim_H2_min,sim_H2_max)
		sim_H2_random<-c(sim_H2_struct,sim_H2_unstruct)
	}else{
	    set.seed(1)
	    sim_H2_random<-runif(sim_cluster_len,sim_H2_min,sim_H2_max)
	}
    
	## get the apportion of PG [random]
	
    sim_apportion_random<-1/(sqrt(1/sim_H2_random-1)+1)

    ### dense sampling
    
	## setting H2 [dense] approximate to H2 [random]
	
    sim_H2_dense<-lapply(1:sim_cluster_len,function(i){
		sim_H2_dense_min<-ifelse(sim_H2_random[i]-0.2<sim_H2_min,sim_H2_min,sim_H2_random[i]-0.2)
		sim_H2_dense_max<-ifelse(sim_H2_random[i]+0.2>sim_H2_max,sim_H2_max,sim_H2_random[i]+0.2)
	    set.seed(1)		
		sim_H2_dense_temp<-runif(sim_cluster[i]-1,sim_H2_dense_min,sim_H2_dense_max)
        return(sim_H2_dense_temp)
	})
	
	## get the apportion of PG [dense]
	
    sim_apportion_dense<-lapply(sim_H2_dense,function(xx){
        1/(sqrt(1/xx-1)+1)
    })
	
 
	### combine
	
	sim_H2<-c(sim_H2_random,unlist(sim_H2_dense))
	
	sim_apportion<-c(sim_apportion_random,unlist(sim_apportion_dense))
	
	sim_index<-c(1:sim_cluster_len,unlist(lapply(1:sim_cluster_len,function(xxx){rep(xxx,sim_cluster[xxx]-1)})))

    ### obtain sim trait, PG and PNG
	
	sim_PG<-foreach(i=1:sim_trait_num,.combine=cbind)%do%{sim_apportion[i]*sim_PG_random[,sim_index[i]]}
	sim_PNG<-foreach(i=1:sim_trait_num,.combine=cbind)%do%{(1-sim_apportion[i])*sim_PNG_random[,sim_index[i]]}
	sim_trait<-sim_PG+sim_PNG
	
	### output
	
    sim_trait_mean<-apply(sim_trait,2,mean)
    sim_trait_sd<-apply(sim_trait,2,sd) 	
	
    ## sim_trait_scale	
	sim_trait_scale<-scale(sim_trait)	
	rownames(sim_trait_scale)<-paste0("S",1:sim_sam_size)
    colnames(sim_trait_scale)<-paste0("T",1:sim_trait_num)

    ## sim_PG_component and sim_PNG_component  
    sim_PG_component<-foreach(i=1:sim_trait_num,.combine=cbind)%do%{(sim_PG[,i]-sim_trait_mean[i])/sim_trait_sd[i]}
    rownames(sim_PG_component)<-paste0("S",1:sim_sam_size)
    colnames(sim_PG_component)<-paste0("T",1:sim_trait_num)

    sim_PNG_component<-foreach(i=1:sim_trait_num,.combine=cbind)%do%{(sim_PNG[,i]-sim_trait_mean[i])/sim_trait_sd[i]}
    rownames(sim_PNG_component)<-paste0("S",1:sim_sam_size)
    colnames(sim_PNG_component)<-paste0("T",1:sim_trait_num) 
    
	## dim usage
	sim_PG_coef_all<-sim_PG_coef[,sim_index]
	sim_PNG_coef_all<-sim_PNG_coef[,sim_index]
		
    ## recalculate H2

    sim_H2_recal1<-sim_apportion^2/(sim_apportion^2+(1-sim_apportion)^2)
	sim_H2_recal2<-apply(sim_PG_component,2,var)
		
    names(sim_H2)<-paste0("T",1:sim_trait_num)
	names(sim_H2_recal1)<-paste0("T",1:sim_trait_num)
	names(sim_H2_recal2)<-paste0("T",1:sim_trait_num)
	

    return(list(sim_trait_scale,sim_PG_component,sim_PNG_component,sim_H2,sim_H2_recal1,sim_H2_recal2,sim_PG_coef_all,sim_PNG_coef_all))

}


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



### linear_learn for matrix
UBHDD<-function(X,Y,unCorCutoff=1,coreNum=2,alpha_seq1=1,nlambda1=200,testSize1=100){

    corMatTemp<-abs(cor(Y,X))
				
	perfMatTemp<-rep(NA,dim(Y)[2])
	coefMatTemp<-matrix(0,nrow=dim(X)[2]+1,ncol=dim(Y)[2])
	testIndexMatTemp<-matrix(NA,nrow=testSize1,ncol=dim(Y)[2])

    clTemp <- makeCluster(coreNum,type="FORK")
		
	tempPARA4<-parLapply(clTemp,1:dim(Y)[2],function(tt){
	
        #		
	    unCorIndTemp<-which(corMatTemp[tt,]<unCorCutoff)
			
	    cvglmTemp4<-linear_learn(X[,unCorIndTemp],Y[,tt],alpha_ll=alpha_seq1,nlambda_ll=nlambda1,testSize_ll=testSize1)
	    perfMatTemp[tt]<-cvglmTemp4[[1]][1]
		coefMatTemp[c(1,unCorIndTemp+1),tt]<-cvglmTemp4[[2]][[1]][,1]
		testIndexMatTemp[,tt]<-cvglmTemp4[[3]][[1]]
			
		return(list(perfMatTemp[tt],coefMatTemp[,tt],testIndexMatTemp[,tt]))
			
	})
		
	stopCluster(clTemp)		
		
	perfMatTemp<-foreach(tt=1:dim(Y)[2],.combine=c)%do%{
	    tempPARA4[[tt]][[1]]
	}
	coefMatTemp<-foreach(tt=1:dim(Y)[2],.combine=cbind)%do%{
	    tempPARA4[[tt]][[2]]
	}
	testIndexMatTemp<-foreach(tt=1:dim(Y)[2],.combine=cbind)%do%{
	    tempPARA4[[tt]][[3]]
	}	
		
    names(perfMatTemp)<-colnames(Y)
	colnames(coefMatTemp)<-colnames(Y)
	rownames(coefMatTemp)<-c("intercept",colnames(X))
	colnames(testIndexMatTemp)<-colnames(Y)
	
    return(list(perfMatTemp,coefMatTemp,testIndexMatTemp))

}


########################## simulation ##########################################





# select from this can reduce time
#set.seed(1)
#mvnormPre<-mvrnorm(1000,mu=rep(0,11000),Sigma=diag(11000),tol=1)
#save(mvnormPre,file="data/mvnormPre.RData")

load("data/mvnormPre.RData")

############################################################################

####################### main || format:  [sim_PG_dim,sim_PNG_dim,sim_cluster[,]]

#### 

### case 1 : unstructural phenotype space

## parameter

sim_PG_dim1=10
sim_PNG_dim1=10000
sim_PG_dim_per_trait1=sim_PG_dim1/2
sim_cluster1=rep(10,100)

## get simulated phenotype space

sim_10_10000_10_100<-get_sim_trait(
sim_PG_dim=sim_PG_dim1,
sim_PNG_dim=sim_PNG_dim1,
sim_PG_dim_per_trait=sim_PG_dim_per_trait1,
sim_cluster=sim_cluster1)

## set uncorrelation cutoff

sim_unCorCutoff<-calUnCorCutoff(
samSize=1000,
traitNum=sum(sim_cluster1),
alpha=0.01)

## UBHDD modelling

sim_UBHDD_10_10000_10_100<-UBHDD(
X=sim_10_10000_10_100[[1]],
Y=sim_10_10000_10_100[[1]],
unCorCutoff=sim_unCorCutoff,
coreNum=50)

sim_UBHDD_10_10000_10_100_perfMat<-sim_UBHDD_10_10000_10_100[[1]]
sim_UBHDD_10_10000_10_100_coefMat<-sim_UBHDD_10_10000_10_100[[2]]
sim_UBHDD_10_10000_10_100_testIndex<-sim_UBHDD_10_10000_10_100[[3]]

save.image("data/sim_UBHDD_10_10000_10_100.RData")


### case 2 : structured phenotype space 

## parameter

sim_PG_dim1=10
sim_PNG_dim1=10000
sim_PG_dim_per_trait1=sim_PG_dim1/2
sim_cluster1=c(300,200,rep(10,50))

## get simulated phenotype space

sim_10_10000_struct<-get_sim_trait(
sim_PG_dim=sim_PG_dim1,
sim_PNG_dim=sim_PNG_dim1,
sim_PG_dim_per_trait=sim_PG_dim_per_trait1,
sim_cluster=sim_cluster1)

## set uncorrelation cutoff

sim_unCorCutoff<-calUnCorCutoff(
samSize=1000,
traitNum=sum(sim_cluster1),
alpha=0.01)

## UBHDD modelling

sim_UBHDD_10_10000_struct<-UBHDD(
X=sim_10_10000_struct[[1]],
Y=sim_10_10000_struct[[1]],
unCorCutoff=sim_unCorCutoff,
coreNum=50)

sim_UBHDD_10_10000_struct_perfMat<-sim_UBHDD_10_10000_struct[[1]]
sim_UBHDD_10_10000_struct_coefMat<-sim_UBHDD_10_10000_struct[[2]]
sim_UBHDD_10_10000_struct_testIndex<-sim_UBHDD_10_10000_struct[[3]]

save.image("data/sim_UBHDD_10_10000_struct.RData")



#####################################################################################

############## the dimensionality of PG (N1=20, 50, 100) ###############################################################

### N1=20

## parameter

sim_PG_dim1=20
sim_PNG_dim1=10000
sim_PG_dim_per_trait1=sim_PG_dim1/2
sim_cluster1=rep(10,100)

## get simulated phenotype space

sim_20_10000_10_100<-get_sim_trait(
sim_PG_dim=sim_PG_dim1,
sim_PNG_dim=sim_PNG_dim1,
sim_PG_dim_per_trait=sim_PG_dim_per_trait1,
sim_cluster=sim_cluster1)

## set uncorrelation cutoff

sim_unCorCutoff<-calUnCorCutoff(
samSize=1000,
traitNum=sum(sim_cluster1),
alpha=0.01)

## UBHDD modelling

sim_UBHDD_20_10000_10_100<-UBHDD(
X=sim_20_10000_10_100[[1]],
Y=sim_20_10000_10_100[[1]],
unCorCutoff=sim_unCorCutoff,
coreNum=50)

sim_UBHDD_20_10000_10_100_perfMat<-sim_UBHDD_20_10000_10_100[[1]]
sim_UBHDD_20_10000_10_100_coefMat<-sim_UBHDD_20_10000_10_100[[2]]
sim_UBHDD_20_10000_10_100_testIndex<-sim_UBHDD_20_10000_10_100[[3]]

save.image("data/sim_UBHDD_20_10000_10_100.RData")



### N1=50

## parameter

sim_PG_dim1=50
sim_PNG_dim1=10000
sim_PG_dim_per_trait1=sim_PG_dim1/2
sim_cluster1=rep(10,100)

## get simulated phenotype space

sim_50_10000_10_100<-get_sim_trait(
sim_PG_dim=sim_PG_dim1,
sim_PNG_dim=sim_PNG_dim1,
sim_PG_dim_per_trait=sim_PG_dim_per_trait1,
sim_cluster=sim_cluster1)

## set uncorrelation cutoff

sim_unCorCutoff<-calUnCorCutoff(
samSize=1000,
traitNum=sum(sim_cluster1),
alpha=0.01)

## UBHDD modelling

sim_UBHDD_50_10000_10_100<-UBHDD(
X=sim_50_10000_10_100[[1]],
Y=sim_50_10000_10_100[[1]],
unCorCutoff=sim_unCorCutoff,
coreNum=50)

sim_UBHDD_50_10000_10_100_perfMat<-sim_UBHDD_50_10000_10_100[[1]]
sim_UBHDD_50_10000_10_100_coefMat<-sim_UBHDD_50_10000_10_100[[2]]
sim_UBHDD_50_10000_10_100_testIndex<-sim_UBHDD_50_10000_10_100[[3]]

save.image("data/sim_UBHDD_50_10000_10_100.RData")



### N1=100

## parameter

sim_PG_dim1=100
sim_PNG_dim1=10000
sim_PG_dim_per_trait1=sim_PG_dim1/2
sim_cluster1=rep(10,100)

## get simulated phenotype space

sim_100_10000_10_100<-get_sim_trait(
sim_PG_dim=sim_PG_dim1,
sim_PNG_dim=sim_PNG_dim1,
sim_PG_dim_per_trait=sim_PG_dim_per_trait1,
sim_cluster=sim_cluster1)

## set uncorrelation cutoff

sim_unCorCutoff<-calUnCorCutoff(
samSize=1000,
traitNum=sum(sim_cluster1),
alpha=0.01)

## UBHDD modelling

sim_UBHDD_100_10000_10_100<-UBHDD(
X=sim_100_10000_10_100[[1]],
Y=sim_100_10000_10_100[[1]],
unCorCutoff=sim_unCorCutoff,
coreNum=50)

sim_UBHDD_100_10000_10_100_perfMat<-sim_UBHDD_100_10000_10_100[[1]]
sim_UBHDD_100_10000_10_100_coefMat<-sim_UBHDD_100_10000_10_100[[2]]
sim_UBHDD_100_10000_10_100_testIndex<-sim_UBHDD_100_10000_10_100[[3]]

save.image("data/sim_UBHDD_100_10000_10_100.RData")



### N1=100 ; sim_cluster1=rep(10,200)

## parameter

sim_PG_dim1=100
sim_PNG_dim1=10000
sim_PG_dim_per_trait1=sim_PG_dim1/2
sim_cluster1=rep(10,200)

## get simulated phenotype space

sim_100_10000_10_200<-get_sim_trait(
sim_PG_dim=sim_PG_dim1,
sim_PNG_dim=sim_PNG_dim1,
sim_PG_dim_per_trait=sim_PG_dim_per_trait1,
sim_cluster=sim_cluster1)

## set uncorrelation cutoff

sim_unCorCutoff<-calUnCorCutoff(
samSize=1000,
traitNum=sum(sim_cluster1),
alpha=0.01)

## UBHDD modelling

sim_UBHDD_100_10000_10_200<-UBHDD(
X=sim_100_10000_10_200[[1]],
Y=sim_100_10000_10_200[[1]],
unCorCutoff=sim_unCorCutoff,
coreNum=50)

sim_UBHDD_100_10000_10_200_perfMat<-sim_UBHDD_100_10000_10_200[[1]]
sim_UBHDD_100_10000_10_200_coefMat<-sim_UBHDD_100_10000_10_200[[2]]
sim_UBHDD_100_10000_10_200_testIndex<-sim_UBHDD_100_10000_10_200[[3]]

save.image("data/sim_UBHDD_100_10000_10_200.RData")



### N1=100 ; sim_cluster1=rep(10,200); samSize=2000

#set.seed(1)
#mvnormPre<-mvrnorm(2000,mu=rep(0,11000),Sigma=diag(11000),tol=1)

## parameter

sim_PG_dim1=100
sim_PNG_dim1=10000
sim_PG_dim_per_trait1=sim_PG_dim1/2
sim_cluster1=rep(10,200)
sim_sam_size1=2000
## get simulated phenotype space

sim_100_10000_10_200_samSize2000<-get_sim_trait(
sim_PG_dim=sim_PG_dim1,
sim_PNG_dim=sim_PNG_dim1,
sim_PG_dim_per_trait=sim_PG_dim_per_trait1,
sim_cluster=sim_cluster1,
sim_sam_size=sim_sam_size1)

## set uncorrelation cutoff

sim_unCorCutoff<-calUnCorCutoff(
samSize=sim_sam_size1,
traitNum=sum(sim_cluster1),
alpha=0.01)

## UBHDD modelling

sim_UBHDD_100_10000_10_200_samSize2000<-UBHDD(
X=sim_100_10000_10_200_samSize2000[[1]],
Y=sim_100_10000_10_200_samSize2000[[1]],
unCorCutoff=sim_unCorCutoff,
coreNum=100)

sim_UBHDD_100_10000_10_200_samSize2000_perfMat<-sim_UBHDD_100_10000_10_200_samSize2000[[1]]
sim_UBHDD_100_10000_10_200_samSize2000_coefMat<-sim_UBHDD_100_10000_10_200_samSize2000[[2]]
sim_UBHDD_100_10000_10_200_samSize2000_testIndex<-sim_UBHDD_100_10000_10_200_samSize2000[[3]]

save.image("data/sim_UBHDD_100_10000_10_200_samSize2000.RData")








