
# 02/13/22
# R script in:
# The trait coding rule in phenotype space
# 
# Jianguo Wang (wangjg36@mail.sysu.edu.cn)
# 


############################################ library ####################################

library(bestNormalize)
library(parallel)
library(foreach)
library(doParallel)
library(apcluster)




############################################ functions #################################################

#format transformation

factor2matrix<-function(X){
    X.m<-foreach(i=1:dim(X)[2],.combine=cbind)%do%{
	    focalX<-X[,i]
		focalX.table<-table(focalX)
		focalX.table.name<-names(table(focalX))
		focalX.m.temp<-matrix(0,nrow=length(focalX),ncol=length(focalX.table))
	
		focalX.m<-foreach(j=1:length(focalX.table),.combine=cbind)%do%{focalX.m.temp[which(focalX==focalX.table.name[j]),j]=1
			return(focalX.m.temp[,j])}
		rownames(focalX.m)<-rownames(X)
	    colnames(focalX.m)<-paste0(colnames(X)[i],"_",focalX.table.name)			
		return(focalX.m)
	}
	
	return(X.m)

}



### normalTransform transform a continuous variable into a normal distributed variable at the basis of function bestNormalize in R package "bestNormalize"

### however, bestNormalize fails for normally distributed variable;

### in addition, the raw trait is selected when the raw trait has a larger correlation with the expected normal distribution generated by rnorm()

### parallel computation is perfered

### the missing values is not dealed with because of potential different study aims



normalTransform<-function(X,coreNum){
   
   
   cl <- makeCluster(coreNum,type="FORK")
   
   #temp<-parLapply(cl,1:dim(X)[2],function(i){
   #    return(bestNormalize(X[,i])$x.t)
   
   #})
   
   #X_1<-X
   #for(i in 1:dim(X)[2]){
   #    X_1[,i]<-temp[[i]]
   #}
   
   #windows
   #registerDoSNOW(cl) 
   #linux
   registerDoParallel(cl)
    
   #errorhandling to be considered
   #X_1 <- foreach(i = 1:dim(X)[2],.combine=cbind,.errorhandling="pass",.options.snow = "dynamic")%dopar%(bestNormalize(X[,i])$x.t)
   #windows
   #X_1 <- foreach(i = 1:dim(X)[2],.combine=cbind,.options.snow = list(preschedule=TRUE))%dopar%(bestNormalize(X[,i])$x.t)
   #linux
   #X_1 <- foreach(i = 1:dim(X)[2],.combine=cbind,.options.multicore=list(preschedule=TRUE))%dopar%(bestNormalize(X[,i])$x.t)
   X_1 <- foreach(i = 1:dim(X)[2],.combine=cbind)%dopar%(bestNormalize(X[,i])$x.t)
   
    
   stopCluster(cl)
   #stopImplicitCluster()
   
   rownames(X_1)<-rownames(X)
   colnames(X_1)<-colnames(X)
   
   ##sd
   
   X_1_sd<-apply(X_1,2,sd,na.rm=T)

   X_1[,which(X_1_sd<0.99|X_1_sd>1)]<-scale(X[,which(X_1_sd<0.99|X_1_sd>1)])
   
   ##cor with normal 

   set.seed(100)
   tempNorm<-matrix(rnorm(dim(X_1)[1]*dim(X_1)[2]),nrow=dim(X_1)[1],ncol=dim(X_1)[2])
   tempNorm[which(is.na(X_1))]<-NA

   normalCor_X_1<-sapply(1:dim(X_1)[2],function(x){cor(sort(na.omit(X_1[,x])),sort(na.omit(tempNorm[,x])))})
   normalCor_X<-sapply(1:dim(X_1)[2],function(x){cor(sort(na.omit(X[,x])),sort(na.omit(tempNorm[,x])))})

   X_1[,which(normalCor_X>normalCor_X_1)]<-scale(X[,which(normalCor_X>normalCor_X_1)])

   
   return(X_1)
   
}


############################### UKB phenome preprocessing ######################################

setwd("/home/wjg/UBHDD/UKB/") 
 
phenoPath<-"/storage/resource/database/UKBiobank_20200304"

#this is the raw phenome data from UK Biobank
uk<-read.csv(paste0(phenoPath,"/ukb40795_1.csv"))
traitAnnotAll<-read.csv("data/basic/traitAnnotation_all18108.csv")

#this is the kinship of 
#genealogy<-read.table("data/genomePreprocessing/sample_MAF_INFO/ukb54148_rel_s488239.dat",header=T)
load("data/basic/genealogy.RData")


##traitName

traitName<-sapply(names(uk)[-1],function(x){unlist(strsplit(x,"X"))[2]})

traitNameField<-sapply(traitName,function(x){unlist(strsplit(x,"\\."))[1]})

write.csv(data.frame(traitName,traitNameField),"data/basic/ukTraitNameField.csv")


##Ethnicity  

British<-uk[which(uk$X21000.0.0==1001|uk$X21000.1.0==1001|uk$X21000.2.0==1001),1]

write.csv(British,"data/basic/British.csv")


##uk_British

uk_British<-uk[match(British,uk[,1]),]

rm(uk)


##brain_trait related subjects


brainRelatedTraits<-traitAnnotAll[which(as.character(traitAnnotAll[,17])=="YES"),1]

uk_British_brain<-uk_British[,c(1,match(brainRelatedTraits,traitName)+1)]

save(brainRelatedTraits,file="data/basic/brainRelatedTraits.RData")

##genealogy


###select kinship>0;
genealogy<-genealogy[which(genealogy$Kinship>0),]

###keep those pairs both members of which are in the uk_British 
genealogy_within<-genealogy[intersect(which(!is.na(match(genealogy[,1],uk_British[,1]))),which(!is.na(match(genealogy[,2],uk_British[,1])))),]

###subject with more relatives will be first deleted until no subjects with more than one relatives and all deleted subjects are recorded in the variable "deleteIDRecord"
IDcomb<-c(genealogy_within[,1],genealogy_within[,2])
IDcombCount<-table(IDcomb)

deleteIDRecord<-c()

count=1

deleteIDInd<-which(IDcombCount==max(IDcombCount))[1]
deleteID<-names(IDcombCount)[deleteIDInd]
deleteIDRecord[count]<-deleteID

deletePairID<-c(genealogy_within[which(genealogy_within[,2]==deleteID),1],genealogy_within[which(genealogy_within[,1]==deleteID),2])
deleteIDPairInd<-match(deletePairID,names(IDcombCount))

IDcombCount_new<-IDcombCount
IDcombCount_new[deleteIDInd]<-0
IDcombCount_new[deleteIDPairInd]<-IDcombCount_new[deleteIDPairInd]-1

IDcombCount<-IDcombCount_new


while(max(IDcombCount)>1){
    
	count=count+1
	print(count)
	
	deleteIDInd<-which(IDcombCount==max(IDcombCount))[1]
    deleteID<-names(IDcombCount)[deleteIDInd]
	deleteIDRecord[count]<-deleteID
		
    deletePairID<-setdiff(c(genealogy_within[which(genealogy_within[,2]==deleteID),1],genealogy_within[which(genealogy_within[,1]==deleteID),2]),deleteIDRecord)
    deleteIDPairInd<-match(deletePairID,names(IDcombCount))
	
    IDcombCount_new<-IDcombCount
    IDcombCount_new[deleteIDInd]<-0
    IDcombCount_new[deleteIDPairInd]<-IDcombCount_new[deleteIDPairInd]-1
    
	IDcombCount<-IDcombCount_new
	
}

deleteIDRecord<-as.numeric(deleteIDRecord)


####

genealogyRemain<-genealogy_within[-union(which(!is.na(match(genealogy_within[,1],deleteIDRecord))),which(!is.na(match(genealogy_within[,2],deleteIDRecord)))),] 

###the number of NA in each subject in uk_British

subjectNonNACount<-dim(uk_British_brain)[2]-apply(is.na(uk_British_brain),1,sum)


#subjectNonNACount<-c()
#for(i in 1:dim(uk_British_brain)[1]){
#    print(i)
#	subjectNonNACount[i]<-length(na.omit(uk_British_brain[i,]))
#}




###in genealogyRemain, those with more NAs will be first deleted

deleteIDRecord2<-c()

pairNonNACount<-cbind(subjectNonNACount[match(genealogyRemain[,1],uk_British[,1])],subjectNonNACount[match(genealogyRemain[,2],uk_British[,1])])

for(i in 1:dim(genealogyRemain)[1]){
    deleteIDRecord2[i]<-genealogyRemain[i,(pairNonNACount[i,1]>pairNonNACount[i,2])+1]

}


###combine two parts and obtain the total subjects to be deleted based on genealogy
deleteIDRecordComb<-c(deleteIDRecord,deleteIDRecord2)

write.csv(deleteIDRecordComb,"data/basic/deleteIDRecordComb.csv")


###

uk_British_genealogy<-uk_British[-match(deleteIDRecordComb,uk_British[,1]),]

rm(uk_British)


######################################################

### select subjects with brain phenotype > 500

uk_British_genealogy_brain<-uk_British_genealogy[,c(1,match(brainRelatedTraits,traitName)+1)]


subjectNonNACount2<-dim(uk_British_genealogy_brain)[2]-apply(is.na(uk_British_genealogy_brain),1,sum)

write.csv(subjectNonNACount2,"data/basic/uk_British_genealogy_brain_subjectNonNACount2.csv")


#
uk_British_genealogy_brain500<-uk_British_genealogy[which(subjectNonNACount2>500),]

traitNonNACount<-dim(uk_British_genealogy_brain500)[1]-apply(is.na(uk_British_genealogy_brain500),2,sum)

write.csv(traitNonNACount,"data/basic/uk_British_genealogy_brain500_traitNonNACount.csv")

#
subjectNonNACount3<-dim(uk_British_genealogy_brain500)[2]-apply(is.na(uk_British_genealogy_brain500),1,sum)

write.csv(subjectNonNACount3,"data/basic/uk_British_genealogy_brain_subjectNonNACount3.csv")

rm(uk_British_brain)
rm(subjectNonNACount)
rm(genealogy)
rm(genealogyRemain)
rm(genealogy_within)

rm(deleteIDRecord)
rm(deleteIDRecord2)
rm(i)
rm(pairNonNACount)

rm(British)
rm(deleteIDRecordComb)

rm(uk_British_genealogy)


save.image("data/basic/uk_British_genealogy_brain500.RData")


###

###############################################################################




load("data/basic/uk_British_genealogy_brain500.RData")
load("data/basic/brainRelatedTraits.RData")

traitAnnotAll<-read.csv("data/basic/traitAnnotation_all18108.csv")


######################### trait and subjects choice ######################

#uk_British_genealogy_brain
###注意数据类型是否为数值型还是类别型

###define brainPhneo, brainCovariate and totalCovariate

load("data/basic/brainPhenoName_b2o.RData")
load("data/basic/brainCovariateName_b2o.RData")
load("data/basic/totalCovariateName_b2o.RData")

##totalCovariateName_brain2other<-c(sexCovariate,ageCovariate,weightCovariate,locationCovariate,centreCovariate)

###define subjects with complete observations for each brainPheno

temp<-uk_British_genealogy_brain500[,c(brainPhenoName_brain2other,totalCovariateName_brain2other)]
brainNACount<-apply(is.na(temp),1,sum)

subjectSelect_brain2other<-uk_British_genealogy_brain500[which(brainNACount==0),1]


###obtain brainPheno, brainCovariate and totalCovariate and otherPheno trait matrix (raw trait)

uk_British_genealogy_brain_brain2other<-uk_British_genealogy_brain500[match(subjectSelect_brain2other,uk_British_genealogy_brain500[,1]),brainPhenoName_brain2other]
rownames(uk_British_genealogy_brain_brain2other)<-subjectSelect_brain2other
uk_British_genealogy_brain_brain2other<-as.matrix(uk_British_genealogy_brain_brain2other)

uk_British_genealogy_brain_brain2other_brainCovariate<-uk_British_genealogy_brain500[match(subjectSelect_brain2other,uk_British_genealogy_brain500[,1]),brainCovariateName_brain2other]
rownames(uk_British_genealogy_brain_brain2other_brainCovariate)<-subjectSelect_brain2other
uk_British_genealogy_brain_brain2other_brainCovariate<-as.matrix(uk_British_genealogy_brain_brain2other_brainCovariate)

uk_British_genealogy_brain_brain2other_totalCovariate<-uk_British_genealogy_brain500[match(subjectSelect_brain2other,uk_British_genealogy_brain500[,1]),totalCovariateName_brain2other]
rownames(uk_British_genealogy_brain_brain2other_totalCovariate)<-subjectSelect_brain2other
uk_British_genealogy_brain_brain2other_totalCovariate<-as.matrix(uk_British_genealogy_brain_brain2other_totalCovariate)

uk_British_genealogy_brain_brain2other_totalCovariate<-cbind(uk_British_genealogy_brain_brain2other_totalCovariate[,1:17],factor2matrix(uk_British_genealogy_brain_brain2other_totalCovariate[,18:19]))


############################## imputation for brainCovariate ########################

###imputation for brainCovariate

which(is.na(uk_British_genealogy_brain_brain2other))
which(is.na(uk_British_genealogy_brain_brain2other_totalCovariate))

which(is.na(uk_British_genealogy_brain_brain2other_brainCovariate))

length(which(is.na(uk_British_genealogy_brain_brain2other_brainCovariate)))/dim(uk_British_genealogy_brain_brain2other_brainCovariate)[1]/dim(uk_British_genealogy_brain_brain2other_brainCovariate)[2]

brainCovariateNonNACount_brain2other<-apply(uk_British_genealogy_brain_brain2other_brainCovariate,2,function(x){length(which(!is.na(x)))})

table(brainCovariateNonNACount_brain2other)


##step 1: get imputation model

cl<-makeCluster(dim(uk_British_genealogy_brain_brain2other_brainCovariate)[2],type="FORK")

brainCovariate_imp_model_brain2other<-parLapply(cl,1:dim(uk_British_genealogy_brain_brain2other_brainCovariate)[2],function(i){
    if(brainCovariateNonNACount_brain2other[i]==length(subjectSelect_brain2other)){
	    return(NA)
	}else{
	    focalSamInd<-which(!is.na(uk_British_genealogy_brain_brain2other_brainCovariate[,i]))
	    focalX<-uk_British_genealogy_brain_brain2other[focalSamInd,]
		focalY<-uk_British_genealogy_brain_brain2other_brainCovariate[focalSamInd,i]
		
		set.seed(100)
	    return(cv.glmnet(focalX,focalY,alpha=1,nlambda=200))
		
	}

})

stopCluster(cl)




##step 2: get imputed brainCovariate

uk_British_genealogy_brain_brain2other_brainCovariate_imp<-uk_British_genealogy_brain_brain2other_brainCovariate



for(i in 1:dim(uk_British_genealogy_brain_brain2other_brainCovariate)[2]){
    print(i)
    if(brainCovariateNonNACount_brain2other[i]!=length(subjectSelect_brain2other)){

	    focalSamInd<-which(is.na(uk_British_genealogy_brain_brain2other_brainCovariate[,i]))
	    focalX<-uk_British_genealogy_brain_brain2other[focalSamInd,]
		if(length(focalSamInd)>1){
		    uk_British_genealogy_brain_brain2other_brainCovariate_imp[focalSamInd,i]<-as.numeric(cbind(rep(1,dim(focalX)[1]),focalX)%*%coef(brainCovariate_imp_model_brain2other[[i]],s="lambda.min")[,1])
		}
		
		if(length(focalSamInd)==1){
		    uk_British_genealogy_brain_brain2other_brainCovariate_imp[focalSamInd,i]<-as.numeric(t(c(1,focalX))%*%coef(brainCovariate_imp_model_brain2other[[i]],s="lambda.min")[,1])
		    
		}		

	}

}


which(is.na(uk_British_genealogy_brain_brain2other_brainCovariate_imp))

#################### save data #########################

save(subjectSelect_brain2other,file="data/basic/subjectSelect_b2o.RData")

save(uk_British_genealogy_brain_brain2other,file="data/basic/brainPheno_b2o.RData")

save(uk_British_genealogy_brain_brain2other_brainCovariate,file="data/basic/brainCovariate_b2o.RData")
save(uk_British_genealogy_brain_brain2other_brainCovariate_imp,file="data/basic/brainCovariate_imp_b2o.RData")

save(uk_British_genealogy_brain_brain2other_totalCovariate,file="data/basic/totalCovariate_b2o.RData")



################## nomal transformation for brainPheno & otherPheno ("R package : bestNormalize" adjusted)##############################

##uk_British_genealogy_brain_brain2other_bestNorm

rm(list=ls())

load("data/basic/brainPheno_b2o.RData")

uk_British_genealogy_brain_brain2other_bestNorm<-normalTransform(uk_British_genealogy_brain_brain2other,55)

date()

quantile(sapply(1:dim(uk_British_genealogy_brain_brain2other)[2],function(x){cor(uk_British_genealogy_brain_brain2other[,x],uk_British_genealogy_brain_brain2other_bestNorm[,x],use="pairwise.complete.obs")}))

save(uk_British_genealogy_brain_brain2other_bestNorm,file="data/basic/brainPheno_bestNorm_b2o.RData")


##############################################################################################################



load("data/basic/brainPheno_bestNorm_b2o.RData")
load("data/basic/brainCovariate_imp_b2o.RData")
load("data/basic/totalCovariate_b2o.RData")


subject_with_genome<-read.table("data/genomePreprocessing/genome_plink_maf_indep_geno_mind_hwe_pruned/ukb_imp_chr1_v3_plink_maf_indep_geno_mind_hwe_pruned.fam")

brain_with_genome<-uk_British_genealogy_brain_brain2other_bestNorm[match(subject_with_genome[,2],rownames(uk_British_genealogy_brain_brain2other_bestNorm)),]
brainCovariate_with_genome<-uk_British_genealogy_brain_brain2other_brainCovariate_imp[match(subject_with_genome[,2],rownames(uk_British_genealogy_brain_brain2other_bestNorm)),]
totalCovariate_with_genome<-uk_British_genealogy_brain_brain2other_totalCovariate[match(subject_with_genome[,2],rownames(uk_British_genealogy_brain_brain2other_bestNorm)),]

###################################



cl<-makeCluster(20,type="FORK")
registerDoParallel(cl)

brain_with_genome_exclude_covariate<-foreach(i=1:dim(brain_with_genome)[2],.combine=cbind)%dopar%{print(i);residuals(lm(brain_with_genome[,i]~cbind(brainCovariate_with_genome,totalCovariate_with_genome)))}

stopCluster(cl)

colnames(brain_with_genome_exclude_covariate)<-colnames(uk_British_genealogy_brain_brain2other_bestNorm)



############### dMRI_trait (water)


result_b2b_water_015_DF<-read.csv("data/b2b/result_b2b_water_015_DF.csv")

load("data/basic/dMRI_trait.RData")


brain_with_genome_exclude_covariate_water_T<-brain_with_genome_exclude_covariate[,match(dMRI_trait,colnames(brain_with_genome_exclude_covariate))]
brain_with_genome_exclude_covariate_water_Tg<-cbind(rep(1,25957),as.matrix(brain_with_genome_exclude_covariate_water_T))%*%t(result_b2b_water_015_DF[,4:dim(result_b2b_water_015_DF)[2]])
colnames(brain_with_genome_exclude_covariate_water_Tg)<-colnames(brain_with_genome_exclude_covariate_water_T)
brain_with_genome_exclude_covariate_water_Tng<-brain_with_genome_exclude_covariate_water_T-brain_with_genome_exclude_covariate_water_Tg


save(brain_with_genome_exclude_covariate_water_T,file="data/b2b/brain_with_genome_exclude_covariate_water_T.RData")
save(brain_with_genome_exclude_covariate_water_Tg,file="data/b2b/brain_with_genome_exclude_covariate_water_Tg.RData")
save(brain_with_genome_exclude_covariate_water_Tng,file="data/b2b/brain_with_genome_exclude_covariate_water_Tng.RData")

write.table(cbind(subject_with_genome[,1:2],brain_with_genome_exclude_covariate_water_T),"data/b2b/brain_with_genome_exclude_covariate_water_T_2GCTA.txt",col.names=F,row.names=F,quote=F)
write.table(cbind(subject_with_genome[,1:2],brain_with_genome_exclude_covariate_water_Tg),"data/b2b/brain_with_genome_exclude_covariate_water_Tg_2GCTA.txt",col.names=F,row.names=F,quote=F)
write.table(cbind(subject_with_genome[,1:2],brain_with_genome_exclude_covariate_water_Tng),"data/b2b/brain_with_genome_exclude_covariate_water_Tng_2GCTA.txt",col.names=F,row.names=F,quote=F)



############# exemplar trait


brain_with_genome_exclude_covariate_water_T_corMat<-cor(brain_with_genome_exclude_covariate_water_T)
save(brain_with_genome_exclude_covariate_water_T_corMat,file="brain_with_genome_exclude_covariate_water_T_corMat.RData")

#
temp<-apcluster(s=abs(brain_with_genome_exclude_covariate_water_T_corMat))
clusterName<-temp@exemplars
save(clusterName,file="data/basic/clusterName.RData")

#

load("data/basic/clusterName.RData")

brain_with_genome_exclude_covariate_water_exemplar_T<-brain_with_genome_exclude_covariate_water_T[,clusterName]
brain_with_genome_exclude_covariate_water_exemplar_Tg<-brain_with_genome_exclude_covariate_water_Tg[,clusterName]
brain_with_genome_exclude_covariate_water_exemplar_Tng<-brain_with_genome_exclude_covariate_water_Tng[,clusterName]

save(brain_with_genome_exclude_covariate_water_exemplar_T,file="data/b2b/brain_with_genome_exclude_covariate_water_exemplar_T.RData")
save(brain_with_genome_exclude_covariate_water_exemplar_Tg,file="data/b2b/brain_with_genome_exclude_covariate_water_exemplar_Tg.RData")
save(brain_with_genome_exclude_covariate_water_exemplar_Tng,file="data/b2b/brain_with_genome_exclude_covariate_water_exemplar_Tng.RData")

write.table(cbind(subject_with_genome[,1:2],brain_with_genome_exclude_covariate_water_exemplar_T),"data/b2b/brain_with_genome_exclude_covariate_water_exemplar_T_2GCTA.txt",col.names=F,row.names=F,quote=F)
write.table(cbind(subject_with_genome[,1:2],brain_with_genome_exclude_covariate_water_exemplar_Tg),"data/b2b/brain_with_genome_exclude_covariate_water_exemplar_Tg_2GCTA.txt",col.names=F,row.names=F,quote=F)
write.table(cbind(subject_with_genome[,1:2],brain_with_genome_exclude_covariate_water_exemplar_Tng),"data/b2b/brain_with_genome_exclude_covariate_water_exemplar_Tng_2GCTA.txt",col.names=F,row.names=F,quote=F)


################################################## UKB genome preprocessing ############################################

setwd("/home/wjg/UBHDD/UKB/data") 


softwarePath<-"/home/wjg/software"
genomePath<-"/data/resource/project/UKBioBankExomeSeqAndRetina.20210306/ukb/"

###qctool_filter

num=1
chr<-read.table(paste0("genomePreprocessing/sample_MAF_INFO/ukb_imp_mfi/ukb_mfi_chr",num,"_v3.txt"))
chr<-na.omit(chr)

#snp filtering: MAF>0.01 & imputation score > 0.7 & non-multiple allele
SNPID<-chr[which(!(duplicated(chr[,1])|duplicated(chr[,2])|duplicated(chr[,3]))&chr[,6]>0.01&chr[,8]>0.7),1]
length(SNPID)
write.table(t(SNPID),paste0("genomePreprocessing/ukb_imp_mfi_qctool_filter/ukb_mfi_chr",num,"_v3_qctool_filter.txt"),col.names=F,row.names=F,quote=F,eol = "")


### Threshhold genotype call probabilities = 0.9; 

qctool_filter<-c()


for(i in 1:22){
    qctool_filter[i]<-paste0(softwarePath,"/qctool_v2.0.6-Ubuntu16.04-x86_64/qctool -g ",genomePath,"ukb_imp_chr",i,"_v3.bgen -s genomePreprocessing/sample_MAF_INFO/chr_sample/ukb22828_c",i,"_b0_v3_s487275.sample -ofiletype binary_ped -threshold 0.9 -og genomePreprocessing/genome_qctool_filter/ukb_imp_chr",i,"_v3_qctool_filter -incl-snpids genomePreprocessing/ukb_imp_mfi_qctool_filter/ukb_mfi_chr",i,"_v3_qctool_filter.txt")
}


write.table(qctool_filter,"genomePreprocessing/qctool_filter.txt",col.names=F,row.names=F,quote=F)



### plink_filter

subjects_b2o<-read.csv("b2o/subjects_b2o.csv")

#

for(i in 1:22){

	chr_i<-read.table(paste0("genomePreprocessing/genome_qctool_filter/ukb_imp_chr",i,"_v3_qctool_filter.fam"))
	write.table(chr_i[na.omit(match(subject_b2o[,2],chr_i[,2])),1:2],paste0("genomePreprocessing/ukb_imp_fam_plink_filter/ukb_imp_chr",i,"_v3_fam_plink_filter.txt"),col.names=F,row.names=F,quote=F)

}


###

plink_filter<-c()

for(i in 1:22){
	plink_filter[i]<-paste0(softwarePath,"/plink_linux_x86_64_20210606/plink --bfile genomePreprocessing/genome_qctool_filter/ukb_imp_chr",i,"_v3_qctool_filter --make-bed --out genomePreprocessing/genome_plink_filter/ukb_imp_chr",i,"_v3_plink_filter --keep genomePreprocessing/ukb_imp_fam_plink_filter/ukb_imp_chr",i,"_v3_fam_plink_filter.txt --threads 80")
	
}

write.table(plink_filter,"genomePreprocessing/plink_filter.txt",col.names=F,row.names=F,quote=F)


###

genome_plink_maf_indep_geno_mind_hwe<-c()

for(i in 1:22){
	genome_plink_maf_indep_geno_mind_hwe_filter[i]<-paste0(softwarePath,"/plink_linux_x86_64_20210606/plink --bfile genomePreprocessing/genome_plink_filter/ukb_imp_chr",i,"_v3_plink_filter --maf 0.01 --geno 0.1 --mind 0.1 --hwe 1e-9 --indep 50 5 5 --make-bed --out genomePreprocessing/genome_plink_maf_indep_geno_mind_hwe/ukb_imp_chr",i,"_v3_plink_maf_indep_geno_mind_hwe --threads 30")	
}

write.table(genome_plink_maf_indep_geno_mind_hwe,"genomePreprocessing/genome_plink_maf_indep_geno_mind_hwe.txt",col.names=F,row.names=F,quote=F)


###

genome_plink_maf_indep_geno_mind_hwe_pruned<-c()

for(i in 1:22){
	genome_plink_maf_indep_geno_mind_hwe_pruned[i]<-paste0(softwarePath,"/plink_linux_x86_64_20210606/plink --bfile genomePreprocessing/genome_plink_maf_indep_geno_mind_hwe/ukb_imp_chr",i,"_v3_plink_maf_indep_geno_mind_hwe --extract genomePreprocessing/genome_plink_maf_indep_geno_mind_hwe/ukb_imp_chr",i,"_v3_plink_maf_indep_geno_mind_hwe.prune.in --make-bed --out genomePreprocessing/genome_plink_maf_indep_geno_mind_hwe_pruned/ukb_imp_chr",i,"_v3_plink_maf_indep_geno_mind_hwe_pruned")
	
}

write.table(genome_plink_maf_indep_geno_mind_hwe_pruned,"genomePreprocessing/genome_plink_maf_indep_geno_mind_hwe_pruned.txt",col.names=F,row.names=F,quote=F)


