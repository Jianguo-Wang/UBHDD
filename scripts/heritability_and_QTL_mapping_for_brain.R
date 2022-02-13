
# 02/13/22
# R script in:
# The trait coding rule in phenotype space
# 
# Jianguo Wang (wangjg36@mail.sysu.edu.cn)
# 




################################# narrow-sense heritability ####################################################

setwd("/home/wjg/UBHDD/UKB/data") 


softwarePath<-"/home/wjg/software"

mpheno_Vg_GREML_water_exempar_T<-c()
for(i in 1:107){
	mpheno_Vg_GREML_water_exempar_T[i]<-paste0(softwarePath,"/gcta_1.93.2beta/gcta64 --grm genomePreprocessing/GRM/GRM --reml --reml-maxit 100 --pheno b2b/brain_with_genome_exclude_covariate_water_exemplar_T_2GCTA.txt --mpheno ",i," --out genomePreprocessing/Vg_GREML_water_exemplar/Vg_GREML_water_exemplar_T_",i," --thread-num 30")
}

mpheno_Vg_GREML_water_exempar_Tg<-c()
for(i in 1:107){
	mpheno_Vg_GREML_water_exempar_Tg[i]<-paste0(softwarePath,"/gcta_1.93.2beta/gcta64 --grm genomePreprocessing/GRM/GRM --reml --reml-maxit 100 --pheno b2b/brain_with_genome_exclude_covariate_water_exemplar_Tg_2GCTA.txt --mpheno ",i," --out genomePreprocessing/Vg_GREML_water_exemplar/Vg_GREML_water_exemplar_Tg_",i," --thread-num 30")
}

mpheno_Vg_GREML_water_exempar_Tng<-c()
for(i in 1:107){
	mpheno_Vg_GREML_water_exempar_Tng[i]<-paste0(softwarePath,"/gcta_1.93.2beta/gcta64 --grm genomePreprocessing/GRM/GRM --reml --reml-maxit 100 --pheno b2b/brain_with_genome_exclude_covariate_water_exemplar_Tng_2GCTA.txt --mpheno ",i," --out genomePreprocessing/Vg_GREML_water_exemplar/Vg_GREML_water_exemplar_Tng_",i," --thread-num 30")
}

write.table(mpheno_Vg_GREML_water_exempar_T,"genomePreprocessing/mpheno_Vg_GREML_water_exempar_T.txt",col.names=F,row.names=F,quote=F)
write.table(mpheno_Vg_GREML_water_exempar_Tg,"genomePreprocessing/mpheno_Vg_GREML_water_exempar_Tg.txt",col.names=F,row.names=F,quote=F)
write.table(mpheno_Vg_GREML_water_exempar_Tng,"genomePreprocessing/mpheno_Vg_GREML_water_exempar_Tng.txt",col.names=F,row.names=F,quote=F)



######################################### QTL mapping ########################################################




### GRM by GCTA

ukb_imp_v3_plink_maf_indep_geno_mind_hwe_pruned<-c()
for(i in 1:22){
	ukb_imp_v3_plink_maf_indep_geno_mind_hwe_pruned[i]<-paste0("ukb_imp_chr",i,"_v3_plink_maf_indep_geno_mind_hwe_pruned")
}

write.table(ukb_imp_v3_plink_maf_indep_geno_mind_hwe_pruned,"genomePreprocessing/ukb_imp_v3_plink_maf_indep_geno_mind_hwe_pruned.txt",col.names=F,row.names=F,quote=F)


###

/home/wjg/software/gcta_1.93.2beta/gcta64 --mbfile genomePreprocessing/ukb_imp_v3_plink_maf_indep_geno_mind_hwe_pruned.txt --make-grm --thread-num 80 --out genomePreprocessing/GRM/GRM

### sparse GRM

/home/wjg/software/gcta_1.93.2beta/gcta64 --grm genomePreprocessing/GRM/GRM --make-bK-sparse 0.05 --out genomePreprocessing/sp_GRM/sp_GRM





### QTL mapping



# fastGWA mixed model (based on the sparse GRM generated above)



## brain trait related to water


mpheno_fastGWA_mlm_water_Tg<-c()
for(i in 1:675){
    mpheno_fastGWA_mlm_water_Tg[i]<-paste0(softwarePath,"/gcta_1.93.2beta/gcta64 --mbfile genomePreprocessing/ukb_imp_v3_plink_maf_indep_geno_mind_hwe_pruned.txt --grm-sparse genomePreprocessing/sp_GRM/sp_GRM --fastGWA-mlm --est-vg REML --pheno b2b/brain_with_genome_exclude_covariate_water_Tg_2GCTA.txt --mpheno ",i," --out genomePreprocessing/geno_assoc_fastGWA_mlm_water/geno_assoc_fastGWA_mlm_water_Tg_",i," --seed 1 --threads 40")

}

write.table(mpheno_fastGWA_mlm_water_Tg,"genomePreprocessing/mpheno_fastGWA_mlm_water_Tg.txt",col.names=F,row.names=F,quote=F)



mpheno_fastGWA_mlm_water_Tng<-c()
for(i in 1:675){
    mpheno_fastGWA_mlm_water_Tng[i]<-paste0(softwarePath,"/gcta_1.93.2beta/gcta64 --mbfile genomePreprocessing/ukb_imp_v3_plink_maf_indep_geno_mind_hwe_pruned.txt --grm-sparse genomePreprocessing/sp_GRM/sp_GRM --fastGWA-mlm --est-vg REML --pheno b2b/brain_with_genome_exclude_covariate_water_Tng_2GCTA.txt --mpheno ",i," --out genomePreprocessing/geno_assoc_fastGWA_mlm_water/geno_assoc_fastGWA_mlm_water_Tng_",i," --seed 1 --threads 40")

}

write.table(mpheno_fastGWA_mlm_water_Tng,"genomePreprocessing/mpheno_fastGWA_mlm_water_Tng.txt",col.names=F,row.names=F,quote=F)


######################################### clumping analysis ##########################################################################


softwarePath<-"/home/wjg/software"


###


filenames_water <- dir("genomePreprocessing/geno_assoc_fastGWA_mlm_water")

Tg_fastGWA<-filenames_water[which(grepl("Tg",filenames_water)&!grepl("log",filenames_water))]

Tng_fastGWA<-filenames_water[which(grepl("Tng",filenames_water)&!grepl("log",filenames_water))]



load("b2b/clusterName.RData")
exemplarTrait<-names(clusterName)

load("b2b/brain_with_genome_exclude_covariate_water_T_corMat.RData")
exeInd_in_water<-match(exemplarTrait,colnames(brain_with_genome_exclude_covariate_water_T_corMat))


##

Tg_fastGWA_No<-sapply(Tg_fastGWA,function(x){as.numeric(unlist(strsplit(unlist(strsplit(x,"_"))[7],"\\."))[1])})
Tg_fastGWA_exemplar<-sapply(1:length(exeInd_in_water),function(x){ifelse(!is.na(match(exeInd_in_water[x],Tg_fastGWA_No)),Tg_fastGWA[match(exeInd_in_water[x],Tg_fastGWA_No)],NA)})
Tng_fastGWA_No<-sapply(Tng_fastGWA,function(x){as.numeric(unlist(strsplit(unlist(strsplit(x,"_"))[7],"\\."))[1])})
Tng_fastGWA_exemplar<-sapply(1:length(exeInd_in_water),function(x){ifelse(!is.na(match(exeInd_in_water[x],Tng_fastGWA_No)),Tng_fastGWA[match(exeInd_in_water[x],Tng_fastGWA_No)],NA)})


###


#

cl<-makeCluster(20,type="FORK")
registerDoParallel(cl)

Tg_fastGWA_exemplar_QTL_result<-foreach(i=1:107)%dopar%{
    if(is.na(Tg_fastGWA_exemplar[i])){
	    return(NA)
	}else{
        dataTemp<-read.table(paste0("genomePreprocessing/geno_assoc_fastGWA_mlm_water/",Tg_fastGWA_exemplar[i]),header=T)
		if(length(which(dataTemp$P<2.68e-8))==0){
		    return(NA)
		}else{
            return(dataTemp[dataTemp$P<2.68e-8,])	    		    
		}
	}    

}

stopCluster(cl)

#

cl<-makeCluster(20,type="FORK")
registerDoParallel(cl)

Tng_fastGWA_exemplar_QTL_result<-foreach(i=1:107)%dopar%{
    if(is.na(Tng_fastGWA_exemplar[i])){
	    return(NA)
	}else{
        dataTemp<-read.table(paste0("genomePreprocessing/geno_assoc_fastGWA_mlm_water/",Tng_fastGWA_exemplar[i]),header=T)
		if(length(which(dataTemp$P<2.68e-8))==0){
		    return(NA)
		}else{
            return(dataTemp[dataTemp$P<2.68e-8,])	    		    
		}	    
	}    

}

stopCluster(cl)

###

#

for(i in 1:length(Tg_fastGWA_exemplar_QTL_result)){
    if(!is.na(Tg_fastGWA_exemplar_QTL_result[[i]])){
	    write.table(Tg_fastGWA_exemplar_QTL_result[[i]],paste0("genomePreprocessing/clumping_analysis/QTL_mapping_result/",Tg_fastGWA_exemplar[i]),quote=F,col.names=T,row.names=F,sep="\t")
	}
}

#

for(i in 1:length(Tng_fastGWA_exemplar_QTL_result)){
    if(!is.na(Tng_fastGWA_exemplar_QTL_result[[i]])){
	    write.table(Tng_fastGWA_exemplar_QTL_result[[i]],paste0("genomePreprocessing/clumping_analysis/QTL_mapping_result/",Tng_fastGWA_exemplar[i]),quote=F,col.names=T,row.names=F,sep="\t")
	}
}



###

sum(sapply(Tg_fastGWA_exemplar_QTL_result,function(x){
    if(is.na(x)){
	    return(0)
	}else{
        return(length(table(x[duplicated(x[,1]),1])))
	}
}))

sum(sapply(Tng_fastGWA_exemplar_QTL_result,function(x){
    if(is.na(x)){
	    return(0)
	}else{
        return(length(table(x[duplicated(x[,1]),1])))
	}
}))


###


#

clump_obj_Tg<-c()
count=0
for(i in 1:length(Tg_fastGWA_exemplar)){
    if(!is.na(Tg_fastGWA_exemplar_QTL_result[[i]])){
	    if(length(which(duplicated(Tg_fastGWA_exemplar_QTL_result[[i]][,1])))>0){
		    focalChr<-as.numeric(names(table(Tg_fastGWA_exemplar_QTL_result[[i]][duplicated(Tg_fastGWA_exemplar_QTL_result[[i]][,1]),1])))
		    for(j in 1:length(focalChr)){
		        count=count+1
                clump_obj_Tg[count]<-paste0(softwarePath,"/plink-1.07-x86_64/plink --bfile genomePreprocessing/genome_plink_maf_indep_geno_mind_hwe_pruned/ukb_imp_chr",focalChr[j],"_v3_plink_maf_indep_geno_mind_hwe_pruned --clump genomePreprocessing/clumping_analysis/QTL_mapping_result/",Tg_fastGWA_exemplar[i]," --allow-no-sex --clump-p1 2.68e-8 --clump-p2 2.68e-8 --noweb --out genomePreprocessing/clumping_analysis/Tg/",Tg_fastGWA_exemplar[i],"_chr",focalChr[j])
		    
		    }	    
		}

	}

}


#

clump_obj_Tng<-c()

count=0
for(i in 1:length(Tng_fastGWA_exemplar)){
    if(!is.na(Tng_fastGWA_exemplar_QTL_result[[i]])){
	    if(length(which(duplicated(Tng_fastGWA_exemplar_QTL_result[[i]][,1])))>0){
		    focalChr<-as.numeric(names(table(Tng_fastGWA_exemplar_QTL_result[[i]][duplicated(Tng_fastGWA_exemplar_QTL_result[[i]][,1]),1])))
		    for(j in 1:length(focalChr)){
		        count=count+1
                clump_obj_Tng[count]<-paste0(softwarePath,"/plink-1.07-x86_64/plink --bfile genomePreprocessing/genome_plink_maf_indep_geno_mind_hwe_pruned/ukb_imp_chr",focalChr[j],"_v3_plink_maf_indep_geno_mind_hwe_pruned --clump genomePreprocessing/clumping_analysis/QTL_mapping_result/",Tng_fastGWA_exemplar[i]," --allow-no-sex --clump-p1 2.68e-8 --clump-p2 2.68e-8 --noweb --out genomePreprocessing/clumping_analysis/Tng/",Tng_fastGWA_exemplar[i],"_chr",focalChr[j])
		    
		    }	    
		}

	}

}


##

write.table(clump_obj_Tg,"genomePreprocessing/clump_obj_Tg.txt",quote=F,col.names=F,row.names=F)
write.table(clump_obj_Tng,"genomePreprocessing/clump_obj_Tng.txt",quote=F,col.names=F,row.names=F)


##################################### clumping QTL result



#
Tg_fastGWA_exemplar_QTL_result_clump<-Tg_fastGWA_exemplar_QTL_result

for(i in 1:length(Tg_fastGWA_exemplar)){
    if(!is.na(Tg_fastGWA_exemplar_QTL_result[[i]])){
	    if(length(which(duplicated(Tg_fastGWA_exemplar_QTL_result[[i]][,1])))>0){
		    focalChr<-as.numeric(names(table(Tg_fastGWA_exemplar_QTL_result[[i]][duplicated(Tg_fastGWA_exemplar_QTL_result[[i]][,1]),1])))
			tempList<-list()
		    for(j in 1:length(focalChr)){
                tempList[[j]]<-read.table(paste0("genomePreprocessing/clumping_analysis/Tg/",Tg_fastGWA_exemplar[i],"_chr",focalChr[j],".clumped"),header=T)
		    
		    }
            dataTemp<-na.omit(do.call("rbind",tempList))
			clumpInd<-match(dataTemp[,3],Tg_fastGWA_exemplar_QTL_result[[i]][,2])
			
			chr_undup<-setdiff(as.numeric(names(table(Tg_fastGWA_exemplar_QTL_result[[i]][,1]))),focalChr)			
			unClumpInd<-match(chr_undup,Tg_fastGWA_exemplar_QTL_result[[i]][,1])
            Tg_fastGWA_exemplar_QTL_result_clump[[i]]<-Tg_fastGWA_exemplar_QTL_result[[i]][sort(c(clumpInd,unClumpInd)),]
            			
		}

	}

}


#

Tng_fastGWA_exemplar_QTL_result_clump<-Tng_fastGWA_exemplar_QTL_result

for(i in 1:length(Tng_fastGWA_exemplar)){
    if(!is.na(Tng_fastGWA_exemplar_QTL_result[[i]])){
	    if(length(which(duplicated(Tng_fastGWA_exemplar_QTL_result[[i]][,1])))>0){
		    focalChr<-as.numeric(names(table(Tng_fastGWA_exemplar_QTL_result[[i]][duplicated(Tng_fastGWA_exemplar_QTL_result[[i]][,1]),1])))
			tempList<-list()
		    for(j in 1:length(focalChr)){
                tempList[[j]]<-read.table(paste0("genomePreprocessing/clumping_analysis/Tng/",Tng_fastGWA_exemplar[i],"_chr",focalChr[j],".clumped"),header=T)
		    
		    }
            dataTemp<-na.omit(do.call("rbind",tempList))
			clumpInd<-match(dataTemp[,3],Tng_fastGWA_exemplar_QTL_result[[i]][,2])
			
			chr_undup<-setdiff(as.numeric(names(table(Tng_fastGWA_exemplar_QTL_result[[i]][,1]))),focalChr)			
			unClumpInd<-match(chr_undup,Tng_fastGWA_exemplar_QTL_result[[i]][,1])
            Tng_fastGWA_exemplar_QTL_result_clump[[i]]<-Tng_fastGWA_exemplar_QTL_result[[i]][sort(c(clumpInd,unClumpInd)),]
            			
		}

	}

}

########################## QTL_num based on clumped QTL results

Tg_QTL_num_clump<-sapply(Tg_fastGWA_exemplar_QTL_result_clump,function(x){if(is.na(x)){return(0)}else{return(dim(x)[1])}})
Tng_QTL_num_clump<-sapply(Tng_fastGWA_exemplar_QTL_result_clump,function(x){if(is.na(x)){return(0)}else{return(dim(x)[1])}})


save(Tg_QTL_num_clump,Tng_QTL_num_clump,file="genomePreprocessing/clumping_analysis/QTL_mapping_result/QTL_num_Tg_Tng.RData")


######################## h2 of Tg, Tng


#cd /data/user/wjg/biobank/genomePreprocessing/Vg_GREML_water_exemplar

filename<-dir("genomePreprocessing/Vg_GREML_water_exemplar")

filename_hsq<-filename[which(grepl("hsq",filename))]
traitNo<-sapply(filename_hsq,function(x){as.numeric(unlist(strsplit(unlist(strsplit(x,"_"))[6],"\\."))[1])})
traitType<-sapply(filename_hsq,function(x){unlist(strsplit(x,"_"))[5]})

filename_hsq_path<-paste0("genomePreprocessing/Vg_GREML_water_exemplar/",filename_hsq)
h2<-foreach(i=1:length(filename_hsq),.combine=c)%do%{ as.numeric(unlist(strsplit(unlist(strsplit(readChar(con = filename_hsq_path[i], nchars = file.info(filename_hsq_path[i])$size, useByte=F),"\n"))[5],"\t"))[2]) }

result<-data.frame(traitNo,traitType,h2)
result<-result[order(result[,1],result[,2]),]
result<-cbind(c(1:dim(result)[1]),result)

h2_Tg<-result[c(1:107)*2-1,4]
h2_Tng<-result[c(1:107)*2,4]


save(h2_Tg,h2_Tng,file="genomePreprocessing/clumping_analysis/QTL_mapping_result/h2_Tg_Tng.RData")













