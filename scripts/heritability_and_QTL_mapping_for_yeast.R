
# 02/13/22
# R script in:
# The trait coding rule in phenotype space
# 
# Jianguo Wang (wangjg36@mail.sysu.edu.cn)
# 


###################### library ##########################################

setwd("E:/papers/UBHDD/yeastSegregants") 


######

library(qtl)
library(abind)
library(rrBLUP)
library(lme4)

library(foreach)
library(doParallel)
library(parallel)




################## H2 and h2 ##############################################


# location of QTL_mappingFx.R
source((paste("functions/QTL_mappingFx.R", sep="")))

# Loads 'pheno_raw' list (contains phenotype measurements for each trait)
load("data/basic/cross.RData")
load("data/basic/pheno_raw_405traits.RData")

##basic
gdata     = extractGenotype(cross)
n.pheno   = countStrainsPerTrait(cross$pheno) 
pdata.01  = as.matrix(cross$pheno)
rownames(pdata.01)<-totalStrain
#first element of each replicate that isn't NA
pdata.02 =  getOnePhenoValue(pheno_raw, pdata.01)
n.pheno.02 = apply(pdata.02, 2, function(x) {sum(!is.na(x))})

pdata.02.2<-pdata.01*2-pdata.02

pdata.01<-na.omit(pdata.01)
pdata.02<-na.omit(pdata.02)
pdata.02.2<-na.omit(pdata.02.2)

pdata.01<-scale(pdata.01)
pdata.02<-scale(pdata.02)
pdata.02.2<-scale(pdata.02.2)

save(pdata.01,file="data/basic/pdata.01.RData")
save(pdata.02,file="data/basic/pdata.02.RData")

################### heritability #########################


#Calculate broad-sense heritability ####################################################################################
broadH2      = sapply(pheno_raw, calc.BroadH2, jackknife=FALSE)
# Calculate SE for broad-sense heritability
broadH2.jk<-matrix(NA,815,405)
for(s in 1:815){
    for(tt in 1:405){
        broadH2.jk[s,tt]<-calcH2(pheno_raw[[tt]][strainSelected[-s]])
	}
}


broadH2.jk.se<-sapply(1:405,function(x){sqrt(sum((broadH2.jk[,x]-broadH2[x])^2)*(length(which(!is.na(broadH2.jk[,x])))-1)/length(which(!is.na(broadH2.jk[,x]))))})

save(broadH2,broadH2.jk,broadH2.jk.se,file="data/basic/broadH2_broadH2.jk.RData")

#########################################################################################################################

# use rrBLUP or regress package
# mixed.solve is in rrBLUP ... constraint is that there can only be one random effect term and random effect term covariance matrix
# regress allows multiple random effects term (slower?) but doesn't compute BLUP estimates
# can also use emma.MLE
# Calculate segregant relatedness
A = A.mat(gdata[match(strainSelected,totalStrain),], shrink=FALSE)/2

# narrow-sense heritability on same scale as broad-sense heritability
narrowh2.one = foreach(i=1:ncol(pdata.02)) %dopar% {print(i); mixed.solve(pdata.02[,i], K=A, method='REML') }
narrowh2.one = sapply(narrowh2.one, function(x) {x$Vu/(x$Ve+x$Vu) })

#calc se of narrow.one
narrow.02.jk<-matrix(NA,815,405)
for(s in 1:815){
    for(tt in 1:405){
	    temp<-mixed.solve(pdata.02[-s,tt], K=A[-s,-s], method='REML')
		narrow.02.jk[s,tt]<-temp$Vu/(temp$Ve+temp$Vu)
	}
}

narrow.02.jk.se<-sapply(1:405,function(x){sqrt(sum((narrow.02.jk[,x]-narrowh2.one[x])^2,na.rm=T)*(length(which(!is.na(narrow.02.jk[,x])))-1)/length(which(!is.na(narrow.02.jk[,x]))))})


save(narrowh2.one,narrow.02.jk,narrow.02.jk.se,file="data/basic/narrow.02.jk.RData")


###

pdata.01_ori=as.matrix(cross$pheno)
rownames(pdata.01_ori)<-totalStrain
pdata.02_ori =  getOnePhenoValue(pheno_raw, pdata.01_ori)

save(pdata.01_ori,file="data/basic/pdata.01_ori.RData")
save(pdata.02_ori,file="data/basic/pdata.02_ori.RData")

#############################################################################

##############################################################################

###

load("data/basic/pdata.02.RData")
load("data/oneRep/perfMat_unCorr_one_rep_20cycles.RData")



### function rewrite

###### get Peak based FDR #############################################################################
getPeakFDR = function(chromosome.peaks.lod, pdata, gdata, perms=100,coreNum=1) {
 
   n.pheno = countStrainsPerTrait(pdata)
   cl<-makeCluster(coreNum,type="FORK")
   registerDoParallel(cl)   
   permpeakLODs = foreach( i = 1:perms ) %dopar% {
        print(i)
        pLODS = get.LOD.by.COR(n.pheno, pdata[sample(1:nrow(pdata)),], gdata)
        sapply(mindex.split, function(markers) { apply(pLODS[,markers], 1, max) })

   }
    
	stopCluster(cl)
	
    permpeakLODs= abind(permpeakLODs, along=c(3))


    ###### CHROMOSOME and PEAK BASED FDR #################################################################
    obsPcnt = sapply(seq(2, 5, .01), function(thresh) { sum(chromosome.peaks.lod>thresh) }   )
    names(obsPcnt) = seq(2,5, .01)

    # expected number of QTL peaks with LOD greater than threshold
    expPcnt = sapply(seq(2, 5, .01),  
             function(thresh) { 
                    print(thresh); 
                    mean(apply(permpeakLODs, 3, function(ll) {sum(ll>thresh) }) )
                } )
    names(expPcnt) = seq(2, 5, .01)

    pFDR = expPcnt/obsPcnt

    return(pFDR)
}
########################################################################################################
###### fix QTLs and get residuals phenotypes ###########################################################
getPhenoResids = function(pdata,gdata, peakArray, intercept=FALSE) {
    presids = pdata
    for( i in 1:ncol(pdata) ) {
        spA = peakArray[peakArray$trait==colnames(pdata)[i],]
        if(nrow(spA)>0){
            if(intercept) {
                  rr = residuals(lm(pdata[,i]~gdata[,spA$markerIndex]))
            }else{
                 rr = residuals(lm(pdata[,i]~gdata[,spA$markerIndex]-1))
            }
            presids[which(!is.na(pdata[,i])),i]=rr
        }
    }
    return(presids)
} 

##############################################################################

################################# QTL mapping for Tg ##########################################

### genetic component of pheno

pdata.02_Tg<-sapply(1:405,function(i){apply(cbind(rep(1,815),pdata.02)%*%coefList_unCorr_one_rep[[i]],1,mean,na.rm=T)})
colnames(pdata.02_Tg)<-colnames(pdata.02)

##########

mindex.split = getMarkerIndexSplit(cross)
# get chromosome offsets  (to convert from marker index per chromosome to marker index across genome)
chr.mindex.offset = sapply(mindex.split, min)-1

######extract phenotypes, genotypes, and number of individuals phenotyped per trait 
gdata     = extractGenotype(cross)[match(strainSelected,totalStrain),]
n.pheno   = countStrainsPerTrait(cross$pheno) 

##

LODS.01       = get.LOD.by.COR(n.pheno, pdata.02_Tg, gdata)
LODS.01s      = LODmatrix.2.scanone(LODS.01, cross)
peaklist.01   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.01)
#LODS.01.FDR   = getPeakFDR(peaklist.01$chr.peaks.lod, pdata.02_Tg, gdata, 1000,coreNum=100)
#LODS.01.FDR_old<-LODS.01.FDR
#FDR.01=min(as.numeric(names(LODS.01.FDR)[which(LODS.01.FDR<0.05)]))
peakArray.01  = getPeakArray(peaklist.01, 2.97)
#FDR.01=2.97

##

pdata.02_Tg_round2      = getPhenoResids(pdata.02_Tg[,names(table(peakArray.01[,1]))], gdata, peakArray.01)
LODS.02       = get.LOD.by.COR(n.pheno[match(colnames(pdata.02_Tg_round2),colnames(pdata.02))], pdata.02_Tg_round2, gdata)
LODS.02s      = LODmatrix.2.scanone(LODS.02, cross, LODS.01s)
peaklist.02   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.02) 

#LODS.02.FDR   = getPeakFDR(peaklist.02$chr.peaks.lod, pdata.02_Tg_round2, gdata, 1000,coreNum=100)
#LODS.02.FDR_old<-LODS.02.FDR
#FDR.02=min(as.numeric(names(LODS.02.FDR)[which(LODS.02.FDR<0.05)]))

peakArray.02  = getPeakArray(peaklist.02, 4.52)
peakArray.02  = rbind(peakArray.01, peakArray.02)
#FDR.02=4.52

##

#pdata.02_Tg_round3      = getPhenoResids(pdata.02_Tg_round2[,names(table(peakArray.02[,1]))], gdata, peakArray.02)
#LODS.03       = get.LOD.by.COR(n.pheno[match(colnames(pdata.02_Tg_round3),colnames(pdata.02))], pdata.02_Tg_round3, gdata)
#LODS.03s      = LODmatrix.2.scanone(LODS.03, cross, LODS.01s)
#peaklist.03   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.03) 

#LODS.03.FDR   = getPeakFDR(peaklist.03$chr.peaks.lod, pdata.02_Tg_round3, gdata, 1000,coreNum=100)
#LODS.03.FDR_old<-LODS.03.FDR
#FDR.03=min(as.numeric(names(LODS.03.FDR)[which(LODS.03.FDR<0.05)]))

#peakArray.03  = getPeakArray(peaklist.03, )
#peakArray.03  = rbind(peakArray.01, peakArray.02,peakArray.03)
#FDR.03=





## combine


peak.index = data.frame(peakArray.02,jump=c(rep("J1",dim(peakArray.01)[1]),rep("J2",dim(peakArray.02)[1]-dim(peakArray.01)[1])))
peak.index = peak.index[order(peak.index$trait, peak.index$markerIndex),]
peak.index.list = split(peak.index, peak.index$trait)

peak.index_one_rep_Tg<-peak.index
peak.index.list_one_rep_Tg<-peak.index.list


save(peak.index_one_rep_Tg,peak.index.list_one_rep_Tg,file="data/QTLmapping/peak.index.list_one_rep_Tg.RData")


################################# QTL mapping for Tng ##########################################


######################## QTL mapping for all traits (mean value of uncorrelated prediction, pdata.01_genetic) ######
### genetic component of pheno

pdata.02_Tng<-sapply(1:405,function(i){apply(pdata.02[,i]-cbind(rep(1,815),pdata.02)%*%coefList_unCorr_one_rep[[i]],1,mean,na.rm=T)})
colnames(pdata.02_Tng)<-colnames(pdata.02)

##########

mindex.split = getMarkerIndexSplit(cross)
# get chromosome offsets  (to convert from marker index per chromosome to marker index across genome)
chr.mindex.offset = sapply(mindex.split, min)-1

######extract phenotypes, genotypes, and number of individuals phenotyped per trait 
gdata     = extractGenotype(cross)[match(strainSelected,totalStrain),]
n.pheno   = countStrainsPerTrait(cross$pheno) 

##

LODS.01       = get.LOD.by.COR(n.pheno, pdata.02_Tng, gdata)
LODS.01s      = LODmatrix.2.scanone(LODS.01, cross)
peaklist.01   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.01)
#LODS.01.FDR   = getPeakFDR(peaklist.01$chr.peaks.lod, pdata.02_Tng, gdata, 1000,coreNum=100)
#LODS.01.FDR_old<-LODS.01.FDR
#FDR.01=min(as.numeric(names(LODS.01.FDR)[which(LODS.01.FDR<0.05)]))
peakArray.01  = getPeakArray(peaklist.01, 4)
#FDR.01=4

##

pdata.02_Tng_round2      = getPhenoResids(pdata.02_Tng[,names(table(peakArray.01[,1]))], gdata, peakArray.01)
LODS.02       = get.LOD.by.COR(n.pheno[match(colnames(pdata.02_Tng_round2),colnames(pdata.02))], pdata.02_Tng_round2, gdata)
LODS.02s      = LODmatrix.2.scanone(LODS.02, cross, LODS.01s)
peaklist.02   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.02) 

#LODS.02.FDR   = getPeakFDR(peaklist.02$chr.peaks.lod, pdata.02_Tng_round2, gdata, 1000,coreNum=100)
#LODS.02.FDR_old<-LODS.02.FDR
#FDR.02=min(as.numeric(names(LODS.02.FDR)[which(LODS.02.FDR<0.05)]))

peakArray.02  = getPeakArray(peaklist.02, 4.57)
peakArray.02  = rbind(peakArray.01, peakArray.02)
#FDR.02=4.57

##

#pdata.02_Tng_round3      = getPhenoResids(pdata.02_Tng_round2[,names(table(peakArray.02[,1]))], gdata, peakArray.02)
#LODS.03       = get.LOD.by.COR(n.pheno[match(colnames(pdata.02_Tng_round3),colnames(pdata.02))], pdata.02_Tng_round3, gdata)
#LODS.03s      = LODmatrix.2.scanone(LODS.03, cross, LODS.01s)
#peaklist.03   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.03) 

#LODS.03.FDR   = getPeakFDR(peaklist.03$chr.peaks.lod, pdata.02_Tng_round3, gdata, 1000,coreNum=100)
#LODS.03.FDR_old<-LODS.03.FDR
#FDR.03=min(as.numeric(names(LODS.03.FDR)[which(LODS.03.FDR<0.05)]))

#peakArray.03  = getPeakArray(peaklist.03, )
#peakArray.03  = rbind(peakArray.01, peakArray.02,peakArray.03)
#FDR.03=





## combine


peak.index = data.frame(peakArray.02,jump=c(rep("J1",dim(peakArray.01)[1]),rep("J2",dim(peakArray.02)[1]-dim(peakArray.01)[1])))
peak.index = peak.index[order(peak.index$trait, peak.index$markerIndex),]
peak.index.list = split(peak.index, peak.index$trait)

peak.index_one_rep_Tng<-peak.index
peak.index.list_one_rep_Tng<-peak.index.list


save(peak.index_one_rep_Tng,peak.index.list_one_rep_Tng,file="data/QTLmapping/peak.index.list_one_rep_Tng.RData")







