
# 02/13/22
# R script in:
# The trait coding rule in phenotype space
# 
# Jianguo Wang (wangjg36@mail.sysu.edu.cn)
# 




##################### features of yeast traits ################################



yeastFeature<-data.frame(
traitName=colnames(pdata.02),
H2=broadH2,
H2_se=broadH2.jk.se,
h2=narrowh2.one,
h2_se=narrow.02.jk.se,
R2_UBHDD=apply(perfMat_unCorr_one_rep,1,mean_replace),
R2_UBHDD_se=apply(perfMat_unCorr_one_rep,1,sd_replace),
R2_UBHDD_shuffle=apply(perfMat_unCorr_one_rep_shuffle,1,mean_replace),
R2_UBHDD_shuffle_se=apply(perfMat_unCorr_one_rep_shuffle,1,sd_replace),
R2_TG_Tg,
R2_TG_Tng,
h2_Tg=h2_one_rep_Tg,
h2_Tng=h2_one_rep_Tng,
QTL_num_Tg=QTL_num_one_rep_Tg,
QTL_num_Tng=QTL_num_one_rep_Tng,
identical_score=idenScore2_Tg_in_del
)

rownames(yeastFeature)<-c(1:405)


save(yeastFeature,file="yeastFeature.RData")
write.csv(yeastFeature,"yeastFeature.csv")



#################### features of brain traits ######################################


brainFeature<-data.frame(
traitName=prediction[,1],
R2_UBHDD=prediction[,2],
R2_UBHDD_se=prediction[,3],
R2_UBHDD_shuffle=prediction_shuffle[,2],
R2_UBHDD_shuffle_se=prediction_shuffle[,3],
if_exemplar=rep(NA,dim(prediction)[1]),
h2_Tg=rep(NA,dim(prediction)[1]),
h2_Tng=rep(NA,dim(prediction)[1]),
QTL_num_Tg=rep(NA,dim(prediction)[1]),
QTL_num_Tng=rep(NA,dim(prediction)[1])
)

brainFeature$if_exemplar[clusterName]<-TRUE
brainFeature$h2_Tg[clusterName]<-h2_Tg
brainFeature$h2_Tng[clusterName]<-h2_Tng
brainFeature$QTL_num_Tg[clusterName]<-Tg_QTL_num_clump
brainFeature$QTL_num_Tng[clusterName]<-Tng_QTL_num_clump


save(brainFeature,file="brainFeature.RData")
write.csv(brainFeature,"brainFeature.csv")



#################### features of symmetrical brain region traits ############################

brainPairFeature<-data.frame(
traitName_left=prediction[na.omit(match(brain_left_traitName,prediction[,1])),1],
traitName_right=prediction[na.omit(match(brain_right_traitName,prediction[,1])),1],
R2_Tg_LR,
R2_Tng_LR,
brain_region=traitAnnot$brainRegion[match(names(R2_T_LR),paste0("X",traitAnnot[,1]))],
measure_type=sapply(traitAnnot$dataType[match(names(R2_T_LR),paste0("X",traitAnnot[,1]))],function(x){unlist(strsplit(x," "))[2]})
)

rownames(brainPairFeature)<-c(1:297)

save(brainPairFeature,file="brainPairFeature.RData")
write.csv(brainPairFeature,"brainPairFeature.csv")



