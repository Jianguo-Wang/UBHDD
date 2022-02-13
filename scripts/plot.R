
# 02/13/22
# R script in:
# The trait coding rule in phenotype space
# 
# Jianguo Wang (wangjg36@mail.sysu.edu.cn)
# 


####################################### library #####################################

setwd("E:/papers/UBHDD/figs")

library(zipfR)
library(ggplot2)

library(Hmisc)

library(parallel)
library(foreach)
library(doParallel)

library(lme4)
library(rrBLUP)




####################################### functions #########################################################



cal_identical_score2<-function(x,y){


    if(class(x)=="numeric"&class(y)=="numeric"){
	    idenScore<-1-sum((x-y)^2)/2/sum(x^2+y^2)
	}
    if(class(x)=="matrix"&class(y)=="matrix"){
	    idenScore<-matrix(NA,nrow=dim(x)[2],ncol=dim(y)[2])
	    for(i in 1:dim(x)[2]){
		    for(j in 1:dim(y)[2]){
				idenScore[i,j]<-1-sum((x[,i]-y[,j])^2)/2/sum(x[,i]^2+y[,j]^2)
			}
		}
	}
    if(class(x)=="numeric"&class(y)=="matrix"){
	    idenScore<-rep(NA,dim(y)[2])
		for(j in 1:dim(y)[2]){
		    idenScore[j]<-1-sum((x-y[,j])^2)/2/sum(x^2+y[,j]^2)
		}
	}
    if(class(x)=="matrix"&class(y)=="numeric"){
	    idenScore<-rep(NA,dim(x)[2])
		for(i in 1:dim(x)[2]){
		    idenScore[i]<-1-sum((x[,i]-y)^2)/2/sum(x[,i]^2+y^2)
		}
	}   	
    
	return(idenScore)
}

###

get_ll_list_mean<-function(x){

    return(sapply(x,function(y){
	    if(length(which(!is.na(y)))==0){
		    return(0)
		}else{
	        return(mean_replace(sapply(y,function(z){
			    if(is.na(z)){
				    return(0)
				}else{
				    z[[1]][1]
				}
			})))  
		}
    }))
    
}

# 

get_ll_list_sd<-function(x){

    return(sapply(x,function(y){
	    if(length(which(!is.na(y)))==0){
		    return(0)
		}else{
	        return(sd_replace(sapply(y,function(z){
			    if(is.na(z)){
				    return(0)
				}else{
				    z[[1]][1]
				}
			})))  
		}
    }))
    
}

# 

get_llm_mean<-function(x){
    return(apply(x[[1]],1,mean_replace))
}

# 

get_llm_sd<-function(x){
    return(apply(x[[1]],1,sd_replace))
}

# 


mean_replace<-function(x){
    x[is.na(x)]<-0
	return(mean(x))
}

# 

NA_replace<-function(x){
	return(sapply(x,function(y){ifelse(is.na(y),0,y)}))
}

# 


sd_replace<-function(x){
    x[is.na(x)]<-0
	return(sd(x))
}

# 


theme_set(theme_classic())
theme_update(panel.border=element_blank(),
             plot.title=element_text(size=rel(8/11),hjust=0.5),
             plot.title.position="plot",
             axis.text=element_text(size=rel(8/11)),
             axis.text.y=element_text(angle=90,hjust=0.5),
             axis.title=element_text(size=rel(9/11)),
             axis.title.x=element_text(hjust=0.5),
             axis.title.y=element_text(hjust=0.5),
             plot.margin=unit(c(0.05,0.15,0.15,0.15),"in"),
             legend.title=element_blank(),
             legend.background = element_rect(colour = NA,fill = NA),
             legend.key.size = unit(8, "pt"),
             legend.key.height = unit(8, "pt"),
             legend.text=element_text(size=rel(7/11))
)






####################################### plot ###############################################################




############################## Fig.1_plot_theory



### independent probability of uncorrelated traits

R=0.15
phi=acos(R)
#phi=acos(0.15)
k=3
n<-c(4:100)

#P<-1/(1+k/(n-k)*(1-Rbeta(sin(phi)^2,(k-1)/2,1/2)))
P=1/(1+3/(n-3)*R)


##ggpot2

dataTemp<-data.frame(n,P)
ggplot(dataTemp,aes(x=n,y=P))+geom_point(show.legend=FALSE,cex=0.2)+geom_line(show.legend=FALSE)+
labs(x=expression(paste("Dimension (",italic(n),")")),y=paste("Probability of independence","\n","of uncorrelated traits"),title="")

#ggsave("theory_probability_of_independence_of_uncorrelated_traits.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### probability of uncorrelated traits sharing common traits


#R<-c(runif(1000,0,0.95),runif(1000,0.94,1))
R<-c(runif(1000,0.1,0.95),runif(1000,0.94,1))
phi=acos(R)
#phi=acos(0.15)
k=2

n10<-10
P10<-3*pi/(3*pi+2*(n10-2)*(1-cos(phi)))

n100<-100
P100<-3*pi/(3*pi+2*(n100-2)*(1-cos(phi)))

n1000<-1000
P1000<-3*pi/(3*pi+2*(n1000-2)*(1-cos(phi)))


##ggplot2

dataTemp<-data.frame(R2=c(R^2,R^2,R^2),P=c(P10,P100,P1000),dims=c(rep("a",length(P10)),rep("b",length(P100)),rep("c",length(P1000))))
ggplot(dataTemp,aes(x=R2,y=P,col=dims))+geom_point(show.legend=FALSE,cex=0.1,stroke=0.1)+geom_line(lwd=0.5,show.legend=FALSE)+
scale_color_manual(labels = c("a","b","c"),values=c("#F8766D","black","#00BFC4"))+
scale_x_continuous(breaks=c(0,0.5,1))+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(x=expression(italic(R^2)),y=paste("Probability of sharing","\n","the same dimensions  "),title="")

#ggsave("theory_probability_of_sharing_the_same_dimensions.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))




############################## Fig.1_Fig.S1_Fig.S2_simulation


### 


load("E:/papers/UBHDD/theory/data/sim_UBHDD_10_10000_10_100.RData")
load("E:/papers/UBHDD/theory/data/sim_UBHDD_10_10000_struct.RData")

load("E:/papers/UBHDD/theory/data/sim_UBHDD_20_10000_10_100.RData")
load("E:/papers/UBHDD/theory/data/sim_UBHDD_50_10000_10_100.RData")
load("E:/papers/UBHDD/theory/data/sim_UBHDD_100_10000_10_100.RData")
load("E:/papers/UBHDD/theory/data/sim_UBHDD_100_10000_10_200.RData")
load("E:/papers/UBHDD/theory/data/sim_UBHDD_100_10000_10_200_samSize2000.RData")




############### dim estimation of PG and PNG with trait sampling #################################################################

###


traitRandomInd<-foreach(i=1:100,.combine=cbind)%do%{
    set.seed(i)
    return(sample(1000,1000))
}

#
	
dim_est_TG<-foreach(j=1:100,.combine=cbind)%do%{
    print(j)
	collect_dim_TG<-foreach(i=c(1:50),.combine=c)%do%{
	    focalUse<-abs(sim_10_10000_10_100[[7]][,traitRandomInd[1:(i*20),j]])
	    temp<-rep(0,10)
	    for(k in 1:dim(focalUse)[2]){
		    temp<-temp+focalUse[,k]
		}
        return(length(which(temp>0)))
   	}
}

#

dim_est_TNG<-foreach(j=1:100,.combine=cbind)%do%{
    print(j)
	collect_dim_TNG<-foreach(i=1:50,.combine=c)%do%{
	    focalUse<-abs(sim_10_10000_10_100[[8]][,traitRandomInd[1:(i*20),j]])
	    temp<-rep(0,10000)
	    for(k in 1:dim(focalUse)[2]){
		    temp<-temp+focalUse[,k]
		}
        return(length(which(temp>0)))
   	}
}

#

dim_est_TG_small<-foreach(j=1:100,.combine=cbind)%do%{
    print(j)
	collect_dim_TG<-foreach(i=c(1:20),.combine=c)%do%{
	    focalUse<-as.matrix(abs(sim_10_10000_10_100[[7]][,traitRandomInd[1:i,j]]))
	    temp<-rep(0,10)
	    for(k in 1:dim(focalUse)[2]){
		    temp<-temp+focalUse[,k]
		}
        return(length(which(temp>0)))
   	}
}


#save(dim_est_TG,dim_est_TNG,dim_est_TG_small,file="dim_est_TG_TNG.RData")

###

load("dim_est_TG_TNG.RData")


dataTemp<-data.frame(dim_est=c(apply(dim_est_TG,1,mean),apply(dim_est_TNG,1,mean)),
lower=c(apply(dim_est_TG,1,quantile,0.025),apply(dim_est_TNG,1,quantile,0.025)),
upper=c(apply(dim_est_TG,1,quantile,0.975),apply(dim_est_TNG,1,quantile,0.975)),
traitNo=rep(c(1:50)*20,2),Subspace=c(rep("TG",50),rep("TNG",50)))

dataTemp2<-data.frame(dim_est=c(dim_est_TG[,1],dim_est_TNG[,1]),
traitNo=rep(c(1:50)*20,2),Subspace=c(rep("TG",50),rep("TNG",50)))

ggplot(dataTemp2,aes(x=traitNo,y=dim_est,col=Subspace))+geom_point(cex=0.2,stroke=0.2)+geom_line(lwd=0.2)+geom_linerange(data=dataTemp,aes(ymin=lower,ymax=upper),lwd=0.2)+
scale_color_manual(labels = c(expression(italic(P^G)),expression(italic(P^NG))),values=c("#F8766D","#00BFC4"))+
scale_x_continuous(breaks=c(0,500,1000),limits=c(0,1000))+
scale_y_continuous(breaks=c(0,250,500),limits=c(0,500))+
labs(x="No. of traits",y="No. of dimensions",title="")+
theme(
legend.position=c(0.8,0.8),
legend.title=element_blank(),
legend.text=element_text(size=8),
legend.key.height=unit(0,"line"),
legend.key.width=unit(0.8,"line"),
legend.text.align=0,
legend.spacing.x=unit(0.2,"line"),
legend.background=element_blank())

#ggsave("simulation_dim_est_PG_PNG_with_trait_sampling.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#

dataTemp<-data.frame(dim_est=apply(dim_est_TG_small,1,mean),
lower=apply(dim_est_TG_small,1,quantile,0.025),
upper=apply(dim_est_TG_small,1,quantile,0.975),
traitNo=c(1:20))

dataTemp2<-data.frame(dim_est=dim_est_TG_small[,1],
traitNo=c(1:20))

ggplot(dataTemp2,aes(x=traitNo,y=dim_est))+geom_point(col="#F8766D",cex=0.2,stroke=0.2)+geom_line(col="#F8766D",lwd=0.2)+geom_linerange(data=dataTemp,aes(ymin=lower,ymax=upper),lwd=0.2,col="#F8766D")+
scale_color_manual(labels = c(expression(italic(P^G)), expression(italic(P^NG))),values=c("#F8766D","#00BFC4"))+
scale_x_continuous(breaks=c(0,20),limits=c(0,20))+
scale_y_continuous(breaks=c(5,10),limits=c(5,10))+
labs(x="No. of traits",y="No. of PC dimensions",title="")+
theme(
legend.position=c(0.8,0.8),
legend.title=element_blank(),
legend.text=element_text(size=8),
legend.key.height=unit(0,"line"),
legend.key.width=unit(0.8,"line"),
legend.text.align=0,
legend.spacing.x=unit(0.2,"line"),
legend.background=element_blank())


#ggsave("simulation_dim_est_PG_with_trait_sampling_inner_panel.pdf",width=unit(0.9,"in"),height=unit(1,"in"))


############### dimension sharing with correlation level #################################################################

###


sim_10_10000_10_100_corMat<-cor(sim_10_10000_10_100[[1]])

sim_10_10000_10_100_PGshare<-matrix(nrow=1000,ncol=1000)
sim_10_10000_10_100_PNGshare<-matrix(nrow=1000,ncol=1000)

for(i in 1:999){
    for(j in (i+1):1000){
	    print(i)
	    sim_10_10000_10_100_PGshare[i,j]<-length(which(sim_10_10000_10_100[[7]][,i]*sim_10_10000_10_100[[7]][,j]!=0))
		sim_10_10000_10_100_PGshare[j,i]<-sim_10_10000_10_100_PGshare[i,j]
	    sim_10_10000_10_100_PNGshare[i,j]<-length(which(sim_10_10000_10_100[[8]][,i]*sim_10_10000_10_100[[8]][,j]!=0))
		sim_10_10000_10_100_PNGshare[j,i]<-sim_10_10000_10_100_PNGshare[i,j]
	}
}

M<-sim_10_10000_10_100_corMat
dimShare<-na.omit(data.frame(row=c(1:dim(M)[1])[row(M)], col=c(1:dim(M)[2])[col(M)], R=c(M),PGShare=c(sim_10_10000_10_100_PGshare),PNGShare=c(sim_10_10000_10_100_PNGshare)))

R_level<-rep(NA,dim(dimShare)[1])
R_level[which(abs(dimShare[,3])<sim_unCorCutoff)]="Uncorrelated"
R_level[which(abs(dimShare[,3])>=sim_unCorCutoff&abs(dimShare[,3])<0.4)]="[0.14,0.4)"
R_level[which(abs(dimShare[,3])>=0.4&abs(dimShare[,3])<0.6)]="[0.4,0.6)"
R_level[which(abs(dimShare[,3])>=0.6&abs(dimShare[,3])<0.8)]="[0.6,0.8)"
R_level[which(abs(dimShare[,3])>=0.8&abs(dimShare[,3])<=1)]="[0.8,1]"


dimShareComb<-data.frame(R_level=c("Uncorrelated","[0.14,0.4)","[0.4,0.6)","[0.6,0.8)","[0.8,1]"),
PGShare_mean=sapply(c("Uncorrelated","[0.14,0.4)","[0.4,0.6)","[0.6,0.8)","[0.8,1]"),function(x){mean(dimShare[R_level==x,4]/5*100)}),
PGShare_sd=sapply(c("Uncorrelated","[0.14,0.4)","[0.4,0.6)","[0.6,0.8)","[0.8,1]"),function(x){sd(dimShare[R_level==x,4]/5*100)}),
PNGShare_mean=sapply(c("Uncorrelated","[0.14,0.4)","[0.4,0.6)","[0.6,0.8)","[0.8,1]"),function(x){mean(dimShare[R_level==x,5]/5*100)}),
PNGShare_sd=sapply(c("Uncorrelated","[0.14,0.4)","[0.4,0.6)","[0.6,0.8)","[0.8,1]"),function(x){sd(dimShare[R_level==x,5]/5*100)}))


#save(dimShareComb,file="dimShareComb.RData")


### x-axis ticks are adjusted by drawing software to transform into square

load("dimShareComb.RData")

dataTemp<-data.frame(Rlevel=c(0.07,0.27,0.5,0.7,0.9,0.07,0.27,0.5,0.7,0.9),Mean=c(dimShareComb[,2],dimShareComb[,4]),lower=c(dimShareComb[,2],dimShareComb[,4])-c(dimShareComb[,3],dimShareComb[,5]),upper=c(dimShareComb[,2],dimShareComb[,4])+c(dimShareComb[,3],dimShareComb[,5]),Subspace=c(rep("PG",dim(dimShareComb)[1]),rep("PNG",dim(dimShareComb)[1])))
ggplot(dataTemp,aes(x=Rlevel,y=Mean,col=Subspace))+geom_point(cex=0.5)+geom_line(lwd=0.5)+geom_linerange(aes(ymin=lower,ymax=upper))+
scale_color_manual(labels = c(expression(italic(P^G)), expression(italic(P^NG))),values=c("#F8766D","#00BFC4"))+
scale_x_continuous(breaks=c(0,0.5,1),limit=c(0,1))+
scale_y_continuous(breaks=c(0,50,100))+
labs(x=expression(paste("Correlation level (",italic(R^2),")")),y=paste("% of common dimensions"),title="")+
theme(
legend.position=c(0.3,0.9),
legend.title=element_blank(),
legend.text=element_text(size=8),
legend.key.height=unit(0,"line"),
legend.key.width=unit(0.8,"line"),
legend.text.align=0,
legend.spacing.x=unit(0.2,"line"),
legend.background=element_blank())

#ggsave("simulation_dimension_sharing_Vs_correlation_level.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

############### correlation with TG of Tg and Tng #################################################################

###

PG<-sim_10_10000_10_100[[2]]
PNG<-sim_10_10000_10_100[[3]]
Pg<-cbind(rep(1,1000),sim_10_10000_10_100[[1]])%*%sim_UBHDD_10_10000_10_100_coefMat
Png<-sim_10_10000_10_100[[1]]-Pg

#

R2_TG_Tg<-sapply(1:1000,function(x){cor(Pg[,x],PG[,x])})^2
R2_TG_Tng<-sapply(1:1000,function(x){cor(Png[,x],PG[,x])})^2
R2_TNG_Tg<-sapply(1:1000,function(x){cor(Pg[,x],PNG[,x])})^2
R2_TNG_Tng<-sapply(1:1000,function(x){cor(Png[,x],PNG[,x])})^2


# 

dataTemp<-data.frame(R2=c(R2_TG_Tg,R2_TG_Tng),subspace=c(rep("a",1000),rep("b",1000)))
ggplot(dataTemp,aes(y=R2,x=subspace))+geom_boxplot(width=0.3, color="black", alpha=1,show.legend=FALSE,outlier.size=0.1,outlier.stroke = 0.3,size=0.3)+
scale_color_discrete(labels = c(expression(italic(T^g)),expression(italic(T^ng))))+
scale_fill_discrete(labels = c(expression(italic(T^g)),expression(italic(T^ng))))+
scale_x_discrete(labels = c(expression(italic(T^g)),expression(italic(T^ng))))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1.2))+
labs(x="",y=expression(paste("Correlation with ",italic(T^G),"(",italic(R^2),")")),title="")

t.test(R2_TG_Tg,R2_TG_Tng,pairs=T)

#ggsave("simulation_correlation_with_TG_of_Tg_Vs_Tng.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


############### comparison_between_H2_and_R2_by_UBHDD #################################################################


###

dataTemp<-data.frame(H2=sim_10_10000_10_100[[4]],R2=NA_replace(sim_UBHDD_10_10000_10_100_perfMat))
ggplot(dataTemp,aes(x=H2,y=R2))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,stroke=0.2,alpha=1)+
scale_x_continuous(breaks=c(0,0.5,1))+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(x=expression(paste("Variance of ",italic(T^G)," (",italic(H^2),")")),y=expression(paste("Variance of ",italic(T^g))),title="")

cor.test(dataTemp[,1],dataTemp[,2])

#ggsave("simulation_H2_Vs_R2_by_UBHDD_10_10000_10_100.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


###

dataTemp<-data.frame(H2=sim_20_10000_10_100[[4]],R2=NA_replace(sim_UBHDD_20_10000_10_100_perfMat))
ggplot(dataTemp,aes(x=H2,y=R2))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,stroke=0.2,alpha=1)+
scale_x_continuous(breaks=c(0,0.5,1))+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(x=expression(paste("Variance of ",italic(T^G)," (",italic(H^2),")")),y=expression(paste("Variance of ",italic(T^g))),title=expression(paste(italic(N[1]),"=20")))

cor.test(dataTemp[,1],dataTemp[,2])

#ggsave("simulation_H2_Vs_R2_by_UBHDD_20_10000_10_100.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


###

dataTemp<-data.frame(H2=sim_50_10000_10_100[[4]],R2=NA_replace(sim_UBHDD_50_10000_10_100_perfMat))
ggplot(dataTemp,aes(x=H2,y=R2))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,stroke=0.2,alpha=1)+
scale_x_continuous(breaks=c(0,0.5,1))+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(x=expression(paste("Variance of ",italic(T^G)," (",italic(H^2),")")),y=expression(paste("Variance of ",italic(T^g))),title=expression(paste(italic(N[1]),"=50")))

cor.test(dataTemp[,1],dataTemp[,2])

#ggsave("simulation_H2_Vs_R2_by_UBHDD_50_10000_10_100.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


###

dataTemp<-data.frame(H2=sim_100_10000_10_100[[4]],R2=NA_replace(sim_UBHDD_100_10000_10_100_perfMat))
ggplot(dataTemp,aes(x=H2,y=R2))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,stroke=0.2,alpha=1)+
scale_x_continuous(breaks=c(0,0.5,1))+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(x=expression(paste("Variance of ",italic(T^G)," (",italic(H^2),")")),y=expression(paste("Variance of ",italic(T^g))),title=expression(paste(italic(N[1]),"=100")))

cor.test(dataTemp[,1],dataTemp[,2])

#ggsave("simulation_H2_Vs_R2_by_UBHDD_100_10000_10_100.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


###

dataTemp<-data.frame(H2=sim_100_10000_10_200[[4]],R2=NA_replace(sim_UBHDD_100_10000_10_200_perfMat))
ggplot(dataTemp,aes(x=H2,y=R2))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,stroke=0.2,alpha=1)+
scale_x_continuous(breaks=c(0,0.5,1))+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(x=expression(paste("Variance of ",italic(T^G)," (",italic(H^2),")")),y=expression(paste("Variance of ",italic(T^g))),title=expression(paste(italic(N[1]),"=100")))

cor.test(dataTemp[,1],dataTemp[,2])

#ggsave("simulation_H2_Vs_R2_by_UBHDD_100_10000_10_200.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


###

dataTemp<-data.frame(H2=sim_100_10000_10_200_samSize2000[[4]],R2=NA_replace(sim_UBHDD_100_10000_10_200_samSize2000_perfMat))
ggplot(dataTemp,aes(x=H2,y=R2))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,stroke=0.2,alpha=1)+
scale_x_continuous(breaks=c(0,0.5,1))+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(x=expression(paste("Variance of ",italic(T^G)," (",italic(H^2),")")),y=expression(paste("Variance of ",italic(T^g))),title=expression(paste(italic(N[1]),"=100")))

cor.test(dataTemp[,1],dataTemp[,2])

#ggsave("simulation_H2_Vs_R2_by_UBHDD_100_10000_10_200_samSize2000.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



############### structured population cluster #################################################################

###

sim_10_10000_struct_PCA<-princomp(sim_10_10000_struct[[1]])

##

quantile(sim_10_10000_struct[[4]])
quantile(sim_UBHDD_10_10000_struct_perfMat)

cumsum(sim_10_10000_struct_PCA$sdev^2)[1:5]/sum(sim_10_10000_struct_PCA$sdev^2)

Var_PCA<-apply(sim_10_10000_struct_PCA$scores[,1:2]%*%t(sim_10_10000_struct_PCA$loadings)[1:2,],2,function(x){var(x)})

#save(Var_PCA,file="Var_PCA_sim_10_10000_struct.RData")

### ggplot2

load("Var_PCA_sim_10_10000_struct.RData")

#

dataTemp<-data.frame(PC1=jitter(t(sim_10_10000_struct_PCA$loadings)[1,],factor=1500),PC2=jitter(t(sim_10_10000_struct_PCA$loadings)[2,],factor=300))
dataTemp<-dataTemp[c(3:52,551:1000,1,53:351,2,352:550),]
col_cluster=c(rep("black",500),rep("red",300),rep("blue",200))

ggplot(dataTemp,aes(x=PC1,y=PC2))+geom_point(cex=0.2,stroke=0.2,alpha=1,col=col_cluster)+
scale_x_continuous(breaks=c(-0.02,0,0.05))+
scale_y_continuous(breaks=c(-0.05,0,0.05))+
labs(x="PC1",y="PC2",title="Structured space")

#ggsave("simulation_structured_population_cluster.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### 


dataTemp<-data.frame(H2=sim_10_10000_struct[[4]],R2=sim_UBHDD_10_10000_struct_perfMat)
col_cluster=c(rep("black",500),rep("red",300),rep("blue",200))
dataTemp<-dataTemp[c(3:52,551:1000,1,53:351,2,352:550),]

ggplot(dataTemp,aes(x=H2,y=R2))+geom_point(cex=0.2,stroke=0.2,alpha=1,col=col_cluster)+geom_abline(col="black",lwd=0.5,intercept=0,slope=1)+
scale_x_continuous(breaks=c(0,0.5,1))+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(x=expression(paste("Variance of ",italic(T^G))),y=expression(paste("Variance of ",italic(T^g))),title="UBHDD")

cor.test(dataTemp[,1],dataTemp[,2])

#ggsave("simulation_structured_population_H2_Vs_R2_UBHDD.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


###PCA


dataTemp<-data.frame(H2=sim_10_10000_struct[[4]],R2=Var_PCA)
dataTemp<-dataTemp[c(3:52,551:1000,1,53:351,2,352:550),]
col_cluster=c(rep("black",500),rep("red",300),rep("blue",200))

ggplot(dataTemp,aes(x=H2,y=R2))+geom_point(cex=0.2,alpha=1,stroke=0.2,col=col_cluster)+geom_abline(col="black",lwd=0.5,intercept=0,slope=1)+
scale_x_continuous(breaks=c(0,0.5,1))+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(x=expression(paste("Variance of",italic(T^G))),y=expression(paste("Variance of ",italic(T^pc))),title="PCA")


#ggsave("simulation_structured_population_H2_Vs_Var_PCA.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



############################## Fig.2_Fig.S3_plot_TG_Tg_Tng_yeastSegregants




load("E:/papers/UBHDD/yeastSegregants/data/basic/pdata.02.RData")
load("E:/papers/UBHDD/yeastSegregants/data/basic/pdata.02_ori.RData")
load("E:/papers/UBHDD/yeastSegregants/data/basic/cross.RData")
load("E:/papers/UBHDD/yeastSegregants/data/basic/pheno_raw_405traits.RData")

load("E:/papers/UBHDD/yeastSegregants/data/basic/broadH2_broadH2.jk.RData")
load("E:/papers/UBHDD/yeastSegregants/data/basic/narrow.02.jk.RData")

load("E:/papers/UBHDD/yeastSegregants/data/oneRep/perfMat_unCorr_one_rep_20cycles.RData")
load("E:/papers/UBHDD/yeastSegregants/data/oneRep/perfMat_unCorr_one_rep_shuffle_20cycles.RData")



###########################################################

### UBHDD

Pg<-sapply(1:405,function(i){apply(cbind(rep(1,815),pdata.02)%*%coefList_unCorr_one_rep[[i]],1,mean)})
colnames(Pg)<-colnames(pdata.02)
Png<-pdata.02-Pg



### LMM

PG<-foreach(i=1:405,.combine=cbind)%do%{
    print(i)
    pdata<-pheno_raw[[i]]
	indName<-names(pdata)
	
	dataTemp<-matrix(NA,815,1)
	rownames(dataTemp)<-indName
	
    pdata.notNAcnt = sapply(pdata, function(x){sum(!is.na(x))})
    pdata[pdata.notNAcnt<2]=NULL

    sp<-stack(pdata)

    reffMod = lmer(values ~ 1 + (1|ind),data=sp)
    fittedValue<-fitted(reffMod)[c(1:(dim(sp)[1]/2))*2-1]
	fittedValue_scale<-(fittedValue-mean(pdata.02_ori[,i],na.rm=T))/sd(pdata.02_ori[,i],na.rm=T)
	
	dataTemp[match(sp[c(1:(dim(sp)[1]/2))*2-1,2],indName),1]<-fittedValue_scale
    return(dataTemp)    
}


colnames(PG)<-colnames(pdata.02)

PNG<-pdata.02-PG


#save(Pg,Png,PG,PNG,file="Pg_Png_PG_PNG.RData")

###

R2_TG_Tg<-sapply(1:405,function(x){cor(PG[,x],Pg[,x],use="pairwise.complete.obs")^2})
R2_TG_Tng<-sapply(1:405,function(x){cor(PG[,x],Png[,x],use="pairwise.complete.obs")^2})
R2_TNG_Tg<-sapply(1:405,function(x){cor(PNG[,x],Pg[,x],use="pairwise.complete.obs")^2})
R2_TNG_Tng<-sapply(1:405,function(x){cor(PNG[,x],Png[,x],use="pairwise.complete.obs")^2})

#save(R2_TG_Tg,R2_TG_Tng,R2_TNG_Tg,R2_TNG_Tng,file="R2_TG_Tg_TNG_Tng.RData")


load("Pg_Png_PG_PNG.RData")
load("R2_TG_Tg_TNG_Tng.RData")


### 

dataTemp<-data.frame(R2=c(R2_TG_Tg,R2_TG_Tng),subspace=c(rep("a",405),rep("b",405)))
ggplot(dataTemp,aes(y=R2,x=subspace))+geom_boxplot(width=0.3, color="black", alpha=1,show.legend=FALSE,outlier.size=0.1,outlier.stroke = 0.3,size=0.3)+
scale_color_discrete(labels = c(expression(italic(T^g)),expression(italic(T^ng))))+
scale_fill_discrete(labels = c(expression(italic(T^g)),expression(italic(T^ng))))+
scale_x_discrete(labels = c(expression(italic(T^g)),expression(italic(T^ng))))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1.2))+
labs(x="",y=expression(paste("Correlation with ",italic(T^G),"(",italic(R^2),")")),title="")

t.test(R2_TG_Tg,R2_TG_Tng,pairs=T)

#ggsave("yeastSegregants_correlation_with_TG_of_Tg_Vs_Tng.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



### a case for the correlation with TG of Tg Vs Tng


dataTemp<-data.frame(TG1=PG[,1],Tg1=Pg[,1],Tng1=Png[,1])

ggplot(dataTemp,aes(x=Tg1,y=TG1))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.3)+
scale_x_continuous(breaks=c(-2,0,2,4,6),limits=c(-2.5,6.5))+
scale_y_continuous(breaks=c(-2,0,2,4,6),limits=c(-2.5,6.5))+
labs(x=expression(italic(T^g)),y=expression(italic(T^G)),title="C11.1_A")

cor.test(dataTemp[,1],dataTemp[,2])

#ggsave("yeastSegregants_R2_TG_Tg_C11.1_A.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#

dataTemp<-data.frame(TG1=PG[,1],Tg1=Pg[,1],Tng1=Png[,1])

ggplot(dataTemp,aes(x=Tng1,y=TG1))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.3)+
scale_x_continuous(breaks=c(-2,0,2,4,6),limits=c(-2.5,6.5))+
scale_y_continuous(breaks=c(-2,0,2,4,6),limits=c(-2.5,6.5))+
labs(x=expression(italic(T^ng)),y=expression(italic(T^G)),title="C11.1_A")


cor.test(dataTemp[,1],dataTemp[,3])

#ggsave("yeastSegregants_R2_TG_Tng_C11.1_A.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



###

dataTemp<-data.frame(H2=broadH2,R2=apply(perfMat_unCorr_one_rep,1,mean_replace))

ggplot(dataTemp,aes(x=R2,y=H2))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.3)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste("Variance of ",italic(T^g))),y=expression(paste("Variance of ",italic(T^G)," (",italic(H^2),")")),title="")

cor.test(dataTemp[,1],dataTemp[,2])

#ggsave("yeastSegregants_variance_of_TG_Vs_Tg.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### R2 of UBHDD

dataTemp<-data.frame(R2=apply(perfMat_unCorr_one_rep,1,mean_replace))

ggplot(dataTemp,aes(x=R2))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.15,alpha=0.4)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste("Variance of ",italic(T^g))),y="No. of traits",title=" ")


#ggsave("yeastSegregants_variance_of_Tg.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### R2 of UBHDD (shuffling)


dataTemp<-data.frame(R2=apply(perfMat_unCorr_one_rep_shuffle,1,mean_replace))

ggplot(dataTemp,aes(x=R2))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.15,alpha=0.4)+
scale_x_continuous(breaks=c(0,0.013),limits=c(0,0.013))+
scale_y_continuous(breaks=c(0,60))+
labs(x=expression(paste("Variance of ",italic(T^g))),y="No. of traits",title=" ")


#ggsave("yeastSegregants_variance_of_Tg_shuffle.pdf",width=unit(1,"in"),height=unit(1,"in"))


### H2 and h2


dataTemp<-data.frame(H2=broadH2,h2=narrowh2.one)

ggplot(dataTemp,aes(x=H2,y=h2))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.3)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(italic(H^2)),y=expression(italic(h^2)),title="")


#ggsave("yeastSegregants_H2_Vs_h2.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))




############################## Fig.2_plot_QTL_mapping_yeastSegregants




load("E:/papers/UBHDD/yeastSegregants/data/QTLmapping/peak.index.list_one_rep_T.RData")
load("E:/papers/UBHDD/yeastSegregants/data/QTLmapping/peak.index.list_one_rep_Tg.RData")
load("E:/papers/UBHDD/yeastSegregants/data/QTLmapping/peak.index.list_one_rep_Tng.RData")

load("E:/papers/UBHDD/yeastSegregants/data/basic/gdata_comm.RData")
load("E:/papers/UBHDD/yeastSegregants/data/basic/pdata.02.RData")
load("E:/papers/UBHDD/yeastSegregants/data/oneRep/perfMat_unCorr_one_rep_20cycles.RData")



###

pdata.02_Tg<-sapply(1:405,function(i){apply(cbind(rep(1,815),pdata.02)%*%coefList_unCorr_one_rep[[i]],1,mean,na.rm=T)})
colnames(pdata.02_Tg)<-colnames(pdata.02)

pdata.02_Tng<-pdata.02-pdata.02_Tg
colnames(pdata.02_Tng)<-colnames(pdata.02)

#save(pdata.02_Tg,pdata.02_Tng,file="pdata.02_Tg_Tng.RData")

load("pdata.02_Tg_Tng.RData")

### h2


A = A.mat(gdata_comm, shrink=FALSE)/2

#

h2_one_rep_T = foreach(i=1:ncol(pdata.02),.combine=c) %do% {
    print(i)
	modelTemp<-mixed.solve(pdata.02[,i], K=A, method='REML')
    modelTemp$Vu/(modelTemp$Ve+modelTemp$Vu)	
}



#

h2_one_rep_Tg = foreach(i=1:ncol(pdata.02_Tg),.combine=c) %do% {
    print(i)
	modelTemp<-mixed.solve(pdata.02_Tg[,i], K=A, method='REML')
    modelTemp$Vu/(modelTemp$Ve+modelTemp$Vu)	
}


#

h2_one_rep_Tng = foreach(i=1:ncol(pdata.02_Tng),.combine=c) %do% {
    print(i)
	modelTemp<-mixed.solve(pdata.02_Tng[,i], K=A, method='REML')
    modelTemp$Vu/(modelTemp$Ve+modelTemp$Vu)	
}


#save(h2_one_rep_T,h2_one_rep_Tg,h2_one_rep_Tng,file="h2_one_rep_T_Tg_Tng.RData")

load("h2_one_rep_T_Tg_Tng.RData")

###


QTL_num_one_rep_T<-foreach(i=colnames(pdata.02),.combine=c)%do%{
    focalInd<-match(i,names(peak.index.list_one_rep_T))
	if(is.na(focalInd)){
	    return(0)
	}else{
	    return(dim(peak.index.list_one_rep_T[[focalInd]])[1])
	
	}
}


QTL_num_one_rep_Tg<-foreach(i=colnames(pdata.02),.combine=c)%do%{
    focalInd<-match(i,names(peak.index.list_one_rep_Tg))
	if(is.na(focalInd)){
	    return(0)
	}else{
	    return(dim(peak.index.list_one_rep_Tg[[focalInd]])[1])
	
	}
}

QTL_num_one_rep_Tng<-foreach(i=colnames(pdata.02),.combine=c)%do%{
    focalInd<-match(i,names(peak.index.list_one_rep_Tng))
	if(is.na(focalInd)){
	    return(0)
	}else{
	    return(dim(peak.index.list_one_rep_Tng[[focalInd]])[1])
	
	}
}


#save(QTL_num_one_rep_T,QTL_num_one_rep_Tg,QTL_num_one_rep_Tng,file="QTL_num_one_rep_T_Tg_Tng.RData")

load("QTL_num_one_rep_T_Tg_Tng.RData")

### h2 


dataTemp<-data.frame(h2_Tg=h2_one_rep_Tg,h2_Tng=h2_one_rep_Tng)
ggplot(dataTemp,aes(x=h2_Tg,y=h2_Tng))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.1)+
scale_x_continuous(breaks=c(0,0.3,0.6),limits=c(0,0.6))+
scale_y_continuous(breaks=c(0,0.3,0.6),limits=c(0,0.6))+
labs(x=expression(italic(T^g)),y=expression(italic(T^ng)),title=expression(italic(h^2)))

length(which(h2_one_rep_Tg>h2_one_rep_Tng))
length(which(h2_one_rep_Tg<h2_one_rep_Tng))

#ggsave("yeastSegregants_h2_Tg_Tng.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### QTL_num - Tg Vs Tng

dataTemp<-data.frame(QTL_num_Tg=jitter(QTL_num_one_rep_Tg),QTL_num_Tng=jitter(QTL_num_one_rep_Tng))
ggplot(dataTemp,aes(x=QTL_num_Tg,y=QTL_num_Tng))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.05,alpha=1,stroke=0.3)+
scale_x_continuous(breaks=c(0,6,12),limits=c(-0.5,12.5))+
scale_y_continuous(breaks=c(0,6,12),limits=c(-0.5,12.5))+
labs(x=expression(italic(T^g)),y=expression(italic(T^ng)),title="No. of QTLs")

length(which(QTL_num_one_rep_Tg>QTL_num_one_rep_Tng))
length(which(QTL_num_one_rep_Tg<QTL_num_one_rep_Tng))


#ggsave("yeastSegregants_QTL_num_Tg_Tng.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))




############################## Fig.2_plot_UBHDD_seg_del




load("E:/papers/UBHDD/yeastSegregants/data/onerep/perfMat_unCorr_one_rep_20cycles.RData")
load("E:/papers/UBHDD/yeastDeletome/data/perfMat_unCorr_mt4718_common_scale_20cycles.RData")


#############################################################


## deletome by deletome

Tg_deletome_by_deletome<-foreach(i=1:405)%do%{cbind(rep(1,4718),mt4718_common_scale)%*%coefList_unCorr_mt4718_common_scale[[i]]}
Tng_deletome_by_deletome<-foreach(i=1:405)%do%{mt4718_common_scale[,i]-Tg_deletome_by_deletome[[i]]}

Tg_deletome_by_deletome_mean<-foreach(i=1:405,.combine=cbind)%do%{
    apply(Tg_deletome_by_deletome[[i]],1,mean_replace)
}

Tng_deletome_by_deletome_mean<-foreach(i=1:405,.combine=cbind)%do%{
    apply(Tng_deletome_by_deletome[[i]],1,mean_replace)
}



## deletome by one_rep

Tg_deletome_by_one_rep<-foreach(i=1:405)%do%{cbind(rep(1,4718),mt4718_common_scale)%*%coefList_unCorr_one_rep[[i]]}
Tng_deletome_by_one_rep<-foreach(i=1:405)%do%{mt4718_common_scale[,i]-Tg_deletome_by_one_rep[[i]]}

Tg_deletome_by_one_rep_mean<-foreach(i=1:405,.combine=cbind)%do%{
    apply(Tg_deletome_by_one_rep[[i]],1,mean_replace)
}

Tng_deletome_by_one_rep_mean<-foreach(i=1:405,.combine=cbind)%do%{
    apply(Tng_deletome_by_one_rep[[i]],1,mean_replace)
}




### comparison between Tg_in_del_by_formula_in_seg and Tg_in_del_by_formula_in_del


idenScore2_Tg_in_del<-sapply(1:405,function(x){mean(cal_identical_score2(Tg_deletome_by_deletome[[x]],Tg_deletome_by_one_rep[[x]]),na.rm=T)})


#save(Tg_deletome_by_deletome,Tng_deletome_by_deletome,Tg_deletome_by_deletome_mean,Tng_deletome_by_deletome_mean,
#Tg_deletome_by_one_rep,Tng_deletome_by_one_rep,Tg_deletome_by_one_rep_mean,Tng_deletome_by_one_rep_mean,idenScore2_Tg_in_del,file="plot_UBHDD_seg_del.RData")

load("plot_UBHDD_seg_del.RData")

### case


dataTemp<-data.frame(Tg_in_del_by_del=Tg_deletome_by_deletome[[1]][,1],Tg_in_del_by_seg=Tg_deletome_by_one_rep[[1]][,1])

ggplot(dataTemp,aes(x=Tg_in_del_by_del,y=Tg_in_del_by_seg))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.1)+
scale_x_continuous(breaks=c(0,5,10),limits=c(-4,14))+
scale_y_continuous(breaks=c(0,5,10),limits=c(-4,14))+
labs(x=expression(paste(italic(phi),"(",italic(T[1]^del),",...,",italic(T[n]^del),")")),y=expression(paste(italic(f),"(",italic(T[1]^del),",...,",italic(T[n]^del),")")),title=paste("C11.1_A","\n","(Identity score=0.88)"))


#ggsave("yeastDeletome_seg_del_idenScore_C11.1_A.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### distribution of identical score 

dataTemp<-data.frame(idenScore2_Tg_in_del)

ggplot(dataTemp,aes(x=idenScore2_Tg_in_del))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.25,alpha=0.4)+
scale_x_continuous(breaks=c(0.25,0.5,0.75,1),limits=c(0.25,1))+
scale_y_continuous(breaks=c(0,30,60),limits=c(0,60))+
labs(x="Identity score",y="No. of traits",title=paste("All traits","\n","(n=405)"))


#ggsave("yeastDeletome_seg_del_idenScore_distribution.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))




############################## Fig.3_plot_QTL_mapping_brain




prediction<-read.csv("E:/papers/UBHDD/UKB/data/b2b/result_b2b_water_015_DF.csv")
prediction_shuffle<-read.csv("E:/papers/UBHDD/UKB/data/b2b/result_b2b_water_015_DF_shuffle.csv")


load("E:/papers/UBHDD/UKB/data/genomePreprocessing/clumping_analysis/QTL_mapping_result/h2_Tg_Tng.RData")
load("E:/papers/UBHDD/UKB/data/genomePreprocessing/clumping_analysis/QTL_mapping_result/QTL_num_Tg_Tng.RData")


### R2 of UBHDD in brain


ggplot(prediction,aes(x=Pred))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.15,alpha=0.4)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste("Variance of ",italic(T^g))),y="No. of traits",title=" ")


#ggsave("brain_variance_of_Tg.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### R2 of UBHDD in brain (shuffling)


dataTemp<-prediction_shuffle[,2]
dataTemp[which(dataTemp<0)]=0
ggplot(prediction_shuffle,aes(x=dataTemp))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.15,alpha=0.4)+
scale_x_continuous(breaks=c(0,4e-4),limits=c(0,4e-4))+
scale_y_continuous(breaks=c(0,600))+
labs(x=expression(paste("Variance of ",italic(T^g))),y="No. of traits",title="Brain phenotype (n=675)")


#ggsave("brain_variance_of_Tg_shuffle.pdf",width=unit(1,"in"),height=unit(1,"in"))




### h2 of Tg and Tng

dataTemp<-data.frame(h2_Tg,h2_Tng)

ggplot(dataTemp,aes(x=h2_Tg,y=h2_Tng))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.5)+
scale_x_continuous(breaks=c(0,0.25,0.5),limits=c(0,0.5))+
scale_y_continuous(breaks=c(0,0.25,0.5),limits=c(0,0.5))+
labs(x=expression(italic(T)^g),y=expression(italic(T)^ng),title=expression(italic(h^2)))

length(which(dataTemp[,1]>dataTemp[,2]))
length(which(dataTemp[,2]>dataTemp[,1]))


#ggsave("brain_h2_Tg_Tng.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### QTL num of Tg and Tng

dataTemp<-data.frame(Tg_QTL_num_clump=jitter(Tg_QTL_num_clump),Tng_QTL_num_clump=jitter(Tng_QTL_num_clump))

ggplot(dataTemp,aes(x=Tg_QTL_num_clump,y=Tng_QTL_num_clump))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.5,alpha=1,stroke=0.4,pch=1)+
scale_x_continuous(breaks=c(0,20,40),limits=c(-1,45))+
scale_y_continuous(breaks=c(0,20,40),limits=c(-1,45))+
labs(x=expression(italic(T)^g),y=expression(italic(T)^ng),title="No. of QTLs")

length(which(Tg_QTL_num_clump>Tng_QTL_num_clump))
length(which(Tng_QTL_num_clump>Tg_QTL_num_clump))


#ggsave("brain_QTL_num_Tg_Tng.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))




############################## Fig.3_plot_brain_symmetry_Tg_Tng




load("E:/papers/UBHDD/UKB/data/basic/brain_left_right_traitName.RData")
traitAnnot<-read.csv("E:/papers/UBHDD/UKB/data/basic/traitAnnotation_all18108.csv")


###

load("E:/papers/UBHDD/UKB/data/b2b/brain_with_genome_exclude_covariate_water_Tg.RData")
load("E:/papers/UBHDD/UKB/data/b2b/brain_with_genome_exclude_covariate_water_Tng.RData")


###

Tg_left<-brain_with_genome_exclude_covariate_water_Tg[,na.omit(match(brain_left_traitName,colnames(brain_with_genome_exclude_covariate_water_Tg)))]
Tg_right<-brain_with_genome_exclude_covariate_water_Tg[,na.omit(match(brain_right_traitName,colnames(brain_with_genome_exclude_covariate_water_Tg)))]

Tng_left<-brain_with_genome_exclude_covariate_water_Tng[,na.omit(match(brain_left_traitName,colnames(brain_with_genome_exclude_covariate_water_Tng)))]
Tng_right<-brain_with_genome_exclude_covariate_water_Tng[,na.omit(match(brain_right_traitName,colnames(brain_with_genome_exclude_covariate_water_Tng)))]

T_left<-Tg_left+Tng_left
T_right<-Tg_right+Tng_right


###

R2_T_LR<-sapply(1:dim(T_left)[2],function(x){cor(T_left[,x],T_right[,x])^2})
R2_Tg_LR<-sapply(1:dim(T_left)[2],function(x){cor(Tg_left[,x],Tg_right[,x])^2})
R2_Tng_LR<-sapply(1:dim(T_left)[2],function(x){cor(Tng_left[,x],Tng_right[,x])^2})

names(R2_T_LR)<-names(R2_T_LR)
names(R2_Tg_LR)<-colnames(Tg_left)
names(R2_Tng_LR)<-colnames(Tng_left)


save(R2_T_LR,R2_Tg_LR,R2_Tng_LR,file="R2_T_Tg_Tng_LR.RData")

load("R2_T_Tg_Tng_LR.RData")

###

asyResults<-data.frame(traitName=rep(names(R2_T_LR),2),R2=c(R2_Tg_LR,R2_Tng_LR),subspace=c(rep("Tg",297),rep("Tng",297)),
brainRegion=rep(traitAnnot$brainRegion[match(names(R2_T_LR),paste0("X",traitAnnot[,1]))],2),
dataType=rep(traitAnnot$dataType[match(names(R2_T_LR),paste0("X",traitAnnot[,1]))],2),
dataType2=sapply(rep(traitAnnot$dataType[match(names(R2_T_LR),paste0("X",traitAnnot[,1]))],2),function(x){unlist(strsplit(x," "))[2]}))


p<-ggplot(asyResults,aes(y=brainRegion,x=R2,group=subspace,col=subspace,fill=subspace))+geom_bar(stat="identity",position="dodge",width=0.5)
p+facet_grid(cols = vars(dataType2))+
scale_x_continuous(breaks=c(0,1),limits=c(0,1))+
scale_color_discrete(labels = c(expression(italic(T)^g),expression(italic(T)^ng)))+
scale_fill_discrete(labels = c(expression(italic(T)^g),expression(italic(T)^ng)))+
labs(x=expression(italic(R^2)),y="Brain region",title="")+
theme(axis.text.y=element_text(angle=0,hjust=1))


#ggsave("brain_asymmetry_Tg_Tng_region_measure.pdf",width=unit(7.2,"in"),height=unit(4.7,"in"))



###

#

dataTemp<-data.frame(R2_Tg_LR,R2_Tng_LR)

ggplot(dataTemp,aes(x=R2_Tg_LR,y=R2_Tng_LR))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.2)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(italic(T^g)),y=expression(italic(T^ng)),title=expression(paste(italic(R^2)," of symmetrical regions")))

#ggsave("brain_asymmetry_Tg_Tng.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

#

dataTemp<-data.frame(R2_Tg_LR,R2_T_LR)

ggplot(dataTemp,aes(x=R2_Tg_LR,y=R2_T_LR))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.2)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(italic(T^g)),y=expression(italic(T)),title=expression(paste(italic(R^2)," of symmetrical regions")))

#ggsave("brain_symmetry_Tg_T.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))




############################## Fig.4_plot_dim_estimation_yeast




load("Pg_Png_PG_PNG.RData")


########### dim estimation of Tg and Tng


Pg_scale<-scale(Pg)
Png_scale<-scale(Png)

traitRandomInd<-foreach(i=1:100,.combine=cbind)%do%{
    set.seed(i)
    return(sample(405,405))
}

cl<-makeCluster(27,type="FORK")
registerDoParallel(cl)
	
dim_est_Tg<-sapply(1:100,function(j){
	print(j)
    foreach(i=1:27,.combine=c)%dopar%{
        focal_Tg<-Pg_scale[,traitRandomInd[1:(i*15),j]]
        focal_Tg_PCA<-princomp(focal_Tg)	
	    return(min(which(cumsum(focal_Tg_PCA$sdev^2)/sum(focal_Tg_PCA$sdev^2)>0.85)))
	}
})


dim_est_Tng<-sapply(1:100,function(j){

    foreach(i=1:27,.combine=c)%dopar%{
	    focal_Tng<-Png_scale[,traitRandomInd[1:(i*15),j]]
	    focal_Tng_PCA<-princomp(focal_Tng)	
	    return(min(which(cumsum(focal_Tng_PCA$sdev^2)/sum(focal_Tng_PCA$sdev^2)>0.85)))	    
	}
})

stopCluster(cl)


#save(dim_est_Tg,dim_est_Tng,file="dim_est_Tg_Vs_Tng.RData")

load("dim_est_Tg_Vs_Tng.RData")

###



dataTemp<-data.frame(dim_est=c(apply(dim_est_Tg,1,mean),apply(dim_est_Tng,1,mean)),
lower=c(apply(dim_est_Tg,1,mean)-apply(dim_est_Tg,1,sd)*1.96,apply(dim_est_Tng,1,mean)-1.96*apply(dim_est_Tng,1,sd)),
upper=c(apply(dim_est_Tg,1,mean)+apply(dim_est_Tg,1,sd)*1.96,apply(dim_est_Tng,1,mean)+1.96*apply(dim_est_Tng,1,sd)),
traitNo=rep(c(1:27)*15,2),Subspace=c(rep("Tg",27),rep("Tng",27)))

dataTemp2<-data.frame(dim_est=c(dim_est_Tg[,1],dim_est_Tng[,1]),
traitNo=rep(c(1:27)*15,2),Subspace=c(rep("Tg",27),rep("Tng",27)))

ggplot(dataTemp2,aes(x=traitNo,y=dim_est,col=Subspace))+geom_point(cex=0.2,stroke=0.2)+geom_line(lwd=0.2)+geom_linerange(data=dataTemp,aes(ymin=lower,ymax=upper),lwd=0.2)+
scale_color_manual(labels = c(expression(italic(P^g)), expression(italic(P^ng))),values=c("#F8766D","#00BFC4"))+
scale_x_continuous(breaks=c(0,200,400),limits=c(0,405))+
scale_y_continuous(breaks=c(0,50,100),limits=c(0,100))+
labs(x="No. of traits",y=paste("No. of dimensions"),title="Yeast")+
theme(
legend.position=c(0.25,0.9),
legend.title=element_blank(),
legend.text=element_text(size=8),
legend.key.height=unit(0,"line"),
legend.key.width=unit(0.8,"line"),
legend.text.align=0,
legend.spacing.x=unit(0.2,"line"),
legend.background=element_blank())

#ggsave("yeasstSegregants_dim_est_Pg_Vs_Png.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


###


dataTemp<-data.frame(dim_est=c(apply(dim_est_Tg,1,mean)[1:22],apply(dim_est_Tng,1,mean)[1:22]),
deltaTrait=c((diff(c(1:27)*15)/diff(apply(dim_est_Tg,1,mean)))[1:22],(diff(c(1:27)*15)/diff(apply(dim_est_Tng,1,mean)))[1:22]),
subspace=c(rep("Pg",22),rep("Png",22)))

ggplot(dataTemp,aes(x=dim_est,y=deltaTrait,col=subspace))+geom_point(cex=0.5,stroke=0.5)+
scale_color_manual(labels = c(expression(italic(P^g)), expression(italic(P^ng))),values=c("#F8766D","#00BFC4"))+
scale_x_continuous(breaks=c(0,50,100),limits=c(0,100))+
scale_y_continuous(breaks=c(0,125,250))+
labs(x="No. of dimensions",y=paste("No. of traits per dimension"),title="Yeast")+
theme(
legend.position=c(0.7,0.9),
legend.title=element_blank(),
legend.text=element_text(size=8),
legend.key.height=unit(0,"line"),
legend.key.width=unit(0.8,"line"),
legend.text.align=0,
legend.spacing.x=unit(0.2,"line"),
legend.background=element_blank())

#ggsave("yeastSegregants_trait_per_dim.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


###

Tg_corMat<-cor(Pg)
Tng_corMat<-cor(Png)

#Pg

dataTemp<-data.frame(R=c(Tg_corMat[lower.tri(Tg_corMat)]))
ggplot(dataTemp,aes(x=R))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.15,alpha=0.4)+geom_vline(xintercept=c(-0.1,0.1),colour="#F8766D",lwd=0.5,stroke=0.3)+
scale_x_continuous(breaks=c(-1,0,1),limits=c(-1,1))+
labs(x=expression(paste("Pariwise correlation (",italic(R),")")),y="No. of trait pairs",title=expression(paste("Yeast ",italic(P^g))))

length(which(abs(dataTemp)>0.1))/(length(which(abs(dataTemp)>0.1))+length(which(abs(dataTemp)<0.1)))

#ggsave("yeastSegregants_pairwise_correlation_Pg.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

#Png

dataTemp<-data.frame(R=c(Tng_corMat[lower.tri(Tng_corMat)]))
ggplot(dataTemp,aes(x=R))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.15,alpha=0.4)+geom_vline(xintercept=c(-0.1,0.1),colour="#F8766D",lwd=0.5,stroke=0.3)+
scale_x_continuous(breaks=c(-1,0,1),limits=c(-1,1))+
labs(x=expression(paste("Pariwise correlation (",italic(R),")")),y="No. of trait pairs",title=expression(paste("Yeast ",italic(P^ng))))

length(which(abs(dataTemp)>0.1))/(length(which(abs(dataTemp)>0.1))+length(which(abs(dataTemp)<0.1)))

#ggsave("yeastSegregants_pairwise_correlation_Png.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



############################## Fig.4_plot_dim_estimation_brain



load("E:/papers/UBHDD/UKB/data/b2b/brain_with_genome_exclude_covariate_water_Tg_corMat.RData")
load("E:/papers/UBHDD/UKB/data/b2b/brain_with_genome_exclude_covariate_water_Tng_corMat.RData")

Tg_corMat<-brain_with_genome_exclude_covariate_water_Tg_corMat
Tng_corMat<-brain_with_genome_exclude_covariate_water_Tng_corMat

### dims of Tg and Tng

load("E:/papers/UBHDD/UKB/data/b2b/brain_with_genome_exclude_covariate_water_Tg_scale.RData")
load("E:/papers/UBHDD/UKB/data/b2b/brain_with_genome_exclude_covariate_water_Tng_scale.RData")

Tg_scale<-brain_with_genome_exclude_covariate_water_Tg_scale
Tng_scale<-brain_with_genome_exclude_covariate_water_Tng_scale


###

traitRandomInd<-foreach(i=1:100,.combine=cbind)%do%{
    set.seed(i)
    return(sample(675,675))
}

cl<-makeCluster(100,type="FORK")
registerDoParallel(cl)
	
dim_est_Tg<-t(sapply(1:45,function(i){

	
    foreach(j=1:100,.combine=c)%dopar%{
        focal_Tg<-Tg_scale[,traitRandomInd[1:(i*15),j]]
        focal_Tg_PCA<-princomp(focal_Tg)	
	    return(min(which(cumsum(focal_Tg_PCA$sdev^2)/sum(focal_Tg_PCA$sdev^2)>0.85)))
	}
}))


dim_est_Tng<-t(sapply(1:45,function(i){

    foreach(j=1:100,.combine=c)%dopar%{
	    focal_Tng<-Tng_scale[,traitRandomInd[1:(i*15),j]]
	    focal_Tng_PCA<-princomp(focal_Tng)	
	    return(min(which(cumsum(focal_Tng_PCA$sdev^2)/sum(focal_Tng_PCA$sdev^2)>0.85)))	    
	}
}))

stopCluster(cl)


#save(dim_est_Tg,dim_est_Tng,file="dim_est_Tg_Tng_brain.RData")

load("dim_est_Tg_Tng_brain.RData")

### dim estimaiton


dataTemp<-data.frame(dim_est=c(apply(dim_est_Tg,1,mean),apply(dim_est_Tng,1,mean)),
lower=c(apply(dim_est_Tg,1,mean)-apply(dim_est_Tg,1,sd)*1.96,apply(dim_est_Tng,1,mean)-1.96*apply(dim_est_Tng,1,sd)),
upper=c(apply(dim_est_Tg,1,mean)+apply(dim_est_Tg,1,sd)*1.96,apply(dim_est_Tng,1,mean)+1.96*apply(dim_est_Tng,1,sd)),
traitNo=rep(c(1:45)*15,2),Subspace=c(rep("Tg",45),rep("Tng",45)))

dataTemp2<-data.frame(dim_est=c(dim_est_Tg[,1],dim_est_Tng[,1]),
traitNo=rep(c(1:45)*15,2),Subspace=c(rep("Tg",45),rep("Tng",45)))

ggplot(dataTemp2,aes(x=traitNo,y=dim_est,col=Subspace))+geom_point(cex=0.2,stroke=0.2)+geom_line(lwd=0.2)+geom_linerange(data=dataTemp,aes(ymin=lower,ymax=upper),lwd=0.2)+
scale_color_manual(labels = c(expression(italic(P^g)), expression(italic(P^ng))),values=c("#F8766D","#00BFC4"))+
scale_x_continuous(breaks=c(0,300,600))+
scale_y_continuous(breaks=c(0,80,160))+
labs(x="No. of traits",y=paste("No. of dimensions"),title="Human brain")+
theme(
legend.position=c(0.25,0.9),
legend.title=element_blank(),
legend.text=element_text(size=8),
legend.key.height=unit(0,"line"),
legend.key.width=unit(0.8,"line"),
legend.text.align=0,
legend.spacing.x=unit(0.2,"line"),
legend.background=element_blank())

#ggsave("brain_dim_est_Pg_Vs_Png.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### Δ # of code  Vs Δ # of trait


dataTemp<-data.frame(dim_est=c(apply(dim_est_Tg,1,mean)[1:43],apply(dim_est_Tng,1,mean)[1:43]),
deltaTrait=c((diff(c(1:45)*15)/diff(apply(dim_est_Tg,1,mean)))[1:43],(diff(c(1:45)*15)/diff(apply(dim_est_Tng,1,mean)))[1:43]),
subspace=c(rep("Pg",43),rep("Png",43)))

ggplot(dataTemp,aes(x=dim_est,y=deltaTrait,col=subspace))+geom_point(cex=0.5,stroke=0.5)+
scale_color_manual(labels = c(expression(italic(P^g)), expression(italic(P^ng))),values=c("#F8766D","#00BFC4"))+
scale_x_continuous(breaks=c(0,50,170))+
scale_y_continuous(breaks=c(0,60,120))+
labs(x="No. of dimensions",y=paste("No. of traits per dimension"),title="Human brain")+
theme(
legend.position=c(0.25,0.9),
legend.title=element_blank(),
legend.text=element_text(size=8),
legend.key.height=unit(0,"line"),
legend.key.width=unit(0.8,"line"),
legend.text.align=0,
legend.spacing.x=unit(0.2,"line"),
legend.background=element_blank())


#ggsave("brain_trait_per_dim.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


###


###


#Pg

dataTemp<-data.frame(R=c(Tg_corMat[lower.tri(Tg_corMat)]))
ggplot(dataTemp,aes(x=R))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.15,alpha=0.4)+geom_vline(xintercept=c(-0.1,0.1),colour="#F8766D",lwd=0.5,stroke=0.3)+
scale_x_continuous(breaks=c(-1,0,1),limits=c(-1,1))+
labs(x=expression(paste("Pariwise correlation (",italic(R),")")),y="No. of trait pairs",title=expression(paste("Brain ",italic(P^g))))

length(which(abs(dataTemp)<0.1))/(length(which(abs(dataTemp)>0.1))+length(which(abs(dataTemp)<0.1)))


#ggsave("brain_pairwise_correlation_Pg.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#Png

dataTemp<-data.frame(R=c(Tng_corMat[lower.tri(Tng_corMat)]))
ggplot(dataTemp,aes(x=R))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.15,alpha=0.4)+geom_vline(xintercept=c(-0.1,0.1),colour="#F8766D",lwd=0.5,stroke=0.3)+
scale_x_continuous(breaks=c(-1,0,1),limits=c(-1,1))+
labs(x=expression(paste("Pariwise correlation (",italic(R),")")),y="No. of trait pairs",title=expression(paste("Brain ",italic(P^ng))))

length(which(abs(dataTemp)<0.1))/(length(which(abs(dataTemp)>0.1))+length(which(abs(dataTemp)<0.1)))


#ggsave("brain_pairwise_correlation_Png.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



############################# Fig.S4_Fig.S7_plot_unCorCutoffRobust_yeast




load("E:/papers/UBHDD/yeastSegregants/data/basic/pdata.02.RData")

load("E:/papers/UBHDD/yeastSegregants/data/oneRep/perfMat_unCorr_one_rep_20cycles.RData")
load("E:/papers/UBHDD/yeastSegregants/data/unCorCutoffRobust/perfMat_unCorr01_one_rep_20cycles.RData")
load("E:/papers/UBHDD/yeastSegregants/data/unCorCutoffRobust/perfMat_unCorr011_one_rep_20cycles.RData")
load("E:/papers/UBHDD/yeastSegregants/data/unCorCutoffRobust/perfMat_unCorr012_one_rep_20cycles.RData")
load("E:/papers/UBHDD/yeastSegregants/data/unCorCutoffRobust/perfMat_unCorr013_one_rep_20cycles.RData")
load("E:/papers/UBHDD/yeastSegregants/data/unCorCutoffRobust/perfMat_unCorr014_one_rep_20cycles.RData")
load("E:/papers/UBHDD/yeastSegregants/data/unCorCutoffRobust/perfMat_unCorr015_one_rep_20cycles.RData")
load("E:/papers/UBHDD/yeastSegregants/data/unCorCutoffRobust/perfMat_unCorr016_one_rep_20cycles.RData")
load("E:/papers/UBHDD/yeastSegregants/data/unCorCutoffRobust/perfMat_unCorr017_one_rep_20cycles.RData")
load("E:/papers/UBHDD/yeastSegregants/data/unCorCutoffRobust/perfMat_unCorr018_one_rep_20cycles.RData")
load("E:/papers/UBHDD/yeastSegregants/data/unCorCutoffRobust/perfMat_unCorr019_one_rep_20cycles.RData")
load("E:/papers/UBHDD/yeastSegregants/data/unCorCutoffRobust/perfMat_unCorr02_one_rep_20cycles.RData")


load("E:/papers/UBHDD/yeastSegregants/data/oneRep/perfMat_total_one_rep_20cycles.RData")


### calculate uncorrelation cutoff

calUnCorCutoff<-function(samSize,traitNum,alpha){
    return(sqrt(qt(1-alpha*0.5/(traitNum-1),samSize-2)^2/(samSize-2+qt(1-alpha*0.5/(traitNum-1),samSize-2)^2)))
}


### robustness of unCorCutoff

perfMat_unCorr_one_rep[which(is.na(perfMat_unCorr_one_rep))]<-0
perfMat_unCorr01_one_rep[which(is.na(perfMat_unCorr01_one_rep))]<-0
perfMat_unCorr011_one_rep[which(is.na(perfMat_unCorr011_one_rep))]<-0
perfMat_unCorr012_one_rep[which(is.na(perfMat_unCorr012_one_rep))]<-0
perfMat_unCorr013_one_rep[which(is.na(perfMat_unCorr013_one_rep))]<-0
perfMat_unCorr014_one_rep[which(is.na(perfMat_unCorr014_one_rep))]<-0
perfMat_unCorr015_one_rep[which(is.na(perfMat_unCorr015_one_rep))]<-0
perfMat_unCorr016_one_rep[which(is.na(perfMat_unCorr016_one_rep))]<-0
perfMat_unCorr017_one_rep[which(is.na(perfMat_unCorr017_one_rep))]<-0
perfMat_unCorr018_one_rep[which(is.na(perfMat_unCorr018_one_rep))]<-0
perfMat_unCorr019_one_rep[which(is.na(perfMat_unCorr019_one_rep))]<-0
perfMat_unCorr02_one_rep[which(is.na(perfMat_unCorr02_one_rep))]<-0


dataTemp<-data.frame(R2_0147=apply(perfMat_unCorr_one_rep,1,mean),
R2_01=apply(perfMat_unCorr01_one_rep,1,mean),
R2_011=apply(perfMat_unCorr011_one_rep,1,mean),
R2_012=apply(perfMat_unCorr012_one_rep,1,mean),
R2_013=apply(perfMat_unCorr013_one_rep,1,mean),
R2_014=apply(perfMat_unCorr014_one_rep,1,mean),
R2_015=apply(perfMat_unCorr015_one_rep,1,mean),
R2_016=apply(perfMat_unCorr016_one_rep,1,mean),
R2_017=apply(perfMat_unCorr017_one_rep,1,mean),
R2_018=apply(perfMat_unCorr018_one_rep,1,mean),
R2_019=apply(perfMat_unCorr019_one_rep,1,mean),
R2_02=apply(perfMat_unCorr02_one_rep,1,mean))


#0.147 Vs 0.1

ggplot(dataTemp,aes(x=R2_0147,y=R2_01))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.3)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.147")),y=expression(paste(italic(R)[u],"=0.1")),title="")+theme_classic()+
theme(panel.border=element_blank(),
plot.title=element_text(size=rel(8/11),hjust=0.5),
plot.title.position="plot",
axis.text=element_text(size=rel(8/11)),
axis.text.y=element_text(angle=90,hjust=0.5),
axis.title=element_text(size=rel(9/11)),
axis.title.x=element_text(hjust=0.5),
axis.title.y=element_text(hjust=0.5),
plot.margin=unit(c(0.05,0.15,0.15,0.15),"in"))

cor(dataTemp$R2_0147,dataTemp$R2_01)

#ggsave("yeastSegregants_Ru_robust_0147_01.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


#0.147 Vs 0.11

ggplot(dataTemp,aes(x=R2_0147,y=R2_011))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.3)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.147")),y=expression(paste(italic(R)[u],"=0.11")),title="")+theme_classic()+
theme(panel.border=element_blank(),
plot.title=element_text(size=rel(8/11),hjust=0.5),
plot.title.position="plot",
axis.text=element_text(size=rel(8/11)),
axis.text.y=element_text(angle=90,hjust=0.5),
axis.title=element_text(size=rel(9/11)),
axis.title.x=element_text(hjust=0.5),
axis.title.y=element_text(hjust=0.5),
plot.margin=unit(c(0.05,0.15,0.15,0.15),"in"))

cor(dataTemp$R2_0147,dataTemp$R2_011)

#ggsave("yeastSegregants_Ru_robust_0147_011.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


#0.147 Vs 0.12

ggplot(dataTemp,aes(x=R2_0147,y=R2_012))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.3)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.147")),y=expression(paste(italic(R)[u],"=0.12")),title="")+theme_classic()+
theme(panel.border=element_blank(),
plot.title=element_text(size=rel(8/11),hjust=0.5),
plot.title.position="plot",
axis.text=element_text(size=rel(8/11)),
axis.text.y=element_text(angle=90,hjust=0.5),
axis.title=element_text(size=rel(9/11)),
axis.title.x=element_text(hjust=0.5),
axis.title.y=element_text(hjust=0.5),
plot.margin=unit(c(0.05,0.15,0.15,0.15),"in"))

cor(dataTemp$R2_0147,dataTemp$R2_012)

#ggsave("yeastSegregants_Ru_robust_0147_012.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


#0.147 Vs 0.13

ggplot(dataTemp,aes(x=R2_0147,y=R2_013))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.3)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.147")),y=expression(paste(italic(R)[u],"=0.13")),title="")+theme_classic()+
theme(panel.border=element_blank(),
plot.title=element_text(size=rel(8/11),hjust=0.5),
plot.title.position="plot",
axis.text=element_text(size=rel(8/11)),
axis.text.y=element_text(angle=90,hjust=0.5),
axis.title=element_text(size=rel(9/11)),
axis.title.x=element_text(hjust=0.5),
axis.title.y=element_text(hjust=0.5),
plot.margin=unit(c(0.05,0.15,0.15,0.15),"in"))

cor(dataTemp$R2_0147,dataTemp$R2_013)

#ggsave("yeastSegregants_Ru_robust_0147_013.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


#0.147 Vs 0.14

ggplot(dataTemp,aes(x=R2_0147,y=R2_014))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.3)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.147")),y=expression(paste(italic(R)[u],"=0.14")),title="")+theme_classic()+
theme(panel.border=element_blank(),
plot.title=element_text(size=rel(8/11),hjust=0.5),
plot.title.position="plot",
axis.text=element_text(size=rel(8/11)),
axis.text.y=element_text(angle=90,hjust=0.5),
axis.title=element_text(size=rel(9/11)),
axis.title.x=element_text(hjust=0.5),
axis.title.y=element_text(hjust=0.5),
plot.margin=unit(c(0.05,0.15,0.15,0.15),"in"))

cor(dataTemp$R2_0147,dataTemp$R2_014)

#ggsave("yeastSegregants_Ru_robust_0147_014.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))



#0.147 Vs 0.15

ggplot(dataTemp,aes(x=R2_0147,y=R2_015))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.3)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.147")),y=expression(paste(italic(R)[u],"=0.15")),title="")+theme_classic()+
theme(panel.border=element_blank(),
plot.title=element_text(size=rel(8/11),hjust=0.5),
plot.title.position="plot",
axis.text=element_text(size=rel(8/11)),
axis.text.y=element_text(angle=90,hjust=0.5),
axis.title=element_text(size=rel(9/11)),
axis.title.x=element_text(hjust=0.5),
axis.title.y=element_text(hjust=0.5),
plot.margin=unit(c(0.05,0.15,0.15,0.15),"in"))

cor(dataTemp$R2_0147,dataTemp$R2_015)

#ggsave("yeastSegregants_Ru_robust_0147_015.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


#0.147 Vs 0.16

ggplot(dataTemp,aes(x=R2_0147,y=R2_016))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.3)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.147")),y=expression(paste(italic(R)[u],"=0.16")),title="")+theme_classic()+
theme(panel.border=element_blank(),
plot.title=element_text(size=rel(8/11),hjust=0.5),
plot.title.position="plot",
axis.text=element_text(size=rel(8/11)),
axis.text.y=element_text(angle=90,hjust=0.5),
axis.title=element_text(size=rel(9/11)),
axis.title.x=element_text(hjust=0.5),
axis.title.y=element_text(hjust=0.5),
plot.margin=unit(c(0.05,0.15,0.15,0.15),"in"))

cor(dataTemp$R2_0147,dataTemp$R2_016)

#ggsave("yeastSegregants_Ru_robust_0147_016.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


#0.147 Vs 0.17

ggplot(dataTemp,aes(x=R2_0147,y=R2_017))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.3)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.147")),y=expression(paste(italic(R)[u],"=0.17")),title="")+theme_classic()+
theme(panel.border=element_blank(),
plot.title=element_text(size=rel(8/11),hjust=0.5),
plot.title.position="plot",
axis.text=element_text(size=rel(8/11)),
axis.text.y=element_text(angle=90,hjust=0.5),
axis.title=element_text(size=rel(9/11)),
axis.title.x=element_text(hjust=0.5),
axis.title.y=element_text(hjust=0.5),
plot.margin=unit(c(0.05,0.15,0.15,0.15),"in"))

cor(dataTemp$R2_0147,dataTemp$R2_017)

#ggsave("yeastSegregants_Ru_robust_0147_017.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


#0.147 Vs 0.18

ggplot(dataTemp,aes(x=R2_0147,y=R2_018))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.3)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.147")),y=expression(paste(italic(R)[u],"=0.18")),title="")+theme_classic()+
theme(panel.border=element_blank(),
plot.title=element_text(size=rel(8/11),hjust=0.5),
plot.title.position="plot",
axis.text=element_text(size=rel(8/11)),
axis.text.y=element_text(angle=90,hjust=0.5),
axis.title=element_text(size=rel(9/11)),
axis.title.x=element_text(hjust=0.5),
axis.title.y=element_text(hjust=0.5),
plot.margin=unit(c(0.05,0.15,0.15,0.15),"in"))

cor(dataTemp$R2_0147,dataTemp$R2_018)

#ggsave("yeastSegregants_Ru_robust_0147_018.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


#0.147 Vs 0.19

ggplot(dataTemp,aes(x=R2_0147,y=R2_019))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.3)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.147")),y=expression(paste(italic(R)[u],"=0.19")),title="")+theme_classic()+
theme(panel.border=element_blank(),
plot.title=element_text(size=rel(8/11),hjust=0.5),
plot.title.position="plot",
axis.text=element_text(size=rel(8/11)),
axis.text.y=element_text(angle=90,hjust=0.5),
axis.title=element_text(size=rel(9/11)),
axis.title.x=element_text(hjust=0.5),
axis.title.y=element_text(hjust=0.5),
plot.margin=unit(c(0.05,0.15,0.15,0.15),"in"))

cor(dataTemp$R2_0147,dataTemp$R2_019)

#ggsave("yeastSegregants_Ru_robust_0147_019.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


#0.147 Vs 0.2

ggplot(dataTemp,aes(x=R2_0147,y=R2_02))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.3)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.147")),y=expression(paste(italic(R)[u],"=0.2")),title="")+theme_classic()+
theme(panel.border=element_blank(),
plot.title=element_text(size=rel(8/11),hjust=0.5),
plot.title.position="plot",
axis.text=element_text(size=rel(8/11)),
axis.text.y=element_text(angle=90,hjust=0.5),
axis.title=element_text(size=rel(9/11)),
axis.title.x=element_text(hjust=0.5),
axis.title.y=element_text(hjust=0.5),
plot.margin=unit(c(0.05,0.15,0.15,0.15),"in"))

cor(dataTemp$R2_0147,dataTemp$R2_02)

#ggsave("yeastSegregants_Ru_robust_0147_02.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))




### UBHDD test between the coef and the marginal R

## UBHDD model

UBHDDTest_unCor_coef2_R2<-matrix(nrow=405,ncol=20)
UBHDDTest_unCor_coef2_R2_p<-matrix(nrow=405,ncol=20)

for(i in 1:405){
	print(i)
    for(j in 1:20){
	    unCorInd<-which(abs(cor(pdata.02[,i],pdata.02))<pcMatCutoff)		
		tempCor2<-cor.test(abs(coefList_unCorr_one_rep[[i]][unCorInd+1,j])^2,abs(cor(pdata.02[,i],pdata.02)[unCorInd])^2)
	    UBHDDTest_unCor_coef2_R2[i,j]<-tempCor2$estimate
	    UBHDDTest_unCor_coef2_R2_p[i,j]<-tempCor2$p.value
		
	}
}


## total model

UBHDDTest_total_coef2_R2<-matrix(nrow=405,ncol=20)
UBHDDTest_total_coef2_R2_p<-matrix(nrow=405,ncol=20)

for(i in 1:405){
	print(i)
    for(j in 1:20){
	    unCorInd<-setdiff(1:405,i)
		tempCor2<-cor.test(abs(coefList_total_one_rep[[i]][unCorInd+1,j])^2,abs(cor(pdata.02[,i],pdata.02)[unCorInd])^2)
	    UBHDDTest_total_coef2_R2[i,j]<-tempCor2$estimate
	    UBHDDTest_total_coef2_R2_p[i,j]<-tempCor2$p.value
		
	}
}


#save(UBHDDTest_unCor_coef2_R2,UBHDDTest_unCor_coef2_R2_p,UBHDDTest_total_coef2_R2,UBHDDTest_total_coef2_R2_p,file="yeast_UBHDDTest_UBHDD_Vs_total.RData")


load("yeast_UBHDDTest_UBHDD_Vs_total.RData")

## ggplot

#unCor

i=3
j=20

unCorInd<-which(abs(cor(pdata.02[,i],pdata.02))<pcMatCutoff)
dataTemp<-data.frame(size_of_coef=abs(coefList_unCorr_one_rep[[i]][unCorInd+1,j])^2,size_of_marginalR=abs(cor(pdata.02[,i],pdata.02)[unCorInd])^2)

ggplot(dataTemp,aes(x=size_of_marginalR,y=size_of_coef))+geom_point(cex=0.2,alpha=1,stroke=0.2)+
scale_x_continuous(breaks=c(0,0.02),limits=c(0,0.022))+
scale_y_continuous(breaks=c(0,2.8),limits=c(0,4))+
labs(x=expression(italic(MC^2)),y=expression(italic(Co^2)),title="UBHDD")

#ggsave("yeastSegregants_UBHDD_Ru_test_case.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#

dataTemp<-data.frame(R2=apply(UBHDDTest_unCor_coef2_R2^2,1,mean_replace))
max(dataTemp)
ggplot(dataTemp,aes(x=R2))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.15,alpha=0.4)+
scale_x_continuous(breaks=c(0,0.08,0.16),limits=c(0,0.16))+
scale_y_continuous(breaks=c(0,40,80),limits=c(0,85))+
labs(x=expression(italic(R^2)),y="No. of traits",title="UBHDD")

#ggsave("yeastSegregants_UBHDD_Ru_test_distribution.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

#

dataTemp<-data.frame(mlogp=apply(-log10(UBHDDTest_unCor_coef2_R2_p),1,mean_replace))
max(dataTemp)
ggplot(dataTemp,aes(x=mlogp))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.15,alpha=0.4)+geom_vline(xintercept=-log10(0.01/405),col="red")+
scale_x_continuous(breaks=c(0,7.5,15),limits=c(0,15.4))+
scale_y_continuous(breaks=c(0,100,200),limits=c(0,200))+
labs(x=expression(paste("-log10(",italic(p),")")),y="No. of traits",title="UBHDD")

#ggsave("yeastSegregants_UBHDD_Ru_test_pvalue.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


# total

i=3
j=20

unCorInd<-setdiff(1:405,i)
dataTemp<-data.frame(size_of_coef=abs(coefList_total_one_rep[[i]][unCorInd+1,j])^2,size_of_marginalR=abs(cor(pdata.02[,i],pdata.02)[unCorInd])^2)
max(dataTemp[,1])
max(dataTemp[,2])

ggplot(dataTemp,aes(x=size_of_marginalR,y=size_of_coef))+geom_point(cex=0.2,alpha=1,stroke=0.2)+
scale_x_continuous(breaks=c(0,0.2,0.4),limits=c(0,0.42))+
scale_y_continuous(breaks=c(0,1.5,0.3),limits=c(0,0.35))+
labs(x=expression(italic(MC^2)),y=expression(italic(Co^2)),title="Total")

#ggsave("yeastSegregants_TOTAL_Ru_test_case.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#

dataTemp<-data.frame(R2=apply(UBHDDTest_total_coef2_R2^2,1,mean_replace))
max(dataTemp)
ggplot(dataTemp,aes(x=R2))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.15,alpha=0.4)+
scale_x_continuous(breaks=c(0,0.4,0.8),limits=c(0,0.81))+
scale_y_continuous(breaks=c(0,35,70),limits=c(0,75))+
labs(x=expression(italic(R^2)),y="No. of traits",title="Total")

#ggsave("yeastSegregants_TOTAL_Ru_test_distribution.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

#

dataTemp<-data.frame(mlogp=apply(-log10(UBHDDTest_total_coef2_R2_p),1,mean_replace))
max(dataTemp)
ggplot(dataTemp,aes(x=mlogp))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.15,alpha=0.4)+geom_vline(xintercept=-log10(0.01/405),col="red")+
scale_x_continuous(breaks=c(0,70,140),limits=c(0,147))+
scale_y_continuous(breaks=c(0,50,100),limits=c(0,115))+
labs(x=expression(paste("-log10(",italic(p),")")),y="No. of traits",title="Total")

#ggsave("yeastSegregants_TOTAL_Ru_test_pvalue.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


## test for case

UBHDDTest_unCor_coef2_R2[3,20]^2
UBHDDTest_unCor_coef2_R2_p[3,20]

UBHDDTest_total_coef2_R2[3,20]^2
UBHDDTest_total_coef2_R2_p[3,20]






############################# Fig.S5_Fig.S7_plot_unCorCutoffRobust_brain



load("E:/papers/UBHDD/UKB/data/b2b/brain_with_genome_exclude_covariate_water_T_corMat.RData")

prediction01<-read.csv("E:/papers/UBHDD/UKB/data/unCorCutoffRobust/result_b2b_water_01_DF.csv")
prediction011<-read.csv("E:/papers/UBHDD/UKB/data/unCorCutoffRobust/result_b2b_water_011_DF.csv")
prediction012<-read.csv("E:/papers/UBHDD/UKB/data/unCorCutoffRobust/result_b2b_water_012_DF.csv")
prediction013<-read.csv("E:/papers/UBHDD/UKB/data/unCorCutoffRobust/result_b2b_water_013_DF.csv")
prediction014<-read.csv("E:/papers/UBHDD/UKB/data/unCorCutoffRobust/result_b2b_water_014_DF.csv")
prediction015<-read.csv("E:/papers/UBHDD/UKB/data/unCorCutoffRobust/result_b2b_water_015_DF.csv")
prediction016<-read.csv("E:/papers/UBHDD/UKB/data/unCorCutoffRobust/result_b2b_water_016_DF.csv")
prediction017<-read.csv("E:/papers/UBHDD/UKB/data/unCorCutoffRobust/result_b2b_water_017_DF.csv")
prediction018<-read.csv("E:/papers/UBHDD/UKB/data/unCorCutoffRobust/result_b2b_water_018_DF.csv")
prediction019<-read.csv("E:/papers/UBHDD/UKB/data/unCorCutoffRobust/result_b2b_water_019_DF.csv")
prediction02<-read.csv("E:/papers/UBHDD/UKB/data/unCorCutoffRobust/result_b2b_water_02_DF.csv")

prediction1<-read.csv("E:/papers/UBHDD/UKB/data/b2b/result_b2b_water_1_DF.csv")






###################################

dataTemp<-data.frame(R2_01=prediction01[,2],
R2_011=prediction011[,2],
R2_012=prediction012[,2],
R2_013=prediction013[,2],
R2_014=prediction014[,2],
R2_015=prediction015[,2],
R2_016=prediction016[,2],
R2_017=prediction017[,2],
R2_018=prediction018[,2],
R2_019=prediction019[,2],
R2_02=prediction02[,2])


###0.15 Vs 0.1

ggplot(dataTemp,aes(x=R2_015,y=R2_01))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.1)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.15")),y=expression(paste(italic(R)[u],"=0.1")),title="")


#ggsave("brain_Ru_robust_015_01.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))

###0.15 Vs 0.11

ggplot(dataTemp,aes(x=R2_015,y=R2_011))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.1)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.15")),y=expression(paste(italic(R)[u],"=0.11")),title="")

#ggsave("brain_Ru_robust_015_011.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))

###0.15 Vs 0.12

ggplot(dataTemp,aes(x=R2_015,y=R2_012))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.1)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.15")),y=expression(paste(italic(R)[u],"=0.12")),title="")

#ggsave("brain_Ru_robust_015_012.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


###0.15 Vs 0.13

ggplot(dataTemp,aes(x=R2_015,y=R2_013))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.1)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.15")),y=expression(paste(italic(R)[u],"=0.13")),title="")

#ggsave("brain_Ru_robust_015_013.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))



###0.15 Vs 0.14

ggplot(dataTemp,aes(x=R2_015,y=R2_014))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.1)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.15")),y=expression(paste(italic(R)[u],"=0.14")),title="")

#ggsave("brain_Ru_robust_015_014.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))



###0.15 Vs 0.16

ggplot(dataTemp,aes(x=R2_015,y=R2_016))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.1)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.15")),y=expression(paste(italic(R)[u],"=0.16")),title="")

#ggsave("brain_Ru_robust_015_016.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))



###0.15 Vs 0.17

ggplot(dataTemp,aes(x=R2_015,y=R2_017))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.1)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.15")),y=expression(paste(italic(R)[u],"=0.17")),title="")

#ggsave("brain_Ru_robust_015_017.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


###0.15 Vs 0.18

ggplot(dataTemp,aes(x=R2_015,y=R2_018))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.1)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.15")),y=expression(paste(italic(R)[u],"=0.18")),title="")

#ggsave("brain_Ru_robust_015_018.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


###0.15 Vs 0.19

ggplot(dataTemp,aes(x=R2_015,y=R2_019))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.1)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.15")),y=expression(paste(italic(R)[u],"=0.19")),title="")

#ggsave("brain_Ru_robust_015_019.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


###0.15 Vs 0.2

ggplot(dataTemp,aes(x=R2_015,y=R2_02))+geom_abline(lwd=0.5,color="black",intercept=0,slope=1)+geom_point(cex=0.2,alpha=1,stroke=0.1)+
scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x=expression(paste(italic(R)[u],"=0.15")),y=expression(paste(italic(R)[u],"=0.2")),title="")

#ggsave("brain_Ru_robust_015_02.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))



####################################
############################# compare UBHDD with Total models ##############################################

## UBHDD model


UBHDDTest_unCor_coef2_R2<-c()
UBHDDTest_unCor_coef2_R2_p<-c()

for(i in 1:675){
	print(i)
    unCorInd<-which(abs(brain_with_genome_exclude_covariate_water_T_corMat[i,])<0.15)
	tempCor2<-cor.test(abs(as.numeric(prediction015[i,5:679])[unCorInd])^2,abs(as.numeric(brain_with_genome_exclude_covariate_water_T_corMat[i,unCorInd]))^2)
    UBHDDTest_unCor_coef2_R2[i]<-tempCor2$estimate
    UBHDDTest_unCor_coef2_R2_p[i]<-tempCor2$p.value
		
}



## total model

UBHDDTest_total_coef2_R2<-c()
UBHDDTest_total_coef2_R2_p<-c()


for(i in 1:675){
	print(i)
    unCorInd<-setdiff(1:675,i)
	tempCor2<-cor.test(abs(as.numeric(prediction1[i,5:679])[unCorInd])^2,abs(as.numeric(brain_with_genome_exclude_covariate_water_T_corMat[i,unCorInd]))^2)
    UBHDDTest_total_coef2_R2[i]<-tempCor2$estimate
    UBHDDTest_total_coef2_R2_p[i]<-tempCor2$p.value
		
}


#save(UBHDDTest_unCor_coef2_R2,UBHDDTest_unCor_coef2_R2_p,UBHDDTest_total_coef2_R2,UBHDDTest_total_coef2_R2_p,file="brain_UBHDDTest_UBHDD_Vs_total.RData")

load("brain_UBHDDTest_UBHDD_Vs_total.RData")

###


## UBHDD case

i=6

prediction015[i,2]

unCorInd<-which(abs(brain_with_genome_exclude_covariate_water_T_corMat[i,])<0.15)
dataTemp<-data.frame(size_of_coef=abs(as.numeric(prediction015[i,5:679])[unCorInd])^2,size_of_marginalR=abs(as.numeric(brain_with_genome_exclude_covariate_water_T_corMat[i,unCorInd]))^2)

ggplot(dataTemp,aes(x=size_of_marginalR,y=size_of_coef))+geom_point(cex=0.2,alpha=1,stroke=0.2)+
scale_x_continuous(breaks=c(0,0.02),limits=c(0,0.025))+
scale_y_continuous(breaks=c(0,0.03),limits=c(0,0.05))+
labs(x=expression(italic(MC^2)),y=expression(italic(C^2)),title="UBHDD")

#ggsave("brain_UBHDDTest_UBHDD_case.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#

dataTemp<-data.frame(R2=UBHDDTest_unCor_coef2_R2^2)
max(dataTemp)
ggplot(dataTemp,aes(x=R2))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.15,alpha=0.4)+
scale_x_continuous(breaks=c(0,0.09,0.18),limits=c(0,0.18))+
scale_y_continuous(breaks=c(0,50,100),limits=c(0,115))+
labs(x=expression(italic(R^2)),y="No. of traits",title="UBHDD")

#ggsave("brain_UBHDDTest_UBHDD_distribution.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

#

dataTemp<-data.frame(mlogp=-log10(UBHDDTest_unCor_coef2_R2_p))
max(dataTemp)
ggplot(dataTemp,aes(x=mlogp))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.15,alpha=0.4)+geom_vline(xintercept=-log10(0.01/675),col="red")+
scale_x_continuous(breaks=c(0,10,20),limits=c(0,20.2))+
scale_y_continuous(breaks=c(0,75,150),limits=c(0,155))+
labs(x=expression(paste("-log10(",italic(p),")")),y="No. of traits",title="UBHDD")

#ggsave("brain_UBHDDTest_UBHDD_pvalue.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


# total

i=6

prediction1[i,2]

unCorInd<-setdiff(1:675,i)
dataTemp<-data.frame(size_of_coef=abs(as.numeric(prediction1[i,5:679])[unCorInd])^2,size_of_marginalR=abs(as.numeric(brain_with_genome_exclude_covariate_water_T_corMat[i,unCorInd]))^2)
max(dataTemp[,1])
max(dataTemp[,2])

ggplot(dataTemp,aes(x=size_of_marginalR,y=size_of_coef))+geom_point(cex=0.2,alpha=1,stroke=0.2)+
scale_x_continuous(breaks=c(0,0.9),limits=c(0,0.9))+
scale_y_continuous(breaks=c(0,1.1),limits=c(0,1.1))+
labs(x=expression(italic(MC^2)),y=expression(italic(C^2)),title="Total")

#ggsave("brain_UBHDDTest_TOTAL_case.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#

dataTemp<-data.frame(R2=UBHDDTest_total_coef2_R2^2)
max(dataTemp)
ggplot(dataTemp,aes(x=R2))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.15,alpha=0.4)+
scale_x_continuous(breaks=c(0,0.3,0.6),limits=c(0,0.66))+
scale_y_continuous(breaks=c(0,60,120),limits=c(0,130))+
labs(x=expression(italic(R^2)),y="No. of traits",title="Total")

#ggsave("brain_UBHDDTest_TOTAL_distribution.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

#

dataTemp<-data.frame(mlogp=-log10(UBHDDTest_total_coef2_R2_p))
max(dataTemp)
ggplot(dataTemp,aes(x=mlogp))+geom_histogram(col="grey",fill="grey",bins=20,boundary=0.15,alpha=0.4)+geom_vline(xintercept=-log10(0.01/675),col="red")+
scale_x_continuous(breaks=c(0,75,150),limits=c(0,157))+
scale_y_continuous(breaks=c(0,100,200),limits=c(0,210))+
labs(x=expression(paste("-log10(",italic(p),")")),y="No. of traits",title="Total")

#ggsave("brain_UBHDDTest_TOTAL_pvalue.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))




############################# Fig.S6_plot_down_sampling_yeast




load("E:/papers/UBHDD/yeastSegregants/data/downSampling/perfList_unCorr_one_rep_down_sampling_1cycle.RData")


###

dataTemp<-data.frame(R2=do.call("rbind",perfList_unCorr_one_rep_down_sampling),
percent=factor(c(rep("10",round(405*0.1)),rep("20",round(405*0.2)),rep("30",round(405*0.3)),rep("40",round(405*0.4)),rep("50",round(405*0.5)),rep("60",round(405*0.6)),rep("70",round(405*0.7)),rep("80",round(405*0.8)),rep("90",round(405*0.9)),rep("100",405)),levels=as.character(c(1:10)*10)))
dataTemp[which(is.na(dataTemp[,1])),1]<-0
ggplot(dataTemp,aes(y=R2,x=percent))+geom_boxplot(width=0.3, color="black", alpha=1,show.legend=FALSE,outlier.size=0.1,outlier.stroke = 0.3,size=0.3)+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x="Percent of traits",y=expression(italic(R^2)),title="Yeast")


#ggsave("yeastSegregants_UBHDD_down_sampling.pdf",width=unit(3,"in"),height=unit(1.7,"in"))





############################ Fig.S6_plot_down_sampling_brain





prediction<-read.csv("E:/papers/UBHDD/UKB/data/downSampling/result_b2b_water_015_DF.csv")
prediction_down01<-read.csv("E:/papers/UBHDD/UKB/data/downSampling/result_b2b_water_015_down01_DF.csv")
prediction_down02<-read.csv("E:/papers/UBHDD/UKB/data/downSampling/result_b2b_water_015_down02_DF.csv")
prediction_down03<-read.csv("E:/papers/UBHDD/UKB/data/downSampling/result_b2b_water_015_down03_DF.csv")
prediction_down04<-read.csv("E:/papers/UBHDD/UKB/data/downSampling/result_b2b_water_015_down04_DF.csv")
prediction_down05<-read.csv("E:/papers/UBHDD/UKB/data/downSampling/result_b2b_water_015_down05_DF.csv")
prediction_down06<-read.csv("E:/papers/UBHDD/UKB/data/downSampling/result_b2b_water_015_down06_DF.csv")
prediction_down07<-read.csv("E:/papers/UBHDD/UKB/data/downSampling/result_b2b_water_015_down07_DF.csv")
prediction_down08<-read.csv("E:/papers/UBHDD/UKB/data/downSampling/result_b2b_water_015_down08_DF.csv")
prediction_down09<-read.csv("E:/papers/UBHDD/UKB/data/downSampling/result_b2b_water_015_down09_DF.csv")


###


perfList_down_sampling<-list(X10=prediction_down01[,2],
X20=prediction_down02[,2],
X30=prediction_down03[,2],
X40=prediction_down04[,2],
X50=prediction_down05[,2],
X60=prediction_down06[,2],
X70=prediction_down07[,2],
X80=prediction_down08[,2],
X90=prediction_down09[,2],
X100=prediction[,2])

boxplot(perfList_down_sampling)


dataTemp<-data.frame(R2=do.call("c",perfList_down_sampling),
percent=factor(c(rep("10",round(675*0.1)),rep("20",round(675*0.2)),rep("30",round(675*0.3)),rep("40",round(675*0.4)),rep("50",round(675*0.5)),rep("60",round(675*0.6)),rep("70",round(675*0.7)),rep("80",round(675*0.8)),rep("90",round(675*0.9)),rep("100",675)),levels=as.character(c(1:10)*10)))

dataTemp[which(is.na(dataTemp[,1])),1]<-0
ggplot(dataTemp,aes(y=R2,x=percent))+geom_boxplot(width=0.3, color="black", alpha=1,show.legend=FALSE,outlier.size=0.1,outlier.stroke = 0.3,size=0.3)+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
labs(x="Percent of traits",y=expression(italic(R^2)),title="Brain")


#ggsave("brain_UBHDD_down_sampling.pdf",width=unit(3,"in"),height=unit(1.7,"in"))










