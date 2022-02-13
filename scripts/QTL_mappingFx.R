
# 02/13/22
# R script in:
# The trait coding rule in phenotype space
# 
# Jianguo Wang (wangjg36@mail.sysu.edu.cn)
# 

The original script is provided by :

# 03/19/13
# R script for heritability and QTL analysis in:
# Finding the sources of missing heritability in a yeast cross
# http://www.nature.com/nature/journal/v494/n7436/full/nature11867.html
# Joshua Bloom (jbloom@princeton.edu)
# This work is licensed under GPL-2.



# functions for QTL Mapping

# Get size of objects in workspace 
getObjS      <- function(y){  sort(sapply(y, function(x) { object.size(get(x))/1e6} ), decreasing=T)}


##### Extract genotype matrix from cross structure and recode as -1,1 #####################################
extractGenotype=function(impcross){ (do.call('cbind', sapply(impcross$geno, function(x) { x$data }))*2)-3 }
###########################################################################################################

##### Count number of strains with data for each phenotype from cross structure ###########################
countStrainsPerTrait = function(pheno) {apply(pheno, 2, function(x){sum(!is.na(x))})}
###########################################################################################################



######  Functions to calculate Broad-Sense Heritabiltiy ####################################################
calcH2 = function(p) {
     sp=stack(p)  
     reffMod = lmer(values ~ 1 + (1|ind),data=sp)
     #Vcomp = as.numeric((attributes(summary(reffMod)))$REmat[,'Variance'])
     Vcomp=as.numeric(data.frame(VarCorr(reffMod))[,4])
     VarG = Vcomp[1]
     VarE = Vcomp[2]
	 ##this needs to be confirmed :H2 = VarG/(VarE/2+VarG)
     H2 = VarG/(VarE+VarG)
     return(H2)
}

calc.BroadH2 = function(pdata, jackknife=FALSE, getcnt=FALSE) {
    # print(names(pd)[counter])
    #     counter <<- counter+1
         tryCatch({
         pdata.notNAcnt = sapply(pdata, function(x){sum(!is.na(x))})
         pdata[pdata.notNAcnt<2]=NULL
         pdata.length=length(pdata)
         if(getcnt==TRUE) {return(pdata.length)}
         if(jackknife==TRUE){  
             jackp= unlist( foreach(i = 1:pdata.length) %dopar% {calcH2(pdata[-i])})
             p = calcH2(pdata)
             return(calc.jkSE(pdata.length, jackp, p))
        } else { return(calcH2(pdata))}},error=function(e){return(NA)})
}

##################################################################################################################



##### Extract phenotype matrix from cross structure and mean center ... optional standardize Variance######
extractScaledPhenotype=function(impcross,scaleVar=FALSE){apply(impcross$pheno, 2, scale, scale=scaleVar)}
###########################################################################################################



#get first element of each replicate that isn't NA ################################################################
getOnePhenoValue=function(pheno_raw, pdata.01) {
    pdata.02 = lapply(pheno_raw, function(x) {
             sapply(x, function(y) { y[!is.na(y)][1] })})
    pdata.holder = matrix(NA,dim(pdata.01)[1], dim(pdata.01)[2])
    colnames(pdata.holder) = colnames(pdata.01)
    rownames(pdata.holder) = totalStrain
    for(i in 1:dim(pdata.01)[2]) {
       pdata.holder[names(pdata.02[[i]]),i]= as.vector(pdata.02[[i]]) }
    pdata.02 = pdata.holder
    return(pdata.02)
}
####################################################################################################################


#### Split marker index by chromosome #################################################################
getMarkerIndexSplit = function(impcross) {
    mr = scanone(impcross, pheno.col=1, method='mr')
    marker.info       = parse.rownames(rownames(mr))
    marker.info$cmPOS = mr$pos
    mindex.split = split(rownames(marker.info), marker.info$NCHR)
    mindex.split = sapply(mindex.split, as.numeric)
    return(mindex.split)
}
#######################################################################################################


# compute standard error of heritability using jackknife ...very very slow -------------------

calc.jkSE = function(l,jackp, p) { sqrt( ((l-1)/l) *sum((p-jackp)^2)) }
calc.jkbias = function(l, jackp, p) (l-1)*(mean(jackp)-p)

calc.mixed.jack.pseudovalues = function(p, A) {
    jackh2s=list()
    for(t in 1:ncol(p)){
        pna = !is.na(p[,t])
        pdata.length=sum(pna)
        y=p[pna,t]
        Ass=A[pna,pna]
        jackh2s[[t]] = foreach(i=1:pdata.length, .combine=c) %dopar% {
            x=mixed.solve(y[-i], K=Ass[-i,-i], method='REML')
            h=x$Vu/(x$Ve+x$Vu)
		    print(paste(colnames(pdata.01)[t],  i, h))
			
            
            return(h)        }  }
    return(jackh2s)
}
# -----------------------------------------------------------------------------------------


###### calculate LOD score given matrices of phenotype and genotype data ##############################
# n.pheno is count of strains with data
# pheno is matrix (strains X phenotypes)
# gdata is matrix (strains X markers)
#.....  be careful about NAs and only using it for single marker scans 

get.LOD.by.COR = function(n.pheno, pheno, gdata) {
    # Lynch and Walsh p. 454
    return( (-n.pheno*log(1-cor(pheno, gdata, use='pairwise.complete.obs')^2))/(2*log(10)) ) }

#######################################################################################################

LODmatrix.2.scanone= function(LODS, cross, LL=NULL) {
    if(is.null(LL)){
        LODSm = t(as.matrix(LODS))
        LODSs = scanone(cross, pheno.col=3, method='mr')
        LODSso = data.frame(LODSs, LODSm)
        LODSso= LODSso[,-3]
        class(LODSso)=class(LODSs)
        return(LODSso)
    } else { 
        LODSso =  data.frame(LL[,c(1,2)],t(as.matrix(LODS)))
        class(LODSso)  = class(LL)
        return(LODSso)
    }
}



#### given reference FASTA file ... get chromosome sizes and convert #########################
# to genome coordinate conversion key if genocoord = TRUE 
getCoordKey = function(ref.fasta, genocoord=TRUE) {
    s288c=read.fasta(ref.fasta)
    s288c.key=sapply(s288c, length)
    if (genocoord==TRUE ) {
        gcoord.key = c(0, cumsum(s288c.key)[-length(s288c.key)])
        names(gcoord.key)=c(names(gcoord.key)[-1], names(s288c.key)[length(s288c.key)])
        return(gcoord.key)    }
    else (return(s288c.key) )  } 
###############################################################################################



####load bioconductor yeast annotation information and extract gene information ###########
getORFs= function(gcoord.key)     {
    orfCommon    =  org.Sc.sgdGENENAME
    orfchrStart  =  org.Sc.sgdCHRLOC
    orfchrEnd    =  org.Sc.sgdCHRLOCEND
    orfchr       =  org.Sc.sgdCHR
    orfDescription = org.Sc.sgdDESCRIPTION

    mapped_genes = mappedkeys(orfchr)
    orf.orf   = mapped_genes
    orf.gn    = as.character(as.list(orfCommon[mapped_genes]))
    orf.chr   = as.numeric(as.list(orfchr[mapped_genes]))
    orf.start = as.numeric(as.list(orfchrStart[mapped_genes]))
    orf.end   = as.numeric(as.list(orfchrEnd[mapped_genes]))
    orf.description = as.character(as.list(orfDescription[mapped_genes]))
    orfs.raw = data.frame(orf.orf, orf.gn, orf.chr, 
                          orf.start, orf.end, orf.description, stringsAsFactors=FALSE)
    names(orfs.raw)=c('ORF','COMMON','CHR', 'START', 'END', 'DESCRIPTION')

    orfs = orfs.raw
    orfs$CHR = names(gcoord.key)[orfs$CHR]
    orfs=orfs[!is.na(orfs$CHR),]
    orfs=orfs[!is.na(orfs$START),]

    # construct gene intervals 
    orfs$START = abs(orfs$START)
    orfs$END   = abs(orfs$END)
    rightorder= orfs$START<orfs$END
    orfs = orfs[rightorder,]
    orfs=orfs[order(orfs$CHR, orfs$START),]
    orfs$GSTART=gcoord.key[orfs$CHR]+orfs$START
    orfs$GEND=gcoord.key[orfs$CHR]+orfs$END
    return(orfs) }
###############################################################################################



#takes care of filling in NAs if missing values for a given strain 
getSegMedianTable = function(phenoData) {
    # without normalization for growth in YPD 10/04/11
    # segm = sapply(phenoData, function(x) { x$segs.median })
    segm = sapply(phenoData, function(x) {x$segs.all.norm.median})
    seg.data = 
        data.frame(t(do.call(rbind, lapply(lapply(segm, unlist), "[",
                    unique(unlist(c(sapply(segm,names))))))))
    return(seg.data) }
getWIMedianTable = function(phenoData) {
    # without normalization for growth in YPD 10/04/11
    # segm = sapply(phenoData, function(x) { x$segs.median })
    segm = sapply(phenoData, function(x) {x$notsegs.all.norm.median})
    seg.data = 
        data.frame(t(do.call(rbind, lapply(lapply(segm, unlist), "[",
                    unique(unlist(c(sapply(segm,names))))))))
    return(seg.data) }
###############################################################################################

 getTable = function(phenoData, value='segs.mean') {
    
    # without normalization for growth in YPD 10/04/11
    # segm = sapply(phenoData, function(x) { x$segs.median })
    
    # want conditions as colnames and segs as rownames
    segm = sapply(phenoData, function(x) {x[value]})
    
    segnames= sort(unique(as.vector(unlist(sapply(segm, names)))))
    
    segd= sapply(segm, function(x) {
        missingStrains = segnames[!segnames %in% names(x)]
        y=x
        if(length(missingStrains)>0) {
            missingD = rep(NA, length(missingStrains)) 
            names(missingD)=missingStraig
            y=c(x,missingD)
        }
        y=y[order(names(y))]
        return(y)  
    })  
    colnames(segd)=gsub(paste('.',value, sep=''),'',colnames(segd))
    return(segd)
}


#Hasman-Elston regression
calch2 = function(pdistM, AmatM) {
    Amat.vec  = as.vector(AmatM)
    pdist.vec = as.vector(pdistM)
    sx = sd(Amat.vec, na.rm=T)
    sy = sd(pdist.vec, na.rm=T)
    h2 = -(cor(Amat.vec/2, pdist.vec, use='complete.obs')*(sy/sx))
    return(h2)
}


peakfinder2D = function(i2d, threshold) {
    #   i2d   = iLODsh[trait,,]
        ipeaks  = which(i2d>threshold, arr.ind=T)
        ipeaks = cbind(ipeaks, i2d[ipeaks])
        ipeaks = cbind(ipeaks, rep(NA, nrow(ipeaks)))
        ipeaks = cbind(ipeaks, ipeaks[,3])
        colnames(ipeaks) = c('x', 'y', 'lod', 'group', 'glod')
        g=1
        while(sum(is.na(ipeaks[,'group']))>0) {
            peak = which.max(ipeaks[,'glod'])
            pdist = sqrt(abs(ipeaks[peak,'x']-ipeaks[,'x'])^2 + abs(ipeaks[peak,'y']-ipeaks[,'y'])^2)
            gpeaks = which(pdist<50 & is.na(ipeaks[,'group']))
            if(length(gpeaks)>2) {
               ipeaks[gpeaks, 'group']=g
               ipeaks[gpeaks, 'glod']=0
               g=g+1
               print(g)
            } else{
               ipeaks[gpeaks, 'group']=0
               ipeaks[gpeaks, 'glod']=0
            }
        }
        ipeaks=ipeaks[ipeaks[,'group']>0,]
        ips = split(data.frame(ipeaks), ipeaks[,'group'])
        ipsp=data.frame(t(sapply(ips, function(x) {
                lmax= which.max(x$lod)
                x[lmax, c(1,2,3,4)]
                })))
        ipsp=        apply(ipsp, 2, as.numeric)
        return(ipsp)
}



# Peak finding for full 2D scan 
getIntPeaks = function(iL, tnames) {
    rownames(iL)=tnames
    adim = dim(iL)
    ipeaks= apply(iL, 1 , peakfinder2D, 3)
    #    x=ipeaks$Magnesium_Sulfate    
    #x=x[order(x[,'lod'],decreasing=T),]
    #plot(x[,'x'], x[,'y'], type='n', xlim=c(0,2800), ylim=c(0,2800))
    #text(x[,'x'], x[,'y'], x[,'group'], pch='.', xlim=c(0,2800), ylim=c(0,2800))
    #abline(0,1)
}





# Scan for interacting QTLs .... restricted 2D search to main effects vs rest of the genome

scanIntwMain = function(fQTLs_FDR05, n.pheno, pcross2, g.pcross.data, 
                            imindex.split, ichr.mindex.offset, doFDR=FALSE, doGWER=FALSE){

intScans = list()
for(t in 1:length(fQTLs_FDR05)){
    print(names(fQTLs_FDR05)[t])
    peaks = fQTLs_FDR05[[t]]$mqtl
    peak.markers=peaks$name
    gms =  gdata[,peak.markers]
    gint =list()
    for(g in 1:length(peak.markers)) {
          gintm=g.pcross.data*gms[,g]   
          # remove 50 cm around each peak
          for (h in g:1) {
              left.m  = find.marker(pcross2, peaks$chr[h], peaks$pos[h] - 30)
              right.m = find.marker(pcross2, peaks$chr[h], peaks$pos[h] + 30)
              print(paste(h, left.m, right.m))
              mNA = match(c(left.m, right.m), colnames(gintm))
              gintm[,mNA[1]:mNA[2]]=NA
          }
          gint[[g]]=gintm
    }
    intLODs.c = matrix(NA, ncol(gms), ncol(g.pcross.data))
    for(g in 1:ncol(gms)){ intLODs.c[g,]=get.LOD.by.COR(n.pheno[t], pcross2$pheno[,t],gint[[g]]) }
    intLODs.c[is.na(intLODs.c)]=0
    intLODs = LODmatrix.2.scanone(intLODs.c, pcross2)
    names(intLODs)[-c(1,2)] = colnames(gms)
   
    if(doFDR) {
        peaklist.01   = getChrPeaks(imindex.split, ichr.mindex.offset, intLODs.c) 
            #LODS.01.FDR   = getPeakFDR(peaklist.01$chr.peaks.lod, pdata.01, gdata, 1000)

            pintLODs.c =  foreach(i = 1:100 ) %dopar% { 
                if(i%%10==0) print(i)
                intLODs.cp = matrix(NA, ncol(gms), ncol(g.pcross.data))
                for(g in 1:ncol(gms)){ 
                    intLODs.cp[g,]=get.LOD.by.COR(n.pheno[t], pdata.05[sample(1:nrow(pdata.05)),t], gint[[g]] ) }
                sapply(imindex.split, function(markers) {apply(intLODs.cp[,markers],1,max,na.rm=T)}) 
            }
            
            permpeakLODs= abind(pintLODs.c, along=c(3))
            obsPcnt = sapply(seq(2, 5, .01), function(thresh) { sum(peaklist.01$chr.peaks.lod>thresh) }   )
            names(obsPcnt) = seq(2,5, .01)
            
            # expected number of QTL peaks with LOD greater than threshold
            expPcnt = sapply(seq(2, 5, .01),  
                     function(thresh) { 
                            print(thresh); 
                            mean(apply(permpeakLODs, 3, function(ll) {sum(ll>thresh) }) )
                        } )
            names(expPcnt) = seq(2, 5, .01)
            pFDR = expPcnt/obsPcnt
            pFDR[!is.finite(pFDR)]=NA
            attr(intLODs, 'pFDR')   = pFDR
        }

        if(doGWER) {
            GWER.0 = quantile( sapply(pintLODs.c, function(x) { sort(x, decreasing=T)[1]}) , .95)
            GWER.1 = quantile( sapply(pintLODs.c, function(x) { 
                     max(apply(x, 1, function(z) {sort(z, decreasing=TRUE)[2] } ))}) , .95)
           attr(intLODs, 'GWER.0') = GWER.0
           attr(intLODs, 'GWER.1') = GWER.1
        }
        if(keepold ) {
               attr(intLODs, 'pFDR') = attr(intScans[[t]],     'pFDR')
               attr(intLODs, 'GWER.0') = attr(intScans[[t]], 'GWER.0')
               attr(intLODs, 'GWER.1') = attr(intScans[[t]], 'GWER.1')
         }
       attr(intLODs, 'Main Peaks')=peak.markers
       intScans[[t]]=intLODs
    }
    names(intScans)=names(fQTLs_FDR05)
    return(intScans)
}




####### Merge data frame of segregant phenotype values with cross data structure ##############
mergePhenosWithCross = function(cross,seg.data) {
    cross$pheno = merge(cross$pheno, seg.data,  by.x='id', by.y='row.names')
    rownames(cross$pheno)=cross$id
    cross$pheno = cross$pheno[,-1]
    return(cross)  }
#################################################################################################


####### for converting pearson R to -log10(p) ###################################################
R.to.T   = function(r,n)  {    return ( r*sqrt((n-2)/(1-r*r))) }
T.to.P   = function(t,df) {  return (  2*(1-pt(abs(t),df)))  }
P.to.NLP = function(p) { 
    nlp=-log10(p)
    nlp[is.infinite(nlp)]=-log10(2.2e-16)
    return(nlp) }
##################################################################################################


######## wrapper for R/QTL read cross function ###################################################
readCross = function(g,p) {
   return( 
       read.cross(format='csvsr', 
          genfile = g,
          phefile = p,
          na.strings='NA', genotypes=c('B','R'), alleles=c('B','R'), 
          sep='\t', estimate.map=FALSE, comment.char='')   ) }
##################################################################################################


#########Add genetic map to cross and remove duplicated segregants ...careful,very specific code##
fixrawCross = function(cross,geneticmap.bin) {
    load(geneticmap.bin)
    
    cross = replace.map(cross, map)
    cross = subset(cross, ind=c(1:1011))     
    # remove individuals that have identical genotypes
    cross = subset(cross, ind=-c(125,170,404))
    return(cross)   }
##################################################################################################


###### convert rownames of cross (marker information) to  a data frame ###########################
parse.rownames = function(rn) {
    rrn = data.frame(do.call('rbind', strsplit(rn, '_')),stringsAsFactors=F)
    rrn[,1] = as.numeric(rrn[,1])
    rrn$nchrom = rrn[,2]
    rrn$nchrom = gsub('chr0', 'chr', rrn$nchrom)
    rrn$nchrom = gsub('chr', '', rrn$nchrom)
    rrn$nchrom = as.numeric(rrn$nchrom)
    rrn[,3] = as.numeric(rrn[,3])
    names(rrn) = c('GPOS', 'CHR', 'POS', 'REF', 'ALT', 'NCHR')
    return(rrn)  }
###################################################################################################


##### make CI data frame ##########################################################################
makeCI.df = function(CIs) {
        CIs.df = do.call('rbind',CIs)
        
        CIs.peak  = CIs.df[seq(2,nrow(CIs.df),3),]
        CIs.peak$Marker = rownames(CIs.peak)
        CIs.peak$cpos  = parse.rownames(CIs.peak$Marker)$POS
        CIs.peak$gpos  = parse.rownames(CIs.peak$Marker)$GPOS
        colnames(CIs.peak)= paste('peak',colnames(CIs.peak), sep='')
        
        CIs.left  = CIs.df[seq(1,nrow(CIs.df),3),c(2,3)]
        CIs.left$Marker = rownames(CIs.left)
        CIs.left$cpos  = parse.rownames(CIs.left$Marker)$POS
        CIs.left$gpos  = parse.rownames(CIs.left$Marker)$GPOS
        colnames(CIs.left)= paste('left',colnames(CIs.left), sep='')
        
        CIs.right = CIs.df[seq(3,nrow(CIs.df),3),c(2,3)]
        CIs.right$Marker = rownames(CIs.right)
        CIs.right$cpos  = parse.rownames(CIs.right$Marker)$POS
        CIs.right$gpos  = parse.rownames(CIs.right$Marker)$GPOS
        colnames(CIs.right)= paste('right', colnames(CIs.right), sep='')
        
        CIs = data.frame(CIs.peak, CIs.left, CIs.right)
        return(CIs)
}

###################################################################################################
imputeMissingGeno  = function(cross) {
    impcross = argmax.geno(cross, error.prob=.00001)
    impcross$geno = lapply(impcross$geno, function(x){x$data = x$argmax; x$argmax=NULL; return(x) })
    dm = findDupMarkers(impcross, exact.only = FALSE , adjacent.only = FALSE )
    names(dm)=NULL
    impcross = drop.markers(impcross, unlist(dm))
    return(impcross) }
#####################################################################################################



#######################################################################################################
getPseudoMarkers = function(fcross, cm.spacing = 1) { 
    impcross.1cm   = argmax.geno(fcross, step=cm.spacing, error.prob=.00001)
    pseudomarkers = sapply(impcross.1cm$geno, function(x) { 
                            g = x$argmax;    
                            return(g[,grep('loc', colnames(g))]) })
    pseudomarker.chrlen = sapply(pseudomarkers, ncol)
    pseudomarker.chrs = rep(names(pseudomarker.chrlen), as.vector(pseudomarker.chrlen))

    pm = t(do.call('cbind', pseudomarkers))
    pm.loc = rownames(pm)
    rownames(pm)=paste(pseudomarker.chrs, pm.loc, sep='_')
    colnames(pm) = rownames(fcross$pheno)

    pm.dup = duplicated(pm)
    pm=pm[!pm.dup,]

    bbrr =apply(pm,2, function(x){y=x; y[x==1]='B'; y[x==2]='R'; return(y)})

    pm.matrix = cbind(rownames(pm), pseudomarker.chrs[!pm.dup], bbrr)
    colnames(pm.matrix)[c(1,2)]=c('id', 'chrom')
    return(pm.matrix)  }
#######################################################################################################


#######################################################################################################
make.pcross = function(pm.matrix, pm.geno.file, pheno.full.file, pm.pheno.file) {
    write.table(pm.matrix, file=pm.geno.file, row.names=FALSE, quote=FALSE, sep='\t')
    phenodf = read.delim(pheno.full.file, header=T ,sep='\t', stringsAsFactors=F)
    phenodf = phenodf[,names(phenodf)%in%colnames(pm.matrix)]
    write.table(file=pm.pheno.file,
                phenodf[c(1,2),], sep='\t', quote=F, row.names=F)
     pcross =readCross(pm.geno.file, pm.pheno.file)
     return(pcross)  }

#######################################################################################################


#######################################################################################################
reduceCross = function(cross, impcross.bin) {
    # first pass ... remove duplicate markers
    dm = findDupMarkers(cross) 
    names(dm)=NULL
    cross = drop.markers(cross, unlist(dm))

    # make genotype calls allowing for some error
    impcross = argmax.geno(cross, error.prob=.0005)

    # replace genotype calls with calls from viterbi algorithm
    impcross$geno = lapply(impcross$geno, function(x){x$data = x$argmax; x$argmax=NULL; return(x) })

    # reiterate removal of duplicate markers
    dm = findDupMarkers(impcross) 
    names(dm)=NULL
    impcross = drop.markers(impcross, unlist(dm))
    save(impcross, file=impcross.bin)
    return(impcross)  }
    
#######################################################################################################


#######################################################################################################
readResults = function(rr, batch=1) {
    # Assuming input is an R binary file with variable called 'Results' from 
    # plate growth script
    load(rr)
    plate.names = as.vector(do.call('rbind', (lapply(Results, function(x) { names(x)} ))))
    R=do.call('rbind', Results)
    names(R)=plate.names
    names(R)= gsub(' ', '', names(R))
    names(R)=gsub('384_', '384-', names(R))
    R = R[!is.na(names(R))]
    names(R)=paste(batch, names(R), sep='_')
    return(R)   }
#######################################################################################################


#######################################################################################################
fix5 = function (Results5) {
    names(Results5)= gsub('_0$', '-0', names(Results5))
    names(Results5)= gsub('_1$', '-1', names(Results5))
    names(Results5)= gsub('_2$', '-2', names(Results5))
    names(Results5)= gsub('_3$', '-3', names(Results5))
    return(Results5)  }
#######################################################################################################

#######################################################################################################
getExpAnnot = function(Results) {
    annot = names(Results)
    annot = data.frame(do.call('rbind', strsplit(annot, '_')))
    annot[,5]=toupper(annot[,5])
    names(annot)=c('SET', 'CONFIG', 'FILE', 'PLATEPOSITION', 'DRUG', 'CONDITION')
    return(annot)  }
#######################################################################################################


#######################################################################################################
plotPhenoDists = function(segdd, allh2s, fileo='~/Desktop/PhenoDists.pdf') {
    pdf(fileo, width=8, height=10)
    for ( i in names(segdd)) { 
            boxplot(segdd[[i]], outline=TRUE ,ylim=quantile(segdd[[i]][[4]],c(.005,.995),na.rm=T), 
                    main=i, sub=paste('H2 = ', round(allh2s[[i]],3)) ) 
            stripchart(segdd[[i]], vertical =T,  jitter=.1, method='jitter', add=T, pch=20, col='blue', cex=.7) 
   }
    dev.off() }
#######################################################################################################

#######################################################################################################
make384 = function(idvec) { 
    nseg384 = matrix(NA, 16,24)
    nseg384[seq(1,16,2), seq(1,24,2)]=idvec[1]
    nseg384[seq(1,16,2), seq(2,24,2)]=idvec[2]
    nseg384[seq(2,16,2), seq(1,24,2)]=idvec[3]
    nseg384[seq(2,16,2), seq(2,24,2)]=idvec[4]
    return(as.vector(t(nseg384))) }
#######################################################################################################



chrPeakFinder = function(x,y,z, LOD.threshold, peak.radius ) {
    peak.lods=c()
    peak.ind=c()
    maxLOD=max(y)
    while(maxLOD>LOD.threshold) {
        maxLOD.cind =which.max(y)
        maxLOD.cpos =x[maxLOD.cind]
        maxLOD.ind  =z[maxLOD.cind]
        l.ind = findInterval(maxLOD.cpos-peak.radius, x, all.inside=T, rightmost.closed=T)
        r.ind = findInterval(maxLOD.cpos+peak.radius, x)
        y[l.ind:r.ind]=0
        peak.lods= c(peak.lods, maxLOD)
        peak.ind = c(peak.ind, maxLOD.ind)
        maxLOD = max(y)
    }
   return( cbind(peak.ind, peak.lods) )     
}

peakFinder = function(LODS, LOD.threshold=NA, peak.radius=NA) { 
    mindex.split      = split(seq(1,nrow(LODS)), LODS$chr)
    genetic.pos.split = split(LODS$pos,LODS$chr)

    allpeaks=apply(LODS[,3:ncol(LODS)], 2, function(x){
          LODS.split = split(x, LODS$chr)
          peaks = mapply(chrPeakFinder, x=genetic.pos.split, y=LODS.split,z=mindex.split, 
                                LOD.threshold, peak.radius, SIMPLIFY=FALSE)
          peaks = do.call('rbind', peaks)
          if(!is.null(peaks)) {  
            peaks = data.matrix(peaks)
            peaks = peaks[order(peaks[,'peak.ind']),];     
          }
      })
    return(allpeaks)
}


getChrPeaks = function(mindex.split, chr.mindex.offset, LODS) {
    chr.peaks.lod    = sapply(mindex.split, function(markers) { apply(LODS[,markers], 1, max) })
    # get marker index of LOD peaks per chromosomes                             
    chr.peaks.index = sapply(mindex.split, function(markers)  { apply(LODS[,markers], 1, which.max) })
    # convert chromosome marker index to genome maker index                             
    chr.peaks.index = t(apply(chr.peaks.index, 1, function(x){x+chr.mindex.offset}))

    return(list(chr.peaks.lod = chr.peaks.lod, chr.peaks.index=chr.peaks.index))
}

###### get Peak based FDR #############################################################################
getPeakFDR = function(chromosome.peaks.lod, pdata, gdata, perms=100) {
 
   n.pheno = countStrainsPerTrait(pdata) 
   permpeakLODs = foreach( i = 1:perms ) %dopar% {
        print(i)
        pLODS = get.LOD.by.COR(n.pheno, pdata[sample(1:nrow(pdata)),], gdata)
        sapply(mindex.split, function(markers) { apply(pLODS[,markers], 1, max) })

   }

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

getPeakGWERk = function(chromosome.peaks.lod, pdata, gdata, k=0, per.trait=TRUE, perms=1000) {
 
   n.pheno = countStrainsPerTrait(pdata) 
   permpeakLODs = foreach( i = 1:perms ) %dopar% {
        print(i)
        pLODS = get.LOD.by.COR(n.pheno, pdata[sample(1:nrow(pdata)),], gdata)
        sapply(mindex.split, function(markers) { apply(pLODS[,markers], 1, max) })

   }

   permpeakLODs= abind(permpeakLODs, along=c(3))

   GWER.stat=   apply(permpeakLODs, 3, function(x) {
               apply(x,1, function(y) {
                    sort(y,decreasing=TRUE)[k+1]}) })
   if(per.trait) {  return( apply(GWER.stat, 1, quantile, .95)   ) } else {return( quantile(GWER.stat, .95) )}
}
########################################################################################################

getPeakArray = function(peaklist, threshold) {
    tryCatch( {
    keepPeaks   = which(peaklist$chr.peaks.lod>threshold, arr.ind=T)
    kP = data.frame(rownames(keepPeaks), peaklist$chr.peaks.index[keepPeaks])
    names(kP)=c('trait', 'markerIndex') 
    kP = kP[order(kP$trait, kP$markerIndex),]
    return(kP)} ,error=function(e) {return(NULL) })
}


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
            presids[as.numeric(names(rr)),i]=rr
        }
    }
    return(presids)
} 
########################################################################################################
getscan2summary = function(fdir, s2files) {
    s2_obs = list()
    for(i in s2files) {
        print(i)
        s2fin=paste(fdir, i, sep='')
        load(s2fin)
        s2_obs[[i]]=summary(s2_results,what='int')
     }
     return(s2_obs)
 }

getscan2summaryLK = function(fdir, s2files,rmarker_peaks) {
    s2_obs = list()
    for(i in s2files) {
        print(i)
        s2fin=paste(fdir, i, sep='')
        load(s2fin)
        
        subs2=rmarker_peaks[[i]]
        print(subs2)
        lod2=matrix(NA,nrow(lod), ncol(lod))
        lod2[subs2,]=s2_results$lod[subs2,]
        lod2[,subs2]=s2_results$lod[,subs2]
        
        s2_results$lod=lod2
        s2_obs[[i]]=summary(s2_results,what='int')
    }
    return(s2_obs)
}

calc.s2fdr = function(s2_obs_int, s2_null_int, threshold=seq(2,8.01)){
    s2_obs_cnt = sapply(seq(2, 8, .01), function(thresh) { sum(s2_obs_int>thresh) }   )
    names(s2_obs_cnt) = seq(2, 8, .01)

    s2_exp_cnt = sapply(seq(2, 8, .01),  
             function(thresh) { 
                    print(thresh); 
                    mean(apply(s2_null_int, 3, function(ll) {sum(ll>thresh) }) )
                } )
    names(s2_exp_cnt) = seq(2, 8, .01)
    return(  s2_exp_cnt/s2_obs_cnt )
          
}


##### Forward search for QTLs .... use each highest peak as regressor for next iteration ###############
#### Running this gives no improvement over peak based FDR, so don't bother 
## NULL function 
forwardSearchQTL = function(impcross,  qtlN=30, forwardScanOutDIR)  {

    # some duplicated code here but easier this way
    # extract genotype data from cross structure 
    # code BY as -1 and RM as 1
    gdata = (do.call('cbind', sapply(impcross$geno, function(x) { x$data }))*2)-3

    # extract phenotype data from cross structure and mean center it 
    pdata = apply(impcross$pheno, 2, scale, scale=FALSE)

    pheno = pdata
    
    # count the number of strains that don't have missing data 
    n.pheno = apply(pheno, 2, function(x){sum(!is.na(x))})
    
    mindex.split = getMarkerIndexSplit(impcross)
    chr.mindex.offset = sapply(mindex.split, min)-1

    # This isn't worth the trouble 
    # Forward search QTL model ... up to qtlN putative QTLs ... pick the best marker at each step
    # ... regress it out and repeat-----------------
    for(qtls in 1:qtlN) {
        # display progress 
        print(qtls)
        
        lod.real = get.LOD.by.COR(n.pheno, pheno,gdata)

        # get LOD peaks per chromosome, return LOD scores
        chromosome.peaks.lod    = sapply(mindex.split, function(markers) { 
                                         apply(lod.real[,markers], 1, max) })
        # get marker index of LOD peaks per chromosomes                             
        chromosome.peaks.index = sapply(mindex.split, function(markers)  { 
                                        apply(lod.real[,markers], 1, which.max) })
        # convert chromosome marker index to genome maker index                             
        chromosome.peaks.index = t(apply(chromosome.peaks.index, 1, function(x){x+chr.mindex.offset}))

        # get peak LOD per trait across genome
        peak.lod  = apply(lod.real, 1, max)

        # get peak LOD marker index
        peak.lod.index = apply(lod.real, 1, which.max)
      
        # save phenotype (residuals) and peaks at each iteration
        fileout= paste(forwardScanOutDIR, 'mQTL_', sprintf('%02d',qtls),sep='')
        save(file=fileout, 
             list=c('chromosome.peaks.lod', 'chromosome.peaks.index','peak.lod', 'peak.lod.index', 'pheno','lod.real'),
             compression_level=2)
        print(paste('Done Saving', qtls))

        # update pheno by regressing out effect of largetst QTL from this scan iteration (no intercept/no problem?)
        npheno = (matrix(NA,nrow=nrow(pheno), ncol=ncol(pheno)))
        colnames(npheno)=colnames(pheno)
        for(p in colnames(pheno)){  npheno[,p]= resid(lm(pheno[,p]~-1+gdata[,peak.lod.index[p]],na.action=na.exclude))  }

        pheno=npheno
    }
    return(NULL)
}
#################################################################################################################




############# Extract and re-save just the residuals from each step of the RET forward search #####################
extractForwardSearchResids = function(fS.data.files , outdir) {
    for(i in 1:length(fS.data.files)) {
       print(i)
       load(fS.data.files[i])
       fileout=paste(outdir, sprintf('%02d', i), sep='')
       save(file=fileout,
            list=c('peak.lod.index', 'pheno'),
            compression_level=1) 
    }
    return(NULL) }
####################################################################################################################



########### Extract chromsome and trait peaks from EXPERIMENTAL/REAL data for each step of forward search ###########
# to do ... add option to extract single chromosome peaks for chromosome based FDR (instead of genome based FDR)
extractForwardPeaks = function(fS.data.files) {
    fS.data = list()
    fS.data$chr.peaks.lod = abind(lapply(fS.data.files, function(x)   { 
                        load(x); return(chromosome.peaks.lod)   }), along=3)
    fS.data$chr.peaks.index = abind(lapply(fS.data.files, function(x) {
                       load(x); return(chromosome.peaks.index) }), along=3)
    fS.data$peaks.lod =  sapply(fS.data.files, function(x) {  load(x); return(peak.lod) })
    fS.data$peaks.index= sapply(fS.data.files, function(x) {  load(x); return(peak.lod.index) })
    return(fS.data)  }
######################################################################################################################



########### Extract chromsome and trait peaks from PERMUTED data for each step of forward search #######################
# to do ... add option to extract single chromosome peaks for chromosome based FDR (instead of genome based FDR)
extractForwardPermPeaks = function(fS.perm.data.files) {
    fS.perm.data = list()
    fS.perm.data$chr.peaks.lod = abind( 
            lapply(fS.perm.data.files, 
                    function(x)  {  load(x); return(results)   }), along=4)
    fS.perm.data$peaks.lod = apply(fS.perm.data$chr.peaks.lod, c(1,3,4), max)
    return(fS.perm.data)
}
######################################################################################################################


findit = function(x,vec){
   y = vec - x
   y[y<=0] = NA
   if(all(is.na(y)))NA else which.min(y)
}


#### Convert Peak List to Peak Array ... fix type problems and expand physical marker information #####################
makePeakArray = function(lodPeaks) {
    lodPeakArray         = peak.2.array(lodPeaks)
    lodPeakArray         = na.omit(lodPeakArray)
    lodPeakArray$chr     = as.numeric(lodPeakArray$chr)
    lodPeakArray$lod     = as.numeric(lodPeakArray$lod)
    lodPeakArray$peak.cM = as.numeric(lodPeakArray$peak.cM)
    lodPeakArray$inf.cM  = as.numeric(lodPeakArray$inf.cM)
    lodPeakArray$sup.cM  = as.numeric(lodPeakArray$sup.cM)

    peak.mi =  parse.rownames(lodPeakArray$mname.peak)[,c(1,3)]
    names(peak.mi)=paste('peak', names(peak.mi), sep='')
    inf.mi  =  parse.rownames(lodPeakArray$mname.inf)[,c(1,3)]
    names(inf.mi)=paste('inf', names(inf.mi), sep='')
    sup.mi  =  parse.rownames(lodPeakArray$mname.sup)[,c(1,3)]
    names(sup.mi)=paste('sup', names(sup.mi), sep='')
    lodPeakArray = data.frame(lodPeakArray, peak.mi, inf.mi, sup.mi, stringsAsFactors=FALSE)
}
########################################################################################################################
doQTLModel = function(i,peak.index, fcross, LODSso, refine=TRUE){
    traitName = names(peak.index)[i]
    print(traitName)
    pheno.colIndex = match(traitName,names(fcross$pheno))

    qtlMarkers = unique(peak.index[[i]][,'markerIndex'])
   
    mqtl=NA
    fqtl=NA
    aqtl=NA
    rqtl=NA
    qtlMarkersCP = NA
    CIs =NA
    nQTL = length(qtlMarkers)
    if(nQTL>0) {
        qtlMarkersNames = rownames(LODSso)[qtlMarkers]
        qtlMarkersCP = LODSso[qtlMarkersNames ,c(1,2)]

        mqtl =  makeqtl(fcross, chr=qtlMarkersCP$chr, pos=qtlMarkersCP$pos, qtl.name=qtlMarkersNames, what='prob')
        if(refine == TRUE) {
            rqtl  = refineqtl(fcross,pheno.col=pheno.colIndex, mqtl, method='hk',
                           model='normal', incl.markers=TRUE,maxit=1)
        } else {
             rqtl = mqtl
        }
        fqtl  = fitqtl(fcross, pheno.col =pheno.colIndex, rqtl, method='hk', get.ests=TRUE, dropone=TRUE)
            
        #scan for interactors amongst QTLs with significant main effects
         if(nQTL==1) { aqtl=NA } 
         if(nQTL>1)  { aqtl = addint(fcross, pheno.col=pheno.colIndex, rqtl, method='hk', qtl.only=TRUE)      }

        # get confidence intervals
        if(refine ==TRUE) {
            CIs = lapply(1:length(qtlMarkers),function(x) {lodint(rqtl,qtl.index=x) })
            CIs = makeCI.df(CIs)
            }
    }
   return(list(qtlMarkers=qtlMarkers, qtlMarkersCP=qtlMarkersCP, mqtl=mqtl, fqtl=fqtl, aqtl=aqtl,rqtl=rqtl,CIs=CIs))
}    

##### Make QTL Model ###################################################################################################
qtlModel  = function(traitPeakArray) {
    traitName =  unique(traitPeakArray$trait)
    pheno.colIndex = match(traitName,names(fcross$pheno))
    
    mqtls = makeqtl(fcross, chr=traitPeakArray$chr, pos=traitPeakArray$peak.cM, what='prob')
    fqtl  = fitqtl(fcross, pheno.col =pheno.colIndex, mqtls, method='hk', get.ests=TRUE, dropone=TRUE)
    
    pheno = pdata[,pheno.colIndex]
    geno  = gdata[,as.character(traitPeakArray$mname.peak)]

    pearsonR = cor(pheno, geno, use='pairwise.complete.obs')
    pearsonR2 =pearsonR^2

    LOD=(-sum(!is.na(pheno))*log(1-pearsonR2))/(2*log(10)) 
    attr(LOD, 'trait')=traitName
    
    if(nrow(traitPeakArray)==1) {aqtl=NA} 
    else { aqtl = addint(fcross, pheno.col=pheno.colIndex, mqtls, method='hk', qtl.only=TRUE)}
    return(list(mqtls=mqtls, fqtl=fqtl,pearsonR=pearsonR, pearsonR2=pearsonR2, LOD=LOD, aqtl=aqtl))
}
########################################################################################################





#########################################################################################################
# x is list of data frames for each experiment set
ols3 <- function (y, x) {
    XtX <- crossprod(x)
    Xty <- crossprod(x, y)
    solve(XtX, Xty)
}
#######################################################################################################



###### Cross - validation ##############################################################################
doCV = function(fcross, LODS.01s, mindex.split, chr.mindex.offset, cv.fold=10) {
    #gsize.notlast = round(1008/cv.fold)
    #gsize.last    = 1008-( gsize.notlast *(cv.fold-1))
    if(cv.fold==10) {
          grps= c(rep('A',101),rep('B',101),rep('C',101),rep('D',101),rep('E',101),
                  rep('F',101),rep('G',101),rep('H',101),rep('I',100),rep('J',100) )
    }
    if(cv.fold==5) {
       grps= c(rep('A',202),rep('B',202),rep('C',201),rep('D',201),rep('E',202))
    }
    sgrps = sample(grps)
    
    


    cv=foreach( gg = unique(grps) ) %dopar% {
        print(gg)
        nAcross = subset(fcross, ind=(sgrps!=gg))
        Across  = subset(fcross, ind=(sgrps==gg))

        Agdata         =  extractGenotype(nAcross)
        Apdata.01      =  extractScaledPhenotype(nAcross, TRUE)
        An.pheno       =  countStrainsPerTrait(nAcross$pheno)
        ALODS.01       =  get.LOD.by.COR(An.pheno, Apdata.01, Agdata)
        ALODS.01s      =  LODmatrix.2.scanone(ALODS.01, nAcross)
        Apeaklist.01   =  getChrPeaks(mindex.split, chr.mindex.offset, ALODS.01) 
        ApeakArray.01  =  getPeakArray(Apeaklist.01, 2.69)
        print('step1')

        Apdata.02      = getPhenoResids(Apdata.01, Agdata, ApeakArray.01) 
        ALODS.02       = get.LOD.by.COR(An.pheno, Apdata.02, Agdata)
        ALODS.02s      = LODmatrix.2.scanone(ALODS.02, nAcross,ALODS.01s)
        Apeaklist.02   = getChrPeaks(mindex.split, chr.mindex.offset, ALODS.02) 
        ApeakArray.02  = getPeakArray(Apeaklist.02, 2.88)
        print('step2')

        Apdata.03      = getPhenoResids(Apdata.02, Agdata, ApeakArray.02) 
        ALODS.03       = get.LOD.by.COR(An.pheno, Apdata.03, Agdata)
        ALODS.03s      = LODmatrix.2.scanone(ALODS.03, nAcross,ALODS.01s)
        Apeaklist.03   = getChrPeaks(mindex.split, chr.mindex.offset, ALODS.03) 
        ApeakArray.03  = getPeakArray(Apeaklist.03, 3.75)
        print('step3')

        Apdata.04      = getPhenoResids(Apdata.03, Agdata, ApeakArray.03) 
        ALODS.04       = get.LOD.by.COR(An.pheno, Apdata.04, Agdata)
        ALODS.04s      = LODmatrix.2.scanone(ALODS.04, nAcross,ALODS.01s)
        Apeaklist.04   = getChrPeaks(mindex.split, chr.mindex.offset, ALODS.04) 
        ApeakArray.04  = getPeakArray(Apeaklist.04, 4.63)
        print('step4')

        ApA1=cbind(ApeakArray.01, 'J1')
        ApA2=cbind(ApeakArray.02, 'J2')
        ApA3=cbind(ApeakArray.03, 'J3')
        if(!is.null(ApeakArray.04)){        ApA4=cbind(ApeakArray.04, 'J4') }

        names(ApA1)[3]='jump'
        names(ApA2)[3]='jump'
        names(ApA3)[3]='jump'
        if(!is.null(ApeakArray.04)){ names(ApA4)[3]='jump' }
        if(!is.null(ApeakArray.04)) { p.df = rbind(ApA1, ApA2, ApA3, ApA4) } else {
             p.df = rbind(ApA1, ApA2, ApA3) }

        Apeak.index = data.frame(p.df)
        Apeak.index = Apeak.index[order(Apeak.index$trait, Apeak.index$markerIndex),]
        Apeak.index = split(Apeak.index, Apeak.index$trait)

        print('Built Peak Index List')
        # Fit QTL Model 
        aQTLS = foreach(i=1:length(Apeak.index)) %dopar% { doQTLModel(i,Apeak.index,nAcross, ALODS.01s,refine=FALSE) }
        names(aQTLS) = names(Apeak.index)
        afQTLs= lapply(1:length(Apeak.index), function(i){ fitqtl(Across, pheno.col=i,aQTLS[[i]]$rqtl, method='hk', get.ests=FALSE, dropone=FALSE) })
        names(afQTLs)=names(aQTLS)

        AVEQA=sapply(afQTLs, function(x) { 
            y=NA
            tryCatch( 
                {
                 y= as.numeric(x$result.full[,'%var'][1]) },
                   error= function(e){y=NA})
           return(y/100)} )
        AVEQA=t(data.frame(t(AVEQA)))
       return(list(peak.index=Apeak.index, QTLs=aQTLS,fQTLs=afQTLs,ve=AVEQA))
    }
    return(cv)
}



l384_96= list(
     '384-01.txt'=make384(c(1,2,3,4)),
     '384-02.txt'=make384(c(5,6,7,8)),
     '384-03.txt'=make384(c(9,10,11,12)),
     '384-04.txt'=make384(c(12,7,6,1)),
     '384-05.txt'=make384(c(8,11,2,5)),
     '384-06.txt'=make384(c(4,3,10,9)))



plotOverlappingHist <- function(a, b, colors=c("white","gray20","gray50"),
                                breaks=NULL, xlim=NULL, ylim=NULL){
 
  ahist=NULL
  bhist=NULL
 
  if(!(is.null(breaks))){
    ahist=hist(a,breaks=breaks,plot=F)
    bhist=hist(b,breaks=breaks,plot=F)
  } else {
    ahist=hist(a,plot=F)
    bhist=hist(b,plot=F)
 
    dist = ahist$breaks[2]-ahist$breaks[1]
    breaks = seq(min(ahist$breaks,bhist$breaks),max(ahist$breaks,bhist$breaks),dist)
 
    ahist=hist(a,breaks=breaks,plot=F)
    bhist=hist(b,breaks=breaks,plot=F)
  }
 
  if(is.null(xlim)){
    xlim = c(min(ahist$breaks,bhist$breaks),max(ahist$breaks,bhist$breaks))
  }
 
  if(is.null(ylim)){
    ylim = c(0,max(ahist$counts,bhist$counts))
  }
 
  overlap = ahist
  for(i in 1:length(overlap$counts)){
    if(ahist$counts[i] > 0 & bhist$counts[i] > 0){
      overlap$counts[i] = min(ahist$counts[i],bhist$counts[i])
    } else {
      overlap$counts[i] = 0
    }
  }
 
  plot(ahist, xlim=xlim, ylim=ylim, col=colors[1])
  plot(bhist, xlim=xlim, ylim=ylim, col=colors[2], add=T)
  plot(overlap, xlim=xlim, ylim=ylim, col=colors[3], add=T)
}

# plot variance components 
scatterHeritability= function(n, b, ne, be, tn,  B2=NULL, p.NA=NULL , 
        xlab='narrow-sense heritability',
        ylab='broad-sense heritability',
        do.diff = FALSE,
        do.name = TRUE,
        do.save = FALSE, pdf.file='~/Desktop/paper_figures/H2normnt.pdf') {

    if(do.save)   pdf(file=pdf.file, width=8, height=8)

    par(xaxs="i", yaxs="i") 
    plotCI(n,  b,  uiw=ne, liw=ne,
            err='x', xlim=c(0,1), ylim=c(0,1), gap=0, 
            xlab= xlab, ylab=ylab,
            sfrac=0, pch=20, lty=2,cex=1.2, barcol='darkgrey')
  
    abline(0,1, col='darkgrey')
    abline(lm(b~n),col='lightblue')

    if(do.name)   {  text(n,b, tn, cex=.7, pos=3, col='blue')}

    if(do.diff) {
        segments(n, b, n , B2)
        segments(n, b, n ,B2, col = (qvalue(p.NA)$qvalues<.1)+1 )
    }
    
    plotCI(n, b,  uiw=be,  liw=be,
            err='y',  xlim=c(0,1), ylim=c(0,1), 
            gap=0, sfrac=.0, pch=20, lty=2,cex=1.2, barcol='darkgrey',add=T)

    if(do.save) dev.off()
}



# modified from e/qtl package
define.peak2=function (scanone, lodcolumn = 1, chr, th = 2.3, si = 1.5, graph = FALSE, 
    window.size = 20, round, save.pict = FALSE, phe.name, ...) 
{
    require(qtl)
    if (!all(attr(scanone, "class", exact = TRUE) %in% c("scanone", 
        "data.frame"))) 
        stop("Input should have class \"scanone\"")
    if (!missing(chr)) {
        if (!is.vector(chr)) 
            stop("Argument chr misspecified: expecting a vector")
        if (!all(chr %in% levels(scanone$chr))) 
            stop("Could not find all chr \"", chr, "\" in scanone$chr")
        if (length(chr) > length(levels(scanone$chr))) 
            stop("Argument chr misspecified: chr could not be longer than the number of chr. in scanone", 
                length(levels(scanone$chr)))
    }
    else chr <- levels(scanone$chr)
    if (!is.numeric(th) | !is.vector(th) | length(th) != 1) 
        stop("Argument th misspecified: expecting a numeric vector of length 1")
    if (!is.numeric(si) | !is.vector(si) | length(si) != 1) 
        stop("Argument si misspecified: expecting a numeric vector of length 1")
    if (!is.logical(graph) & !is.vector(graph) & length(graph) != 
        1) 
        stop("Argument graph misspecified: expecting a logical vector of length 1")
    if (!is.vector(window.size) | !is.numeric(window.size) | 
        length(window.size) > 1) 
        stop("Argument window.size misspecified: expecting a numeric vector of length 1")
    if (window.size < 0) 
        stop("Argument window.size misspecified: wrong value (centiMorgan)")
    if (missing(lodcolumn)) 
        stop("Argument lodcolumn unspecified")
    if (!is.vector(lodcolumn) & (!is.numeric(lodcolumn) | !is.character(lodcolumn))) 
        stop("Argument lodcolumn misspecified: expecting numeric vector or character vector")
    if (length(lodcolumn) > length(names(scanone)) - 2) 
        stop("Argument lodcolumn misspecified: cannot refered to more lodcolumn than scanone contains. The vector is too big.")
    if (length(names(scanone)) == 3) 
        lodcolumn <- "all"
    if (length(lodcolumn) == 1 & lodcolumn[1] == "all") {
        trait <- names(scanone[3:length(names(scanone))])
        if (ncol(scanone) == 3 & names(scanone)[3] == "lod") {
            if (missing(phe.name)) 
                stop("Warning: the lodcolumn name is 'lod', argument phe.name is necessary")
            trait <- phe.name
            if (toupper(phe.name) == "ID") 
                stop("Warning: ID is reserved name for indexing the phenotypes")
        }
        num <- c(1:(length(names(scanone)) - 2))
    }
    else {
        trait <- "NA"
        num <- NA
        for (i in seq(length(lodcolumn))) {
            if (is.numeric(lodcolumn[i])) {
                if (lodcolumn[i] < 1 | lodcolumn[i] >= length(names(scanone)) - 
                  2) 
                  stop("lodcolumn values should be between 1 and the no. of lod columns")
                else trait <- c(trait, names(scanone[lodcolumn[i] + 
                  2]))
            }
            if (is.character(lodcolumn[i])) {
                trait <- c(trait, lodcolumn[i])
                if (!toupper(lodcolumn[i]) %in% toupper(names(scanone))) 
                  stop("Could not identify trait lod column \"", 
                    lodcolumn[i], "\"")
                lodcolumn[i] <- as.numeric(grep(toupper(lodcolumn[i]), 
                  toupper(names(scanone)))) - 2
            }
            num <- c(num, as.numeric(lodcolumn[i]))
        }
        lodcolumn <- as.numeric(lodcolumn)
        trait <- trait[-1]
        num <- num[-1]
    }
    res <- list(init = NA)
    resbytrait <- list(init = NA)
    cat("no. of traits: ", length(trait), "\n", sep = "")
    for (i in seq(trait)) cat("trait: ", trait[i], "\tlodcolumn: ", 
        num[i], "\n", sep = "")
    cat("define.peak in process...\n")
    if (!missing(round) && !is.numeric(round) && round < 0) 
        stop("Argument round misspecified: round should be an integer >= 0")
    define.curve <- function(scanone, num, chr) {
        curve <- list(un = NA)
        for (i in chr) {
            bool <- scanone$chr %in% i
            res <- list(list(mname = row.names(scanone)[bool], 
                pos = scanone$pos[bool], lod = scanone[bool, 
                  num]))
            attributes(res)$names <- paste(i)
            curve <- c(curve, res)
        }
        curve <- curve[-1]
        invisible(curve)
    }
    curve <- define.curve(scanone, num + 2, chr)
    for (z in 1:length(trait)) {
        for (i in seq(length(chr))) {
            c <- grep(paste("^", chr[i], "$", sep = ""), names(curve))
            if (length(trait) > 1) 
                lod <- curve[[c]]$lod[[z]]
            else lod <- curve[[c]]$lod
            pos <- curve[[c]]$pos
            if (max(lod) < th) {
                resbychr <- list(NA)
                attributes(resbychr)$names <- chr[i]
                resbytrait <- c(resbytrait, resbychr)
                next
            }
            filter.peak <- function(pos = pos, lod = lod, dmin = 20, 
                th) {
                if (!any(lod >= th)) 
                  return()
                mx <- NA
                for (i in seq(length(pos))) {
                  if (i == 1 && lod[i] > lod[i + 1] && lod[i] >= 
                    th) 
                    mx <- c(mx, i)
                  if (i == length(lod) && lod[i] > lod[i - 1] && 
                    lod[i] >= th) 
                    mx <- c(mx, i)
                  if (i != 1 && i != length(lod) && lod[i] > 
                    lod[i - 1] && lod[i] > lod[i + 1] && lod[i] >= 
                    th) 
                    mx <- c(mx, i)
                }
                if (length(mx) == 1) {
                  return()
                }
                else mx <- mx[-1]
                if (length(mx) == 1) {
                  return(mx)
                }
                step <- NA
                for (i in seq(length(mx))) {
                  if (i == 1) 
                    next
                  step <- c(step, pos[mx[i]] - pos[mx[i - 1]])
                }
                step <- step[-1]
                step <- min(step)
                del_peak <- NA
                for (p in seq(1, length(pos), by = step)) {
                  if (all(!(pos[mx] >= p & pos[mx] <= p + dmin))) 
                    next
                  cple <- mx[pos[mx] >= p & pos[mx] <= p + dmin]
                  if (all(cple %in% del_peak)) 
                    next
                  cple <- cple[!cple %in% del_peak]
                  del_peak <- c(del_peak, cple[lod[cple] != max(lod[cple])])
                }
                del_peak <- del_peak[-1]
                mx <- mx[!mx %in% del_peak]
                return(mx)
            }
            max <- filter.peak(pos = pos, lod = lod, dmin = window.size, 
                th = th)
            if (length(max) < 1) 
                next
            if (length(trait) > 1) 
                trait_curve <- data.frame(pos = curve[[c]]$pos, 
                  lod = curve[[c]]$lod[[z]])
            else trait_curve <- data.frame(pos = curve[[c]]$pos, 
                lod = curve[[c]]$lod)
            limit.peak <- function(curve, peak, s, m = 50) {
                inf <- NA
                sup <- NA
                infq <- NA
                supq <- NA
                for (i in seq(length(curve$lod))) {
                  for (y in curve$pos[peak]) {
                    if (curve$pos[i] == y) {
                      b <- i
                      a <- i
                      if (i == 1) {
                        for (b in i:(length(curve$lod) - 1)) {
                          if (curve$lod[b] <= 0) {
                            q <- "*"
                            b <- b - 1
                            break
                          }
                          if (curve$pos[b] - curve$pos[i] > m & 
                            m != 0) {
                            q <- "->"
                            break
                          }
                          if (curve$lod[b] < curve$lod[i] - s) {
                            q <- "+"
                            break
                          }
                          b <- b + 1
                        }
                        sup <- c(sup, b)
                        inf <- c(inf, i)
                        infq <- c(infq, "|")
                        supq <- c(supq, q)
                      }
                      if (i == length(curve$lod)) {
                        for (a in i:2) {
                          if (curve$lod[a] <= 0) {
                            q <- "*"
                            a <- a + 1
                            break
                          }
                          if (curve$pos[i] - curve$pos[a] > m & 
                            m != 0) {
                            q <- "<-"
                            break
                          }
                          if (curve$lod[a] < curve$lod[i] - s) {
                            q <- "+"
                            break
                          }
                          a <- a - 1
                        }
                        inf <- c(inf, a)
                        sup <- c(sup, b)
                        supq <- c(supq, "|")
                        infq <- c(infq, q)
                      }
                      if (i != 1 && i != length(curve$lod)) {
                        for (a in i:2) {
                          if (curve$lod[a] <= 0) {
                            q <- "*"
                            a <- a + 1
                            break
                          }
                          if (curve$pos[i] - curve$pos[a] > m & 
                            m != 0) {
                            q <- "<-"
                            break
                          }
                          if (curve$lod[a] < curve$lod[i] - s) {
                            q <- "+"
                            break
                          }
                          q <- "|"
                          a <- a - 1
                        }
                        inf <- c(inf, a)
                        infq <- c(infq, q)
                        for (b in i:(length(curve$lod) - 1)) {
                          if (curve$lod[b] <= 0) {
                            q <- "*"
                            b <- b - 1
                            break
                          }
                          if (curve$pos[b] - curve$pos[i] > m & 
                            m != 0) {
                            q <- "->"
                            break
                          }
                          if (curve$lod[b] < curve$lod[i] - s) {
                            q <- "+"
                            break
                          }
                          q <- "|"
                          b <- b + 1
                        }
                        sup <- c(sup, b)
                        supq <- c(supq, q)
                      }
                      break
                    }
                  }
                }
                inf <- inf[-1]
                sup <- sup[-1]
                supq <- supq[-1]
                infq <- infq[-1]
                qual <- paste(infq, supq, sep = "")
                invisible(data.frame(inf = inf, max = peak, sup = sup, 
                  qual = qual))
            }
            peak <- limit.peak(trait_curve, max, si, ...)
            if (!missing(round)) 
                trait_curve$lod[peak$max] <- round(as.numeric(as.vector(trait_curve$lod[peak$max])), 
                  round)
            resbychr <- list(data.frame(lod = trait_curve$lod[peak$max], 
                mname.peak = curve[[c]]$mname[peak$max], peak.cM = curve[[c]]$pos[peak$max], 
                mname.inf = curve[[c]]$mname[peak$inf], inf.cM = curve[[c]]$pos[peak$inf], 
                mname.sup = curve[[c]]$mname[peak$sup], sup.cM = curve[[c]]$pos[peak$sup], 
                si.quality = peak$qual))
            attributes(resbychr)$names <- chr[i]
            resbytrait <- c(resbytrait, resbychr)
            if (save.pict) 
                png(filename = paste(trait[z], i, z, ".png", 
                  sep = "_"), width = 1280, height = 1024)
            if (graph | save.pict) {
                qtl:::plot.scanone(scanone, lodcolumn = num[z], 
                  chr = chr[i], show.marker.names = TRUE, lwd = 1, 
                  main = paste("chr", chr[i]))
                abline(h = th, col = "pink", lwd = 1)
                abline(v = curve[[c]]$pos[peak$inf], col = "blue")
                abline(v = curve[[c]]$pos[peak$sup], col = "blue")
                abline(v = curve[[c]]$pos[peak$max], col = "red", 
                  lty = 2)
            }
            if (save.pict) 
                dev.off()
        }
        resbytrait <- list(resbytrait[-1])
        if (length(trait) > 1) 
            attributes(resbytrait)$names <- paste(names(curve[[c]]$lod[z]))
        else attributes(resbytrait)$names <- paste(trait)
        res <- c(res, resbytrait)
    }
    res <- res[-1]
    attributes(res)$class <- c("peak", "list")
    attributes(res)$features <- c("lod", "mname.peak", "peak.cM", 
        "mname.inf", "inf.cM", "mname.sup", "sup.cM", "si.quality")
    attributes(res)$scanone <- deparse(substitute(scanone))
    attributes(res)$lod.th <- th
    attributes(res)$si <- si
    attributes(res)$window <- window.size
    return(res)
}
# resolve errors in eQTL package define.peak algorithm by collapsing QTLs if their confidence intervals
# overlap
cleanPeaks = function(lp1){
    lp2=lapply(lp1, function(x){ lapply(x, function(y){
        if(is.null(dim(y))){ return(y)}
        while(nrow(y)>1) {
            z=Intervals(y[,c('inf.cM', 'sup.cM')])
            zz=interval_overlap(z,z)
            zz.length=sapply(zz, length)
            zz.max= max(zz.length)
            if(zz.max>1) {
                zz.max.ind = which.max(zz.length)
                olap.ind=zz[[zz.max.ind]]
                y[olap.ind[1],]=y[olap.ind[which.max(y[olap.ind,'lod'])],]
                y=y[-olap.ind[-1],]
            } else {  return(y) }
        }
        return(y)
    }) }) 
    attr(lp2, "class")   = attr(lp1, "class")
    attr(lp2, "features")=  attr(lp1, "features")
    attr(lp2, "scanone") =  attr(lp1, "scanone") 
    attr(lp2, "lod.th")  =  attr(lp1, "lod.th")  
    attr(lp2, "si")      =  attr(lp1, "si")      
    attr(lp2, "window")  =  attr(lp1, "window") 
    return(lp2)
}



chrPeakFinder = function(x,y,z, LOD.threshold, peak.radius ) {
    peak.lods=c()
    peak.ind=c()
    maxLOD=max(y)
    while(maxLOD>LOD.threshold) {
        maxLOD.cind =which.max(y)
        maxLOD.cpos =x[maxLOD.cind]
        maxLOD.ind  =z[maxLOD.cind]
        l.ind = findInterval(maxLOD.cpos-peak.radius, x, all.inside=T, rightmost.closed=T)
        r.ind = findInterval(maxLOD.cpos+peak.radius, x)
        y[l.ind:r.ind]=0
        peak.lods= c(peak.lods, maxLOD)
        peak.ind = c(peak.ind, maxLOD.ind)
        maxLOD = max(y)
    }
   return( cbind(peak.ind, peak.lods) )     
}








