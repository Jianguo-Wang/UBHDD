# UBHDD


# 02/13/22
# R script in:
# The trait coding rule in phenotype space
# 
# Jianguo Wang (wangjg36@mail.sysu.edu.cn)
# 


############### /UBHDD/scripts description ################################ 

	/brain_genotype_phenotype_preprocessing.R is used to preprocess the genome and phenome data in UK Biobank.
	/data_summary.R is used to generate the three tables in folder /UBHDD/summary.
	/heritability_and_QTL_mapping_for_brain.R is used to conduct narrow-sense heritability and QTL mapping for brain traits.
	/heritability_and_QTL_mapping_for_yeast.R is used to conduct broad-sense, narrow-sense heritability and QTL mapping for yeast traits.
	/plot.R is used to draw the figures.
	/QTL_mappingFx.R is a script previously provided by the paper "Finding the missing heritability in a yeast cross" with a few out-of-date packages updated.
	/ROIEditor_usage.txt is a short description on the usage of 'ROIEditor' to obtain the Fig.3a.
	/simulation.R is used to generate a simulated phenotype space in which UBHDD model is conducted.
	/UBHDD_for_brain.py is used to conduct linear model in a learning framework for brain, the core package is 'glmnet'.
	/UBHDD_for_yeast.R is used to conduct linear model in a learning framework for yeast, the core package is 'glmnet'.


############### /UBHDD/data description ################################ 

	/UBHDD/data/yeast includes the basic data of yeast: 
		/cross.RData contains the geneotype and mean phenotype information for each segregants of yeast cross population
		/pheno_raw_405traits.RData contains the two replicate values of each segregants for each of 405 morphological traits.
		/mt4718.RData contains the morphological traits in yeast deletome population
	
	/UBHDD/data/brain includes the basic data of brain in UK Biobank
		/brain_left_right_traitName.RData contains  the brain trait name in left/right brain regions
		/brain_region_abbreviation.xlsx contains the abbreviation and complete brain region name
		/brainCovariateName_b2o.RData is the MRI measuring sysmtem related covariates according to the suggestions in UK Biobank.
		/brainPhenoName_b2o.RData is the continuous brain traits
		/brainRelatedTraits.RData is various brain-related traits
		/clusterName.RData is the exemplars of brain traits by dMRI and obtained by R package 'apcluster'.
		/dMRI_trait.RData is the brain traits by dMRI used for main analysis
		/totalCovariateName_b2o.RData is the covariates including age, sex, weight, home location and measuring centers
		/traitAnnotation_all18108.csv is an annotation file for the available UK Biobank traits [by hand; traits in UKB are not strictly mutually exclusive across categories]


############### /UBHDD/summary description ################################ 

	/brainFeature.csv is the main features of brain traits obtained in this study.
	/brainPairFeature.csv is the main features of brain traits in symmetrical brain regions in this study.
	/yeastFeature.csv is the main features of yeast in this study.


############## remarks #########################################################

#the UK Biobank data are not fully shared due to the policy of UK Biobank.

#the file/folder paths in these scripts are organized as follows:

	/UBHDD/figs
		/UBHDD/figs/functions

    /UBHDD/theory
		/UBHDD/theory/functions
		/UBHDD/theory/scripts
		/UBHDD/theory/data		
	
	/UBHDD/yeastSegregants
		/UBHDD/yeastSegregants/functions
		/UBHDD/yeastSegregants/scripts
		/UBHDD/yeastSegregants/data
			/UBHDD/yeastSegregants/data/basic
			/UBHDD/yeastSegregants/data/downSampling
			/UBHDD/yeastSegregants/data/oneRep
			/UBHDD/yeastSegregants/data/QTLmapping
			/UBHDD/yeastSegregants/data/unCorCutoffRobust
		
	/UBHDD/yeastDeletome
		/UBHDD/yeastDeletome/functions
		/UBHDD/yeastDeletome/scripts
		/UBHDD/yeastDeletome/data
	
	/UBHDD/UKB
		/UBHDD/UKB/functions
		/UBHDD/UKB/scripts
		/UBHDD/UKB/data
			/UBHDD/UKB/data/b2b
			/UBHDD/UKB/data/basic
			/UBHDD/UKB/data/downSampling
			/UBHDD/UKB/data/unCorCutoffRobust
			/UBHDD/UKB/data/genomeProcessing
				/UBHDD/UKB/data/genomeProcessing/geno_assoc_fastGWA_mlm_water
				/UBHDD/UKB/data/genomeProcessing/genome_plink_filter
				/UBHDD/UKB/data/genomeProcessing/genome_plink_maf_indep_geno_mind_hwe
				/UBHDD/UKB/data/genomeProcessing/genome_plink_maf_indep_geno_mind_hwe_pruned
				/UBHDD/UKB/data/genomeProcessing/genome_qctool_filter
				/UBHDD/UKB/data/genomeProcessing/GRM
				/UBHDD/UKB/data/genomeProcessing/sp_GRM
				/UBHDD/UKB/data/genomeProcessing/ukb_imp_fam_plink_filter
				/UBHDD/UKB/data/genomeProcessing/ukb_imp_mfi_qctool_filter
				/UBHDD/UKB/data/genomeProcessing/Vg_GREML_water_exemplar			
				/UBHDD/UKB/data/genomeProcessing/clumping_analysis
					/UBHDD/UKB/data/genomeProcessing/clumping_analysis/QTL_mapping_result
					/UBHDD/UKB/data/genomeProcessing/clumping_analysis/Tg
					/UBHDD/UKB/data/genomeProcessing/clumping_analysis/Tng
				/UBHDD/UKB/data/genomeProcessing/sample_MAF_INFO
					/UBHDD/UKB/data/genomeProcessing/sample_MAF_INFO/chr_sample
					/UBHDD/UKB/data/genomeProcessing/sample_MAF_INFO/ukb_imp_mfi				

	

#Based on the organization, one can easily reproduce our results. 

#However,it will be welcom to tell us if you have any question.











