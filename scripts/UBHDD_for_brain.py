
# 02/13/22
# python script in:
# The trait coding rule in phenotype space
# 
# Jianguo Wang (wangjg36@mail.sysu.edu.cn)
# 


cd UBHDD/UKB/

import multiprocessing
import os
import time
import numpy

import numpy as np
import pandas as pd
from pandas import  *
from numpy import *

from glmnet import ElasticNet

import random

##

brain = pd.read_csv('data/b2b/brain_with_genome_exclude_covariate_water_T.csv');
brainColName=brain.columns[1:(brain.shape[1]+1)]
brainRowName=brain.iloc[:,0]

brain2=brain.iloc[:,1:(brain.shape[1]+1)];

brain3=brain2.to_numpy()

del brain,brain2


##

corMat = pd.read_csv('data/b2b/brain_with_genome_exclude_covariate_water_T_corMat.csv');
corMat2=corMat.iloc[:,1:(corMat.shape[1]+1)];

corMat3=corMat2.to_numpy()

del corMat,corMat2


##

shuffle = pd.read_csv('data/b2b/brain_with_genome_exclude_covariate_water_T_shuffleIndex.csv');
shuffle2=shuffle.iloc[:,1:(shuffle.shape[1]+1)];

shuffle3=shuffle2.to_numpy()

del shuffle,shuffle2




##################### UBHDD #############

	
def linear_learn_unCor(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	focalUnCorInd=np.nonzero(abs(corMat3[:,i])<0.15)[0]
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###

pool = multiprocessing.Pool(processes=80)
result = pool.map(linear_learn_unCor, np.arange(brain3.shape[1]))

pool.close()
pool.join()

time.localtime()




result_b2b_water_015_DF=pd.DataFrame(result,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName)
result_b2b_water_015_DF.to_csv('data/b2b/result_b2b_water_015_DF.csv')


######################################################################

############### shuffle ######################


#####################

	
def linear_learn_unCor_shuffle(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	focalUnCorInd=np.nonzero(abs(corMat3[:,i])<0.15)[0]
	fit=m.fit(brain3[:,focalUnCorInd],brain3[shuffle3[:,i]-1,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###

pool = multiprocessing.Pool(processes=100)
result = pool.map(linear_learn_unCor_shuffle, np.arange(brain3.shape[1]))

pool.close()
pool.join()

time.localtime()




result_b2b_water_015_DF_shuffle=pd.DataFrame(result,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName)
result_b2b_water_015_DF_shuffle.to_csv('data/b2b/result_b2b_water_015_DF_shuffle.csv')



######################### different uncorrelation cutoff #############################

### threshold=1

del linear_learn_unCor095, result095, result_b2b_water_095_DF

def linear_learn_unCor1(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	focalUnCorInd=np.nonzero(abs(corMat3[:,i])<1)[0]
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###

pool = multiprocessing.Pool(processes=100)
result1 = pool.map(linear_learn_unCor1, np.arange(brain3.shape[1]))

pool.close()
pool.join()

time.localtime()




result_b2b_water_1_DF=pd.DataFrame(result1,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName)
result_b2b_water_1_DF.to_csv('data/b2b/result_b2b_water_1_DF.csv')


######################################################################################


######################### different uncorrelation cutoff #############################

### threshold=0.1

del linear_learn_unCor005, result005, result_b2b_water_005_DF

def linear_learn_unCor01(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	focalUnCorInd=np.nonzero(abs(corMat3[:,i])<0.1)[0]
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###

pool = multiprocessing.Pool(processes=100)
result01 = pool.map(linear_learn_unCor01, np.arange(brain3.shape[1]))

pool.close()
pool.join()

time.localtime()




result_b2b_water_01_DF=pd.DataFrame(result01,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName)
result_b2b_water_01_DF.to_csv('data/unCorCutoffRobust/result_b2b_water_01_DF.csv')



######################### different uncorrelation cutoff #############################

### threshold=0.11

del linear_learn_unCor01, result01, result_b2b_water_01_DF

def linear_learn_unCor011(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	focalUnCorInd=np.nonzero(abs(corMat3[:,i])<0.11)[0]
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###

pool = multiprocessing.Pool(processes=100)
result011 = pool.map(linear_learn_unCor011, np.arange(brain3.shape[1]))

pool.close()
pool.join()

time.localtime()




result_b2b_water_011_DF=pd.DataFrame(result011,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName)
result_b2b_water_011_DF.to_csv('data/unCorCutoffRobust/result_b2b_water_011_DF.csv')




######################### different uncorrelation cutoff #############################

### threshold=0.12

del linear_learn_unCor011, result011, result_b2b_water_011_DF

def linear_learn_unCor012(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	focalUnCorInd=np.nonzero(abs(corMat3[:,i])<0.12)[0]
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###

pool = multiprocessing.Pool(processes=100)
result012 = pool.map(linear_learn_unCor012, np.arange(brain3.shape[1]))

pool.close()
pool.join()

time.localtime()




result_b2b_water_012_DF=pd.DataFrame(result012,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName)
result_b2b_water_012_DF.to_csv('data/unCorCutoffRobust/result_b2b_water_012_DF.csv')




######################### different uncorrelation cutoff #############################

### threshold=0.13

del linear_learn_unCor012, result012, result_b2b_water_012_DF

def linear_learn_unCor013(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	focalUnCorInd=np.nonzero(abs(corMat3[:,i])<0.13)[0]
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###

pool = multiprocessing.Pool(processes=100)
result013 = pool.map(linear_learn_unCor013, np.arange(brain3.shape[1]))

pool.close()
pool.join()

time.localtime()




result_b2b_water_013_DF=pd.DataFrame(result013,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName)
result_b2b_water_013_DF.to_csv('data/unCorCutoffRobust/result_b2b_water_013_DF.csv')




######################### different uncorrelation cutoff #############################

### threshold=0.14

del linear_learn_unCor013, result013, result_b2b_water_013_DF

def linear_learn_unCor014(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	focalUnCorInd=np.nonzero(abs(corMat3[:,i])<0.14)[0]
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###

pool = multiprocessing.Pool(processes=100)
result014 = pool.map(linear_learn_unCor014, np.arange(brain3.shape[1]))

pool.close()
pool.join()

time.localtime()




result_b2b_water_014_DF=pd.DataFrame(result014,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName)
result_b2b_water_014_DF.to_csv('data/unCorCutoffRobust/result_b2b_water_014_DF.csv')




######################### different uncorrelation cutoff #############################

### threshold=0.16

del linear_learn_unCor014, result014, result_b2b_water_014_DF

def linear_learn_unCor016(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	focalUnCorInd=np.nonzero(abs(corMat3[:,i])<0.16)[0]
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###

pool = multiprocessing.Pool(processes=100)
result016 = pool.map(linear_learn_unCor016, np.arange(brain3.shape[1]))

pool.close()
pool.join()

time.localtime()




result_b2b_water_016_DF=pd.DataFrame(result016,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName)
result_b2b_water_016_DF.to_csv('data/unCorCutoffRobust/result_b2b_water_016_DF.csv')




######################### different uncorrelation cutoff #############################

### threshold=0.17

del linear_learn_unCor016, result016, result_b2b_water_016_DF

def linear_learn_unCor017(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	focalUnCorInd=np.nonzero(abs(corMat3[:,i])<0.17)[0]
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###

pool = multiprocessing.Pool(processes=100)
result017 = pool.map(linear_learn_unCor017, np.arange(brain3.shape[1]))

pool.close()
pool.join()

time.localtime()




result_b2b_water_017_DF=pd.DataFrame(result017,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName)
result_b2b_water_017_DF.to_csv('data/unCorCutoffRobust/result_b2b_water_017_DF.csv')



######################### different uncorrelation cutoff #############################

### threshold=0.18

del linear_learn_unCor017, result017, result_b2b_water_017_DF

def linear_learn_unCor018(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	focalUnCorInd=np.nonzero(abs(corMat3[:,i])<0.18)[0]
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###

pool = multiprocessing.Pool(processes=100)
result018 = pool.map(linear_learn_unCor018, np.arange(brain3.shape[1]))

pool.close()
pool.join()

time.localtime()




result_b2b_water_018_DF=pd.DataFrame(result018,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName)
result_b2b_water_018_DF.to_csv('data/unCorCutoffRobust/result_b2b_water_018_DF.csv')



######################### different uncorrelation cutoff #############################

### threshold=0.19

del linear_learn_unCor018, result018, result_b2b_water_018_DF

def linear_learn_unCor019(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	focalUnCorInd=np.nonzero(abs(corMat3[:,i])<0.19)[0]
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###

pool = multiprocessing.Pool(processes=100)
result019 = pool.map(linear_learn_unCor019, np.arange(brain3.shape[1]))

pool.close()
pool.join()

time.localtime()




result_b2b_water_019_DF=pd.DataFrame(result019,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName)
result_b2b_water_019_DF.to_csv('data/unCorCutoffRobust/result_b2b_water_019_DF.csv')




######################### different uncorrelation cutoff #############################

### threshold=0.2

del linear_learn_unCor019, result019, result_b2b_water_019_DF

def linear_learn_unCor02(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	focalUnCorInd=np.nonzero(abs(corMat3[:,i])<0.2)[0]
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###

pool = multiprocessing.Pool(processes=100)
result02 = pool.map(linear_learn_unCor02, np.arange(brain3.shape[1]))

pool.close()
pool.join()

time.localtime()




result_b2b_water_02_DF=pd.DataFrame(result02,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName)
result_b2b_water_02_DF.to_csv('data/unCorCutoffRobust/result_b2b_water_02_DF.csv')


######################################################################################

############## down sampling ###############################################

### down sample ratio = 0.9
	
def linear_learn_unCor015_down09(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	np.random.seed(9)
	downInd=np.random.choice(range(675),round(675*0.9))
	unCorrInd=np.nonzero(abs(corMat3[:,i])<0.15)[0]	
	focalUnCorInd=list(set(unCorrInd).intersection(set(downInd)))
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###
np.random.seed(9)
downInd=np.random.choice(range(675),round(675*0.9))
pool = multiprocessing.Pool(processes=100)
result_down09 = pool.map(linear_learn_unCor015_down09, downInd)

pool.close()
pool.join()

time.localtime()



result_b2b_water_015_down09_DF=pd.DataFrame(result_down09,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName[downInd])
result_b2b_water_015_down09_DF.to_csv('data/downSampling/result_b2b_water_015_down09_DF.csv')



############## down sampling ###############################################

del downInd, result_down09

### down sample ratio = 0.8
	
def linear_learn_unCor015_down08(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	np.random.seed(8)
	downInd=np.random.choice(range(675),round(675*0.8))
	unCorrInd=np.nonzero(abs(corMat3[:,i])<0.15)[0]	
	focalUnCorInd=list(set(unCorrInd).intersection(set(downInd)))
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###
np.random.seed(8)
downInd=np.random.choice(range(675),round(675*0.8))
pool = multiprocessing.Pool(processes=100)
result_down08 = pool.map(linear_learn_unCor015_down08, downInd)

pool.close()
pool.join()

time.localtime()



result_b2b_water_015_down08_DF=pd.DataFrame(result_down08,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName[downInd])
result_b2b_water_015_down08_DF.to_csv('data/downSampling/result_b2b_water_015_down08_DF.csv')





############## down sampling ###############################################
del result_down09
del downInd, result_down08

### down sample ratio = 0.7
	
def linear_learn_unCor015_down07(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	np.random.seed(7)
	downInd=np.random.choice(range(675),round(675*0.7))
	unCorrInd=np.nonzero(abs(corMat3[:,i])<0.15)[0]	
	focalUnCorInd=list(set(unCorrInd).intersection(set(downInd)))
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###
np.random.seed(7)
downInd=np.random.choice(range(675),round(675*0.7))
pool = multiprocessing.Pool(processes=100)
result_down07 = pool.map(linear_learn_unCor015_down07, downInd)

pool.close()
pool.join()

time.localtime()




result_b2b_water_015_down07_DF=pd.DataFrame(result_down07,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName[downInd])
result_b2b_water_015_down07_DF.to_csv('data/downSampling/result_b2b_water_015_down07_DF.csv')




############## down sampling ###############################################

del downInd, result_down07

### down sample ratio = 0.6
	
def linear_learn_unCor015_down06(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	np.random.seed(6)
	downInd=np.random.choice(range(675),round(675*0.6))
	unCorrInd=np.nonzero(abs(corMat3[:,i])<0.15)[0]	
	focalUnCorInd=list(set(unCorrInd).intersection(set(downInd)))
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###
np.random.seed(6)
downInd=np.random.choice(range(675),round(675*0.6))
pool = multiprocessing.Pool(processes=100)
result_down06 = pool.map(linear_learn_unCor015_down06, downInd)

pool.close()
pool.join()

time.localtime()




result_b2b_water_015_down06_DF=pd.DataFrame(result_down06,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName[downInd])
result_b2b_water_015_down06_DF.to_csv('data/downSampling/result_b2b_water_015_down06_DF.csv')




############## down sampling ###############################################

del downInd, result_down06

### down sample ratio = 0.5
	
def linear_learn_unCor015_down05(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	np.random.seed(5)
	downInd=np.random.choice(range(675),round(675*0.5))
	unCorrInd=np.nonzero(abs(corMat3[:,i])<0.15)[0]	
	focalUnCorInd=list(set(unCorrInd).intersection(set(downInd)))
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###
np.random.seed(5)
downInd=np.random.choice(range(675),round(675*0.5))
pool = multiprocessing.Pool(processes=100)
result_down05 = pool.map(linear_learn_unCor015_down05, downInd)

pool.close()
pool.join()

time.localtime()




result_b2b_water_015_down05_DF=pd.DataFrame(result_down05,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName[downInd])
result_b2b_water_015_down05_DF.to_csv('data/downSampling/result_b2b_water_015_down05_DF.csv')




############## down sampling ###############################################

del downInd, result_down05

### down sample ratio = 0.4
	
def linear_learn_unCor015_down04(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	np.random.seed(4)
	downInd=np.random.choice(range(675),round(675*0.4))
	unCorrInd=np.nonzero(abs(corMat3[:,i])<0.15)[0]	
	focalUnCorInd=list(set(unCorrInd).intersection(set(downInd)))
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###
np.random.seed(4)
downInd=np.random.choice(range(675),round(675*0.4))
pool = multiprocessing.Pool(processes=100)
result_down04 = pool.map(linear_learn_unCor015_down04, downInd)

pool.close()
pool.join()

time.localtime()




result_b2b_water_015_down04_DF=pd.DataFrame(result_down04,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName[downInd])
result_b2b_water_015_down04_DF.to_csv('data/downSampling/result_b2b_water_015_down04_DF.csv')




############## down sampling ###############################################

del downInd, result_down04

### down sample ratio = 0.3
	
def linear_learn_unCor015_down03(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	np.random.seed(3)
	downInd=np.random.choice(range(675),round(675*0.3))
	unCorrInd=np.nonzero(abs(corMat3[:,i])<0.15)[0]	
	focalUnCorInd=list(set(unCorrInd).intersection(set(downInd)))
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###
np.random.seed(3)
downInd=np.random.choice(range(675),round(675*0.3))
pool = multiprocessing.Pool(processes=100)
result_down03 = pool.map(linear_learn_unCor015_down03, downInd)

pool.close()
pool.join()

time.localtime()




result_b2b_water_015_down03_DF=pd.DataFrame(result_down03,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName[downInd])
result_b2b_water_015_down03_DF.to_csv('data/downSampling/result_b2b_water_015_down03_DF.csv')




############## down sampling ###############################################

del downInd, result_down03

### down sample ratio = 0.2
	
def linear_learn_unCor015_down02(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	np.random.seed(2)
	downInd=np.random.choice(range(675),round(675*0.2))
	unCorrInd=np.nonzero(abs(corMat3[:,i])<0.15)[0]	
	focalUnCorInd=list(set(unCorrInd).intersection(set(downInd)))
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###
np.random.seed(2)
downInd=np.random.choice(range(675),round(675*0.2))
pool = multiprocessing.Pool(processes=100)
result_down02 = pool.map(linear_learn_unCor015_down02, downInd)

pool.close()
pool.join()

time.localtime()




result_b2b_water_015_down02_DF=pd.DataFrame(result_down02,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName[downInd])
result_b2b_water_015_down02_DF.to_csv('data/downSampling/result_b2b_water_015_down02_DF.csv')




############## down sampling ###############################################

del downInd, result_down02

### down sample ratio = 0.1
	
def linear_learn_unCor015_down01(i):
	print("PID=",os.getpid(),"i=",i)
	m=ElasticNet(n_splits=10,n_lambda=200,scoring='r2')
	np.random.seed(1)
	downInd=np.random.choice(range(675),round(675*0.1))
	unCorrInd=np.nonzero(abs(corMat3[:,i])<0.15)[0]	
	focalUnCorInd=list(set(unCorrInd).intersection(set(downInd)))
	fit=m.fit(brain3[:,focalUnCorInd],brain3[:,i])
	lambda_max_ind=np.argwhere(fit.lambda_path_==fit.lambda_max_)[0,0]
	coef=np.zeros(brain3.shape[1])
	coef[focalUnCorInd]=fit.coef_path_[:,lambda_max_ind]
	return([fit.cv_mean_score_[lambda_max_ind],fit.cv_standard_error_[lambda_max_ind],fit.intercept_path_[lambda_max_ind]]+coef.reshape(1,brain3.shape[1]).tolist()[0])




###
np.random.seed(1)
downInd=np.random.choice(range(675),round(675*0.1))
pool = multiprocessing.Pool(processes=100)
result_down01 = pool.map(linear_learn_unCor015_down01, downInd)

pool.close()
pool.join()

time.localtime()




result_b2b_water_015_down01_DF=pd.DataFrame(result_down01,columns=['Pred','Pred_se']+['Intercept']+list(brainColName),index=brainColName[downInd])
result_b2b_water_015_down01_DF.to_csv('data/downSampling/result_b2b_water_015_down01_DF.csv')

