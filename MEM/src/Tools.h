#ifndef _Tools_h
#define _Tools_h

#include "Permutations.h"

void CombineHypotheses(Permutations Perm1, Permutations Perm2, double *weight_mean, double* weight_log, double *weight_err, float* weight_chi2, float* weight_time, double *weight_avg, double* weight_max, double* weight_logmean, double* weight_kin_log, double* weight_kin_logint, double* weight_kinmax, double* weight_kinmaxint, double* weight_JEC_up, double* weight_JEC_down, double* weight_JER_up, double* weight_JER_down){

  cout << "CombineHypotheses"<<endl;

  double weight = 0;

  *weight_mean = 0;
  *weight_log = 0;
  *weight_err = 0;
  *weight_time = 0;
  *weight_avg = 0;
  *weight_max = 0;
  *weight_logmean = 0;

  double nHypAllowed = (double)(Perm1.nHypAllowed+Perm2.nHypAllowed);
  double nNullResult = (double)(Perm1.nNullResult+Perm2.nNullResult);

  for (unsigned int i=0; i<Perm1.resMEM_all.size(); i++){
    weight += Perm1.resMEM_all.at(i).weight;
    if (Perm1.resMEM_all.at(i).weight>0) *weight_logmean += log(Perm1.resMEM_all.at(i).weight);
    if (Perm1.resMEM_all.at(i).weight>0) *weight_err += Perm1.resMEM_all.at(i).err * Perm1.resMEM_all.at(i).err;
    if (Perm1.resMEM_all.at(i).weight>0) *weight_time += Perm1.resMEM_all.at(i).time;
    if (Perm1.resMEM_all.at(i).weight>0) *weight_chi2 += Perm1.resMEM_all.at(i).chi2;
    if (*weight_max < Perm1.resMEM_all.at(i).weight) *weight_max = Perm1.resMEM_all.at(i).weight;
  }

  for (unsigned int i=0; i<Perm2.resMEM_all.size(); i++){
    weight += Perm2.resMEM_all.at(i).weight;
    if (Perm2.resMEM_all.at(i).weight>0) *weight_logmean += log(Perm2.resMEM_all.at(i).weight);
    if (Perm2.resMEM_all.at(i).weight>0) *weight_err += Perm2.resMEM_all.at(i).err * Perm2.resMEM_all.at(i).err;
    if (Perm2.resMEM_all.at(i).weight>0) *weight_time += Perm2.resMEM_all.at(i).time;
    if (Perm2.resMEM_all.at(i).weight>0) *weight_chi2 += Perm2.resMEM_all.at(i).chi2;
    if (*weight_max < Perm2.resMEM_all.at(i).weight) *weight_max = Perm2.resMEM_all.at(i).weight;
  }

  for (unsigned int i=0; i<Perm1.resMEM_all_JEC_up.size(); i++) *weight_JEC_up += Perm1.resMEM_all_JEC_up.at(i).weight;
  for (unsigned int i=0; i<Perm1.resMEM_all_JEC_down.size(); i++) *weight_JEC_down += Perm1.resMEM_all_JEC_down.at(i).weight;
  for (unsigned int i=0; i<Perm1.resMEM_all_JER_up.size(); i++) *weight_JER_up += Perm1.resMEM_all_JER_up.at(i).weight;
  for (unsigned int i=0; i<Perm1.resMEM_all_JER_down.size(); i++) *weight_JER_down += Perm1.resMEM_all_JER_down.at(i).weight;

  for (unsigned int i=0; i<Perm2.resMEM_all_JEC_up.size(); i++) *weight_JEC_up += Perm2.resMEM_all_JEC_up.at(i).weight;
  for (unsigned int i=0; i<Perm2.resMEM_all_JEC_down.size(); i++) *weight_JEC_down += Perm2.resMEM_all_JEC_down.at(i).weight;
  for (unsigned int i=0; i<Perm2.resMEM_all_JER_up.size(); i++) *weight_JER_up += Perm2.resMEM_all_JER_up.at(i).weight;
  for (unsigned int i=0; i<Perm2.resMEM_all_JER_down.size(); i++) *weight_JER_down += Perm2.resMEM_all_JER_down.at(i).weight;

  *weight_mean = weight / nHypAllowed;
  *weight_avg = weight / (nHypAllowed+nNullResult);
  *weight_logmean = *weight_logmean / nHypAllowed;

  *weight_JEC_up /= nHypAllowed;
  *weight_JEC_down /= nHypAllowed;
  *weight_JER_up /= nHypAllowed;
  *weight_JER_down /= nHypAllowed;

  if (!(*weight_max > 0)) *weight_max = 1e-300;
  if (!(*weight_mean > 0)) *weight_mean = 1e-300;
  if (!(*weight_avg > 0)) *weight_avg = 1e-300;
  if (!(*weight_logmean > 0)) *weight_logmean = log(1e-300);
 
  if (!(*weight_JEC_up > 0)) *weight_JEC_up = 1e-300;
  if (!(*weight_JEC_down > 0)) *weight_JEC_down = 1e-300;
  if (!(*weight_JER_up > 0)) *weight_JER_up = 1e-300;
  if (!(*weight_JER_down > 0)) *weight_JER_down = 1e-300;


  *weight_log = log(*weight_mean); 

  *weight_err = sqrt(*weight_err)/nHypAllowed;
  *weight_chi2 /= nHypAllowed;
  
  *weight_kin_log = (Perm1.resKin_maxKinFit.weight > Perm2.resKin_maxKinFit.weight) ? Perm1.resKin_maxKinFit.weight : Perm2.resKin_maxKinFit.weight;
  *weight_kin_log = log(*weight_kin_log);

  *weight_kin_logint = (Perm1.resKin_maxKinFit_Int.weight > Perm2.resKin_maxKinFit_Int.weight) ? Perm1.resKin_maxKinFit_Int.weight : Perm2.resKin_maxKinFit_Int.weight;
  *weight_kin_logint = log(*weight_kin_logint);

  *weight_kinmax = (Perm1.resKin_maxKinFit.weight > Perm2.resKin_maxKinFit.weight) ? Perm1.resMEM_maxKinFit.weight : Perm2.resMEM_maxKinFit.weight;

  *weight_kinmaxint = (Perm1.resKin_maxKinFit_Int.weight > Perm2.resKin_maxKinFit_Int.weight) ? Perm1.resMEM_maxKinFit_Int.weight : Perm2.resMEM_maxKinFit_Int.weight;

  return;
}

void FillWeightVectors(Permutations Perm, std::vector<double>* MEAllWeights, std::vector<float>* MEAllWeights_log){

  cout << "FillWeightVectors" << endl;

  for (unsigned int i=0; i<Perm.resMEM_all.size(); i++)
  if (Perm.resMEM_all.at(i).weight>0){
    (*MEAllWeights).push_back(Perm.resMEM_all.at(i).weight);
    (*MEAllWeights_log).push_back(log(Perm.resMEM_all.at(i).weight));
  }
  else{
    (*MEAllWeights).push_back(1e-300);
    (*MEAllWeights_log).push_back(log(1e-300));
  }

  return;
}

#endif
