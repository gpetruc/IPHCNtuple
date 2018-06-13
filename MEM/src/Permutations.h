#ifndef PERMUTATIONS_H
#define PERMUTATIONS_H

#include "HypIntegrator.h"
#include "MultiLepton.h"

using namespace std;

class Permutations {

  public:

  MultiLepton multiLepton;

  int Hypothesis;
  string HypothesisName;

  bool computeHyp;
  
  bool doMinimization;

  bool doPermutationLep;
  bool doPermutationJet;
  bool doPermutationBjet;
  int nPermutationJetSyst;

  int stage;
  int iterationNumber;
 
  IntegrationResult resMEM_max;
  IntegrationResult resMEM_maxKinFit;
  IntegrationResult resKin_maxKinFit;
  IntegrationResult resMEM_maxKinFit_Int;
  IntegrationResult resKin_maxKinFit_Int;

  std::vector<IntegrationResult> resMEM_all;
  std::vector<IntegrationResult> resMEM_all_JEC_up; 
  std::vector<IntegrationResult> resMEM_all_JEC_down; 
  std::vector<IntegrationResult> resMEM_all_JER_up;
  std::vector<IntegrationResult> resMEM_all_JER_down;
  IntegrationResult resMEM_logavgExl0;
  IntegrationResult resMEM_avgExl0;
  IntegrationResult resMEM_avgIncl0;
  IntegrationResult resMEM_avgExl0_JEC_up;
  IntegrationResult resMEM_avgExl0_JEC_down;
  IntegrationResult resMEM_avgExl0_JER_up;
  IntegrationResult resMEM_avgExl0_JER_down;

  int nHypAllowed;
  int nNullResult;

  Permutations();
  ~Permutations(); 
  void SetMultiLepton(MultiLepton*, HypIntegrator*);
  int InitializeHyp(HypIntegrator*, int, int, string, bool, string, int);
  void LoopPermutations(HypIntegrator*);

  IntegrationResult GetMEMKinFitResult(HypIntegrator*, int);
  void GetMEMResult_Syst(HypIntegrator*, IntegrationResult**);
  IntegrationResult GetMEMResult(HypIntegrator*);

  void UpdateAverageWeight(IntegrationResult*, std::vector<IntegrationResult>, bool, string);

  private:

  double xs;
  double xsTTH;
  double xsTTLL;
  double xsTTW;
  double xsTTbar;
  double xsTLLJ;

  IntegrationResult res;
  IntegrationResult* res_syst;
  IntegrationResult resMinimized;

};

Permutations::Permutations(){

  nHypAllowed = 0;
  nNullResult = 0;
  computeHyp = false;
}

Permutations::~Permutations(){
}

void Permutations::SetMultiLepton(MultiLepton* multilepton_, HypIntegrator* hypIntegrator){

  multiLepton = *multilepton_;
  cout << "Nleptons="<< multiLepton.Leptons.size()<<", Category: "<<multiLepton.kCatJets<<endl;
  
  (*hypIntegrator).meIntegrator->SetNleptonMode(multiLepton.Leptons.size());

  return;
}

int Permutations::InitializeHyp(HypIntegrator* hypIntegrator, int hyp, int nPointsHyp, string shyp, bool doMinimization_, string JetChoice, int nPermutationJetSyst_){

  cout << "InitializeHyp " << hyp << " (" << shyp << "), nPoints=" <<  nPointsHyp << " doMinimization="<<doMinimization_<<", JetChoice=" << JetChoice << ", nsyst="<<nPermutationJetSyst_<<endl;

  nHypAllowed = 0;
  nNullResult = 0;

  resMEM_max.weight = 0;
  res.weight = 0;
  resKin_maxKinFit.weight = 0;
  resMEM_maxKinFit.weight = 0;
  resKin_maxKinFit_Int.weight = 0;
  resMEM_maxKinFit_Int.weight = 0;
  if (nPermutationJetSyst_>=1) res_syst = new IntegrationResult[nPermutationJetSyst_];

  resMEM_all.clear();
  resMEM_all_JEC_up.clear();
  resMEM_all_JEC_down.clear();
  resMEM_all_JER_up.clear();
  resMEM_all_JER_down.clear();

  resMEM_logavgExl0.weight = 0;
  resMEM_avgExl0.weight = 0;
  resMEM_avgIncl0.weight = 0;
  resMEM_avgExl0_JEC_up.weight = 0;
  resMEM_avgExl0_JEC_down.weight = 0;
  resMEM_avgExl0_JER_up.weight = 0;
  resMEM_avgExl0_JER_down.weight = 0;

  computeHyp = false;

  xsTTH = (*hypIntegrator).meIntegrator->xsTTH * (*hypIntegrator).meIntegrator->brTopHad * (*hypIntegrator).meIntegrator->brTopLep * (*hypIntegrator).meIntegrator->brHiggsFullLep;
  xsTTLL = (*hypIntegrator).meIntegrator->xsTTLL * (*hypIntegrator).meIntegrator->brTopHad * (*hypIntegrator).meIntegrator->brTopLep;
  xsTTW = (*hypIntegrator).meIntegrator->xsTTW * (*hypIntegrator).meIntegrator->brTopLep * (*hypIntegrator).meIntegrator->brTopLep;
  xsTTbar = (*hypIntegrator).meIntegrator->xsTTbar * (*hypIntegrator).meIntegrator->brTopHad * (*hypIntegrator).meIntegrator->brTopLep;
  xsTLLJ = (*hypIntegrator).meIntegrator->xsTLLJ * (*hypIntegrator).meIntegrator->brTopLep;

  if (hyp==kMEM_TTH_TopAntitopHiggsDecay || hyp==kMEM_TTH_TopAntitopHiggsSemiLepDecay) xs = xsTTH;
  if (hyp==kMEM_TTLL_TopAntitopDecay) xs = xsTTLL;
  if (hyp==kMEM_TTW_TopAntitopDecay || hyp==kMEM_TTWJJ_TopAntitopDecay) xs = xsTTW;
  if (hyp==kMEM_TTbar_TopAntitopSemiLepDecay || hyp==kMEM_TTbar_TopAntitopFullyLepDecay) xs = xsTTbar;
  if (hyp==kMEM_TLLJ_TopLepDecay) xs = xsTLLJ;

  cout << "Process "<<hyp<< " (" << shyp<< ")" << " xs="<< xs<< endl;

  Hypothesis = hyp;
  HypothesisName = shyp;

  doMinimization = doMinimization_;

  nPermutationJetSyst = nPermutationJetSyst_;

  //Setup integrator
  (*hypIntegrator).SetupIntegrationHypothesis(hyp, multiLepton.kCatJets, nPointsHyp);
  if (doMinimization) (*hypIntegrator).SetupMinimizerHypothesis(hyp, multiLepton.kCatJets, 0, nPointsHyp);

  //Check if MEM can be computed for a given hypothesis
  if (multiLepton.Leptons.size()==4 && (hyp!=kMEM_TTH_TopAntitopHiggsDecay && hyp!=kMEM_TTLL_TopAntitopDecay && hyp!=kMEM_TTbar_TopAntitopFullyLepDecay)) return 0;
  if (multiLepton.Leptons.size()==2 && (hyp!=kMEM_TTH_TopAntitopHiggsSemiLepDecay && hyp!=kMEM_TTW_TopAntitopDecay && hyp!=kMEM_TTbar_TopAntitopSemiLepDecay)) return 0;
  //if (multiLepton.Leptons.size()==3 && hyp==kMEM_TLLJ_TopLepDecay && multiLepton.kCatJets!=kCat_3l_1b_1j) return 0;

  //Select the non b-jets
     if (multiLepton.kCatJets==kCat_3l_2b_2j || multiLepton.kCatJets==kCat_3l_1b_2j){
       if (JetChoice=="HighestPt") multiLepton.SwitchJetsFromAllJets(kJetPair_HighestPt);
       else if (JetChoice=="MjjLowest") multiLepton.SwitchJetsFromAllJets(kJetPair_MjjLowest);
       else if (JetChoice=="MwClosest") multiLepton.SwitchJetsFromAllJets(kJetPair_MwClosest);
       else if (JetChoice=="Mixed") {
         if (hyp==kMEM_TTH_TopAntitopHiggsSemiLepDecay) multiLepton.SwitchJetsFromAllJets(kJetPair_MjjLowest);
         else if (hyp==kMEM_TTWJJ_TopAntitopDecay) multiLepton.SwitchJetsFromAllJets(kJetPair_HighestPt);
         else multiLepton.SwitchJetsFromAllJets(kJetPair_MwClosest);
       }
       if (hyp==kMEM_TLLJ_TopLepDecay) multiLepton.SwitchJetsFromAllJets(kJetSingle);
     }
     else if (multiLepton.kCatJets==kCat_3l_2b_1j || multiLepton.kCatJets==kCat_3l_1b_1j) multiLepton.SwitchJetsFromAllJets(kJetSingle);
     else if (multiLepton.kCatJets==kCat_2lss_2b_4j || multiLepton.kCatJets==kCat_2lss_1b_4j){
       if (JetChoice=="Mixed"){
         if (hyp==kMEM_TTH_TopAntitopHiggsSemiLepDecay) multiLepton.SwitchJetsFromAllJets(kTwoJetPair_MwClosest_2ndMjjLowest);
         else if (hyp==kMEM_TTW_TopAntitopDecay || hyp==kMEM_TTbar_TopAntitopSemiLepDecay) multiLepton.SwitchJetsFromAllJets(kJetPair_MwClosest);
         else if (hyp==kMEM_TTWJJ_TopAntitopDecay) multiLepton.SwitchJetsFromAllJets(kTwoJetPair_MwClosest_2ndHighestPt);
       }
       else if (JetChoice=="MwClosest") {
         if (hyp==kMEM_TTW_TopAntitopDecay || hyp==kMEM_TTbar_TopAntitopSemiLepDecay) multiLepton.SwitchJetsFromAllJets(kJetPair_MwClosest);
         else multiLepton.SwitchJetsFromAllJets(kTwoJetPair_MwClosest_2ndMwClosest);
       }
     }
     else if (multiLepton.kCatJets==kCat_2lss_2b_3j || multiLepton.kCatJets==kCat_2lss_1b_3j) {
         if (hyp==kMEM_TTH_TopAntitopHiggsSemiLepDecay) multiLepton.SwitchJetsFromAllJets(kThreeJets_MwClosest_HighestPt);
         else if (hyp==kMEM_TTW_TopAntitopDecay || hyp==kMEM_TTbar_TopAntitopSemiLepDecay) multiLepton.SwitchJetsFromAllJets(kJetPair_MwClosest);
     }
     else if (multiLepton.kCatJets==kCat_2lss_2b_2j) {
       multiLepton.SwitchJetsFromAllJets(kJetPair_MwClosest);
     }

  doPermutationLep = true;
  doPermutationJet = true;
  doPermutationBjet = true;

  //Check if hypothesis can do jet permutation
  if ((multiLepton.Leptons.size()==3 && hyp==kMEM_TTW_TopAntitopDecay)
     //|| (hyp==kMEM_TLLJ_TopLepDecay && (multiLepton.kCatJets==kCat_3l_1b_1j || multiLepton.kCatJets==kCat_3l_2b_1j))
     || hyp==kMEM_TLLJ_TopLepDecay
     || hyp==kMEM_TTbar_TopAntitopFullyLepDecay 
     || multiLepton.kCatJets==kCat_3l_2b_0j 
     || multiLepton.Leptons.size()==4)
	doPermutationJet=false;

  if (hyp==kMEM_TLLJ_TopLepDecay && (multiLepton.kCatJets==kCat_3l_1b_1j || multiLepton.kCatJets==kCat_3l_1b_2j)) doPermutationBjet = false;

  //Sort leptons and jets by Pt
  multiLepton.DoSort(&multiLepton.Leptons);
  multiLepton.DoSort(&multiLepton.Bjets);
  multiLepton.DoSort(&multiLepton.Jets);

  computeHyp = true;

  return 1;
}

void Permutations::LoopPermutations(HypIntegrator* hypIntegrator){

  int combcheck = 1;

  int permreslep = 1;
  int permresjet = 1;
  int permresbjet = 1;

  stage = 0;
  iterationNumber = 5;

  int iperm=0;

  cout << "Start looping" << endl;

     multiLepton.DoSort(&multiLepton.Bjets);
     do {
       multiLepton.DoSort(&multiLepton.Jets);
       do {
         multiLepton.DoSort(&multiLepton.Leptons);
         do {

             combcheck = multiLepton.CheckPermutationHyp(Hypothesis);

             if (combcheck) {
               for (unsigned int il=0; il<multiLepton.Leptons.size(); il++) cout << " Lepton"<< il<<"Id="<<multiLepton.Leptons.at(il).Id; 
               cout<<endl;
               for (unsigned int ib=0; ib<multiLepton.Bjets.size(); ib++) cout << " Bjet"<< ib<<"Pt="<<multiLepton.Bjets.at(ib).P4.Pt(); 
               cout<<endl;
               for (unsigned int ij=0; ij<multiLepton.Jets.size(); ij++) cout << " Jet"<< ij<<"Pt="<<multiLepton.Jets.at(ij).P4.Pt();
               cout<<endl;

               multiLepton.FillParticlesHypothesis(Hypothesis, &((*hypIntegrator).meIntegrator));
               multiLepton.SwitchJetSyst(0);

	       //Use MEM as kinematic fit, minimization from random initial value of integration variables
	       if (doMinimization) resMinimized = GetMEMKinFitResult(hypIntegrator, 10);

	       //Compute MEM weight
	       res = GetMEMResult(hypIntegrator);
	       if (resMEM_max.weight < res.weight) resMEM_max = res;
	       resMEM_all.push_back(res);
	
	       //Save MEM weight and Kin Fit result for best Kin Fit (minimization from initial value of integration variables)
	       if (doMinimization) {
                 if (resMinimized.weight/xs < resKin_maxKinFit.weight) {
                   resKin_maxKinFit.weight = resMinimized.weight / xs;
                   resMEM_maxKinFit = res;
                 }
	       }

	       //Save MEM weight and Kin Fit result for best Kin Fit (minimization performed along with MEM integration)
               if ((*hypIntegrator).meIntegrator->weight_max/xs > resKin_maxKinFit_Int.weight) {
                 cout << "better kinweight_maxint"<<endl;
                 resKin_maxKinFit_Int.weight = (*hypIntegrator).meIntegrator->weight_max / xs;
                 resMEM_maxKinFit_Int = res;
               }

	       //MEM Systematics JEC/JER
               if (nPermutationJetSyst>=1) GetMEMResult_Syst(hypIntegrator, &res_syst);

	     }


             iperm++;

           if (doPermutationLep) {
             if (Hypothesis!=kMEM_TTbar_TopAntitopSemiLepDecay) permreslep = multiLepton.DoPermutation(&multiLepton.Leptons);
             else permreslep = multiLepton.DoPermutationLinear(&multiLepton.Leptons);
           }
         } while (doPermutationLep && permreslep);

         if (doPermutationJet) {
           if (multiLepton.kCatJets!=kCat_3l_2b_1j && multiLepton.kCatJets!=kCat_3l_1b_1j && multiLepton.kCatJets!=kCat_2lss_2b_3j && multiLepton.kCatJets!=kCat_2lss_1b_3j && multiLepton.kCatJets!=kCat_2lss_2b_2j) permresjet = multiLepton.DoPermutation(&multiLepton.Jets);
           else permresjet = multiLepton.DoPermutationMissingJet("jet");
         }
       } while (doPermutationJet && permresjet);

       if (doPermutationBjet) {
	 if (Hypothesis==kMEM_TLLJ_TopLepDecay) permresbjet = multiLepton.DoPermutationLinear(&multiLepton.Bjets);
         else if (multiLepton.kCatJets!=kCat_3l_1b_2j && multiLepton.kCatJets!=kCat_3l_1b_1j && multiLepton.kCatJets!=kCat_4l_1b && multiLepton.kCatJets!=kCat_2lss_1b_4j && multiLepton.kCatJets!=kCat_2lss_1b_3j) permresbjet = multiLepton.DoPermutation(&multiLepton.Bjets);
         else permresbjet = multiLepton.DoPermutationMissingJet("bjet");
       }
     } while (doPermutationBjet && permresbjet);


  if (resMEM_max.weight==0) resMEM_max.weight = 1e-300;

  if (resKin_maxKinFit.weight==0) resMEM_maxKinFit.weight = 1e-300;
  if (resKin_maxKinFit_Int.weight==0) resMEM_maxKinFit_Int.weight = 1e-300; 

  if (resKin_maxKinFit.weight==0) resKin_maxKinFit.weight = 1e-300; 
  if (resKin_maxKinFit_Int.weight==0) resKin_maxKinFit_Int.weight = 1e-300;

  UpdateAverageWeight(&resMEM_avgExl0, resMEM_all, true, "avg");
  UpdateAverageWeight(&resMEM_avgIncl0, resMEM_all, false, "avg");
  UpdateAverageWeight(&resMEM_logavgExl0,resMEM_all,true,"logavg");

  if (nPermutationJetSyst>=1) UpdateAverageWeight(&resMEM_avgExl0_JEC_up, resMEM_all_JEC_up, true, "avg");
  if (nPermutationJetSyst>=2) UpdateAverageWeight(&resMEM_avgExl0_JEC_down, resMEM_all_JEC_down, true, "avg");
  if (nPermutationJetSyst>=3) UpdateAverageWeight(&resMEM_avgExl0_JER_up, resMEM_all_JER_up, true, "avg");
  if (nPermutationJetSyst>=4) UpdateAverageWeight(&resMEM_avgExl0_JER_down, resMEM_all_JER_down, true, "avg");

}

IntegrationResult Permutations::GetMEMKinFitResult(HypIntegrator* hypIntegrator, int ntry){

  IntegrationResult resMin;

  double kinresweightmin = 1000;
  double * var;
  
  (*hypIntegrator).meIntegrator->SetMinimization(1);
             
  (*hypIntegrator).meIntegrator->weight_max = 0;

  for (int itry=0; itry<ntry; itry++){

    var = (*hypIntegrator).FindMinimizationiInitialValues(multiLepton.xL, multiLepton.xU);
    resMin = (*hypIntegrator).DoMinimization(multiLepton.xL, multiLepton.xU, var);

    cout << "KIN TRY Hyp "<< HypothesisName <<" Vegas Ncall="<<(*hypIntegrator).nPointsCatHyp <<" -log(max)=" << resMin.weight<<" Time(s)="<<resMin.time<<endl;
 
    if ( resMin.weight < kinresweightmin) kinresweightmin = resMin.weight;
  }

  (*hypIntegrator).meIntegrator->SetMinimization(0);

  resMin.weight = exp(-kinresweightmin);

  cout << "KIN FINAL Hyp "<< HypothesisName <<" Vegas Ncall="<<(*hypIntegrator).nPointsCatHyp <<" max=" << resMin.weight<<" Time(s)="<<resMin.time<<endl;

  if (!(resMin.weight>0)) {
     resMin.weight = 0;
  }

  return resMin;
}

IntegrationResult Permutations::GetMEMResult(HypIntegrator* hypIntegrator){

  IntegrationResult resMEM;

  stage = 0;
  iterationNumber = 5;

  (*hypIntegrator).meIntegrator->weight_max = 0;

  resMEM = (*hypIntegrator).DoIntegration(multiLepton.xL, multiLepton.xU, stage, iterationNumber);

  cout << "MEM Hyp "<< HypothesisName<<" Vegas Ncall="<<(*hypIntegrator).nPointsCatHyp <<" Cross section (pb) : " << resMEM.weight<< " +/- "<< resMEM.err<<" chi2/ndof="<< resMEM.chi2<<" Time(s)="<<resMEM.time<<endl;
  cout << "KIN from MEM Hyp "<< HypothesisName<<" Vegas Ncall="<<(*hypIntegrator).nPointsCatHyp <<" max=" << (*hypIntegrator).meIntegrator->weight_max<<" Time(s)="<<resMEM.time<<endl;

  if (resMEM.weight>0){
    resMEM.weight /= xs;
    resMEM.err /= xs;
    nHypAllowed++;
  }
  else {
    resMEM.weight = 0;
    resMEM.err = 0;
    resMEM.chi2 = 0;
    cout << "Setting result weight to almost 0:"<< resMEM.weight << endl;
    nNullResult++;
  }

  cout << "After norm, Cross section (pb) : " << resMEM.weight<< " +/- "<< resMEM.err << endl;

  return resMEM;
}

void Permutations::GetMEMResult_Syst(HypIntegrator* hypIntegrator, IntegrationResult** resMEM_syst){

    stage = 1;
    iterationNumber = 1;

    for (int iSyst=1; iSyst<=nPermutationJetSyst; iSyst++){

      multiLepton.SwitchJetSyst(iSyst);

      multiLepton.FillParticlesHypothesis(Hypothesis, &((*hypIntegrator).meIntegrator));

      (*resMEM_syst)[iSyst-1] = (*hypIntegrator).DoIntegration(multiLepton.xL, multiLepton.xU, stage, iterationNumber);
      (*resMEM_syst)[iSyst-1].weight /= xs;
      (*resMEM_syst)[iSyst-1].err /= xs;

      cout << "MEM Hyp, SYST"<< iSyst-1<<" "<< HypothesisName<<" Vegas Ncall="<<(*hypIntegrator).nPointsCatHyp <<" Cross section (pb) : " << (*resMEM_syst)[iSyst-1].weight<< " +/- "<< (*resMEM_syst)[iSyst-1].err<<" chi2/ndof="<< (*resMEM_syst)[iSyst-1].chi2<<" Time(s)="<<(*resMEM_syst)[iSyst-1].time<<endl;

      if ((*resMEM_syst)[iSyst-1].weight>0);
      else {
        (*resMEM_syst)[iSyst-1].weight = 0;
        (*resMEM_syst)[iSyst-1].err = 0;
        (*resMEM_syst)[iSyst-1].chi2 = 0;
      }
    }

  return;
}

void Permutations::UpdateAverageWeight(IntegrationResult* ResAvg, std::vector<IntegrationResult> Res_all, bool doExclude0weight, string mode){

  (*ResAvg).weight = 0;
  (*ResAvg).err = 0;
  (*ResAvg).time = 0;
  (*ResAvg).chi2 = 0;

  for (unsigned int i=0; i<Res_all.size(); i++){
    if (!(Res_all.at(i).weight>0) && doExclude0weight) continue;
    if (mode=="avg") (*ResAvg).weight += Res_all.at(i).weight;
    if (mode=="logavg" && Res_all.at(i).weight>0) (*ResAvg).weight += log(Res_all.at(i).weight);
    (*ResAvg).err += Res_all.at(i).err * Res_all.at(i).err;
    (*ResAvg).time += Res_all.at(i).time;
    (*ResAvg).chi2 += Res_all.at(i).chi2;
  }

  double nWeightsTot = 0;
  if (doExclude0weight) nWeightsTot = (double)nHypAllowed;
  else nWeightsTot = (double)(nHypAllowed+nNullResult);
 
  (*ResAvg).weight /= nWeightsTot;
  (*ResAvg).err = sqrt((*ResAvg).err)/nWeightsTot;
  //(*ResAvg).time /= nWeightsTot;
  (*ResAvg).chi2 /= nWeightsTot; 

  if (!((*ResAvg).weight>0)) (*ResAvg).weight = 1e-300;
  if (!((*ResAvg).weight>0) && mode=="logavg") (*ResAvg).weight = log(1e-300);

  return;
}

#endif
