#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TInterpreter.h"

#include "Math/Integrator.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/IOptions.h"
#include "Math/IntegratorOptions.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/AdaptiveIntegratorMultiDim.h"

#include "MEPhaseSpace.h"
#include "HypIntegrator.h"
#include "MultiLepton.h"
#include "ReadGenFlatTree.h"
#include "ConfigParser.h"
#include "Tools.h"
#include "Permutations.h"

#include <iostream>
#include <ctime>

using namespace std;

const float mZ = 91.18;
const float mW = 80.38;

int main(int argc, char *argv[])
{

  gInterpreter->EnableAutoLoading();

  string InputFileName = string(argv[1]);
  ReadGenFlatTree tree;

  int evMax = atoi(argv[3]);
  cout << "Will run on "<<argv[1]<<" from event "<<argv[2]<<" to "<<evMax<<" with option "<< argv[4]<<endl;

  int nhyp;
  string* shyp;
  int* hyp;
  int* nPointsHyp;
  int* index_hyp;

  ConfigParser cfgParser;
  cfgParser.GetConfigFromFile(string(argv[5]));
  cfgParser.LoadHypotheses(&nhyp, &shyp, &hyp, &nPointsHyp, &index_hyp);

  MultiLepton multiLepton;
  cfgParser.LoadIntegrationRange(&multiLepton.JetTFfracmin, &multiLepton.JetTFfracmax, &multiLepton.NeutMaxE);


  if (strcmp(argv[4],"--dryRun")==0){

    tree.InitializeDryRun(InputFileName);
    if (atoi(argv[3])==-1 || evMax>tree.tInput->GetEntries()) evMax = tree.tInput->GetEntries();

    cout << "Start running on events with --dryRun option"<<endl;
    for (Long64_t iEvent=atoi(argv[2]); iEvent<evMax; iEvent++){

      if (iEvent % 1000 == 0) cout << "Event "<<iEvent<<endl;

      tree.mc_event = iEvent;
      tree.FillGenMultilepton(iEvent, &multiLepton);

      bool sel = tree.ApplyGenSelection(iEvent, &multiLepton);
      if (sel==0) continue;

      tree.WriteMultilepton(&multiLepton);
   
      tree.tOutput->Fill();
    }

    tree.tOutput->Write();
    tree.fOutput->Close();
  }
  else if (strcmp(argv[4],"--MEMRun")==0){ 
 
    cout << "Test MEM weight computation" <<endl;

    int index[7];
    for (int ih=0; ih<nhyp; ih++){
      if (shyp[ih]=="TTLL") index[ih] = 0;
      if (shyp[ih]=="TTHfl") index[ih] = 1;
      if (shyp[ih]=="TTHsl") index[ih] = 2;
      if (shyp[ih]=="TTW") index[ih] = 3;
      if (shyp[ih]=="TTWJJ") index[ih] = 4;
      if (shyp[ih]=="TTbarfl") index[ih] = 5;
      if (shyp[ih]=="TTbarsl") index[ih] = 6;
    }
 
    int doOptim;
    cfgParser.LoadOptim(&doOptim);
    int doOptimTopHad, doOptimTopLep, doOptimHiggs, doOptimW;
    cfgParser.LoadOptim(&doOptimTopHad, &doOptimTopLep, &doOptimHiggs, &doOptimW);

    HypIntegrator hypIntegrator;
    hypIntegrator.InitializeIntegrator(&cfgParser);

    int doMinimization = cfgParser.valDoMinimization;

    string JetChoice;
    cfgParser.LoadJetChoice(&JetChoice);

    int nPermutationJetSyst = cfgParser.nJetSyst;

    double xsTTH = hypIntegrator.meIntegrator->xsTTH * hypIntegrator.meIntegrator->brTopHad * hypIntegrator.meIntegrator->brTopLep * hypIntegrator.meIntegrator->brHiggsFullLep;
    double xsTTLL = hypIntegrator.meIntegrator->xsTTLL * hypIntegrator.meIntegrator->brTopHad * hypIntegrator.meIntegrator->brTopLep;
    double xsTTW = hypIntegrator.meIntegrator->xsTTW * hypIntegrator.meIntegrator->brTopLep * hypIntegrator.meIntegrator->brTopLep;
    double xsTTbar = hypIntegrator.meIntegrator->xsTTbar * hypIntegrator.meIntegrator->brTopHad * hypIntegrator.meIntegrator->brTopLep;

    tree.InitializeMEMRun(InputFileName);
    if (atoi(argv[3])==-1 || evMax>tree.tInput->GetEntries()) evMax = tree.tInput->GetEntries();

  int nEvents = evMax - atoi(argv[2]);
  int nErrHist = nEvents*12*nhyp;
  tree.hMEPhaseSpace_Error = new TH1F*[nErrHist];
  tree.hMEPhaseSpace_ErrorTot = new TH1F*[nhyp];
  tree.hMEPhaseSpace_ErrorTot_Pass = new TH1F*[nhyp];
  tree.hMEPhaseSpace_ErrorTot_Fail = new TH1F*[nhyp];
  string snum1[12];
  stringstream ss1[12];
  for (int i=0; i<12; i++) {
    ss1[i] << i;
    ss1[i] >> snum1[i];
  }
  int i0=0;
  string snum2[nEvents];
  stringstream ss2[nEvents];
  for (int i=atoi(argv[2]); i<evMax; i++) {
    ss2[i0] << i0;
    ss2[i0] >> snum2[i0];
    i0++;
  }
  string sname="";
  int it=0, ie0=0;
  for (int ie=atoi(argv[2]); ie<evMax; ie++){
    for (int ip=0; ip<12; ip++){
      for (int ih=0; ih<nhyp; ih++){
        sname = "hErrHist_ev" + snum2[ie0] + "_perm" + snum1[ip] + "_hyp" + shyp[ih];
        tree.hMEPhaseSpace_Error[it] = new TH1F(sname.c_str(), sname.c_str(), 10, 0, 10);
        cout << "Creating error histo, event "<<ie<<", permutation "<<ip<<", hyp "<<ih<<" "<<tree.hMEPhaseSpace_Error[it]->GetName() <<endl;
        it++;
      }
    }
    ie0++;
  }
  for (int ih=0; ih<nhyp; ih++){
    sname = "hErrHist_Tot_hyp" + shyp[ih];
    tree.hMEPhaseSpace_ErrorTot[ih] = new TH1F(sname.c_str(), sname.c_str(), 10, 0, 10);
    sname = "hErrHist_Tot_hyp" + shyp[ih] + "_Pass";
    tree.hMEPhaseSpace_ErrorTot_Pass[ih] = new TH1F(sname.c_str(), sname.c_str(), 10, 0, 10);
    sname = "hErrHist_Tot_hyp" + shyp[ih] + "_Fail";
    tree.hMEPhaseSpace_ErrorTot_Fail[ih] = new TH1F(sname.c_str(), sname.c_str(), 10, 0, 10);
  }
  
  string index_CatJets[12];
  index_CatJets[0] = "3l_2b_2j";
  index_CatJets[1] = "3l_1b_2j";
  index_CatJets[2] = "3l_2b_1j";
  index_CatJets[3] = "3l_1b_1j";
  index_CatJets[4] = "3l_2b_0j";
  index_CatJets[5] = "4l_2b";
  index_CatJets[6] = "4l_1b";
  index_CatJets[7] = "2lss_2b_4j";
  index_CatJets[8] = "2lss_1b_4j";
  index_CatJets[9] = "2lss_2b_3j";
  index_CatJets[10] = "2lss_1b_3j";
  index_CatJets[11] = "2lss_2b_2j";

  Double_t index_XS[7];
  for(int ih=0; ih<nhyp; ih++){
     if(shyp[ih]=="TTLL") index_XS[ih]=xsTTLL;
     if(shyp[ih]=="TTW" || shyp[ih]=="TTWJJ") index_XS[ih]=xsTTW;
     if(shyp[ih]=="TTH" || shyp[ih]=="TTHsl" || shyp[ih]=="TTHfl") index_XS[ih]=xsTTH;
     if(shyp[ih]=="TTbar" || shyp[ih]=="TTbarsl" || shyp[ih]=="TTbarfl") index_XS[ih]=xsTTbar;
  }

  bool sel_tmp=0;
  
  for(int ih=0; ih<nhyp; ih++){
     for(int index_cat=0; index_cat<12; index_cat++){
   	  sname = "hMEAllWeights_" + shyp[ih] + "_" + index_CatJets[index_cat];
   	  tree.hMEAllWeights[ih][index_cat] = new TH1D(sname.c_str(),sname.c_str(),1000,0,1e-9);
   	  sname = "hMEAllWeights_nlog_" + shyp[ih] + "_" + index_CatJets[index_cat];
   	  tree.hMEAllWeights_nlog[ih][index_cat] = new TH1D(sname.c_str(),sname.c_str(),700,0,700);
     }
  }
  
    int nHypAllowed_TTH = 0, nHypAllowed_TTbar = 0;

    Permutations* MEMpermutations = new Permutations[nhyp];
                                                                                                                   
    int ievent=0;
    for (Long64_t iEvent=atoi(argv[2]); iEvent<evMax; iEvent++){

      cout << "iEvent="<<iEvent<<endl;

      tree.ReadMultilepton(iEvent, &multiLepton);

      tree.multilepton_Bjet1_P4 = *tree.multilepton_Bjet1_P4_ptr;
      tree.multilepton_Bjet2_P4 = *tree.multilepton_Bjet2_P4_ptr;
      tree.multilepton_Lepton1_P4 = *tree.multilepton_Lepton1_P4_ptr;
      tree.multilepton_Lepton2_P4 = *tree.multilepton_Lepton2_P4_ptr;
      tree.multilepton_Lepton3_P4 = *tree.multilepton_Lepton3_P4_ptr;
      tree.multilepton_Lepton4_P4 = *tree.multilepton_Lepton4_P4_ptr;
      tree.multilepton_JetHighestPt1_P4 = *tree.multilepton_JetHighestPt1_P4_ptr;
      tree.multilepton_JetHighestPt2_P4 = *tree.multilepton_JetHighestPt2_P4_ptr;
      tree.multilepton_JetClosestMw1_P4 = *tree.multilepton_JetClosestMw1_P4_ptr;
      tree.multilepton_JetClosestMw2_P4 = *tree.multilepton_JetClosestMw2_P4_ptr;
      tree.multilepton_JetLowestMjj1_P4 = *tree.multilepton_JetLowestMjj1_P4_ptr;
      tree.multilepton_JetLowestMjj2_P4 = *tree.multilepton_JetLowestMjj2_P4_ptr;
      tree.multilepton_JetHighestPt1_2ndPair_P4 = *tree.multilepton_JetHighestPt1_2ndPair_P4_ptr;
      tree.multilepton_JetHighestPt2_2ndPair_P4 = *tree.multilepton_JetHighestPt2_2ndPair_P4_ptr;
      tree.multilepton_JetClosestMw1_2ndPair_P4 = *tree.multilepton_JetClosestMw1_2ndPair_P4_ptr;
      tree.multilepton_JetClosestMw2_2ndPair_P4 = *tree.multilepton_JetClosestMw2_2ndPair_P4_ptr;
      tree.multilepton_JetLowestMjj1_2ndPair_P4 = *tree.multilepton_JetLowestMjj1_2ndPair_P4_ptr;
      tree.multilepton_JetLowestMjj2_2ndPair_P4 = *tree.multilepton_JetLowestMjj2_2ndPair_P4_ptr;

      tree.multilepton_mET = *tree.multilepton_mET_ptr;
      tree.multilepton_Ptot = *tree.multilepton_Ptot_ptr;

      tree.mc_mem_ttz_weight_logmean = 0;
      tree.mc_mem_ttz_weight_log = 0;
      tree.mc_mem_ttz_weight = 0;
      tree.mc_mem_ttz_weight_time = 0;
      tree.mc_mem_ttz_weight_err = 0;
      tree.mc_mem_ttz_weight_chi2 = 0;
      tree.mc_mem_ttz_weight_max = 0;
      tree.mc_mem_ttz_weight_avg = 0;
      tree.mc_kin_ttz_weight_logmax = 0;
      tree.mc_kin_ttz_weight_logmaxint = 0;
      tree.mc_mem_ttz_weight_kinmax = 0;
      tree.mc_mem_ttz_weight_kinmaxint = 0;
      tree.mc_mem_ttz_weight_JEC_up = 0;
      tree.mc_mem_ttz_weight_JEC_down = 0;
      tree.mc_mem_ttz_weight_JER_up = 0;
      tree.mc_mem_ttz_weight_JER_down = 0;

      tree.mc_mem_ttw_weight_logmean = 0;
      tree.mc_mem_ttw_weight_log = 0;
      tree.mc_mem_ttw_weight = 0;
      tree.mc_mem_ttw_weight_time = 0;
      tree.mc_mem_ttw_weight_err = 0;
      tree.mc_mem_ttw_weight_chi2 = 0;
      tree.mc_mem_ttw_weight_max = 0;
      tree.mc_mem_ttw_weight_avg = 0;
      tree.mc_kin_ttw_weight_logmax = 0;
      tree.mc_kin_ttw_weight_logmaxint = 0;
      tree.mc_mem_ttw_weight_kinmax = 0;
      tree.mc_mem_ttw_weight_kinmaxint = 0;
      tree.mc_mem_ttw_weight_JEC_up = 0;
      tree.mc_mem_ttw_weight_JEC_down = 0;
      tree.mc_mem_ttw_weight_JER_up = 0;
      tree.mc_mem_ttw_weight_JER_down = 0;

      tree.mc_mem_ttwjj_weight_logmean = 0;
      tree.mc_mem_ttwjj_weight_log = 0;
      tree.mc_mem_ttwjj_weight = 0;
      tree.mc_mem_ttwjj_weight_time = 0;
      tree.mc_mem_ttwjj_weight_err = 0;
      tree.mc_mem_ttwjj_weight_chi2 = 0;
      tree.mc_mem_ttwjj_weight_max = 0;
      tree.mc_mem_ttwjj_weight_avg = 0;
      tree.mc_kin_ttwjj_weight_logmax = 0;
      tree.mc_kin_ttwjj_weight_logmaxint = 0;
      tree.mc_mem_ttwjj_weight_kinmax = 0;
      tree.mc_mem_ttwjj_weight_kinmaxint = 0;
      tree.mc_mem_ttwjj_weight_JEC_up = 0;
      tree.mc_mem_ttwjj_weight_JEC_down = 0;
      tree.mc_mem_ttwjj_weight_JER_up = 0;
      tree.mc_mem_ttwjj_weight_JER_down = 0;

      tree.mc_mem_tthfl_weight_logmean = 0;
      tree.mc_mem_tthfl_weight_log = 0;
      tree.mc_mem_tthfl_weight = 0;
      tree.mc_mem_tthfl_weight_time = 0;
      tree.mc_mem_tthfl_weight_err = 0;
      tree.mc_mem_tthfl_weight_chi2 = 0;
      tree.mc_mem_tthfl_weight_max = 0;
      tree.mc_mem_tthfl_weight_avg = 0;   
      tree.mc_kin_tthfl_weight_logmax = 0;
      tree.mc_kin_tthfl_weight_logmaxint = 0;
      tree.mc_mem_tthfl_weight_kinmax = 0;
      tree.mc_mem_tthfl_weight_kinmaxint = 0;
      tree.mc_mem_tthfl_weight_JEC_up = 0;
      tree.mc_mem_tthfl_weight_JEC_down = 0;
      tree.mc_mem_tthfl_weight_JER_up = 0;
      tree.mc_mem_tthfl_weight_JER_down = 0;

      tree.mc_mem_tthsl_weight_logmean = 0;
      tree.mc_mem_tthsl_weight_log = 0;
      tree.mc_mem_tthsl_weight = 0;
      tree.mc_mem_tthsl_weight_time = 0;
      tree.mc_mem_tthsl_weight_err = 0;
      tree.mc_mem_tthsl_weight_chi2 = 0;
      tree.mc_mem_tthsl_weight_max = 0;
      tree.mc_mem_tthsl_weight_avg = 0;
      tree.mc_kin_tthsl_weight_logmax = 0;
      tree.mc_kin_tthsl_weight_logmaxint = 0;
      tree.mc_mem_tthsl_weight_kinmax = 0;
      tree.mc_mem_tthsl_weight_kinmaxint = 0;
      tree.mc_mem_tthsl_weight_JEC_up = 0;
      tree.mc_mem_tthsl_weight_JEC_down = 0;
      tree.mc_mem_tthsl_weight_JER_up = 0;
      tree.mc_mem_tthsl_weight_JER_down = 0;

      tree.mc_mem_tth_weight_logmean = 0;
      tree.mc_mem_tth_weight_log = 0;
      tree.mc_mem_tth_weight = 0;
      tree.mc_mem_tth_weight_time = 0;
      tree.mc_mem_tth_weight_err = 0;
      tree.mc_mem_tth_weight_chi2 = 0;
      tree.mc_mem_tth_weight_max = 0;
      tree.mc_mem_tth_weight_avg = 0;
      tree.mc_kin_tth_weight_logmax = 0;
      tree.mc_kin_tth_weight_logmaxint = 0;
      tree.mc_mem_tth_weight_kinmax = 0;
      tree.mc_mem_tth_weight_kinmaxint = 0;
      tree.mc_mem_tth_weight_JEC_up = 0;
      tree.mc_mem_tth_weight_JEC_down = 0;
      tree.mc_mem_tth_weight_JER_up = 0;
      tree.mc_mem_tth_weight_JER_down = 0;

      tree.mc_mem_ttbarfl_weight_logmean = 0;
      tree.mc_mem_ttbarfl_weight_log = 0;
      tree.mc_mem_ttbarfl_weight = 0;
      tree.mc_mem_ttbarfl_weight_time = 0;
      tree.mc_mem_ttbarfl_weight_err = 0;
      tree.mc_mem_ttbarfl_weight_chi2 = 0;
      tree.mc_mem_ttbarfl_weight_max = 0;
      tree.mc_mem_ttbarfl_weight_avg = 0;
      tree.mc_kin_ttbarfl_weight_logmax = 0;
      tree.mc_kin_ttbarfl_weight_logmaxint = 0;
      tree.mc_mem_ttbarfl_weight_kinmax = 0;
      tree.mc_mem_ttbarfl_weight_kinmaxint = 0;
      tree.mc_mem_ttbarfl_weight_JEC_up = 0;
      tree.mc_mem_ttbarfl_weight_JEC_down = 0;
      tree.mc_mem_ttbarfl_weight_JER_up = 0;
      tree.mc_mem_ttbarfl_weight_JER_down = 0;

      tree.mc_mem_ttbarsl_weight_logmean = 0;
      tree.mc_mem_ttbarsl_weight_log = 0;
      tree.mc_mem_ttbarsl_weight = 0;
      tree.mc_mem_ttbarsl_weight_time = 0;
      tree.mc_mem_ttbarsl_weight_err = 0;
      tree.mc_mem_ttbarsl_weight_chi2 = 0;
      tree.mc_mem_ttbarsl_weight_max = 0;
      tree.mc_mem_ttbarsl_weight_avg = 0;
      tree.mc_kin_ttbarsl_weight_logmax = 0;
      tree.mc_kin_ttbarsl_weight_logmaxint = 0;
      tree.mc_mem_ttbarsl_weight_kinmax = 0;
      tree.mc_mem_ttbarsl_weight_kinmaxint = 0;
      tree.mc_mem_ttbarsl_weight_JEC_up = 0;
      tree.mc_mem_ttbarsl_weight_JEC_down = 0;
      tree.mc_mem_ttbarsl_weight_JER_up = 0;
      tree.mc_mem_ttbarsl_weight_JER_down = 0;

      tree.mc_mem_ttbar_weight_logmean = 0;
      tree.mc_mem_ttbar_weight_log = 0;
      tree.mc_mem_ttbar_weight = 0;
      tree.mc_mem_ttbar_weight_time = 0;
      tree.mc_mem_ttbar_weight_err = 0;
      tree.mc_mem_ttbar_weight_chi2 = 0;
      tree.mc_mem_ttbar_weight_max = 0;
      tree.mc_mem_ttbar_weight_avg = 0;
      tree.mc_kin_ttbar_weight_logmax = 0;
      tree.mc_kin_ttbar_weight_logmaxint = 0;
      tree.mc_mem_ttbar_weight_kinmax = 0;
      tree.mc_mem_ttbar_weight_kinmaxint = 0;
      tree.mc_mem_ttbar_weight_JEC_up = 0;
      tree.mc_mem_ttbar_weight_JEC_down = 0;
      tree.mc_mem_ttbar_weight_JER_up = 0;
      tree.mc_mem_ttbar_weight_JER_down = 0;

      tree.MEAllWeights_TTLL->clear();
      tree.MEAllWeights_TTHfl->clear();
      tree.MEAllWeights_TTHsl->clear();
      tree.MEAllWeights_TTH->clear();
      tree.MEAllWeights_TTW->clear();
      tree.MEAllWeights_TTWJJ->clear();
      tree.MEAllWeights_TTbarfl->clear();
      tree.MEAllWeights_TTbarsl->clear();
      tree.MEAllWeights_TTbar->clear();

      tree.MEAllWeights_TTLL_log->clear();
      tree.MEAllWeights_TTHfl_log->clear();
      tree.MEAllWeights_TTHsl_log->clear();
      tree.MEAllWeights_TTH_log->clear();
      tree.MEAllWeights_TTW_log->clear();
      tree.MEAllWeights_TTWJJ_log->clear();
      tree.MEAllWeights_TTbarfl_log->clear();
      tree.MEAllWeights_TTbarsl_log->clear();
      tree.MEAllWeights_TTbar_log->clear();

    nHypAllowed_TTH = 0;
    nHypAllowed_TTbar = 0;

   int initresult = 0;
   for (int ih=0; ih<nhyp; ih++){

     MEMpermutations[ih].SetMultiLepton(&multiLepton, &hypIntegrator);

     initresult = MEMpermutations[ih].InitializeHyp(&hypIntegrator, hyp[ih], nPointsHyp[ih], shyp[ih], doMinimization, JetChoice, nPermutationJetSyst);

     if (initresult==1) MEMpermutations[ih].LoopPermutations(&hypIntegrator);


               for(int index_cat=0; index_cat<12; index_cat++){
                  if(index_cat<=6)
                     sel_tmp = (tree.is_3l_TTH_SR)*(tree.catJets==index_cat);
                  if(index_cat>=7 && index_cat<=11)
                     sel_tmp = (tree.is_2lss_TTH_SR)*(tree.catJets==index_cat);
                  if(shyp[ih]=="TTLL" && index_cat<=6) 
                     sel_tmp = (tree.is_3l_TTH_SR)*(tree.catJets==index_cat)*(tree.mc_ttZhypAllowed==1);
                  if(shyp[ih]=="TTLL" && index_cat>=7 && index_cat<=11) 
                     sel_tmp = (tree.is_2lss_TTH_SR)*(tree.catJets==index_cat)*(tree.mc_ttZhypAllowed==1);
                  if(sel_tmp){
		    for (unsigned int ip=0; ip<MEMpermutations[ih].resMEM_all.size() ; ip++){
                      if(MEMpermutations[ih].resMEM_all.at(ip).weight>0){
                        tree.hMEAllWeights[ih][index_cat]->Fill(MEMpermutations[ih].resMEM_all.at(ip).weight,tree.weight);
                        tree.hMEAllWeights_nlog[ih][index_cat]->Fill(-log(MEMpermutations[ih].resMEM_all.at(ip).weight),tree.weight);
                      }
                      else{
                        tree.hMEAllWeights[ih][index_cat]->Fill(1e-300,tree.weight);
                        tree.hMEAllWeights_nlog[ih][index_cat]->Fill(-log(1e-300),tree.weight);
                      }
		    }
                  }
               }


     TLorentzVector Pnull(0,0,0,0);
     if (shyp[ih]=="TTLL"){
          tree.mc_kin_ttz_tophad_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Top_P4 : Pnull;
          tree.mc_kin_ttz_tophad_W_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.W_P4 : Pnull;
          tree.mc_kin_ttz_tophad_Bjet_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Bjet_P4 : Pnull;
          tree.mc_kin_ttz_tophad_Jet1_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Jet1_P4 : Pnull;
          tree.mc_kin_ttz_tophad_Jet2_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Jet2_P4 : Pnull;
	  tree.mc_kin_ttz_tophad_Pt = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Top_Pt : -99;
          tree.mc_kin_ttz_tophad_Wmass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.W_Mass : -99;
          tree.mc_kin_ttz_tophad_Benergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.B_E : -99;
          tree.mc_kin_ttz_tophad_Jet1energy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Jet1_E : -99;
          tree.mc_kin_ttz_tophad_Jet2energy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Jet2_E : -99;
          tree.mc_kin_ttz_toplep_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Top_P4 : Pnull;
          tree.mc_kin_ttz_toplep_W_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.W_P4 : Pnull;
          tree.mc_kin_ttz_toplep_Bjet_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Bjet_P4 : Pnull;
          tree.mc_kin_ttz_toplep_Lep_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Lep_P4 : Pnull;
          tree.mc_kin_ttz_toplep_Neut_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Neut_P4 : Pnull;
	  tree.mc_kin_ttz_toplep_Pt = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Top_Pt : -99;
          tree.mc_kin_ttz_toplep_Wmass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.W_Mass : -99;
          tree.mc_kin_ttz_toplep_Benergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.B_E : -99;
          tree.mc_kin_ttz_toplep_Neutenergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Neut_E : -99;
          tree.mc_kin_ttz_toplep2_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Top_P4 : Pnull;
          tree.mc_kin_ttz_toplep2_W_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.W_P4 : Pnull;
          tree.mc_kin_ttz_toplep2_Bjet_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Bjet_P4 : Pnull;
          tree.mc_kin_ttz_toplep2_Lep_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Lep_P4 : Pnull;
          tree.mc_kin_ttz_toplep2_Neut_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Neut_P4 : Pnull;
          tree.mc_kin_ttz_toplep2_Pt = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Top_Pt : -99;
          tree.mc_kin_ttz_toplep2_Wmass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.W_Mass : -99;
          tree.mc_kin_ttz_toplep2_Benergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.B_E : -99;
          tree.mc_kin_ttz_toplep2_Neutenergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Neut_E : -99;
          tree.mc_kin_ttz_zll_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Zll.Z_P4 : Pnull;
          tree.mc_kin_ttz_zll_Lep1_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Zll.Lep1_P4 : Pnull;
          tree.mc_kin_ttz_zll_Lep2_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Zll.Lep2_P4 : Pnull;
	  tree.mc_kin_ttz_zll_Pt = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Zll.Z_Pt : -99;
          tree.mc_kin_ttz_zll_Zmass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Zll.Z_Mass : -99;
     }
     if (shyp[ih]=="TTHfl") {
	  tree.mc_kin_tthfl_tophad_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Top_P4 : Pnull;
          tree.mc_kin_tthfl_tophad_W_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.W_P4 : Pnull;
          tree.mc_kin_tthfl_tophad_Bjet_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Bjet_P4 : Pnull;
          tree.mc_kin_tthfl_tophad_Jet1_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Jet1_P4 : Pnull;
          tree.mc_kin_tthfl_tophad_Jet2_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Jet2_P4 : Pnull;
          tree.mc_kin_tthfl_tophad_Pt = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Top_Pt : -99;
          tree.mc_kin_tthfl_tophad_Wmass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.W_Mass : -99;
          tree.mc_kin_tthfl_tophad_Benergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.B_E : -99;
          tree.mc_kin_tthfl_tophad_Jet1energy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Jet1_E : -99;
          tree.mc_kin_tthfl_tophad_Jet2energy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Jet2_E : -99;
	  tree.mc_kin_tthfl_toplep_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Top_P4 : Pnull;
          tree.mc_kin_tthfl_toplep_W_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.W_P4 : Pnull;
          tree.mc_kin_tthfl_toplep_Bjet_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Bjet_P4 : Pnull;
          tree.mc_kin_tthfl_toplep_Lep_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Lep_P4 : Pnull;
          tree.mc_kin_tthfl_toplep_Neut_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Neut_P4 : Pnull;
          tree.mc_kin_tthfl_toplep_Pt = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Top_Pt : -99;
          tree.mc_kin_tthfl_toplep_Wmass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.W_Mass : -99;
          tree.mc_kin_tthfl_toplep_Benergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.B_E : -99;
          tree.mc_kin_tthfl_toplep_Neutenergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==4) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Neut_E : -99;
          tree.mc_kin_tthfl_toplep2_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Top_P4 : Pnull;
          tree.mc_kin_tthfl_toplep2_W_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.W_P4 : Pnull;
          tree.mc_kin_tthfl_toplep2_Bjet_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Bjet_P4 : Pnull;
          tree.mc_kin_tthfl_toplep2_Lep_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Lep_P4 : Pnull;
          tree.mc_kin_tthfl_toplep2_Neut_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Neut_P4 : Pnull;
          tree.mc_kin_tthfl_toplep2_Pt = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Top_Pt : -99;
          tree.mc_kin_tthfl_toplep2_Wmass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.W_Mass : -99;
          tree.mc_kin_tthfl_toplep2_Benergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.B_E : -99;
          tree.mc_kin_tthfl_toplep2_Neutenergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Neut_E : -99;
          tree.mc_kin_tthfl_h2l2nu_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_H2l2nu.Higgs_P4 : Pnull;
          tree.mc_kin_tthfl_h2l2nu_W1_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_H2l2nu.W1_P4 : Pnull;
          tree.mc_kin_tthfl_h2l2nu_W2_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_H2l2nu.W2_P4 : Pnull;
          tree.mc_kin_tthfl_h2l2nu_Lep1_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_H2l2nu.Lep1_P4 : Pnull;
          tree.mc_kin_tthfl_h2l2nu_Lep2_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_H2l2nu.Lep2_P4 : Pnull;
          tree.mc_kin_tthfl_h2l2nu_Neut1_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_H2l2nu.Neut1_P4 : Pnull;
          tree.mc_kin_tthfl_h2l2nu_Neut2_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_H2l2nu.Neut2_P4 : Pnull;
          tree.mc_kin_tthfl_h2l2nu_Pt = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_H2l2nu.Higgs_Pt : -99;
          tree.mc_kin_tthfl_h2l2nu_W1mass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_H2l2nu.W1_Mass : -99;
          tree.mc_kin_tthfl_h2l2nu_Neut1energy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_H2l2nu.Neut1_E : -99;
          tree.mc_kin_tthfl_h2l2nu_W2mass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_H2l2nu.W2_Mass : -99;
          tree.mc_kin_tthfl_h2l2nu_Neut2energy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_H2l2nu.Neut2_E : -99;
       //}
     }
     if (shyp[ih]=="TTHsl"){
          tree.mc_kin_tthsl_tophad_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Top_P4 : Pnull;
          tree.mc_kin_tthsl_tophad_W_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.W_P4 : Pnull;
          tree.mc_kin_tthsl_tophad_Bjet_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Bjet_P4 : Pnull;
          tree.mc_kin_tthsl_tophad_Jet1_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Jet1_P4 : Pnull;
          tree.mc_kin_tthsl_tophad_Jet2_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Jet2_P4 : Pnull;
          tree.mc_kin_tthsl_tophad_Pt = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Top_Pt : -99;
          tree.mc_kin_tthsl_tophad_Wmass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.W_Mass : -99;
          tree.mc_kin_tthsl_tophad_Benergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.B_E : -99;
          tree.mc_kin_tthsl_tophad_Jet1energy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Jet1_E : -99;
          tree.mc_kin_tthsl_tophad_Jet2energy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Jet2_E : -99;
          tree.mc_kin_tthsl_toplep_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Top_P4 : Pnull;
          tree.mc_kin_tthsl_toplep_W_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.W_P4 : Pnull;
          tree.mc_kin_tthsl_toplep_Bjet_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Bjet_P4 : Pnull;
          tree.mc_kin_tthsl_toplep_Lep_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Lep_P4 : Pnull;
          tree.mc_kin_tthsl_toplep_Neut_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Neut_P4 : Pnull;
          tree.mc_kin_tthsl_toplep_Pt = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Top_Pt : -99;
          tree.mc_kin_tthsl_toplep_Wmass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.W_Mass : -99;
          tree.mc_kin_tthsl_toplep_Benergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.B_E : -99;
          tree.mc_kin_tthsl_toplep_Neutenergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Neut_E : -99;
          tree.mc_kin_tthsl_toplep2_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Top_P4 : Pnull;
          tree.mc_kin_tthsl_toplep2_W_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.W_P4 : Pnull;
          tree.mc_kin_tthsl_toplep2_Bjet_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Bjet_P4 : Pnull;
          tree.mc_kin_tthsl_toplep2_Lep_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Lep_P4 : Pnull;
          tree.mc_kin_tthsl_toplep2_Neut_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Neut_P4 : Pnull;
          tree.mc_kin_tthsl_toplep2_Pt = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Top_Pt : -99;
          tree.mc_kin_tthsl_toplep2_Wmass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.W_Mass : -99;
          tree.mc_kin_tthsl_toplep2_Benergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.B_E : -99;
          tree.mc_kin_tthsl_toplep2_Neutenergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Neut_E : -99;
          tree.mc_kin_tthsl_hlnujj_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Hlnujj.Higgs_P4 : Pnull;
          tree.mc_kin_tthsl_hlnujj_W1_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Hlnujj.W1_P4 : Pnull;
          tree.mc_kin_tthsl_hlnujj_W2_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Hlnujj.W2_P4 : Pnull;
          tree.mc_kin_tthsl_hlnujj_Lep_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Hlnujj.Lep_P4 : Pnull;
          tree.mc_kin_tthsl_hlnujj_Neut_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Hlnujj.Neut_P4 : Pnull;
          tree.mc_kin_tthsl_hlnujj_Jet1_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Hlnujj.Jet1_P4 : Pnull;
          tree.mc_kin_tthsl_hlnujj_Jet2_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Hlnujj.Jet2_P4 : Pnull;
          tree.mc_kin_tthsl_hlnujj_Pt = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Hlnujj.Higgs_Pt : -99;
          tree.mc_kin_tthsl_hlnujj_W1mass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Hlnujj.Wlnu_Mass : -99;
          tree.mc_kin_tthsl_hlnujj_Neut1energy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Hlnujj.Neut_E : -99;
          tree.mc_kin_tthsl_hlnujj_W2mass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Hlnujj.Wjj_Mass : -99;
          tree.mc_kin_tthsl_hlnujj_Jet1energy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Hlnujj.Jet1_E : -99;
          tree.mc_kin_tthsl_hlnujj_Jet2energy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Hlnujj.Jet2_E : -99;
     }
     if (shyp[ih]=="TTW"){
          tree.mc_kin_ttw_tophad_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Top_P4 : Pnull;
          tree.mc_kin_ttw_tophad_W_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.W_P4 : Pnull;
          tree.mc_kin_ttw_tophad_Bjet_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Bjet_P4 : Pnull;
          tree.mc_kin_ttw_tophad_Jet1_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Jet1_P4 : Pnull;
          tree.mc_kin_ttw_tophad_Jet2_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0  && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Jet2_P4 : Pnull;
          tree.mc_kin_ttw_tophad_Pt = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Top_Pt : -99;
          tree.mc_kin_ttw_tophad_Wmass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.W_Mass : -99;
          tree.mc_kin_ttw_tophad_Benergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.B_E : -99;
          tree.mc_kin_ttw_tophad_Jet1energy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Jet1_E : -99;
          tree.mc_kin_ttw_tophad_Jet2energy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==2) ? hypIntegrator.meIntegrator->MEMKin_TopHad.Jet2_E : -99;
          tree.mc_kin_ttw_toplep_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Top_P4 : Pnull;
          tree.mc_kin_ttw_toplep_W_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.W_P4 : Pnull;
          tree.mc_kin_ttw_toplep_Bjet_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Bjet_P4 : Pnull;
          tree.mc_kin_ttw_toplep_Lep_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Lep_P4 : Pnull;
          tree.mc_kin_ttw_toplep_Neut_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Neut_P4 : Pnull;
          tree.mc_kin_ttw_toplep_Pt = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Top_Pt : -99;
          tree.mc_kin_ttw_toplep_Wmass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.W_Mass : -99;
          tree.mc_kin_ttw_toplep_Benergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.B_E : -99;
          tree.mc_kin_ttw_toplep_Neutenergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0 && hypIntegrator.meIntegrator->iNleptons==3) ? hypIntegrator.meIntegrator->MEMKin_TopLep1.Neut_E : -99;
          tree.mc_kin_ttw_toplep2_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Top_P4 : Pnull;
          tree.mc_kin_ttw_toplep2_W_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.W_P4 : Pnull;
          tree.mc_kin_ttw_toplep2_Bjet_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Bjet_P4 : Pnull;
          tree.mc_kin_ttw_toplep2_Lep_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Lep_P4 : Pnull;
          tree.mc_kin_ttw_toplep2_Neut_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Neut_P4 : Pnull;
          tree.mc_kin_ttw_toplep2_Pt = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Top_Pt : -99;
          tree.mc_kin_ttw_toplep2_Wmass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.W_Mass : -99;
          tree.mc_kin_ttw_toplep2_Benergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.B_E : -99;
          tree.mc_kin_ttw_toplep2_Neutenergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_TopLep2.Neut_E : -99;
          tree.mc_kin_ttw_wlnu_W_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Wlnu.W_P4 : Pnull;
          tree.mc_kin_ttw_wlnu_Lep_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Wlnu.Lep_P4 : Pnull;
          tree.mc_kin_ttw_wlnu_Neut_P4 = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Wlnu.Neut_P4 : Pnull;
          tree.mc_kin_ttw_wlnu_Pt = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Wlnu.W_Pt : -99;
          tree.mc_kin_ttw_wlnu_Wmass = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Wlnu.W_Mass : -99;
          tree.mc_kin_ttw_wlnu_Neutenergy = (MEMpermutations[ih].resKin_maxKinFit_Int.weight>0) ? hypIntegrator.meIntegrator->MEMKin_Wlnu.Neut_E : -99;
     }

   } //end hypothesis for loop
/*
   for (int ih=0; ih<nhyp; ih++) {
     cout<< "hyp "<< shyp[ih]<<" nPermAllowed="<<nHypAllowed[index[ih]]<<" (permutation or kinematics), nPermForbidden="<<nNullResult[index[ih]]<<endl;
     if (nHypAllowed[ih]+nNullResult[ih]>0)
	tree.hMEPhaseSpace_ErrorTot[ih]->Scale(1./(double)(nHypAllowed[ih]+nNullResult[ih]));
     if (nHypAllowed[ih]>0)
	tree.hMEPhaseSpace_ErrorTot_Pass[ih]->Scale(1./(double)(nHypAllowed[ih]));
     if (nNullResult[ih]>0)
        tree.hMEPhaseSpace_ErrorTot_Fail[ih]->Scale(1./(double)(nNullResult[ih]));
   }
*/
   ievent++;
 

   if (nhyp>=0){

     if (index_hyp[0]!=-1){
       if (MEMpermutations[index_hyp[0]].computeHyp){
       tree.mc_mem_ttz_weight = MEMpermutations[index_hyp[0]].resMEM_avgExl0.weight;
       tree.mc_mem_ttz_weight_err = MEMpermutations[index_hyp[0]].resMEM_avgExl0.err;
       tree.mc_mem_ttz_weight_chi2 = MEMpermutations[index_hyp[0]].resMEM_avgExl0.chi2;
       tree.mc_mem_ttz_weight_time = MEMpermutations[index_hyp[0]].resMEM_avgExl0.time;
       tree.mc_mem_ttz_weight_log = log(MEMpermutations[index_hyp[0]].resMEM_avgExl0.weight);
       tree.mc_mem_ttz_weight_logmean = MEMpermutations[index_hyp[0]].resMEM_logavgExl0.weight;
       tree.mc_mem_ttz_weight_avg = MEMpermutations[index_hyp[0]].resMEM_avgIncl0.weight;
       tree.mc_mem_ttz_weight_max = MEMpermutations[index_hyp[0]].resMEM_max.weight;
       tree.mc_mem_ttz_weight_JEC_up = MEMpermutations[index_hyp[0]].resMEM_avgExl0_JEC_up.weight;
       tree.mc_mem_ttz_weight_JEC_down = MEMpermutations[index_hyp[0]].resMEM_avgExl0_JEC_down.weight;
       tree.mc_mem_ttz_weight_JER_up = MEMpermutations[index_hyp[0]].resMEM_avgExl0_JER_up.weight;
       tree.mc_mem_ttz_weight_JER_down = MEMpermutations[index_hyp[0]].resMEM_avgExl0_JER_down.weight;
       tree.mc_kin_ttz_weight_logmax = log(MEMpermutations[index_hyp[0]].resKin_maxKinFit.weight);
       tree.mc_kin_ttz_weight_logmaxint = log(MEMpermutations[index_hyp[0]].resKin_maxKinFit_Int.weight);
       tree.mc_mem_ttz_weight_kinmax = MEMpermutations[index_hyp[0]].resMEM_maxKinFit.weight;
       tree.mc_mem_ttz_weight_kinmaxint = MEMpermutations[index_hyp[0]].resMEM_maxKinFit_Int.weight;
       FillWeightVectors(MEMpermutations[index_hyp[0]], tree.MEAllWeights_TTLL, tree.MEAllWeights_TTLL_log);
      }
    }
    if (index_hyp[1]!=-1){
      if (MEMpermutations[index_hyp[1]].computeHyp){
        nHypAllowed_TTH += MEMpermutations[index_hyp[1]].nHypAllowed;
        tree.mc_mem_tthsl_weight = MEMpermutations[index_hyp[1]].resMEM_avgExl0.weight;
        tree.mc_mem_tthsl_weight_err = MEMpermutations[index_hyp[1]].resMEM_avgExl0.err;
        tree.mc_mem_tthsl_weight_chi2 = MEMpermutations[index_hyp[1]].resMEM_avgExl0.chi2;
        tree.mc_mem_tthsl_weight_time = MEMpermutations[index_hyp[1]].resMEM_avgExl0.time;
        tree.mc_mem_tthsl_weight_log = log(MEMpermutations[index_hyp[1]].resMEM_avgExl0.weight);
        tree.mc_mem_tthsl_weight_logmean = MEMpermutations[index_hyp[1]].resMEM_logavgExl0.weight;
        tree.mc_mem_tthsl_weight_avg = MEMpermutations[index_hyp[1]].resMEM_avgIncl0.weight;
        tree.mc_mem_tthsl_weight_max = MEMpermutations[index_hyp[1]].resMEM_max.weight;
        tree.mc_mem_tthsl_weight_JEC_up = MEMpermutations[index_hyp[1]].resMEM_avgExl0_JEC_up.weight;
        tree.mc_mem_tthsl_weight_JEC_down = MEMpermutations[index_hyp[1]].resMEM_avgExl0_JEC_down.weight;
        tree.mc_mem_tthsl_weight_JER_up = MEMpermutations[index_hyp[1]].resMEM_avgExl0_JER_up.weight;
        tree.mc_mem_tthsl_weight_JER_down = MEMpermutations[index_hyp[1]].resMEM_avgExl0_JER_down.weight;
        tree.mc_kin_tthsl_weight_logmax = log(MEMpermutations[index_hyp[1]].resKin_maxKinFit.weight);
        tree.mc_kin_tthsl_weight_logmaxint = log(MEMpermutations[index_hyp[1]].resKin_maxKinFit_Int.weight);
        tree.mc_mem_tthsl_weight_kinmax = MEMpermutations[index_hyp[1]].resMEM_maxKinFit.weight;
        tree.mc_mem_tthsl_weight_kinmaxint = MEMpermutations[index_hyp[1]].resMEM_maxKinFit_Int.weight;
        FillWeightVectors(MEMpermutations[index_hyp[1]], tree.MEAllWeights_TTHfl, tree.MEAllWeights_TTHfl_log);
      }
   }
   if (index_hyp[2]!=-1){
     if (MEMpermutations[index_hyp[2]].computeHyp){
       nHypAllowed_TTH += MEMpermutations[index_hyp[2]].nHypAllowed;
       tree.mc_mem_tthfl_weight = MEMpermutations[index_hyp[2]].resMEM_avgExl0.weight;
       tree.mc_mem_tthfl_weight_err = MEMpermutations[index_hyp[2]].resMEM_avgExl0.err;
       tree.mc_mem_tthfl_weight_chi2 = MEMpermutations[index_hyp[2]].resMEM_avgExl0.chi2;
       tree.mc_mem_tthfl_weight_time = MEMpermutations[index_hyp[2]].resMEM_avgExl0.time;
       tree.mc_mem_tthfl_weight_log = log(MEMpermutations[index_hyp[2]].resMEM_avgExl0.weight);
       tree.mc_mem_tthfl_weight_logmean = MEMpermutations[index_hyp[2]].resMEM_logavgExl0.weight;
       tree.mc_mem_tthfl_weight_avg = MEMpermutations[index_hyp[2]].resMEM_avgIncl0.weight;
       tree.mc_mem_tthfl_weight_max = MEMpermutations[index_hyp[2]].resMEM_max.weight;
       tree.mc_mem_tthfl_weight_JEC_up = MEMpermutations[index_hyp[2]].resMEM_avgExl0_JEC_up.weight;
       tree.mc_mem_tthfl_weight_JEC_down = MEMpermutations[index_hyp[2]].resMEM_avgExl0_JEC_down.weight;
       tree.mc_mem_tthfl_weight_JER_up = MEMpermutations[index_hyp[2]].resMEM_avgExl0_JER_up.weight;
       tree.mc_mem_tthfl_weight_JER_down = MEMpermutations[index_hyp[2]].resMEM_avgExl0_JER_down.weight;
       tree.mc_kin_tthfl_weight_logmax = log(MEMpermutations[index_hyp[2]].resKin_maxKinFit.weight);
       tree.mc_kin_tthfl_weight_logmaxint = log(MEMpermutations[index_hyp[2]].resKin_maxKinFit_Int.weight);
       tree.mc_mem_tthfl_weight_kinmax = MEMpermutations[index_hyp[2]].resMEM_maxKinFit.weight;
       tree.mc_mem_tthfl_weight_kinmaxint = MEMpermutations[index_hyp[2]].resMEM_maxKinFit_Int.weight;
       FillWeightVectors(MEMpermutations[index_hyp[2]], tree.MEAllWeights_TTHfl, tree.MEAllWeights_TTHfl_log);
     }
   }
   if (index_hyp[3]!=-1){ 
     if (MEMpermutations[index_hyp[3]].computeHyp){
       tree.mc_mem_ttw_weight = MEMpermutations[index_hyp[3]].resMEM_avgExl0.weight;
       tree.mc_mem_ttw_weight_err = MEMpermutations[index_hyp[3]].resMEM_avgExl0.err;
       tree.mc_mem_ttw_weight_chi2 = MEMpermutations[index_hyp[3]].resMEM_avgExl0.chi2;
       tree.mc_mem_ttw_weight_time = MEMpermutations[index_hyp[3]].resMEM_avgExl0.time;
       tree.mc_mem_ttw_weight_log = log(MEMpermutations[index_hyp[3]].resMEM_avgExl0.weight);
       tree.mc_mem_ttw_weight_logmean = MEMpermutations[index_hyp[3]].resMEM_logavgExl0.weight;
       tree.mc_mem_ttw_weight_avg = MEMpermutations[index_hyp[3]].resMEM_avgIncl0.weight;
       tree.mc_mem_ttw_weight_max = MEMpermutations[index_hyp[3]].resMEM_max.weight;
       tree.mc_mem_ttw_weight_JEC_up = MEMpermutations[index_hyp[3]].resMEM_avgExl0_JEC_up.weight;
       tree.mc_mem_ttw_weight_JEC_down = MEMpermutations[index_hyp[3]].resMEM_avgExl0_JEC_down.weight;
       tree.mc_mem_ttw_weight_JER_up = MEMpermutations[index_hyp[3]].resMEM_avgExl0_JER_up.weight;
       tree.mc_mem_ttw_weight_JER_down = MEMpermutations[index_hyp[3]].resMEM_avgExl0_JER_down.weight;
       tree.mc_kin_ttw_weight_logmax = log(MEMpermutations[index_hyp[3]].resKin_maxKinFit.weight);
       tree.mc_kin_ttw_weight_logmaxint = log(MEMpermutations[index_hyp[3]].resKin_maxKinFit_Int.weight);
       tree.mc_mem_ttw_weight_kinmax = MEMpermutations[index_hyp[3]].resMEM_maxKinFit.weight;
       tree.mc_mem_ttw_weight_kinmaxint = MEMpermutations[index_hyp[3]].resMEM_maxKinFit_Int.weight;
       FillWeightVectors(MEMpermutations[index_hyp[3]], tree.MEAllWeights_TTW, tree.MEAllWeights_TTW_log);
     }
   }
   if (index_hyp[4]!=-1){
     if (MEMpermutations[index_hyp[4]].computeHyp){
       tree.mc_mem_ttwjj_weight = MEMpermutations[index_hyp[4]].resMEM_avgExl0.weight;
       tree.mc_mem_ttwjj_weight_err = MEMpermutations[index_hyp[4]].resMEM_avgExl0.err;
       tree.mc_mem_ttwjj_weight_chi2 = MEMpermutations[index_hyp[4]].resMEM_avgExl0.chi2;
       tree.mc_mem_ttwjj_weight_time = MEMpermutations[index_hyp[4]].resMEM_avgExl0.time;
       tree.mc_mem_ttwjj_weight_log = log(MEMpermutations[index_hyp[4]].resMEM_avgExl0.weight);
       tree.mc_mem_ttwjj_weight_logmean = MEMpermutations[index_hyp[4]].resMEM_logavgExl0.weight;
       tree.mc_mem_ttwjj_weight_avg = MEMpermutations[index_hyp[4]].resMEM_avgIncl0.weight;
       tree.mc_mem_ttwjj_weight_max = MEMpermutations[index_hyp[4]].resMEM_max.weight;
       tree.mc_mem_ttwjj_weight_JEC_up = MEMpermutations[index_hyp[4]].resMEM_avgExl0_JEC_up.weight;
       tree.mc_mem_ttwjj_weight_JEC_down = MEMpermutations[index_hyp[4]].resMEM_avgExl0_JEC_down.weight;
       tree.mc_mem_ttwjj_weight_JER_up = MEMpermutations[index_hyp[4]].resMEM_avgExl0_JER_up.weight;
       tree.mc_mem_ttwjj_weight_JER_down = MEMpermutations[index_hyp[4]].resMEM_avgExl0_JER_down.weight;
       tree.mc_kin_ttwjj_weight_logmax = log(MEMpermutations[index_hyp[4]].resKin_maxKinFit.weight);
       tree.mc_kin_ttwjj_weight_logmaxint = log(MEMpermutations[index_hyp[4]].resKin_maxKinFit_Int.weight);
       tree.mc_mem_ttwjj_weight_kinmax = MEMpermutations[index_hyp[4]].resMEM_maxKinFit.weight;
       tree.mc_mem_ttwjj_weight_kinmaxint = MEMpermutations[index_hyp[4]].resMEM_maxKinFit_Int.weight;
       FillWeightVectors(MEMpermutations[index_hyp[4]], tree.MEAllWeights_TTWJJ, tree.MEAllWeights_TTWJJ_log);
     }
   }
   if (index_hyp[5]!=-1){
     if (MEMpermutations[index_hyp[5]].computeHyp){
        nHypAllowed_TTbar += MEMpermutations[index_hyp[5]].nHypAllowed;
       tree.mc_mem_ttbarfl_weight = MEMpermutations[index_hyp[5]].resMEM_avgExl0.weight;
       tree.mc_mem_ttbarfl_weight_err = MEMpermutations[index_hyp[5]].resMEM_avgExl0.err;
       tree.mc_mem_ttbarfl_weight_chi2 = MEMpermutations[index_hyp[5]].resMEM_avgExl0.chi2;
       tree.mc_mem_ttbarfl_weight_time = MEMpermutations[index_hyp[5]].resMEM_avgExl0.time;
       tree.mc_mem_ttbarfl_weight_log = log(MEMpermutations[index_hyp[5]].resMEM_avgExl0.weight);
       tree.mc_mem_ttbarfl_weight_logmean = MEMpermutations[index_hyp[5]].resMEM_logavgExl0.weight;
       tree.mc_mem_ttbarfl_weight_avg = MEMpermutations[index_hyp[5]].resMEM_avgIncl0.weight;
       tree.mc_mem_ttbarfl_weight_max = MEMpermutations[index_hyp[5]].resMEM_max.weight;
       tree.mc_mem_ttbarfl_weight_JEC_up = MEMpermutations[index_hyp[5]].resMEM_avgExl0_JEC_up.weight;
       tree.mc_mem_ttbarfl_weight_JEC_down = MEMpermutations[index_hyp[5]].resMEM_avgExl0_JEC_down.weight;
       tree.mc_mem_ttbarfl_weight_JER_up = MEMpermutations[index_hyp[5]].resMEM_avgExl0_JER_up.weight;
       tree.mc_mem_ttbarfl_weight_JER_down = MEMpermutations[index_hyp[5]].resMEM_avgExl0_JER_down.weight;
       tree.mc_kin_ttbarfl_weight_logmax = log(MEMpermutations[index_hyp[5]].resKin_maxKinFit.weight);
       tree.mc_kin_ttbarfl_weight_logmaxint = log(MEMpermutations[index_hyp[5]].resKin_maxKinFit_Int.weight);
       tree.mc_mem_ttbarfl_weight_kinmax = MEMpermutations[index_hyp[5]].resMEM_maxKinFit.weight;
       tree.mc_mem_ttbarfl_weight_kinmaxint = MEMpermutations[index_hyp[5]].resMEM_maxKinFit_Int.weight;
       FillWeightVectors(MEMpermutations[index_hyp[5]], tree.MEAllWeights_TTbarfl, tree.MEAllWeights_TTbarfl_log);
     }
   }
   if (index_hyp[6]!=-1){
     if (MEMpermutations[index_hyp[6]].computeHyp){
       nHypAllowed_TTbar += MEMpermutations[index_hyp[6]].nHypAllowed;
       tree.mc_mem_ttbarsl_weight = MEMpermutations[index_hyp[6]].resMEM_avgExl0.weight;
       tree.mc_mem_ttbarsl_weight_err = MEMpermutations[index_hyp[6]].resMEM_avgExl0.err;
       tree.mc_mem_ttbarsl_weight_chi2 = MEMpermutations[index_hyp[6]].resMEM_avgExl0.chi2;
       tree.mc_mem_ttbarsl_weight_time = MEMpermutations[index_hyp[6]].resMEM_avgExl0.time;
       tree.mc_mem_ttbarsl_weight_log = log(MEMpermutations[index_hyp[6]].resMEM_avgExl0.weight);
       tree.mc_mem_ttbarsl_weight_avg = MEMpermutations[index_hyp[6]].resMEM_avgIncl0.weight;
       tree.mc_mem_ttbarsl_weight_logmean = MEMpermutations[index_hyp[6]].resMEM_logavgExl0.weight;
       tree.mc_mem_ttbarsl_weight_max = MEMpermutations[index_hyp[6]].resMEM_max.weight;
       tree.mc_mem_ttbarsl_weight_JEC_up = MEMpermutations[index_hyp[6]].resMEM_avgExl0_JEC_up.weight;
       tree.mc_mem_ttbarsl_weight_JEC_down = MEMpermutations[index_hyp[6]].resMEM_avgExl0_JEC_down.weight;
       tree.mc_mem_ttbarsl_weight_JER_up = MEMpermutations[index_hyp[6]].resMEM_avgExl0_JER_up.weight;
       tree.mc_mem_ttbarsl_weight_JER_down = MEMpermutations[index_hyp[6]].resMEM_avgExl0_JER_down.weight;
       tree.mc_kin_ttbarsl_weight_logmax = log(MEMpermutations[index_hyp[6]].resKin_maxKinFit.weight);
       tree.mc_kin_ttbarsl_weight_logmaxint = log(MEMpermutations[index_hyp[6]].resKin_maxKinFit_Int.weight);
       tree.mc_mem_ttbarsl_weight_kinmax = MEMpermutations[index_hyp[6]].resMEM_maxKinFit.weight;
       tree.mc_mem_ttbarsl_weight_kinmaxint = MEMpermutations[index_hyp[6]].resMEM_maxKinFit_Int.weight;
       FillWeightVectors(MEMpermutations[index_hyp[6]], tree.MEAllWeights_TTbarsl, tree.MEAllWeights_TTbarsl_log);
      }
    }

    if (index_hyp[1]!=-1 && index_hyp[2]!=-1) { 
      if (MEMpermutations[index_hyp[1]].computeHyp || MEMpermutations[index_hyp[2]].computeHyp)
        CombineHypotheses(MEMpermutations[index_hyp[1]], MEMpermutations[index_hyp[2]], &tree.mc_mem_tth_weight, &tree.mc_mem_tth_weight_log, &tree.mc_mem_tth_weight_err, &tree.mc_mem_tth_weight_chi2, &tree.mc_mem_tth_weight_time, &tree.mc_mem_tth_weight_avg, &tree.mc_mem_tth_weight_max, &tree.mc_mem_tth_weight_logmean, &tree.mc_kin_tth_weight_logmax, &tree.mc_kin_tth_weight_logmaxint, &tree.mc_mem_tth_weight_kinmax, &tree.mc_mem_tth_weight_kinmaxint, &tree.mc_mem_tth_weight_JEC_up, &tree.mc_mem_tth_weight_JEC_down, &tree.mc_mem_tth_weight_JER_up, &tree.mc_mem_tth_weight_JER_down);
    }
    else if (index_hyp[1]!=-1 && index_hyp[2]==-1){
      tree.mc_mem_tth_weight = tree.mc_mem_tthsl_weight;
      tree.mc_mem_tth_weight_log = tree.mc_mem_tthsl_weight_log;
      tree.mc_mem_tth_weight_err = tree.mc_mem_tthsl_weight_err;
      tree.mc_mem_tth_weight_chi2 = tree.mc_mem_tthsl_weight_chi2;
      tree.mc_mem_tth_weight_time = tree.mc_mem_tthsl_weight_time;
      tree.mc_mem_tth_weight_avg = tree.mc_mem_tthsl_weight_avg;
      tree.mc_mem_tth_weight_max = tree.mc_mem_tthsl_weight_max;
      tree.mc_mem_tth_weight_logmean = tree.mc_mem_tthsl_weight_logmean;
      tree.mc_kin_tth_weight_logmax = tree.mc_kin_tthsl_weight_logmax;
      tree.mc_kin_tth_weight_logmaxint = tree.mc_kin_tthsl_weight_logmaxint;
      tree.mc_mem_tth_weight_kinmax = tree.mc_mem_tthsl_weight_kinmax;
      tree.mc_mem_tth_weight_kinmaxint = tree.mc_mem_tthsl_weight_kinmaxint;
    }
    else if (index_hyp[1]==-1 && index_hyp[2]!=-1){
      tree.mc_mem_tth_weight = tree.mc_mem_tthfl_weight;
      tree.mc_mem_tth_weight_log = tree.mc_mem_tthfl_weight_log;
      tree.mc_mem_tth_weight_err = tree.mc_mem_tthfl_weight_err;
      tree.mc_mem_tth_weight_chi2 = tree.mc_mem_tthfl_weight_chi2;
      tree.mc_mem_tth_weight_time = tree.mc_mem_tthfl_weight_time;
      tree.mc_mem_tth_weight_avg = tree.mc_mem_tthfl_weight_avg;
      tree.mc_mem_tth_weight_max = tree.mc_mem_tthfl_weight_max;
      tree.mc_mem_tth_weight_logmean = tree.mc_mem_tthfl_weight_logmean;
      tree.mc_kin_tth_weight_logmax = tree.mc_kin_tthfl_weight_logmax;
      tree.mc_kin_tth_weight_logmaxint = tree.mc_kin_tthfl_weight_logmaxint;
      tree.mc_mem_tth_weight_kinmax = tree.mc_mem_tthfl_weight_kinmax;
      tree.mc_mem_tth_weight_kinmaxint = tree.mc_mem_tthfl_weight_kinmaxint;
    }

    if (index_hyp[5]!=-1 && index_hyp[6]!=-1) { 
      if (MEMpermutations[index_hyp[5]].computeHyp || MEMpermutations[index_hyp[6]].computeHyp)
        CombineHypotheses(MEMpermutations[index_hyp[5]], MEMpermutations[index_hyp[6]], &tree.mc_mem_ttbar_weight, &tree.mc_mem_ttbar_weight_log, &tree.mc_mem_ttbar_weight_err, &tree.mc_mem_ttbar_weight_chi2, &tree.mc_mem_ttbar_weight_time, &tree.mc_mem_ttbar_weight_avg, &tree.mc_mem_ttbar_weight_max, &tree.mc_mem_ttbar_weight_logmean, &tree.mc_kin_ttbar_weight_logmax, &tree.mc_kin_ttbar_weight_logmaxint, &tree.mc_mem_ttbar_weight_kinmax, &tree.mc_mem_ttbar_weight_kinmaxint, &tree.mc_mem_ttbar_weight_JEC_up, &tree.mc_mem_ttbar_weight_JEC_down, &tree.mc_mem_ttbar_weight_JER_up, &tree.mc_mem_ttbar_weight_JER_down);
    }
    else if (index_hyp[5]!=-1 && index_hyp[6]==-1){
      tree.mc_mem_ttbar_weight = tree.mc_mem_ttbarfl_weight;
      tree.mc_mem_ttbar_weight_log = tree.mc_mem_ttbarfl_weight_log;
      tree.mc_mem_ttbar_weight_err = tree.mc_mem_ttbarfl_weight_err;
      tree.mc_mem_ttbar_weight_chi2 = tree.mc_mem_ttbarfl_weight_chi2;
      tree.mc_mem_ttbar_weight_time = tree.mc_mem_ttbarfl_weight_time;
      tree.mc_mem_ttbar_weight_avg = tree.mc_mem_ttbarfl_weight_avg;
      tree.mc_mem_ttbar_weight_max = tree.mc_mem_ttbarfl_weight_max;
      tree.mc_mem_ttbar_weight_logmean = tree.mc_mem_ttbarfl_weight_logmean;
      tree.mc_kin_ttbar_weight_logmax = tree.mc_kin_ttbarfl_weight_logmax;
      tree.mc_kin_ttbar_weight_logmaxint = tree.mc_kin_ttbarfl_weight_logmaxint;
      tree.mc_mem_ttbar_weight_kinmax = tree.mc_mem_ttbarfl_weight_kinmax;
      tree.mc_mem_ttbar_weight_kinmaxint = tree.mc_mem_ttbarfl_weight_kinmaxint;
    }
    else if (index_hyp[5]==-1 && index_hyp[6]!=-1){
      tree.mc_mem_ttbar_weight = tree.mc_mem_ttbarsl_weight;
      tree.mc_mem_ttbar_weight_log = tree.mc_mem_ttbarsl_weight_log;
      tree.mc_mem_ttbar_weight_err = tree.mc_mem_ttbarsl_weight_err;
      tree.mc_mem_ttbar_weight_chi2 = tree.mc_mem_ttbarsl_weight_chi2;
      tree.mc_mem_ttbar_weight_time = tree.mc_mem_ttbarsl_weight_time;
      tree.mc_mem_ttbar_weight_avg = tree.mc_mem_ttbarsl_weight_avg;
      tree.mc_mem_ttbar_weight_max = tree.mc_mem_ttbarsl_weight_max;
      tree.mc_mem_ttbar_weight_logmean = tree.mc_mem_ttbarsl_weight_logmean;
      tree.mc_kin_ttbar_weight_logmax = tree.mc_kin_ttbarsl_weight_logmax;
      tree.mc_kin_ttbar_weight_logmaxint = tree.mc_kin_ttbarsl_weight_logmaxint;
      tree.mc_mem_ttbar_weight_kinmax = tree.mc_mem_ttbarsl_weight_kinmax;
      tree.mc_mem_ttbar_weight_kinmaxint = tree.mc_mem_ttbarsl_weight_kinmaxint;
    }


     if (tree.mc_mem_ttz_weight>0 && tree.mc_mem_tth_weight>0){
       tree.mc_mem_ttz_tth_likelihood = xsTTLL*tree.mc_mem_ttz_weight / (xsTTLL*tree.mc_mem_ttz_weight + xsTTH*tree.mc_mem_tth_weight);
       tree.mc_mem_ttz_tth_likelihood_nlog = -log(xsTTLL*tree.mc_mem_ttz_weight / (xsTTLL*tree.mc_mem_ttz_weight + xsTTH*tree.mc_mem_tth_weight));
       tree.mc_mem_ttz_tth_likelihood_max = xsTTLL*tree.mc_mem_ttz_weight_max / (xsTTLL*tree.mc_mem_ttz_weight_max + xsTTH*tree.mc_mem_tth_weight_max);
       tree.mc_mem_ttz_tth_likelihood_avg = xsTTLL*tree.mc_mem_ttz_weight_avg / (xsTTLL*tree.mc_mem_ttz_weight_avg + xsTTH*tree.mc_mem_tth_weight_avg);
     }
     if (tree.mc_mem_ttw_weight>0 && tree.mc_mem_tth_weight>0){
       tree.mc_mem_ttw_tth_likelihood = xsTTW*tree.mc_mem_ttw_weight / (xsTTW*tree.mc_mem_ttw_weight + xsTTH*tree.mc_mem_tth_weight);
       tree.mc_mem_ttw_tth_likelihood_nlog = -log(xsTTW*tree.mc_mem_ttw_weight / (xsTTW*tree.mc_mem_ttw_weight + xsTTH*tree.mc_mem_tth_weight));
       tree.mc_mem_ttw_tth_likelihood_max = xsTTW*tree.mc_mem_ttw_weight_max / (xsTTW*tree.mc_mem_ttw_weight_max + xsTTH*tree.mc_mem_tth_weight_max);
       tree.mc_mem_ttw_tth_likelihood_avg = xsTTW*tree.mc_mem_ttw_weight_avg / (xsTTW*tree.mc_mem_ttw_weight_avg + xsTTH*tree.mc_mem_tth_weight_avg);
     }
     if (tree.mc_mem_ttwjj_weight>0 && tree.mc_mem_tth_weight>0){
       tree.mc_mem_ttwjj_tth_likelihood = xsTTW*tree.mc_mem_ttwjj_weight / (xsTTW*tree.mc_mem_ttwjj_weight + xsTTH*tree.mc_mem_tth_weight);
       tree.mc_mem_ttwjj_tth_likelihood_nlog = -log(xsTTW*tree.mc_mem_ttwjj_weight / (xsTTW*tree.mc_mem_ttwjj_weight + xsTTH*tree.mc_mem_tth_weight));
       tree.mc_mem_ttwjj_tth_likelihood_max = xsTTW*tree.mc_mem_ttwjj_weight_max / (xsTTW*tree.mc_mem_ttwjj_weight_max + xsTTH*tree.mc_mem_tth_weight_max);
       tree.mc_mem_ttwjj_tth_likelihood_avg = xsTTW*tree.mc_mem_ttwjj_weight_avg / (xsTTW*tree.mc_mem_ttwjj_weight_avg + xsTTH*tree.mc_mem_tth_weight_avg);
     }
     if (tree.mc_mem_ttbar_weight>0 && tree.mc_mem_tth_weight>0){
       tree.mc_mem_ttbar_tth_likelihood = xsTTbar*tree.mc_mem_ttbar_weight / (xsTTbar*tree.mc_mem_ttbar_weight + xsTTH*tree.mc_mem_tth_weight);
       tree.mc_mem_ttbar_tth_likelihood_nlog = -log(xsTTbar*tree.mc_mem_ttbar_weight / (xsTTbar*tree.mc_mem_ttbar_weight + xsTTH*tree.mc_mem_tth_weight));
       tree.mc_mem_ttbar_tth_likelihood_max = xsTTbar*tree.mc_mem_ttbar_weight_max / (xsTTbar*tree.mc_mem_ttbar_weight_max + xsTTH*tree.mc_mem_tth_weight_max);
       tree.mc_mem_ttbar_tth_likelihood_avg = xsTTbar*tree.mc_mem_ttbar_weight_avg / (xsTTbar*tree.mc_mem_ttbar_weight_avg + xsTTH*tree.mc_mem_tth_weight_avg);
     }

     if (tree.mc_mem_ttz_weight>0 && tree.mc_mem_ttw_weight>0 && tree.mc_mem_tth_weight>0){
       tree.mc_mem_ttv_tth_likelihood = (xsTTLL*tree.mc_mem_ttz_weight + xsTTW*tree.mc_mem_ttw_weight) / (xsTTLL*tree.mc_mem_ttz_weight + xsTTW*tree.mc_mem_ttw_weight + xsTTH*tree.mc_mem_tth_weight);
       tree.mc_mem_ttv_tth_likelihood_nlog = -log((xsTTLL*tree.mc_mem_ttz_weight + xsTTW*tree.mc_mem_ttw_weight) / (xsTTLL*tree.mc_mem_ttz_weight + xsTTW*tree.mc_mem_ttw_weight + xsTTH*tree.mc_mem_tth_weight));
       tree.mc_mem_ttv_tth_likelihood_max = (xsTTLL*tree.mc_mem_ttz_weight_max + xsTTW*tree.mc_mem_ttw_weight_max) / (xsTTLL*tree.mc_mem_ttz_weight_max + xsTTW*tree.mc_mem_ttw_weight_max + xsTTH*tree.mc_mem_tth_weight_max);    
       tree.mc_mem_ttv_tth_likelihood_avg = (xsTTLL*tree.mc_mem_ttz_weight_avg + xsTTW*tree.mc_mem_ttw_weight_avg) / (xsTTLL*tree.mc_mem_ttz_weight_avg + xsTTW*tree.mc_mem_ttw_weight_avg + xsTTH*tree.mc_mem_tth_weight_avg);
     }
     if (tree.mc_mem_ttz_weight>0 && tree.mc_mem_ttwjj_weight>0 && tree.mc_mem_tth_weight>0){
       tree.mc_mem_ttvjj_tth_likelihood = (xsTTLL*tree.mc_mem_ttz_weight + xsTTW*tree.mc_mem_ttwjj_weight) / (xsTTLL*tree.mc_mem_ttz_weight + xsTTW*tree.mc_mem_ttwjj_weight + xsTTH*tree.mc_mem_tth_weight);
       tree.mc_mem_ttvjj_tth_likelihood_nlog = -log((xsTTLL*tree.mc_mem_ttz_weight + xsTTW*tree.mc_mem_ttwjj_weight) / (xsTTLL*tree.mc_mem_ttz_weight + xsTTW*tree.mc_mem_ttwjj_weight + xsTTH*tree.mc_mem_tth_weight));
       tree.mc_mem_ttvjj_tth_likelihood_max = (xsTTLL*tree.mc_mem_ttz_weight_max + xsTTW*tree.mc_mem_ttwjj_weight_max) / (xsTTLL*tree.mc_mem_ttz_weight_max + xsTTW*tree.mc_mem_ttwjj_weight_max + xsTTH*tree.mc_mem_tth_weight_max);
       tree.mc_mem_ttvjj_tth_likelihood_avg = (xsTTLL*tree.mc_mem_ttz_weight_avg + xsTTW*tree.mc_mem_ttwjj_weight_avg) / (xsTTLL*tree.mc_mem_ttz_weight_avg + xsTTW*tree.mc_mem_ttwjj_weight_avg + xsTTH*tree.mc_mem_tth_weight_avg);
     }

     cout << "PRINT DISCRIMINATION"<<endl;
     if (tree.mc_mem_tth_weight>0){
       cout << "TTH nHypAllowed="<<nHypAllowed_TTH<<endl;
       cout << "TTH weight="<<tree.mc_mem_tth_weight<<"; log(weight)="<<tree.mc_mem_tth_weight_log<<endl;
     }
     if (tree.mc_mem_ttz_weight>0){
       cout << "TTZ nHypAllowed="<<MEMpermutations[index_hyp[0]].nHypAllowed<<endl;
       cout << "TTZ weight="<<tree.mc_mem_ttz_weight<<" ; log(weight)="<<tree.mc_mem_ttz_weight_log<<endl;
     }
     if (tree.mc_mem_ttw_weight>0){
       cout << "TTW nHypAllowed="<<MEMpermutations[index_hyp[3]].nHypAllowed<<endl;
       cout << "TTW weight="<<tree.mc_mem_ttw_weight<<" log(weight)="<<tree.mc_mem_ttw_weight_log<<endl;
     }
     if (tree.mc_mem_ttwjj_weight>0){
       cout << "TTWJJ nHypAllowed="<<MEMpermutations[index_hyp[4]].nHypAllowed<<endl;
       cout << "TTWJJ weight="<<tree.mc_mem_ttwjj_weight<<" log(weight)="<<tree.mc_mem_ttwjj_weight_log<<endl;
     }
     if (tree.mc_mem_ttbar_weight>0){
       cout << "TTbar nHypAllowed="<<nHypAllowed_TTbar<<endl;
       cout << "TTbar weight="<<tree.mc_mem_ttbar_weight<<" ; log(weight)="<<tree.mc_mem_ttbar_weight_log<<endl;
     }
     if (tree.mc_mem_tth_weight>0 && tree.mc_mem_ttz_weight>0)
       cout << "TTHvsTTZ discrim="<<tree.mc_mem_ttz_tth_likelihood<<"; -log(discrim)="<< tree.mc_mem_ttz_tth_likelihood_nlog<<endl;
     if (tree.mc_mem_ttw_weight>0 && tree.mc_mem_tth_weight>0)
       cout << "TTHvsTTW discrim="<<tree.mc_mem_ttw_tth_likelihood<<"; -log(discrim)="<< tree.mc_mem_ttw_tth_likelihood_nlog<<endl;
     if (tree.mc_mem_ttwjj_tth_likelihood>0 && tree.mc_mem_tth_weight>0)
       cout << "TTHvsTTWJJ discrim="<<tree.mc_mem_ttwjj_tth_likelihood<<"; -log(discrim)="<< tree.mc_mem_ttwjj_tth_likelihood_nlog<<endl;
     if (tree.mc_mem_tth_weight>0 && tree.mc_mem_ttz_weight>0 && tree.mc_mem_ttw_weight>0)
       cout << "TTHvsTTV discrim="<<tree.mc_mem_ttv_tth_likelihood<<"; -log(discrim)="<< tree.mc_mem_ttv_tth_likelihood_nlog<<endl;
     if (tree.mc_mem_tth_weight>0 && tree.mc_mem_ttz_weight>0 && tree.mc_mem_ttwjj_weight>0)
       cout << "TTHvsTTVjj discrim="<<tree.mc_mem_ttvjj_tth_likelihood<<"; -log(discrim)="<< tree.mc_mem_ttvjj_tth_likelihood_nlog<<endl;
     
     cout << "End of event ---"<<endl;
     cout << endl;

 
   }

   tree.tOutput->Fill();

  }

  tree.fOutput->cd();

  for (int ihist=0; ihist<nErrHist; ihist++) tree.hMEPhaseSpace_Error[ihist]->Write();
  for (int ihyp=0; ihyp<nhyp; ihyp++) {
    tree.hMEPhaseSpace_ErrorTot[ihyp]->Write();    
    tree.hMEPhaseSpace_ErrorTot_Pass[ihyp]->Write();
    tree.hMEPhaseSpace_ErrorTot_Fail[ihyp]->Write();
  }
  for(int ih=0; ih<nhyp; ih++){
     for(int index_cat=0; index_cat<12; index_cat++){
        tree.hMEAllWeights[ih][index_cat]->Write();
        tree.hMEAllWeights_nlog[ih][index_cat]->Write();
     }
  }

  tree.tOutput->Write();
  tree.fOutput->Close();

 }

}
