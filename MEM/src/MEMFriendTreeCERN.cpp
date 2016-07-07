#include "../interface/MEMFriendTreeCERN.h"

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
#include "ConfigParser.h"
#include "Tools.h"
#include "Permutations.h"

#include <iostream>
#include <ctime>

using namespace std;

const float mZ = 91.18;
const float mW = 80.38;

MEMFriendTreeCERN::~MEMFriendTreeCERN() {
    delete multiLepton;
    delete cfgParser;
    delete hypIntegrator;
    for (auto & p : MEMpermutations) delete p;
}

void MEMFriendTreeCERN::init(const std::string &config) {
    multiLepton = new MultiLepton();
    cfgParser = new ConfigParser();
    hypIntegrator = new HypIntegrator(); 

    cfgParser->GetConfigFromFile(config);
    cfgParser->LoadHypotheses(&nhyp, &shyp, &hyp, &nPointsHyp, &index_hyp);

    cfgParser->LoadIntegrationRange(&multiLepton->JetTFfracmin, &multiLepton->JetTFfracmax, &multiLepton->NeutMaxE);

    cout << "MEM weight computation: initialization" << endl;

    for (int ih=0; ih<nhyp; ih++){
        if (shyp[ih]=="TTLL") index[ih] = 0;
        if (shyp[ih]=="TTHfl") index[ih] = 1;
        if (shyp[ih]=="TTHsl") index[ih] = 2;
        if (shyp[ih]=="TTW") index[ih] = 3;
        if (shyp[ih]=="TTWJJ") index[ih] = 4;
        if (shyp[ih]=="TTbarfl") index[ih] = 5;
        if (shyp[ih]=="TTbarsl") index[ih] = 6;
    }

    cfgParser->LoadOptim(&doOptim);

    cfgParser->LoadOptim(&doOptimTopHad, &doOptimTopLep, &doOptimHiggs, &doOptimW);

    hypIntegrator->InitializeIntegrator(cfgParser);
    doMinimization = cfgParser->valDoMinimization;

    cfgParser->LoadJetChoice(&JetChoice);

    nPermutationJetSyst = cfgParser->nJetSyst;
    xsTTH = hypIntegrator->meIntegrator->xsTTH * hypIntegrator->meIntegrator->brTopHad * hypIntegrator->meIntegrator->brTopLep * hypIntegrator->meIntegrator->brHiggsFullLep;
    xsTTLL = hypIntegrator->meIntegrator->xsTTLL * hypIntegrator->meIntegrator->brTopHad * hypIntegrator->meIntegrator->brTopLep;
    xsTTW = hypIntegrator->meIntegrator->xsTTW * hypIntegrator->meIntegrator->brTopLep * hypIntegrator->meIntegrator->brTopLep;
    xsTTbar = hypIntegrator->meIntegrator->xsTTbar * hypIntegrator->meIntegrator->brTopHad * hypIntegrator->meIntegrator->brTopLep;

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

    for(int ih=0; ih<nhyp; ih++){
        if(shyp[ih]=="TTLL") index_XS[ih]=xsTTLL;
        if(shyp[ih]=="TTW" || shyp[ih]=="TTWJJ") index_XS[ih]=xsTTW;
        if(shyp[ih]=="TTH" || shyp[ih]=="TTHsl" || shyp[ih]=="TTHfl") index_XS[ih]=xsTTH;
        if(shyp[ih]=="TTbar" || shyp[ih]=="TTbarsl" || shyp[ih]=="TTbarfl") index_XS[ih]=xsTTbar;
    }

    MEMpermutations.reserve(nhyp);
    for (int i = 0; i < nhyp; ++i) MEMpermutations.push_back(new Permutations());
}

void MEMFriendTreeCERN::clear() {
  multiLepton->Leptons.clear();
  multiLepton->Jets.clear();
  multiLepton->Bjets.clear();
  multiLepton->AllJets.clear();
  multiLepton->JetsHighestPt.clear();
  multiLepton->JetsClosestMw.clear();
  multiLepton->JetsLowestMjj.clear();
}

void MEMFriendTreeCERN::addLepton(TLorentzVector p4, int pdgId) {
  multiLepton->FillParticle("lepton", pdgId, p4);
}

void MEMFriendTreeCERN::addJet(const std::string &what,  TLorentzVector p4, float CSV) {
    multiLepton->FillParticle(what, 0, CSV, 0,0,0,0, p4);
}

void MEMFriendTreeCERN::addBJet(TLorentzVector p4, float CSV) {
    multiLepton->FillParticle("bjet", 5, CSV, 0,0,0,0, p4);
}

void MEMFriendTreeCERN::setMET(TLorentzVector p4, double cov00, double cov01, double cov10, double cov11, double mHT) {
  multiLepton->mET = p4;
  multiLepton->mET_cov00 = cov00;
  multiLepton->mET_cov01 = cov01;
  multiLepton->mET_cov10 = cov10;
  multiLepton->mET_cov11 = cov11;
  multiLepton->mHT = mHT;
}

bool MEMFriendTreeCERN::setCategory(const std::string &cat) {
    bool found = false;
    for (int catJets = 0; catJets < 12; ++catJets) {
        if (index_CatJets[catJets] == cat) {
            multiLepton->kCatJets = catJets;
            found = true;
            break;
        }
    }
    return found;
}

std::map<std::string,float> MEMFriendTreeCERN::compute() {
    std::map<std::string,float> ret;
    for (int ih=0; ih<nhyp; ih++){
        MEMpermutations[ih]->SetMultiLepton(multiLepton, hypIntegrator);
        int initresult = MEMpermutations[ih]->InitializeHyp(hypIntegrator, hyp[ih], nPointsHyp[ih], shyp[ih], doMinimization, JetChoice, nPermutationJetSyst);
        if (initresult==1) {
            MEMpermutations[ih]->LoopPermutations(hypIntegrator);
            ret[shyp[ih]] = MEMpermutations[ih]->resMEM_avgExl0.weight;
            ret[shyp[ih]+"_err"] = MEMpermutations[ih]->resMEM_avgExl0.err;
            ret[shyp[ih]+"_time"] = MEMpermutations[ih]->resMEM_avgExl0.time;
            ret[shyp[ih]+"_chi2"] = MEMpermutations[ih]->resMEM_avgExl0.chi2;
            ret[shyp[ih]+"_nHypAll"] = MEMpermutations[ih]->nHypAllowed;
            ret[shyp[ih]+"_nNull"] = MEMpermutations[ih]->nNullResult;
        }
    }
    if (index[1]!=-1 && index[2]!=-1) {
        if (MEMpermutations[index[1]]->computeHyp || MEMpermutations[index[2]]->computeHyp) {
            double mem_tth_weight, mem_tth_weight_log, mem_tth_weight_err, mem_tth_weight_avg, mem_tth_weight_max, mem_tth_weight_logmean, kin_tth_weight_logmax, kin_tth_weight_logmaxint, mem_tth_weight_kinmax, mem_tth_weight_kinmaxint, mem_tth_weight_JEC_up, mem_tth_weight_JEC_down, mem_tth_weight_JER_up, mem_tth_weight_JER_down;
            float mem_tth_weight_chi2, mem_tth_weight_time;
            CombineHypotheses(*MEMpermutations[index[1]], *MEMpermutations[index[2]], &mem_tth_weight, &mem_tth_weight_log, &mem_tth_weight_err, &mem_tth_weight_chi2, &mem_tth_weight_time, &mem_tth_weight_avg, &mem_tth_weight_max, &mem_tth_weight_logmean, &kin_tth_weight_logmax, &kin_tth_weight_logmaxint, &mem_tth_weight_kinmax, &mem_tth_weight_kinmaxint, &mem_tth_weight_JEC_up, &mem_tth_weight_JEC_down, &mem_tth_weight_JER_up, &mem_tth_weight_JER_down);
            ret["TTH"] = mem_tth_weight;
            ret["TTH_err"] = mem_tth_weight_err;
            ret["TTH_chi2"] = mem_tth_weight_chi2;
            ret["TTH_time"] = mem_tth_weight_time;
        }
    }
    return ret;
}
