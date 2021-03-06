#include "../include/TTbarHiggsMultileptonAnalysis.h"
#include "TSystem.h"
#include "SignalExtractionMVA.cxx"
#include "Helper.cxx"
#include "BTagging.cxx"
#include "FakeRate.cxx"
#include "ChargeFlip.cxx"

#define kCat_3l_2b_2j   0
#define kCat_3l_1b_2j   1
#define kCat_3l_2b_1j   2
#define kCat_3l_1b_1j   3
#define kCat_3l_2b_0j   4
#define kCat_4l_2b      5
#define kCat_4l_1b      6
#define kCat_2lss_2b_4j 7
#define kCat_2lss_1b_4j 8
#define kCat_2lss_2b_3j 9
#define kCat_2lss_1b_3j 10
#define kCat_2lss_2b_2j 11

#define DEBUG           false

using namespace std;

TTbarHiggsMultileptonAnalysis::TTbarHiggsMultileptonAnalysis() 
{
    _printLHCO_MC = false;
    _processLHCO_MC = -1;

    _printLHCO_RECO = false;
    _processLHCO_RECO = -1;
}

TTbarHiggsMultileptonAnalysis::~TTbarHiggsMultileptonAnalysis() 
{
    delete theHistoManager;

    if (_doSystCombine) 
    {
        delete  histoManager_2lss_mm_0tau_bl_neg;
        delete  histoManager_2lss_mm_0tau_bt_neg;
        delete  histoManager_2lss_ee_0tau_bl_neg;
        delete  histoManager_2lss_ee_0tau_bt_neg;
        delete  histoManager_2lss_em_0tau_bl_neg;
        delete  histoManager_2lss_em_0tau_bt_neg;

        delete  histoManager_2lss_mm_0tau_bl_pos;
        delete  histoManager_2lss_mm_0tau_bt_pos;
        delete  histoManager_2lss_ee_0tau_bl_pos;
        delete  histoManager_2lss_ee_0tau_bt_pos;
        delete  histoManager_2lss_em_0tau_bl_pos;
        delete  histoManager_2lss_em_0tau_bt_pos;

        delete  histoManager_3l_bl_neg;
        delete  histoManager_3l_bt_neg;

        delete  histoManager_3l_bl_pos;
        delete  histoManager_3l_bt_pos;
    }
}

TTbarHiggsMultileptonAnalysis::TTbarHiggsMultileptonAnalysis(TString inputFileName, TChain *tree,
        TString sampleName, TString treeName, TString outputFileName, bool isdata, bool doSystCombine, float xsec, float lumi, int nowe, int nmax)
{    

    //
    _isdata = isdata;
    _xsec = xsec;
    _lumi = lumi;
    _nowe = nowe;
    _nmax = nmax;
    _outputFileName = outputFileName;
    _sampleName = sampleName;
    _process = "toto";
    _doSystCombine = doSystCombine;
    if (_isdata) _doSystCombine = false;

    //
    _printLHCO_MC = false;
    _processLHCO_MC = -1;

    _printLHCO_RECO = false;
    _processLHCO_RECO = -1;

    //
    //_file_PVreweighting = TFile::Open("/home-pbs/lebihan/someone/ttH_070116/ttH/NtupleAnalyzer/test/PUweight.root");
    //_h_PV = (TH1F*)_file_PVreweighting->Get("PU_reweighting");

    //
    tree = new TChain(treeName.Data());

    std::ifstream infile;
    infile.open(inputFileName.Data());
    std::string ifile = "";
    while( getline(infile, ifile) )
    {
        std::string fnameStr = std::string(ifile);
        tree->Add(fnameStr.c_str());
    }   
    infile.close();
    Init(tree);

    theHistoManager = new HistoManager();

    if (_doSystCombine)
    { 
        std::cout <<"doSystCombine" << std::endl;
        histoManager_2lss_mm_0tau_bl_neg = new HistoManager();
        histoManager_2lss_mm_0tau_bt_neg = new HistoManager();
        histoManager_2lss_ee_0tau_bl_neg = new HistoManager();
        histoManager_2lss_ee_0tau_bt_neg = new HistoManager();
        histoManager_2lss_em_0tau_bl_neg = new HistoManager();
        histoManager_2lss_em_0tau_bt_neg = new HistoManager();

        histoManager_2lss_mm_0tau_bl_pos = new HistoManager();
        histoManager_2lss_mm_0tau_bt_pos = new HistoManager();
        histoManager_2lss_ee_0tau_bl_pos = new HistoManager();
        histoManager_2lss_ee_0tau_bt_pos = new HistoManager();
        histoManager_2lss_em_0tau_bl_pos = new HistoManager();
        histoManager_2lss_em_0tau_bt_pos = new HistoManager();

        histoManager_3l_bl_neg = new HistoManager();
        histoManager_3l_bt_neg = new HistoManager();
        histoManager_3l_bl_pos = new HistoManager();
        histoManager_3l_bt_pos = new HistoManager();

        histoManager_2lss_mm_0tau_bl_neg = new HistoManager();
        histoManager_2lss_mm_0tau_bt_neg = new HistoManager();
        histoManager_2lss_ee_0tau_bl_neg = new HistoManager();
        histoManager_2lss_ee_0tau_bt_neg = new HistoManager();
        histoManager_2lss_em_0tau_bl_neg = new HistoManager();
        histoManager_2lss_em_0tau_bt_neg = new HistoManager();

        histoManager_2lss_mm_0tau_bl_pos = new HistoManager();
        histoManager_2lss_mm_0tau_bt_pos = new HistoManager();
        histoManager_2lss_ee_0tau_bl_pos = new HistoManager();
        histoManager_2lss_ee_0tau_bt_pos = new HistoManager();
        histoManager_2lss_em_0tau_bl_pos = new HistoManager();
        histoManager_2lss_em_0tau_bt_pos = new HistoManager();

        histoManager_3l_bl_neg = new HistoManager();
        histoManager_3l_bt_neg = new HistoManager();
        histoManager_3l_bl_pos = new HistoManager();
        histoManager_3l_bt_pos = new HistoManager();
    }

    TString outputfileNameRoot = _outputFileName+".root";
    outputfile = new TFile(outputfileNameRoot.Data(), "recreate");  

}

void TTbarHiggsMultileptonAnalysis::InitLHCO(int process_MC, int process_RECO) 
{  
    if (!_isdata)
    {    
        _printLHCO_MC = true;
        _processLHCO_MC = process_MC;     
        TString fout_MC_path(_outputFileName+"_LHCO_MC.txt"); 
        fout_MC.open(fout_MC_path.Data());
    }

    _printLHCO_RECO = true;
    _processLHCO_RECO = process_RECO;
    TString fout_RECO_path(_outputFileName+"_LHCO_RECO.txt");
    fout_RECO.open(fout_RECO_path.Data());

    fline00 = "#   typ	  eta	 phi	   pt  jmass  ntrk  btag   had/em  dummy dummy";
    del = "    ";
    trig = "8";   
}

void TTbarHiggsMultileptonAnalysis::createHistograms()
{    
    outputfile->cd();
    initializeOutputTree();

    // General
    theHistoManager->addHisto("CutFlow",        "noSel",        "noChannel",   "",  10,   0,     10);

    // Preselection variables

    theHistoManager->addHisto("MuonPt",         "noSel",        "noChannel",   "",   25,   0,   200);
    theHistoManager->addHisto("MuonEta",        "noSel",        "noChannel",   "",   25,  -3,     3);
    theHistoManager->addHisto("MuonMVA",        "noSel",        "noChannel",   "",   20,   0,     1);
    theHistoManager->addHisto("ElectronPt",     "noSel",        "noChannel",   "",   25,   0,   200);
    theHistoManager->addHisto("ElectronEta",    "noSel",        "noChannel",   "",   25,  -3,     3);
    theHistoManager->addHisto("ElectronMVA",    "noSel",        "noChannel",   "",   20,   0,     1);
    theHistoManager->addHisto("TauPt",          "noSel",        "noChannel",   "",   25,   0,   200);
    theHistoManager->addHisto("TauEta",         "noSel",        "noChannel",   "",   25,  -3,     3);
    theHistoManager->addHisto("JetPt",          "noSel",        "noChannel",   "",   25,   0,   200);
    theHistoManager->addHisto("JetEta",         "noSel",        "noChannel",   "",   25,  -3,     3);
    theHistoManager->addHisto("JetCSVv2",       "noSel",        "noChannel",   "",   20,   0,     1);
    theHistoManager->addHisto("MET",            "noSel",        "noChannel",   "",   50,   0,   200);

    // Signal Region Two Leptons

    theHistoManager->addHisto2D("SelectedLeptonsVsJets",            "noSel",        "ttH_2lss",   "",    8,    0,    7,    8,    0,    7);
    theHistoManager->addHisto2D("SelectedLeptonsVsBJets",           "noSel",        "ttH_2lss",   "",    8,    0,    7,    8,    0,    7);

    theHistoManager->addHisto("WeightCSV_min",                      "finalSel",     "ttH_2lss",   "",    30,    0,    3);
    theHistoManager->addHisto("WeightCSV_max",                      "finalSel",     "ttH_2lss",   "",    30,    0,    3);

    theHistoManager->addHisto("CutFlow",                            "finalSel",     "ttH_2lss",   "",    10,   -1,    9);

    theHistoManager->addHisto("CutFlow",                            "finalSel",     "ttH_2lss_ee",   "",    10,   -1,    9);
    theHistoManager->addHisto("LeadingLeptonPt",                    "finalSel",     "ttH_2lss_ee",   "",    20,    0,  200);
    theHistoManager->addHisto("SubLeadingLeptonPt",                 "finalSel",     "ttH_2lss_ee",   "",    20,    0,  200);
    theHistoManager->addHisto("MET",                                "finalSel",     "ttH_2lss_ee",   "",    20,    0,  200);
    theHistoManager->addHisto("MHT",                                "finalSel",     "ttH_2lss_ee",   "",    20,    0,  200);
    theHistoManager->addHisto("MetLD",                              "finalSel",     "ttH_2lss_ee",   "",    15, -0.2,  1.4);
    theHistoManager->addHisto("TauMultiplicity",                    "finalSel",     "ttH_2lss_ee",   "",    10,   -1,    9);
    theHistoManager->addHisto("JetMultiplicity",                    "finalSel",     "ttH_2lss_ee",   "",    10,   -1,    9);
    theHistoManager->addHisto("LooseBJetMultiplicity",              "finalSel",     "ttH_2lss_ee",   "",    10,   -1,    9);
    theHistoManager->addHisto("MediumBJetMultiplicity",             "finalSel",     "ttH_2lss_ee",   "",    10,   -1,    9);
    theHistoManager->addHisto("NumberOfSelectedLeptons",            "finalSel",     "ttH_2lss_ee",   "",    10,   -1,    9);

    theHistoManager->addHisto("CutFlow",                            "finalSel",     "ttH_2lss_em",   "",    10,   -1,    9);
    theHistoManager->addHisto("LeadingLeptonPt",                    "finalSel",     "ttH_2lss_em",   "",    20,    0,  200);
    theHistoManager->addHisto("SubLeadingLeptonPt",                 "finalSel",     "ttH_2lss_em",   "",    20,    0,  200);
    theHistoManager->addHisto("MET",                                "finalSel",     "ttH_2lss_em",   "",    20,    0,  200);
    theHistoManager->addHisto("MHT",                                "finalSel",     "ttH_2lss_em",   "",    20,    0,  200);
    theHistoManager->addHisto("MetLD",                              "finalSel",     "ttH_2lss_em",   "",    15, -0.2,  1.4);
    theHistoManager->addHisto("TauMultiplicity",                    "finalSel",     "ttH_2lss_em",   "",    10,   -1,    9);
    theHistoManager->addHisto("JetMultiplicity",                    "finalSel",     "ttH_2lss_em",   "",    10,   -1,    9);
    theHistoManager->addHisto("LooseBJetMultiplicity",              "finalSel",     "ttH_2lss_em",   "",    10,   -1,    9);
    theHistoManager->addHisto("MediumBJetMultiplicity",             "finalSel",     "ttH_2lss_em",   "",    10,   -1,    9);
    theHistoManager->addHisto("NumberOfSelectedLeptons",            "finalSel",     "ttH_2lss_em",   "",    10,   -1,    9);

    theHistoManager->addHisto("CutFlow",                            "finalSel",     "ttH_2lss_mm",   "",    10,   -1,    9);
    theHistoManager->addHisto("LeadingLeptonPt",                    "finalSel",     "ttH_2lss_mm",   "",    20,    0,  200);
    theHistoManager->addHisto("SubLeadingLeptonPt",                 "finalSel",     "ttH_2lss_mm",   "",    20,    0,  200);
    theHistoManager->addHisto("MET",                                "finalSel",     "ttH_2lss_mm",   "",    20,    0,  200);
    theHistoManager->addHisto("MHT",                                "finalSel",     "ttH_2lss_mm",   "",    20,    0,  200);
    theHistoManager->addHisto("MetLD",                              "finalSel",     "ttH_2lss_mm",   "",    15, -0.2,  1.4);
    theHistoManager->addHisto("TauMultiplicity",                    "finalSel",     "ttH_2lss_mm",   "",    10,   -1,    9);
    theHistoManager->addHisto("JetMultiplicity",                    "finalSel",     "ttH_2lss_mm",   "",    10,   -1,    9);
    theHistoManager->addHisto("LooseBJetMultiplicity",              "finalSel",     "ttH_2lss_mm",   "",    10,   -1,    9);
    theHistoManager->addHisto("MediumBJetMultiplicity",             "finalSel",     "ttH_2lss_mm",   "",    10,   -1,    9);
    theHistoManager->addHisto("NumberOfSelectedLeptons",            "finalSel",     "ttH_2lss_mm",   "",    10,   -1,    9);

     // Signal Region Two Leptons - Fake Rate

    theHistoManager->addHisto2D("SelectedLeptonsVsJets",            "noSel",        "ttH_2lss_fr",   "",    8,    0,    7,    8,    0,    7);
    theHistoManager->addHisto2D("SelectedLeptonsVsBJets",           "noSel",        "ttH_2lss_fr",   "",    8,    0,    7,    8,    0,    7);

    theHistoManager->addHisto("WeightCSV_min",                      "noSel",        "ttH_2lss_fr",   "",    30,    0,    3);
    theHistoManager->addHisto("WeightCSV_max",                      "noSel",        "ttH_2lss_fr",   "",    30,    0,    3);

    theHistoManager->addHisto("CutFlow",                            "finalSel",     "ttH_2lss_fr",   "",    10,   -1,    9);

    theHistoManager->addHisto("CutFlow",                            "finalSel",     "ttH_2lss_ee_fr",   "",    10,   -1,    9);
    theHistoManager->addHisto("LeadingLeptonPt",                    "finalSel",     "ttH_2lss_ee_fr",   "",    20,    0,  200);
    theHistoManager->addHisto("SubLeadingLeptonPt",                 "finalSel",     "ttH_2lss_ee_fr",   "",    20,    0,  200);
    theHistoManager->addHisto("MET",                                "finalSel",     "ttH_2lss_ee_fr",   "",    20,    0,  200);
    theHistoManager->addHisto("MHT",                                "finalSel",     "ttH_2lss_ee_fr",   "",    20,    0,  200);
    theHistoManager->addHisto("MetLD",                              "finalSel",     "ttH_2lss_ee_fr",   "",    15, -0.2,  1.4);
    theHistoManager->addHisto("TauMultiplicity",                    "finalSel",     "ttH_2lss_ee_fr",   "",    10,   -1,    9);
    theHistoManager->addHisto("JetMultiplicity",                    "finalSel",     "ttH_2lss_ee_fr",   "",    10,   -1,    9);
    theHistoManager->addHisto("LooseBJetMultiplicity",              "finalSel",     "ttH_2lss_ee_fr",   "",    10,   -1,    9);
    theHistoManager->addHisto("MediumBJetMultiplicity",             "finalSel",     "ttH_2lss_ee_fr",   "",    10,   -1,    9);
    theHistoManager->addHisto("NumberOfSelectedLeptons",            "finalSel",     "ttH_2lss_ee_fr",   "",    10,   -1,    9);

    theHistoManager->addHisto("CutFlow",                            "finalSel",     "ttH_2lss_em_fr",   "",    10,   -1,    9);
    theHistoManager->addHisto("LeadingLeptonPt",                    "finalSel",     "ttH_2lss_em_fr",   "",    20,    0,  200);
    theHistoManager->addHisto("SubLeadingLeptonPt",                 "finalSel",     "ttH_2lss_em_fr",   "",    20,    0,  200);
    theHistoManager->addHisto("MET",                                "finalSel",     "ttH_2lss_em_fr",   "",    20,    0,  200);
    theHistoManager->addHisto("MHT",                                "finalSel",     "ttH_2lss_em_fr",   "",    20,    0,  200);
    theHistoManager->addHisto("MetLD",                              "finalSel",     "ttH_2lss_em_fr",   "",    15, -0.2,  1.4);
    theHistoManager->addHisto("TauMultiplicity",                    "finalSel",     "ttH_2lss_em_fr",   "",    10,   -1,    9);
    theHistoManager->addHisto("JetMultiplicity",                    "finalSel",     "ttH_2lss_em_fr",   "",    10,   -1,    9);
    theHistoManager->addHisto("LooseBJetMultiplicity",              "finalSel",     "ttH_2lss_em_fr",   "",    10,   -1,    9);
    theHistoManager->addHisto("MediumBJetMultiplicity",             "finalSel",     "ttH_2lss_em_fr",   "",    10,   -1,    9);
    theHistoManager->addHisto("NumberOfSelectedLeptons",            "finalSel",     "ttH_2lss_em_fr",   "",    10,   -1,    9);

    theHistoManager->addHisto("CutFlow",                            "finalSel",     "ttH_2lss_mm_fr",   "",    10,   -1,    9);
    theHistoManager->addHisto("LeadingLeptonPt",                    "finalSel",     "ttH_2lss_mm_fr",   "",    20,    0,  200);
    theHistoManager->addHisto("SubLeadingLeptonPt",                 "finalSel",     "ttH_2lss_mm_fr",   "",    20,    0,  200);
    theHistoManager->addHisto("MET",                                "finalSel",     "ttH_2lss_mm_fr",   "",    20,    0,  200);
    theHistoManager->addHisto("MHT",                                "finalSel",     "ttH_2lss_mm_fr",   "",    20,    0,  200);
    theHistoManager->addHisto("MetLD",                              "finalSel",     "ttH_2lss_mm_fr",   "",    15, -0.2,  1.4);
    theHistoManager->addHisto("TauMultiplicity",                    "finalSel",     "ttH_2lss_mm_fr",   "",    10,   -1,    9);
    theHistoManager->addHisto("JetMultiplicity",                    "finalSel",     "ttH_2lss_mm_fr",   "",    10,   -1,    9);
    theHistoManager->addHisto("LooseBJetMultiplicity",              "finalSel",     "ttH_2lss_mm_fr",   "",    10,   -1,    9);
    theHistoManager->addHisto("MediumBJetMultiplicity",             "finalSel",     "ttH_2lss_mm_fr",   "",    10,   -1,    9);
    theHistoManager->addHisto("NumberOfSelectedLeptons",            "finalSel",     "ttH_2lss_mm_fr",   "",    10,   -1,    9);

   // Signal Region Two Leptons - Charge Flip

    theHistoManager->addHisto2D("SelectedLeptonsVsJets",            "noSel",        "ttH_2lss_cf",   "",    8,    0,    7,    8,    0,    7);
    theHistoManager->addHisto2D("SelectedLeptonsVsBJets",           "noSel",        "ttH_2lss_cf",   "",    8,    0,    7,    8,    0,    7);

    theHistoManager->addHisto("WeightCSV_min",                      "noSel",        "ttH_2lss_cf",   "",    30,    0,    3);
    theHistoManager->addHisto("WeightCSV_max",                      "noSel",        "ttH_2lss_cf",   "",    30,    0,    3);

    theHistoManager->addHisto("CutFlow",                            "finalSel",     "ttH_2lss_cf",   "",    10,   -1,    9);

    theHistoManager->addHisto("CutFlow",                            "finalSel",     "ttH_2lss_ee_cf",   "",    10,   -1,    9);
    theHistoManager->addHisto("LeadingLeptonPt",                    "finalSel",     "ttH_2lss_ee_cf",   "",    20,    0,  200);
    theHistoManager->addHisto("SubLeadingLeptonPt",                 "finalSel",     "ttH_2lss_ee_cf",   "",    20,    0,  200);
    theHistoManager->addHisto("MET",                                "finalSel",     "ttH_2lss_ee_cf",   "",    20,    0,  200);
    theHistoManager->addHisto("MHT",                                "finalSel",     "ttH_2lss_ee_cf",   "",    20,    0,  200);
    theHistoManager->addHisto("MetLD",                              "finalSel",     "ttH_2lss_ee_cf",   "",    15, -0.2,  1.4);
    theHistoManager->addHisto("TauMultiplicity",                    "finalSel",     "ttH_2lss_ee_cf",   "",    10,   -1,    9);
    theHistoManager->addHisto("JetMultiplicity",                    "finalSel",     "ttH_2lss_ee_cf",   "",    10,   -1,    9);
    theHistoManager->addHisto("LooseBJetMultiplicity",              "finalSel",     "ttH_2lss_ee_cf",   "",    10,   -1,    9);
    theHistoManager->addHisto("MediumBJetMultiplicity",             "finalSel",     "ttH_2lss_ee_cf",   "",    10,   -1,    9);
    theHistoManager->addHisto("NumberOfSelectedLeptons",            "finalSel",     "ttH_2lss_ee_cf",   "",    10,   -1,    9);

    // Signal Region Three Leptons

    theHistoManager->addHisto("CutFlow",                            "finalSel",     "ttH_3l",   "",   10,   -1,    9);

    theHistoManager->addHisto("LeadingLeptonPt",                    "finalSel",     "ttH_3l",   "",   20,    0,  200);
    theHistoManager->addHisto("SubLeadingLeptonPt",                 "finalSel",     "ttH_3l",   "",   20,    0,  200);
    theHistoManager->addHisto("ThirdLeptonPt",                      "finalSel",     "ttH_3l",   "",   20,    0,  200);
    theHistoManager->addHisto("MET",                                "finalSel",     "ttH_3l",   "",   20,    0,  200);
    theHistoManager->addHisto("MHT",                                "finalSel",     "ttH_3l",   "",   20,    0,  200);
    theHistoManager->addHisto("MetLD",                              "finalSel",     "ttH_3l",   "",   15, -0.2,  1.4);
    theHistoManager->addHisto("TauMultiplicity",                    "finalSel",     "ttH_3l",   "",   10,   -1,    9);
    theHistoManager->addHisto("JetMultiplicity",                    "finalSel",     "ttH_3l",   "",   10,   -1,    9);
    theHistoManager->addHisto("LooseBJetMultiplicity",              "finalSel",     "ttH_3l",   "",   10,   -1,    9);
    theHistoManager->addHisto("MediumBJetMultiplicity",             "finalSel",     "ttH_3l",   "",   10,   -1,    9);
    theHistoManager->addHisto("SumOfLeptonsCharges",                "finalSel",     "ttH_3l",   "",   10,   -5,    5);
    theHistoManager->addHisto("SumOfThreeLeptonsCharges",           "finalSel",     "ttH_3l",   "",   10,   -5,    5);
    theHistoManager->addHisto("NumberOfSelectedLeptons",            "finalSel",     "ttH_3l",   "",   10,   -1,    9);

    // Signal Region Three Leptons - Fake Rate

    theHistoManager->addHisto("CutFlow",                            "finalSel",     "ttH_3l_fr",   "",   10,   -1,    9);

    theHistoManager->addHisto("LeadingLeptonPt",                    "finalSel",     "ttH_3l_fr",   "",   20,    0,  200);
    theHistoManager->addHisto("SubLeadingLeptonPt",                 "finalSel",     "ttH_3l_fr",   "",   20,    0,  200);
    theHistoManager->addHisto("ThirdLeptonPt",                      "finalSel",     "ttH_3l_fr",   "",   20,    0,  200);
    theHistoManager->addHisto("MET",                                "finalSel",     "ttH_3l_fr",   "",   20,    0,  200);
    theHistoManager->addHisto("MHT",                                "finalSel",     "ttH_3l_fr",   "",   20,    0,  200);
    theHistoManager->addHisto("MetLD",                              "finalSel",     "ttH_3l_fr",   "",   15, -0.2,  1.4);
    theHistoManager->addHisto("TauMultiplicity",                    "finalSel",     "ttH_3l_fr",   "",   10,   -1,    9);
    theHistoManager->addHisto("JetMultiplicity",                    "finalSel",     "ttH_3l_fr",   "",   10,   -1,    9);
    theHistoManager->addHisto("LooseBJetMultiplicity",              "finalSel",     "ttH_3l_fr",   "",   10,   -1,    9);
    theHistoManager->addHisto("MediumBJetMultiplicity",             "finalSel",     "ttH_3l_fr",   "",   10,   -1,    9);
    theHistoManager->addHisto("SumOfLeptonsCharges",                "finalSel",     "ttH_3l_fr",   "",   10,   -5,    5);
    theHistoManager->addHisto("SumOfThreeLeptonsCharges",           "finalSel",     "ttH_3l_fr",   "",   10,   -5,    5);
    theHistoManager->addHisto("NumberOfSelectedLeptons",            "finalSel",     "ttH_3l_fr",   "",   10,   -1,    9);

    // WZ variables

    theHistoManager->addHisto("CutFlow",                            "noSel",        "WZ_CR",   "",   10,   -1,    9);

    theHistoManager->addHisto("ZCandidateInvariantMass",            "finalSel",     "WZ_CR",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateTransverseMomentum",       "finalSel",     "WZ_CR",   "",  12,    0,   500);
    theHistoManager->addHisto("MET",                                "finalSel",     "WZ_CR",   "",  25,    0,   200);
    theHistoManager->addHisto("MHT",                                "finalSel",     "WZ_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("MetLD",                              "finalSel",     "WZ_CR",   "",  15, -0.2,   1.4);
    theHistoManager->addHisto("TauMultiplicity",                    "finalSel",     "WZ_CR",   "",   8,    0,     8);
    theHistoManager->addHisto("JetMultiplicity",                    "finalSel",     "WZ_CR",   "",   8,    0,     8);

    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR",        "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_1002",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_1003",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_1004",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_1005",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_1006",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_1007",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_1008",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_1009",   "",  20,    0,   300);
    
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2001",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2002",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2003",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2004",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2005",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2006",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2007",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2008",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2009",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2010",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2011",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2012",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2013",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2014",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2015",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2016",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2017",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2018",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2019",   "",  20,    0,   300);
    
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2020",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2021",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2022",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2023",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2024",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2025",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2026",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2027",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2028",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2029",   "",  20,    0,   300);
   
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2030",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2031",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2032",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2033",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2034",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2035",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2036",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2037",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2038",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2039",   "",  20,    0,   300);
  
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2040",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2041",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2042",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2043",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2044",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2045",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2046",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2047",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2048",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2049",   "",  20,    0,   300);
    
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2050",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2051",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2052",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2053",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2054",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2055",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2056",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2057",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2058",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2059",   "",  20,    0,   300);
    
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2060",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2061",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2062",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2063",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2064",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2065",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2066",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2067",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2068",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2069",   "",  20,    0,   300);
    
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2070",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2071",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2072",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2073",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2074",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2075",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2076",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2077",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2078",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2079",   "",  20,    0,   300);
    
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2080",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2081",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2082",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2083",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2084",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2085",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2086",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2087",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2088",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2089",   "",  20,    0,   300);
    
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2090",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2091",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2092",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2093",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2094",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2095",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2096",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2097",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2098",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2099",   "",  20,    0,   300);
    
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2100",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2101",   "",  20,    0,   300);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_2102",   "",  20,    0,   300);
 
    theHistoManager->addHisto("InvariantMassOfSelectedLeptons",     "finalSel",     "WZ_CR",   "",  15,    0,   600);
    theHistoManager->addHisto("SumOfSelectedLeptonsCharges",        "finalSel",     "WZ_CR",   "",  10,   -5,     5);
    theHistoManager->addHisto("SumVecPtSelectedLeptons",            "finalSel",     "WZ_CR",   "",  12,    0,   500);

    theHistoManager->addHisto("InvMassRemainingLepton",             "finalSel",     "WZ_CR",   "",  15,    0,   300);
    theHistoManager->addHisto("SumOfRemainingLeptonsCharges",       "finalSel",     "WZ_CR",   "",  10,   -5,     5);
    theHistoManager->addHisto("SumVecPtRemainingLeptons",           "finalSel",     "WZ_CR",   "",  25,    0,   500);
   
    theHistoManager->addHisto("JetMultiplicity",                    "finalSel",     "WZ_CR_BSFUp",   "",   8,    0,     8);
    theHistoManager->addHisto("SumOfSelectedLeptonsCharges",        "finalSel",     "WZ_CR_BSFUp",   "",  10,   -5,     5);
    theHistoManager->addHisto("MET",                                "finalSel",     "WZ_CR_BSFUp",   "",  20,    0,   200);
    theHistoManager->addHisto("ZCandidateInvariantMass",            "finalSel",     "WZ_CR_BSFUp",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateTransverseMomentum",       "finalSel",     "WZ_CR_BSFUp",   "",  12,    0,   500);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_BSFUp",   "",  15,    0,   300);
    
    theHistoManager->addHisto("JetMultiplicity",                    "finalSel",     "WZ_CR_BSFDo",   "",   8,    0,     8);
    theHistoManager->addHisto("SumOfSelectedLeptonsCharges",        "finalSel",     "WZ_CR_BSFDo",   "",  10,   -5,     5);
    theHistoManager->addHisto("MET",                                "finalSel",     "WZ_CR_BSFDo",   "",  20,    0,   200);
    theHistoManager->addHisto("ZCandidateInvariantMass",            "finalSel",     "WZ_CR_BSFDo",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateTransverseMomentum",       "finalSel",     "WZ_CR_BSFDo",   "",  12,    0,   500);
    theHistoManager->addHisto("MTW",                                "finalSel",     "WZ_CR_BSFDo",   "",  15,    0,   300);
   
    // WZ variables - relaxed CR

    theHistoManager->addHisto("CutFlow",                            "noSel",        "WZrel_CR",   "",   10,   -1,    9);

    theHistoManager->addHisto("ZCandidateInvariantMass",            "finalSel",     "WZrel_CR",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateTransverseMomentum",       "finalSel",     "WZrel_CR",   "",  12,    0,   500);
    theHistoManager->addHisto("MET",                                "finalSel",     "WZrel_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("MHT",                                "finalSel",     "WZrel_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("MetLD",                              "finalSel",     "WZrel_CR",   "",  15, -0.2,   1.4);
    theHistoManager->addHisto("TauMultiplicity",                    "finalSel",     "WZrel_CR",   "",   8,    0,     8);
    theHistoManager->addHisto("JetMultiplicity",                    "finalSel",     "WZrel_CR",   "",   8,    0,     8);

    theHistoManager->addHisto("InvariantMassOfSelectedLeptons",     "finalSel",     "WZrel_CR",   "",  30,    0,   600);
    theHistoManager->addHisto("SumOfSelectedLeptonsCharges",        "finalSel",     "WZrel_CR",   "",  10,   -5,     5);
    theHistoManager->addHisto("SumVecPtSelectedLeptons",            "finalSel",     "WZrel_CR",   "",  25,    0,   500);

    theHistoManager->addHisto("InvMassRemainingLepton",             "finalSel",     "WZrel_CR",   "",  15,    0,   300);
    theHistoManager->addHisto("SumOfRemainingLeptonsCharges",       "finalSel",     "WZrel_CR",   "",  10,   -5,     5);
    theHistoManager->addHisto("SumVecPtRemainingLeptons",           "finalSel",     "WZrel_CR",   "",  25,    0,   500);

    // TTZ variables - relaxed CR

    theHistoManager->addHisto("CutFlow",                            "noSel",        "TTZ_CR",   "",   10,   -1,    9);

    theHistoManager->addHisto("LeadingLeptonPt",                    "finalSel",     "TTZ_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("SubleadingLeptonPt",                 "finalSel",     "TTZ_CR",   "",  20,    0,   200);

    theHistoManager->addHisto("ZCandidateInvariantMass",            "finalSel",     "TTZ_CR",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateTransverseMomentum",       "finalSel",     "TTZ_CR",   "",  12,    0,   500);
    theHistoManager->addHisto("MET",                                "finalSel",     "TTZ_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("MHT",                                "finalSel",     "TTZ_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("MetLD",                              "finalSel",     "TTZ_CR",   "",  15, -0.2,   1.4);
    theHistoManager->addHisto("TauMultiplicity",                    "finalSel",     "TTZ_CR",   "",   8,    0,     8);
    theHistoManager->addHisto("JetMultiplicity",                    "finalSel",     "TTZ_CR",   "",   8,    0,     8);

    theHistoManager->addHisto("InvariantMassOfSelectedLeptons",     "finalSel",     "TTZ_CR",   "",  30,    0,   600);
    theHistoManager->addHisto("SumOfSelectedLeptonsCharges",        "finalSel",     "TTZ_CR",   "",  10,   -5,     5);
    theHistoManager->addHisto("SumVecPtSelectedLeptons",            "finalSel",     "TTZ_CR",   "",  25,    0,   500);

    theHistoManager->addHisto("InvMassRemainingLepton",             "finalSel",     "TTZ_CR",   "",  15,    0,   300);
    theHistoManager->addHisto("SumOfRemainingLeptonsCharges",       "finalSel",     "TTZ_CR",   "",  10,   -5,     5);
    theHistoManager->addHisto("SumVecPtRemainingLeptons",           "finalSel",     "TTZ_CR",   "",  25,    0,   500);

    theHistoManager->addHisto("LeadingLeptonPt",                    "finalSel",     "TTZ4j_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("SubleadingLeptonPt",                 "finalSel",     "TTZ4j_CR",   "",  20,    0,   200);

    theHistoManager->addHisto("MET",                                "finalSel",     "TTZ4j_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("JetMultiplicity",                    "finalSel",     "TTZ4j_CR",   "",   8,    0,     8);
    theHistoManager->addHisto("ZCandidateInvariantMass",            "finalSel",     "TTZ4j_CR",   "",  15,   60,   120);

    // PU reweighting
    theHistoManager->addHisto("NumberOfPrimaryVertex",              "noSel",        "",   "", 100,    0,    99);

    // 2D histo

    theHistoManager->addHisto2D("SelectedLeptonsVsJets",            "noSel",        "ttH_3l",   "",    8,    0,    7,    8,    0,    7);
    theHistoManager->addHisto2D("SelectedLeptonsVsBJets",           "noSel",        "ttH_3l",   "",    8,    0,    7,    8,    0,    7);

    theHistoManager->addHisto2D("SelectedLeptonsVsJets",            "noSel",        "WZ_CR",   "",    8,    0,    7,    8,    0,    7);
    theHistoManager->addHisto2D("SelectedLeptonsVsBJets",           "noSel",        "WZ_CR",   "",    8,    0,    7,    8,    0,    7);

    theHistoManager->addHisto2D("InvMassLastLeptonVSZMass",         "finalSel",     "WZ_CR",   "",   15,    0,  300,   15,   60,  120);
    theHistoManager->addHisto2D("SumPtLepVSZMass",                  "finalSel",     "WZ_CR",   "",   25,    0,  500,   15,   60,  120);
    theHistoManager->addHisto2D("METLDVSZMass",                     "finalSel",     "WZ_CR",   "",   15, -0.2,  1.4,   15,   60,  120);

    // MVA
    theHistoManager->addHisto("Signal_2lss_TT_MVA",                 "finalSel",     "ttH_2lss",   "",  20,   -1,     1);
    theHistoManager->addHisto("Signal_2lss_TTV_MVA",                "finalSel",     "ttH_2lss",   "",  20,   -1,     1);
    theHistoManager->addHisto("Signal_3l_TT_MVA",                   "finalSel",     "ttH_3l",   "",  20,   -1,     1);
    theHistoManager->addHisto("Signal_3l_TTV_MVA",                  "finalSel",     "ttH_3l",   "",  20,   -1,     1);

    // from Daniel 

    theHistoManager->addHisto("nLep",         "PreSel", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("nLep loose",   "PreSel", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("lep1Pt",       "PreSel", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep2Pt",       "PreSel", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep3Pt",       "PreSel", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1Eta",      "PreSel", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep2Eta",      "PreSel", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep3Eta",      "PreSel", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lepQ",         "PreSel", "", "", 9, -4.5, 4.5);
    theHistoManager->addHisto("Mll"     ,     "PreSel", "", "", 10, 0., 200.);
    theHistoManager->addHisto("Mll"     ,     "PreSel 3l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("Mll"     ,     "PreSel btag", "", "", 10, 0., 200.);
    theHistoManager->addHisto("Mll"     ,     "PreSel btag 3l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("nJets",        "PreSel", "", "", 8, -0.5, 7.5);
    theHistoManager->addHisto("nLooseB",      "PreSel", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("nMediumB",     "PreSel", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("JetPt",        "PreSel", "", "", 28, 20., 300.);
    theHistoManager->addHisto("JetEta",       "PreSel", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("CSVv2" ,       "PreSel", "", "", 102, -0.01, 1.01);
    theHistoManager->addHisto("METpx",        "PreSel", "", "",14,-140.,140.);
    theHistoManager->addHisto("METpy",        "PreSel", "", "",14,-140.,140.);
    theHistoManager->addHisto("MET"  ,        "PreSel", "", "",10,0.,400.);
    theHistoManager->addHisto("METphi",       "PreSel", "", "",16,-3.2,3.2);
    theHistoManager->addHisto("METsum",       "PreSel", "", "",15,0.,3000.);
    theHistoManager->addHisto("MHT",          "PreSel", "", "", 10, 0., 200.);
    theHistoManager->addHisto("MET LD",       "PreSel", "", "", 10, 0., 2.);

    // TTH3l
    theHistoManager->addHisto("nLep",         "TTH3l", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("nFake",        "TTH3l", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("lep1Pt",       "TTH3l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep2Pt",       "TTH3l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep3Pt",       "TTH3l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep4Pt",       "TTH3l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1Eta",      "TTH3l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep2Eta",      "TTH3l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep3Eta",      "TTH3l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep4Eta",      "TTH3l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lepQ",         "TTH3l", "", "", 9, -4.5, 4.5);
    theHistoManager->addHisto("nJets",        "TTH3l", "", "", 8, -0.5, 7.5);
    theHistoManager->addHisto("nLooseB",      "TTH3l", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("nMediumB",     "TTH3l", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("JetPt",        "TTH3l", "", "", 28, 20., 300.);
    theHistoManager->addHisto("JetEta",       "TTH3l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("METpx",        "TTH3l", "", "",14,-140.,140.);
    theHistoManager->addHisto("METpy",        "TTH3l", "", "",14,-140.,140.);
    theHistoManager->addHisto("MET"  ,        "TTH3l", "", "",10,0.,400.);
    theHistoManager->addHisto("METphi",       "TTH3l", "", "",16,-3.2,3.2);
    theHistoManager->addHisto("METsum",       "TTH3l", "", "",15,0.,3000.);
    theHistoManager->addHisto("MHT",          "TTH3l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("MET LD",       "TTH3l", "", "", 10, 0., 2.);
    theHistoManager->addHisto("lep1 Mt",      "TTH3l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1 dRmin",   "TTH3l", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep2 dRmin",   "TTH3l", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep Eta max",  "TTH3l", "", "", 12, 0., 2.4);
    theHistoManager->addHisto("jet dR av",    "TTH3l", "", "", 10, 0., 5.);
    theHistoManager->addHisto("LepId",        "TTH3l", "", "", 4, -0.5, 3.5);
    theHistoManager->addHisto("HT",           "TTH3l", "", "",10, 0., 1000.);
    theHistoManager->addHisto("Mll min",      "TTH3l", "", "",10, 0., 200.);
    theHistoManager->addHisto("best MZ",      "TTH3l", "", "",12,60., 120.);

    // WZ
    theHistoManager->addHisto("nLep",         "WZ", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("nFake",        "WZ", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("lep1Pt",       "WZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep2Pt",       "WZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep3Pt",       "WZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1Eta",      "WZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep2Eta",      "WZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep3Eta",      "WZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lepQ",         "WZ", "", "", 9, -4.5, 4.5);
    theHistoManager->addHisto("nJets",        "WZ", "", "", 8, -0.5, 7.5);
    theHistoManager->addHisto("nLooseB",      "WZ", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("nMediumB",     "WZ", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("JetPt",        "WZ", "", "", 28, 20., 300.);
    theHistoManager->addHisto("JetEta",       "WZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("METpx",        "WZ", "", "",14,-140.,140.);
    theHistoManager->addHisto("METpy",        "WZ", "", "",14,-140.,140.);
    theHistoManager->addHisto("MET"  ,        "WZ", "", "",10,0.,400.);
    theHistoManager->addHisto("METphi",       "WZ", "", "",16,-3.2,3.2);
    theHistoManager->addHisto("METsum",       "WZ", "", "",15,0.,3000.);
    theHistoManager->addHisto("MHT",          "WZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("MET LD",       "WZ", "", "", 10, 0., 2.);
    theHistoManager->addHisto("lep1 Mt",      "WZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1 dRmin",   "WZ", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep2 dRmin",   "WZ", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep Eta max",  "WZ", "", "", 12, 0., 2.4);
    theHistoManager->addHisto("jet dR av",    "WZ", "", "", 10, 0., 5.);
    theHistoManager->addHisto("W Mt",         "WZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("Z Pt",         "WZ", "", "", 20,  0., 500.);
    theHistoManager->addHisto("Pt Sum(l)",    "WZ", "", "", 20,  0., 500.);
    theHistoManager->addHisto("inv.mass(l)",  "WZ", "", "", 20,  0., 600.);
    theHistoManager->addHisto("LepId",        "WZ", "", "", 4, -0.5, 3.5);
    theHistoManager->addHisto("HT",           "WZ", "", "",10, 0., 1000.);
    theHistoManager->addHisto("Mll min",      "WZ", "", "",10, 0., 200.);
    theHistoManager->addHisto("best MZ",      "WZ", "", "",12,60., 120.);

    // WZrelaxed
    theHistoManager->addHisto("nLep",         "WZrelaxed", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("nFake",        "WZrelaxed", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("lep1Pt",       "WZrelaxed", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep2Pt",       "WZrelaxed", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep3Pt",       "WZrelaxed", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1Eta",      "WZrelaxed", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep2Eta",      "WZrelaxed", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep3Eta",      "WZrelaxed", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lepQ",         "WZrelaxed", "", "", 9, -4.5, 4.5);
    theHistoManager->addHisto("nJets",        "WZrelaxed", "", "", 8, -0.5, 7.5);
    theHistoManager->addHisto("nLooseB",      "WZrelaxed", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("nMediumB",     "WZrelaxed", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("JetPt",        "WZrelaxed", "", "", 28, 20., 300.);
    theHistoManager->addHisto("JetEta",       "WZrelaxed", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("METpx",        "WZrelaxed", "", "",14,-140.,140.);
    theHistoManager->addHisto("METpy",        "WZrelaxed", "", "",14,-140.,140.);
    theHistoManager->addHisto("MET"  ,        "WZrelaxed", "", "",10,0.,400.);
    theHistoManager->addHisto("METphi",       "WZrelaxed", "", "",16,-3.2,3.2);
    theHistoManager->addHisto("METsum",       "WZrelaxed", "", "",15,0.,3000.);
    theHistoManager->addHisto("MHT",          "WZrelaxed", "", "", 10, 0., 200.);
    theHistoManager->addHisto("MET LD",       "WZrelaxed", "", "", 10, 0., 2.);
    theHistoManager->addHisto("lep1 Mt",      "WZrelaxed", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1 dRmin",   "WZrelaxed", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep2 dRmin",   "WZrelaxed", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep Eta max",  "WZrelaxed", "", "", 12, 0., 2.4);
    theHistoManager->addHisto("jet dR av",    "WZrelaxed", "", "", 10, 0., 5.);
    theHistoManager->addHisto("W Mt",         "WZrelaxed", "", "", 10, 0., 200.);
    theHistoManager->addHisto("Z Pt",         "WZrelaxed", "", "", 20,  0., 500.);
    theHistoManager->addHisto("Pt Sum(l)",    "WZrelaxed", "", "", 20,  0., 500.);
    theHistoManager->addHisto("inv.mass(l)",  "WZrelaxed", "", "", 20,  0., 600.);
    theHistoManager->addHisto("LepId",        "WZrelaxed", "", "", 4, -0.5, 3.5);
    theHistoManager->addHisto("HT",           "WZrelaxed", "", "",10, 0., 1000.);
    theHistoManager->addHisto("Mll min",      "WZrelaxed", "", "",10, 0., 200.);
    theHistoManager->addHisto("best MZ",      "WZrelaxed", "", "",12,60., 120.);

    // TTZ
    theHistoManager->addHisto("nLep",         "TTZ", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("nFake",        "TTZ", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("lep1Pt",       "TTZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep2Pt",       "TTZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep3Pt",       "TTZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep4Pt",       "TTZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1Eta",      "TTZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep2Eta",      "TTZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep3Eta",      "TTZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep4Eta",      "TTZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lepQ",         "TTZ", "", "", 9, -4.5, 4.5);
    theHistoManager->addHisto("nJets",        "TTZ", "", "", 8, -0.5, 7.5);
    theHistoManager->addHisto("nLooseB",      "TTZ", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("nMediumB",     "TTZ", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("JetPt",        "TTZ", "", "", 28, 20., 300.);
    theHistoManager->addHisto("JetEta",       "TTZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("METpx",        "TTZ", "", "",14,-140.,140.);
    theHistoManager->addHisto("METpy",        "TTZ", "", "",14,-140.,140.);
    theHistoManager->addHisto("MET"  ,        "TTZ", "", "",10,0.,400.);
    theHistoManager->addHisto("METphi",       "TTZ", "", "",16,-3.2,3.2);
    theHistoManager->addHisto("METsum",       "TTZ", "", "",15,0.,3000.);
    theHistoManager->addHisto("MHT",          "TTZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("MET LD",       "TTZ", "", "", 10, 0., 2.);
    theHistoManager->addHisto("lep1 Mt",      "TTZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1 dRmin",   "TTZ", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep2 dRmin",   "TTZ", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep Eta max",  "TTZ", "", "", 12, 0., 2.4);
    theHistoManager->addHisto("jet dR av",    "TTZ", "", "", 10, 0., 5.);
    theHistoManager->addHisto("LepId",        "TTZ", "", "", 4, -0.5, 3.5);
    theHistoManager->addHisto("HT",           "TTZ", "", "",10, 0., 1000.);
    theHistoManager->addHisto("Mll min",      "TTZ", "", "",10, 0., 200.);
    theHistoManager->addHisto("best MZ",      "TTZ", "", "",12,60., 120.);

    // Zl
    theHistoManager->addHisto("nLep",         "Zl", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("nFake",        "Zl", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("lep1Pt",       "Zl", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep2Pt",       "Zl", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep3Pt",       "Zl", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1Eta",      "Zl", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep2Eta",      "Zl", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep3Eta",      "Zl", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lepQ",         "Zl", "", "", 9, -4.5, 4.5);
    theHistoManager->addHisto("nJets",        "Zl", "", "", 8, -0.5, 7.5);
    theHistoManager->addHisto("nLooseB",      "Zl", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("nMediumB",     "Zl", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("JetPt",        "Zl", "", "", 28, 20., 300.);
    theHistoManager->addHisto("JetEta",       "Zl", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("METpx",        "Zl", "", "",14,-140.,140.);
    theHistoManager->addHisto("METpy",        "Zl", "", "",14,-140.,140.);
    theHistoManager->addHisto("MET"  ,        "Zl", "", "",10,0.,400.);
    theHistoManager->addHisto("METphi",       "Zl", "", "",16,-3.2,3.2);
    theHistoManager->addHisto("METsum",       "Zl", "", "",15,0.,3000.);
    theHistoManager->addHisto("MHT",          "Zl", "", "", 10, 0., 200.);
    theHistoManager->addHisto("MET LD",       "Zl", "", "", 10, 0., 2.);
    theHistoManager->addHisto("lep1 Mt",      "Zl", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1 dRmin",   "Zl", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep2 dRmin",   "Zl", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep Eta max",  "Zl", "", "", 12, 0., 2.4);
    theHistoManager->addHisto("jet dR av",    "Zl", "", "", 10, 0., 5.);
    theHistoManager->addHisto("LepId",        "Zl", "", "", 4, -0.5, 3.5);
    theHistoManager->addHisto("HT",           "Zl", "", "",10, 0., 1000.);
    theHistoManager->addHisto("Mll min",      "Zl", "", "",10, 0., 200.);
    theHistoManager->addHisto("best MZ",      "Zl", "", "",12,60., 120.);

    // TTH2l
    theHistoManager->addHisto("nLep",         "TTH2l", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("nFake",        "TTH2l", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("lep1Pt",       "TTH2l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep2Pt",       "TTH2l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep3Pt",       "TTH2l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1Eta",      "TTH2l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep2Eta",      "TTH2l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep3Eta",      "TTH2l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lepQ",         "TTH2l", "", "", 9, -4.5, 4.5);
    theHistoManager->addHisto("nJets",        "TTH2l", "", "", 8, -0.5, 7.5);
    theHistoManager->addHisto("nLooseB",      "TTH2l", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("nMediumB",     "TTH2l", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("JetPt",        "TTH2l", "", "", 28, 20., 300.);
    theHistoManager->addHisto("JetEta",       "TTH2l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("METpx",        "TTH2l", "", "",14,-140.,140.);
    theHistoManager->addHisto("METpy",        "TTH2l", "", "",14,-140.,140.);
    theHistoManager->addHisto("MET"  ,        "TTH2l", "", "",10,0.,400.);
    theHistoManager->addHisto("METphi",       "TTH2l", "", "",16,-3.2,3.2);
    theHistoManager->addHisto("METsum",       "TTH2l", "", "",15,0.,3000.);
    theHistoManager->addHisto("MHT",          "TTH2l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("MET LD",       "TTH2l", "", "", 10, 0., 2.);
    theHistoManager->addHisto("lep1 Mt",      "TTH2l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1 dRmin",   "TTH2l", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep2 dRmin",   "TTH2l", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep Eta max",  "TTH2l", "", "", 12, 0., 2.4);
    theHistoManager->addHisto("jet dR av",    "TTH2l", "", "", 10, 0., 5.);
    theHistoManager->addHisto("LepId",        "TTH2l", "", "", 4, -0.5, 3.5);
    theHistoManager->addHisto("HT",           "TTH2l", "", "",10, 0., 1000.);
    theHistoManager->addHisto("Mll min",      "TTH2l", "", "",10, 0., 200.);
    theHistoManager->addHisto("best MZ",      "TTH2l", "", "",12,60., 120.);

    // TTdilep
    theHistoManager->addHisto("nLep",         "TTemu", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("nFake",        "TTemu", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("lep1Pt",       "TTemu", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep2Pt",       "TTemu", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep3Pt",       "TTemu", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1Eta",      "TTemu", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep2Eta",      "TTemu", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep3Eta",      "TTemu", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lepQ",         "TTemu", "", "", 9, -4.5, 4.5);
    theHistoManager->addHisto("nJets",        "TTemu", "", "", 8, -0.5, 7.5);
    theHistoManager->addHisto("nLooseB",      "TTemu", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("nMediumB",     "TTemu", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("JetPt",        "TTemu", "", "", 28, 20., 300.);
    theHistoManager->addHisto("JetEta",       "TTemu", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("METpx",        "TTemu", "", "",14,-140.,140.);
    theHistoManager->addHisto("METpy",        "TTemu", "", "",14,-140.,140.);
    theHistoManager->addHisto("MET"  ,        "TTemu", "", "",10,0.,400.);
    theHistoManager->addHisto("METphi",       "TTemu", "", "",16,-3.2,3.2);
    theHistoManager->addHisto("METsum",       "TTemu", "", "",15,0.,3000.);
    theHistoManager->addHisto("MHT",          "TTemu", "", "", 10, 0., 200.);
    theHistoManager->addHisto("MET LD",       "TTemu", "", "", 10, 0., 2.);
    theHistoManager->addHisto("lep1 Mt",      "TTemu", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1 dRmin",   "TTemu", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep2 dRmin",   "TTemu", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep Eta max",  "TTemu", "", "", 12, 0., 2.4);
    theHistoManager->addHisto("jet dR av",    "TTemu", "", "", 10, 0., 5.);
    theHistoManager->addHisto("LepId",        "TTemu", "", "", 4, -0.5, 3.5);
    theHistoManager->addHisto("HT",           "TTemu", "", "",10, 0., 1000.);
    theHistoManager->addHisto("Mll min",      "TTemu", "", "",10, 0., 200.);
    theHistoManager->addHisto("best MZ",      "TTemu", "", "",12,60., 120.);

    // "Unrolled" 2D MVAs for Combine   

    if (_doSystCombine)
    {
        histoManager_2lss_mm_0tau_bl_neg->addHisto(_process, "",6,0,7);
        histoManager_2lss_mm_0tau_bt_neg->addHisto(_process, "",6,0,7);
        histoManager_2lss_ee_0tau_bl_neg->addHisto(_process, "",6,0,7);
        histoManager_2lss_ee_0tau_bt_neg->addHisto(_process, "",6,0,7);
        histoManager_2lss_em_0tau_bl_neg->addHisto(_process, "",6,0,7);
        histoManager_2lss_em_0tau_bt_neg->addHisto(_process, "",6,0,7);

        histoManager_2lss_mm_0tau_bl_pos->addHisto(_process, "",6,0,7);
        histoManager_2lss_mm_0tau_bt_pos->addHisto(_process, "",6,0,7);
        histoManager_2lss_ee_0tau_bl_pos->addHisto(_process, "",6,0,7);
        histoManager_2lss_ee_0tau_bt_pos->addHisto(_process, "",6,0,7);
        histoManager_2lss_em_0tau_bl_pos->addHisto(_process, "",6,0,7);
        histoManager_2lss_em_0tau_bt_pos->addHisto(_process, "",6,0,7);

        histoManager_3l_bl_neg->addHisto(_process, "",6,0,7);
        histoManager_3l_bt_neg->addHisto(_process, "",6,0,7);
        histoManager_3l_bl_pos->addHisto(_process, "",6,0,7);
        histoManager_3l_bt_pos->addHisto(_process, "",6,0,7);

        histoManager_2lss_mm_0tau_bl_neg->addHisto(_process, "JES", 6,0,7);
        histoManager_2lss_mm_0tau_bt_neg->addHisto(_process, "JES", 6,0,7);
        histoManager_2lss_ee_0tau_bl_neg->addHisto(_process, "JES", 6,0,7);
        histoManager_2lss_ee_0tau_bt_neg->addHisto(_process, "JES", 6,0,7);
        histoManager_2lss_em_0tau_bl_neg->addHisto(_process, "JES", 6,0,7);
        histoManager_2lss_em_0tau_bt_neg->addHisto(_process, "JES", 6,0,7);

        histoManager_2lss_mm_0tau_bl_pos->addHisto(_process, "JES", 6,0,7);
        histoManager_2lss_mm_0tau_bt_pos->addHisto(_process, "JES", 6,0,7);
        histoManager_2lss_ee_0tau_bl_pos->addHisto(_process, "JES", 6,0,7);
        histoManager_2lss_ee_0tau_bt_pos->addHisto(_process, "JES", 6,0,7);
        histoManager_2lss_em_0tau_bl_pos->addHisto(_process, "JES", 6,0,7);
        histoManager_2lss_em_0tau_bt_pos->addHisto(_process, "JES", 6,0,7);

        histoManager_3l_bl_neg->addHisto(_process, "JES", 6,0,7);
        histoManager_3l_bt_neg->addHisto(_process, "JES", 6,0,7);
        histoManager_3l_bl_pos->addHisto(_process, "JES", 6,0,7);
        histoManager_3l_bt_pos->addHisto(_process, "JES", 6,0,7);
    }

    // Loading weight files and creating corrsponding histograms

    // b-tagging
    std::string inputFileHF = "/opt/sbg/scratch1/cms/TTH/weight/csv_rwt_fit_hf_76x_2016_02_08.root";
    std::string inputFileLF = "/opt/sbg/scratch1/cms/TTH/weight/csv_rwt_fit_lf_76x_2016_02_08.root";
    TFile* f_CSVwgt_HF = new TFile ((inputFileHF).c_str());
    TFile* f_CSVwgt_LF = new TFile ((inputFileLF).c_str());
    fillCSVhistos(f_CSVwgt_HF, f_CSVwgt_LF);

    // new b-tagging (using BTagCalibrationXStandaloneWhatever)
    // setup calibration + reader
    // BTagCalibrationX *      
    calib = BTagCalibrationX("csvv2", "/opt/sbg/scratch1/cms/TTH/weight/CSVv2_ichep.csv");

    // BTagCalibrationXReader *
    reader = BTagCalibrationXReader(  BTagEntryX::OP_LOOSE,  // operating point
                                                            "central",            // central sys type
                                                            {"up", "down"});      // other sys types

    reader.load(   calib,                // calibration instance
                    BTagEntryX::FLAV_B,    // btag flavour
                    "comb");               // measurement type

    reader.load(    calib,                // calibration instance
                    BTagEntryX::FLAV_C,    // btag flavour
                    "comb");              // measurement type

    reader.load(    calib,                // calibration instance
                    BTagEntryX::FLAV_UDSG, // btag flavour
                    "comb");              // measurement type


    // BTag Efficiencies
    std::string inputFileBTagEff = "/opt/sbg/scratch1/cms/TTH/weight/TT_TuneCUETP8M1_13TeV-powheg-pythia8_0.root";
    TFile* f_BTag_eff = new TFile ((inputFileBTagEff).c_str());
    fill_eff_btagging_histos(f_BTag_eff);

    // charge flip
    std::string inputFileQF = "/opt/sbg/scratch1/cms/TTH/weight/QF_data_el.root";
    TFile * f_QFwgt    = new TFile ((inputFileQF).c_str());
    fillQFhistos(f_QFwgt);

    // fake rate
    std::string inputFileFR = "/opt/sbg/scratch1/cms/TTH/weight/FR_data_ttH_mva.root";
    TFile * f_FRwgt    = new TFile ((inputFileFR).c_str());
    fillFRhistos(f_FRwgt);

    // for synchronization
    stat_2lss_SR_ee     = 0;    stat_2lss_SR_eetau  = 0;    stat_2lss_lepMVA_SB_ee  = 0;    stat_2lss_os_SB_ee  = 0;
    stat_2lss_SR_em     = 0;    stat_2lss_SR_emtau  = 0;    stat_2lss_lepMVA_SB_em  = 0;    stat_2lss_os_SB_em  = 0;
    stat_2lss_SR_mm     = 0;    stat_2lss_SR_mmtau  = 0;    stat_2lss_lepMVA_SB_mm  = 0;    stat_2lss_os_SB_mm  = 0;
    stat_2lss_SR_tau    = 0;    stat_2lss_fr_tau    = 0;
    stat_3l_SR          = 0;    stat_3l_lepMVA_SB   = 0;

    // for categorization
    cat_ee      = 0;   cat_em      = 0;   cat_mm      = 0; cat_2ltau   = 0;   cat_3l        = 0;
    cat_ee_fake = 0;   cat_em_fake = 0;   cat_mm_fake = 0;                    cat_3l_fake   = 0;
    cat_ee_flip = 0;   cat_em_flip = 0;
    cat_HtoWW   = 0;   cat_HtoZZ   = 0;   cat_Htott   = 0;
}


void TTbarHiggsMultileptonAnalysis::writeHistograms()
{  
    outputfile->cd();

    std::vector<TH1F*> the1DHisto =  theHistoManager->getHisto1D_list();
    std::vector<TH2F*> the2DHisto =  theHistoManager->getHisto2D_list();

    for(unsigned int i=0; i<the1DHisto.size(); i++)  the1DHisto[i]->Write();
    for(unsigned int i=0; i<the2DHisto.size(); i++)  the2DHisto[i]->Write();

    tOutput->Write();

    //COMBINE
    if (_doSystCombine)
    { 
        //
        outputfile->mkdir("ttH_2lss_mm_0tau_bl_neg");
        outputfile->cd("ttH_2lss_mm_0tau_bl_neg");
        std::vector<TH1F*> h_2lss_mm_0tau_bl_neg = histoManager_2lss_mm_0tau_bl_neg->getHisto1D_list();
        for(unsigned int i=0; i<h_2lss_mm_0tau_bl_neg.size(); i++)  h_2lss_mm_0tau_bl_neg[i]->Write();

        outputfile->cd();
        outputfile->mkdir("ttH_2lss_mm_0tau_bt_neg");
        outputfile->cd("ttH_2lss_mm_0tau_bt_neg");
        std::vector<TH1F*> h_2lss_mm_0tau_bt_neg = histoManager_2lss_mm_0tau_bt_neg->getHisto1D_list();
        for(unsigned int i=0; i<h_2lss_mm_0tau_bt_neg.size(); i++)  h_2lss_mm_0tau_bt_neg[i]->Write();

        outputfile->mkdir("ttH_2lss_ee_0tau_bl_neg");
        outputfile->cd("ttH_2lss_ee_0tau_bl_neg");
        std::vector<TH1F*> h_2lss_ee_0tau_bl_neg = histoManager_2lss_ee_0tau_bl_neg->getHisto1D_list();
        for(unsigned int i=0; i<h_2lss_ee_0tau_bl_neg.size(); i++)  h_2lss_ee_0tau_bl_neg[i]->Write();

        outputfile->cd();
        outputfile->mkdir("ttH_2lss_ee_0tau_bt_neg");
        outputfile->cd("ttH_2lss_ee_0tau_bt_neg");
        std::vector<TH1F*> h_2lss_ee_0tau_bt_neg = histoManager_2lss_ee_0tau_bt_neg->getHisto1D_list();
        for(unsigned int i=0; i<h_2lss_ee_0tau_bt_neg.size(); i++)  h_2lss_ee_0tau_bt_neg[i]->Write();

        outputfile->mkdir("ttH_2lss_em_0tau_bl_neg");
        outputfile->cd("ttH_2lss_em_0tau_bl_neg");
        std::vector<TH1F*> h_2lss_em_0tau_bl_neg = histoManager_2lss_em_0tau_bl_neg->getHisto1D_list();
        for(unsigned int i=0; i<h_2lss_em_0tau_bl_neg.size(); i++)  h_2lss_em_0tau_bl_neg[i]->Write();

        outputfile->cd();
        outputfile->mkdir("ttH_2lss_em_0tau_bt_neg");
        outputfile->cd("ttH_2lss_em_0tau_bt_neg");
        std::vector<TH1F*> h_2lss_em_0tau_bt_neg = histoManager_2lss_em_0tau_bt_neg->getHisto1D_list();
        for(unsigned int i=0; i<h_2lss_em_0tau_bt_neg.size(); i++)  h_2lss_em_0tau_bt_neg[i]->Write();

        //
        outputfile->mkdir("ttH_2lss_mm_0tau_bl_pos");
        outputfile->cd("ttH_2lss_mm_0tau_bl_pos");
        std::vector<TH1F*> h_2lss_mm_0tau_bl_pos = histoManager_2lss_mm_0tau_bl_pos->getHisto1D_list();
        for(unsigned int i=0; i<h_2lss_mm_0tau_bl_pos.size(); i++)  h_2lss_mm_0tau_bl_pos[i]->Write();

        outputfile->cd();
        outputfile->mkdir("ttH_2lss_mm_0tau_bt_pos");
        outputfile->cd("ttH_2lss_mm_0tau_bt_pos");
        std::vector<TH1F*> h_2lss_mm_0tau_bt_pos = histoManager_2lss_mm_0tau_bt_pos->getHisto1D_list();
        for(unsigned int i=0; i<h_2lss_mm_0tau_bt_pos.size(); i++)  h_2lss_mm_0tau_bt_pos[i]->Write();

        outputfile->mkdir("ttH_2lss_ee_0tau_bl_pos");
        outputfile->cd("ttH_2lss_ee_0tau_bl_pos");
        std::vector<TH1F*> h_2lss_ee_0tau_bl_pos = histoManager_2lss_ee_0tau_bl_pos->getHisto1D_list();
        for(unsigned int i=0; i<h_2lss_ee_0tau_bl_pos.size(); i++)  h_2lss_ee_0tau_bl_pos[i]->Write();

        outputfile->cd();
        outputfile->mkdir("ttH_2lss_ee_0tau_bt_pos");
        outputfile->cd("ttH_2lss_ee_0tau_bt_pos");
        std::vector<TH1F*> h_2lss_ee_0tau_bt_pos = histoManager_2lss_ee_0tau_bt_pos->getHisto1D_list();
        for(unsigned int i=0; i<h_2lss_ee_0tau_bt_pos.size(); i++)  h_2lss_ee_0tau_bt_pos[i]->Write();

        outputfile->mkdir("ttH_2lss_em_0tau_bl_pos");
        outputfile->cd("ttH_2lss_em_0tau_bl_pos");
        std::vector<TH1F*> h_2lss_em_0tau_bl_pos = histoManager_2lss_em_0tau_bl_pos->getHisto1D_list();
        for(unsigned int i=0; i<h_2lss_em_0tau_bl_pos.size(); i++)  h_2lss_em_0tau_bl_pos[i]->Write();

        outputfile->cd();
        outputfile->mkdir("ttH_2lss_em_0tau_bt_pos");
        outputfile->cd("ttH_2lss_em_0tau_bt_pos");
        std::vector<TH1F*> h_2lss_em_0tau_bt_pos = histoManager_2lss_em_0tau_bt_pos->getHisto1D_list();
        for(unsigned int i=0; i<h_2lss_em_0tau_bt_pos.size(); i++)  h_2lss_em_0tau_bt_pos[i]->Write();

        //    
        outputfile->cd();
        outputfile->mkdir("ttH_3l_bl_neg");
        outputfile->cd("ttH_3l_bl_neg");
        std::vector<TH1F*> h_3l_bl_neg = histoManager_3l_bl_neg->getHisto1D_list();
        for(unsigned int i=0; i<h_3l_bl_neg.size(); i++)  h_3l_bl_neg[i]->Write();

        outputfile->cd();
        outputfile->mkdir("ttH_3l_bt_neg");
        outputfile->cd("ttH_3l_bt_neg");
        std::vector<TH1F*> h_3l_bt_neg = histoManager_3l_bt_neg->getHisto1D_list();
        for(unsigned int i=0; i<h_3l_bt_neg.size(); i++)  h_3l_bt_neg[i]->Write();

        outputfile->cd();
        outputfile->mkdir("ttH_3l_bl_pos");
        outputfile->cd("ttH_3l_bl_pos");
        std::vector<TH1F*> h_3l_bl_pos = histoManager_3l_bl_pos->getHisto1D_list();
        for(unsigned int i=0; i<h_3l_bl_pos.size(); i++)  h_3l_bl_pos[i]->Write();

        outputfile->cd();
        outputfile->mkdir("ttH_3l_bt_pos");
        outputfile->cd("ttH_3l_bt_pos");
        std::vector<TH1F*> h_3l_bt_pos = histoManager_3l_bt_pos->getHisto1D_list();
        for(unsigned int i=0; i<h_3l_bt_pos.size(); i++)  h_3l_bt_pos[i]->Write();
    }

    outputfile->Close();
}


void TTbarHiggsMultileptonAnalysis::Init(TChain *tree)
{
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;

    vEvent    = new std::vector<Event>();
    vElectron = new std::vector<Electron>();
    vMuon     = new std::vector<Muon>();
    vTau      = new std::vector<Tau>();
    vJet      = new std::vector<Jet>();
    vTruth    = new std::vector<Truth>();

    fChain->SetBranchAddress("Event",    &vEvent   );
    fChain->SetBranchAddress("Electron", &vElectron);
    fChain->SetBranchAddress("Muon",     &vMuon    );
    fChain->SetBranchAddress("Tau",      &vTau     );
    fChain->SetBranchAddress("Jet",      &vJet     );
    fChain->SetBranchAddress("Truth",    &vTruth   );

    Load_MVA();
}


void TTbarHiggsMultileptonAnalysis::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();
    int nentries_max = nentries;
    if ( _nmax != -1 && _nmax < nentries ) nentries_max = _nmax;

    std::cout << "Number of input events = " << nentries << std::endl;
    std::cout << "Number of processed events = " << nentries_max << std::endl;

    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries_max;jentry++) 
    {

        Long64_t ientry = fChain->LoadTree(jentry);
        if (ientry < 0) break;

        if(jentry%10000 == 0) std::cout << "number of processed events " << jentry << std::endl;

        //if(jentry > 100000) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        //
        int pvn = vEvent->at(0).pv_n();
        theHistoManager->fillHisto("NumberOfPrimaryVertex", "noSel", "", "",  pvn, 1);

        mc_event = vEvent->at(0).id();

        weights_pdf.clear();
        ids_pdf.clear();

        if ( !_isdata )
        {
            weight = _lumi*_xsec/_nowe;
            mc_weight = vEvent->at(0).mc_weight();
            //weight_PV = _h_PV->GetBinContent(pvn);
            weight = weight * mc_weight; //*weight_PV;

            weight_scale_muF0p5 = vEvent->at(0).weight_scale_muF0p5();
            weight_scale_muF2 = vEvent->at(0).weight_scale_muF2();
            weight_scale_muR0p5 = vEvent->at(0).weight_scale_muR0p5();
            weight_scale_muR2 = vEvent->at(0).weight_scale_muR2();

            weights_pdf = vEvent->at(0).pdf_weights();
            ids_pdf = vEvent->at(0).pdf_ids();
	    
	    if (weights_pdf.size()<110) std::cout << "problem with size of PDF weights " << weights_pdf.size()<< std::endl;

        /*
        for(unsigned int i=0; i < weights_pdf.size() ; i++)
        {
            std::cout << weights_pdf.at(i) << " " <<ids_pdf.at(i) << std::endl;
        }   
        */

        }
        else 
        {
            weight    = 1.;
            mc_weight = 1.;
            weight_PV = 1.; 

            //
            bool TRIGm   = vEvent->at(0).is_TRIGm()  ;
            bool TRIGe   = vEvent->at(0).is_TRIGe()  ;
            bool TRIGmTk = vEvent->at(0).is_TRIGmTk(); 
            bool TRIGee  = vEvent->at(0).is_TRIGee() ;
            bool TRIGmm  = vEvent->at(0).is_TRIGmm() ; 
            bool TRIGme  = vEvent->at(0).is_TRIGme() ;
            bool TRIGem  = vEvent->at(0).is_TRIGem() ;
            bool TRIGmmTk= vEvent->at(0).is_TRIGmmTk();
            bool TRIGeee = vEvent->at(0).is_TRIGeee();
            bool TRIGmme = vEvent->at(0).is_TRIGmme();
            bool TRIGeem = vEvent->at(0).is_TRIGeem();
            bool TRIGmmm = vEvent->at(0).is_TRIGmmm();

            bool E = false, M = false, EE = false, MM = false, EM = false;
            if ( TRIGme || TRIGem || TRIGeem || TRIGmme ) EM = true;
            if ( TRIGmm || TRIGmmTk || TRIGmmm )          MM = true;
            if ( TRIGee || TRIGeee )	                  EE = true;
            if ( TRIGm  || TRIGmTk )                      M  = true;
            if ( TRIGe  )                                 E  = true;

            bool emdataset = _sampleName.Contains("MuonE");
            bool mmdataset = _sampleName.Contains("DoubleM");
            bool eedataset = _sampleName.Contains("DoubleE");
            bool mdataset  = _sampleName.Contains("SingleM");
            bool edataset  = _sampleName.Contains("SingleE");

            is_trigger = false;
            if ( EM  &&                               (emdataset) ) is_trigger = true;
            if ( !EM && MM  &&                        (mmdataset) ) is_trigger = true;
            if ( !EM && !MM && EE  &&                 (eedataset) ) is_trigger = true;
            if ( !EM && !MM && !EE && M  &&           (mdataset ) ) is_trigger = true;
            if ( !EM && !MM && !EE && !M && E &&      (edataset ) ) is_trigger = true;


            if(is_trigger)
            {
                weight *= 1;
                //std::cout << " EM: " << EM 
                //          << " MM: " << MM 
                //          << " EE: " << EE
                //          << " M:  " << M
                //          << " E:  " << E                    << std::endl;
                //std::cout << " sampleName: " << _sampleName   << std::endl;
                //std::cout << " weight:     " << weight       << std::endl; 
            }
            else
            {
                weight *= 0;
                //std::cout << " EM: " << EM 
                //          << " MM: " << MM  
                //          << " EE: " << EE
                //          << " M:  " << M
                //          << " E:  " << E                    << std::endl;
                //std::cout << " sampleName: " << _sampleName   << std::endl;
                //std::cout << " weight:     " << weight       << std::endl;  
            }
        }


        // ###########################################################
        // #  _       _ _   _       _ _           _   _              #
        // # (_)_ __ (_) |_(_) __ _| (_)___  __ _| |_(_) ___  _ __   #
        // # | | '_ \| | __| |/ _` | | / __|/ _` | __| |/ _ \| '_ \  #
        // # | | | | | | |_| | (_| | | \__ \ (_| | |_| | (_) | | | | #
        // # |_|_| |_|_|\__|_|\__,_|_|_|___/\__,_|\__|_|\___/|_| |_| #
        // #                                                         #
        // ###########################################################
        
        vLeptons.clear();
        vSelectedMuons.clear();
        vSelectedElectrons.clear();
        vSelectedTaus.clear();
        vSelectedLeptons.clear();
        vFakeMuons.clear();
        vFakeElectrons.clear();
        vInclusiveFakeLeptons.clear();
        vFakeLeptons.clear();	
        vSelectedNonBTagJets.clear();
        vSelectedBTagJets.clear();
        vSelectedMediumBTagJets.clear();
        vSelectedJets.clear();
        vSelectedTaus.clear();
        
        weightfake = 0;
        weightflip = 0;

        is_2lss_TTH_SR      = false;
        is_2lss_AppFakes_SR = false;
        is_2lss_AppFlips_SR = false;
        is_2lss_JM_SB       = false;
        is_2lss_LepMVA_SB   = false;
        is_emu_TT_CR        = false;

        is_3l_TTH_SR        = false;
        is_3l_AppFakes_SR   = false;
        is_3l_WZ_CR         = false; 
        is_3l_WZrel_CR      = false;
        is_3l_TTZ_CR        = false;

        is_2bTight          = 0;

        cat_ee_tau      = 0.;   cat_em_tau      = 0.;   cat_mm_tau      = 0.;
        cat_ee_2lss_FR  = 0.;   cat_em_2lss_FR  = 0.;   cat_mm_2lss_FR  = 0.;
        cat_ee_2lss_QF  = 0.;

        cat_ee          = 0;   cat_em           = 0;   cat_mm           = 0;    cat_2ltau   = 0;    cat_3l          = 0;
        cat_ee_fake     = 0;   cat_em_fake      = 0;   cat_mm_fake      = 0;                        cat_3l_fake     = 0;
        cat_ee_flip     = 0;   cat_em_flip      = 0;

        cat_HtoWW       = 0;   cat_HtoZZ       = 0;   cat_Htott     = 0;

        n_tight       = 0;

        // ######################################
        // #  _        _                        #
        // # | |_ _ __(_) __ _  __ _  ___ _ __  #
        // # | __| '__| |/ _` |/ _` |/ _ \ '__| #
        // # | |_| |  | | (_| | (_| |  __/ |    #
        // #  \__|_|  |_|\__, |\__, |\___|_|    #
        // #             |___/ |___/            #
        // #                                    #
        // ######################################

        //AC see above
        //is_trigger = false;
        //if ( vEvent->at(0).ev_trigger_pass_byname_1() >= 1 ) is_trigger = true;

        // #####################################
        // #  _ __ ___  _   _  ___  _ __  ___  #
        // # | '_ ` _ \| | | |/ _ \| '_ \/ __| #
        // # | | | | | | |_| | (_) | | | \__ \ #
        // # |_| |_| |_|\__,_|\___/|_| |_|___/ #
        // #                                   #
        // #####################################

        for(unsigned int imuon=0; imuon < vMuon->size() ; imuon++)
        {   
            Lepton l; l.setLepton(&vMuon->at(imuon),imuon,0,1);

            if ( vMuon->at(imuon).isFakeableTTH() )
            {
                vSelectedLeptons.push_back(l);
                if(vMuon->at(imuon).isTightTTH() ) n_tight = n_tight + 1;
            }


            //if (vMuon->at(imuon).isLooseTTH()) 
                vLeptons.push_back(l);

            theHistoManager->fillHisto("MuonPt",        "noSel",    "noChannel",   "",  vMuon->at(imuon).pt(),      weight);
            theHistoManager->fillHisto("MuonEta",       "noSel",    "noChannel",   "",  vMuon->at(imuon).eta(),     weight);
            theHistoManager->fillHisto("MuonMVA",       "noSel",    "noChannel",   "",  vMuon->at(imuon).lepMVA(),  weight);
        }     

        // ##############################################
        // #       _           _                        #
        // #   ___| | ___  ___| |_ _ __ ___  _ __  ___  #
        // #  / _ \ |/ _ \/ __| __| '__/ _ \| '_ \/ __| #
        // # |  __/ |  __/ (__| |_| | | (_) | | | \__ \ #
        // #  \___|_|\___|\___|\__|_|  \___/|_| |_|___/ #
        // #                                            #
        // ##############################################

        for(unsigned int ielectron=0; ielectron < vElectron->size() ; ielectron++)
        {   
            Lepton l; l.setLepton(&vElectron->at(ielectron),ielectron,1,0);

            if ( vElectron->at(ielectron).isFakeableTTH() )
            {
                vSelectedLeptons.push_back(l);
                if(vElectron->at(ielectron).isTightTTH()) n_tight = n_tight + 1;
            }

            //if (vElectron->at(ielectron).isLooseTTH()) 
                vLeptons.push_back(l);

            theHistoManager->fillHisto("ElectronPt",    "noSel",    "noChannel",    "", vElectron->at(ielectron).pt(),      weight);
            theHistoManager->fillHisto("ElectronEta",   "noSel",    "noChannel",    "", vElectron->at(ielectron).eta(),     weight);
            theHistoManager->fillHisto("ElectronMVA",   "noSel",    "noChannel",    "", vElectron->at(ielectron).lepMVA(),  weight);
        }  

        // ########################
        // #  _                   #
        // # | |_ __ _ _   _ ___  #
        // # | __/ _` | | | / __| #
        // # | || (_| | |_| \__ \ #
        // #  \__\__,_|\__,_|___/ #
        // #                      #
        // ########################            

        for(unsigned int itau=0; itau < vTau->size() ; itau++)
        {
            Lepton l; l.setLepton(&vTau->at(itau),itau,0,0);

            if( vTau->at(itau).isTightTTH() )
            {    
                vSelectedTaus.push_back(vTau->at(itau));
            }

            theHistoManager->fillHisto("TauPt",     "noSel",    "noChannel",    "", vTau->at(itau).pt(),    weight);
            theHistoManager->fillHisto("TauEta",    "noSel",    "noChannel",    "", vTau->at(itau).eta(),   weight);
        }

        // #############################################
        // #                _           _              #
        // #   ___  _ __ __| | ___ _ __(_)_ __   __ _  #
        // #  / _ \| '__/ _` |/ _ \ '__| | '_ \ / _` | #
        // # | (_) | | | (_| |  __/ |  | | | | | (_| | #
        // #  \___/|_|  \__,_|\___|_|  |_|_| |_|\__, | #
        // #                                    |___/  #
        // #                                           #
        // #############################################

        std::sort(vLeptons.begin(), vLeptons.end(), SortingLeptonPt);
        std::sort(vSelectedLeptons.begin(), vSelectedLeptons.end(), SortingLeptonPt);

        if ( vLeptons.size() >= 2) {
            for (unsigned int ilep = 0; ilep<vLeptons.size()-1; ilep++) {
                if ( vLeptons.at(ilep).pt() < vLeptons.at(ilep+1).pt() ) {
                    std::cout << "Run Event " << vEvent->at(0).run() <<  " " << vEvent->at(0).id() 
                        << " nlep/tight/fake " << vLeptons.size() <<  "/" << vSelectedLeptons.size() <<  "/" <<  vFakeLeptons.size() << std::endl; 
                    std::cout << "  all pt[" << ilep << "]: " << vLeptons.at(ilep).pt() 
                        << " pt[" << ilep+1 << "]: " << vLeptons.at(ilep+1).pt() << std::endl;
                }
            }
        }

        if ( vSelectedLeptons.size() >= 2) {
            for (unsigned int ilep = 0; ilep<vSelectedLeptons.size()-1; ilep++) {
                if ( vSelectedLeptons.at(ilep).pt() < vSelectedLeptons.at(ilep+1).pt() ) {
                    std::cout << "Run Event " << vEvent->at(0).run() <<  " " << vEvent->at(0).id() 
                        << " nlep/tight/fake " << vLeptons.size() <<  "/" << vSelectedLeptons.size() <<  "/" <<  vFakeLeptons.size() << std::endl; 
                    std::cout << "tight pt[" << ilep << "]: " << vSelectedLeptons.at(ilep).pt() 
                        << " pt[" << ilep+1 << "]: " << vSelectedLeptons.at(ilep+1).pt() << std::endl;
                }
            }
        }

        // ################################
        // #                              #
        // #  _           _      _        #
        // # | |__       (_) ___| |_ ___  #
        // # | '_ \ _____| |/ _ \ __/ __| #
        // # | |_) |_____| |  __/ |_\__ \ #
        // # |_.__/     _/ |\___|\__|___/ #
        // #           |__/               #
        // #                              #
        // ################################

        nLooseBJets  = 0;
        nMediumBJets = 0;

        for(unsigned int ijet=0; ijet < vJet->size() ; ijet++)
        {

            // updated for 801
            if( vJet->at(ijet).CSVv2() > 0.46  ) nLooseBJets++;
            if( vJet->at(ijet).CSVv2() > 0.80  ) nMediumBJets++;

            if(vJet->at(ijet).CSVv2() >= 0.46 ) vSelectedBTagJets.push_back(vJet->at(ijet));
            else                                vSelectedNonBTagJets.push_back(vJet->at(ijet));
            if(vJet->at(ijet).CSVv2() >= 0.80 ) vSelectedMediumBTagJets.push_back(vJet->at(ijet));

            vSelectedJets.push_back(vJet->at(ijet));

            theHistoManager->fillHisto("JetPt",     "noSel",    "noChannel",    "", vJet->at(ijet).pt(),    weight);
            theHistoManager->fillHisto("JetEta",    "noSel",    "noChannel",    "", vJet->at(ijet).eta(),   weight);
            theHistoManager->fillHisto("JetCSVv2",  "noSel",    "noChannel",    "", vJet->at(ijet).CSVv2(), weight);
        }

        theHistoManager->fillHisto("MET",           "noSel",    "noChannel",    "", vEvent->at(0).metpt(), weight );
        
        theHistoManager->fillHisto("CutFlow",       "noSel",    "noChannel",    "", 1, 1);

        // ###################################
        // #  _____ ____  _   _ _____ _   _  #
        // # |_   _|  _ \| | | |_   _| | | | #
        // #   | | | |_) | | | | | | | |_| | #
        // #   | | |  _ <| |_| | | | |  _  | #
        // #   |_| |_| \_\\___/  |_| |_| |_| #
        // #                                 #
        // ###################################

        //std::cout << "Event..." << std::endl;
        //std::cout << "WW: " << cat_HtoWW << " ZZ: " << cat_HtoZZ << " tau tau: " << cat_Htott << std::endl;
            
        if ( !_isdata )
        {
            for(unsigned int itruth = 0; itruth < vTruth->at(0).mc_truth_label().size() ; itruth++)
            {

                //std::cout << vTruth->at(0).mc_truth_label().at(itruth) << std::endl;
                if( vTruth->at(0).mc_truth_label().at(itruth) == 12 )
                {   cat_HtoWW = true; 
                    //std::cout << "H to WW" << std::endl;
                }
                else if( vTruth->at(0).mc_truth_label().at(itruth) == 14 )
                {   cat_HtoZZ = true; 
                    //std::cout << "H to ZZ" << std::endl;
                }
                else if( vTruth->at(0).mc_truth_label().at(itruth) == 16 )
                {   cat_Htott = true; 
                    //std::cout << "H to tau tau" << std::endl;
                }

                //std::cout << " label: "    << vTruth->at(0).mc_truth_label().at(itruth)  ;
                //std::cout << " id: "       << vTruth->at(0).mc_truth_id().at(itruth)     ; //<< std::endl;
                //std::cout << " pt: "       << vTruth->at(0).mc_truth_pt().at(itruth)     ; //<< std::endl;
                //std::cout << " eta: "      << vTruth->at(0).mc_truth_eta().at(itruth)    ; //<< std::endl;
                //std::cout << " phi: "      << vTruth->at(0).mc_truth_phi().at(itruth)    ; //<< std::endl;
                //std::cout << " E:   "      << vTruth->at(0).mc_truth_E().at(itruth)      << std::endl;
            
            }
        }

        // ################################################################################
        // #  ____  ____    ____  ____ _____                   _       _     _            #
        // # |___ \|  _ \  | __ )|  _ \_   _| __   ____ _ _ __(_) __ _| |__ | | ___  ___  #
        // #   __) | | | | |  _ \| | | || |   \ \ / / _` | '__| |/ _` | '_ \| |/ _ \/ __| #
        // #  / __/| |_| | | |_) | |_| || |    \ V / (_| | |  | | (_| | |_) | |  __/\__ \ #
        // # |_____|____/  |____/|____/ |_|     \_/ \__,_|_|  |_|\__,_|_.__/|_|\___||___/ #
        // #                                                                              #
        // ################################################################################

        max_Lep_eta     = 0. ;
        nJet25_Recl     = 0  ;
        mindr_lep1_jet  = 0. ;
        mindr_lep2_jet  = 0. ;
        met             = 0. ;
        avg_dr_jet      = 0. ;
        MT_met_lep1     = 0. ;
        LepGood_conePt0 = 0. ;
        LepGood_conePt1 = 0. ;
        mhtJet25_Recl   = 0. ;  

        // ############################################
        // #           _           _   _              #
        // #  ___  ___| | ___  ___| |_(_) ___  _ __   #
        // # / __|/ _ \ |/ _ \/ __| __| |/ _ \| '_ \  #
        // # \__ \  __/ |  __/ (__| |_| | (_) | | | | #
        // # |___/\___|_|\___|\___|\__|_|\___/|_| |_| #
        // #                                          #
        // ############################################

        TwoLeptonsSameSignSelection_TTH2l(jentry);
        TwoLeptonsSameSignSelection_ApplicationFakes(jentry);
        TwoLeptonsSameSignSelection_ApplicationFlips(jentry);
        //TwoLeptonsSameSignSelection_LepMVA_sideband(jentry);
        //TwoLeptonsSameSignSelection_JetMultiplicity_sideband(jentry);
        //TwoLeptonsSameSignSelection_TTbar(jentry);

        ThreeLeptonSelection_TTH3l(jentry);
        ThreeLeptonSelection_ApplicationFakes(jentry);
        ThreeLeptonSelection_CR_WZ(jentry);
        ThreeLeptonSelection_CR_WZrelaxed(jentry);
        ThreeLeptonSelection_TTZ(jentry);
        //ThreeLeptonSelection_CR_Zl(jentry);

        //if ( is_2lss_TTH_SR || is_3l_TTH_SR )                                     fillOutputTree();
        //if ( is_2lss_TTH_SR || is_3l_TTH_SR || is_3l_TTZ_CR || is_3l_WZrel_CR )   fillOutputTree();

        // #####################################################################################################
        // #                             .__                        .__                __   .__                #
        // #   _________.__. ____   ____ |  |_________  ____   ____ |__|____________ _/  |_ |__| ____   ____   #
        // #  /  ___<   |  |/    \_/ ___\|  |  \_  __ \/  _ \ /    \|  \___   /\__  \\   __\|  |/  _ \ /    \  #
        // #  \___ \ \___  |   |  \  \___|   Y  \  | \(  <_> )   |  \  |/    /  / __ \|  |  |  (  <_> )   |  \ #
        // # /____  >/ ____|___|  /\___  >___|  /__|   \____/|___|  /__/_____ \(____  /__|  |__|\____/|___|  / #
        // #      \/ \/         \/     \/     \/                  \/         \/     \/                     \/  #
        // #                                                                                                   #
        // #####################################################################################################

        bool print_all                        = false;
        bool test_stat                        = false;
        bool produce_table_2lss_SR_ee         = false;
        bool produce_table_2lss_lepMVA_SB_ee  = false;
        bool produce_table_2lss_os_SB_ee      = false;
        bool produce_table_2lss_SR_em         = false;
        bool produce_table_2lss_lepMVA_SB_em  = false;
        bool produce_table_2lss_SR_mm         = false;
        bool produce_table_2lss_lepMVA_SB_mm  = false;
        bool produce_table_3l_SR              = false;
        bool produce_table_3l_lepMVA_SB       = false;

        if ( print_all )
        {
            std::cout << std::endl;
            std::cout << "==================================================================" << std::endl;
            std::cout << "====  Event  ==== " << std::endl;

            std::cout << " Event:     "         << vEvent->at(0).id()
                      << " n leptons: "         << vSelectedLeptons.size()
                      << " n tight:   "         << n_tight
                      << " n taus:    "         << vSelectedTaus.size()
                      << " n jets:    "         << vSelectedJets.size()
                      << " n b-jets loose: "    << vSelectedBTagJets.size()
                      << " n b-jets medium: "   << vSelectedMediumBTagJets.size() << std::endl;

            std::cout << "==== Category ==== " << std::endl;

            std::cout << "Passing [ee ss SR]:  " << (is_2lss_TTH_SR         && abs(cat_ee)          )
                      <<       "  [ee ss FR]:  " << (is_2lss_AppFakes_SR    && abs(cat_ee_2lss_FR)  )
                      <<       "  [ee os QF]:  " << (is_2lss_AppFlips_SR    && abs(cat_ee_2lss_QF)  )
                      <<       "  [em ss SR]:  " << (is_2lss_TTH_SR         && abs(cat_em)          )
                      <<       "  [em ss FR]:  " << (is_2lss_AppFakes_SR    && abs(cat_em_2lss_FR)  )
                      <<       "  [mm ss SR]:  " << (is_2lss_TTH_SR         && abs(cat_mm)          )
                      <<       "  [mm ss FR]:  " << (is_2lss_AppFakes_SR    && abs(cat_mm_2lss_FR)  )
                      <<       "  [3l    SR]:  " << (is_3l_TTH_SR                                   )
                      <<       "  [3l    FR]:  " << (is_3l_AppFakes_SR                              )                     
                      << std::endl;

            std::cout << "==== Leptons ==== " << std::endl;

            for(int i=0; i<vSelectedLeptons.size() ; i++)
            {
                std::cout << "Lepton[" << i << "] : pt "        << vSelectedLeptons.at(i).pt()
                                            << "  eta  "        << vSelectedLeptons.at(i).eta()
                                            << "  phi  "        << vSelectedLeptons.at(i).phi()
                                            << "  isE  "        << vSelectedLeptons.at(i).isElectron()
                                            << "  isM  "        << vSelectedLeptons.at(i).isMuon()
                                            << "  istight  "    << vSelectedLeptons.at(i).isTightTTH()
                                            << "  lepMVA   "    << vSelectedLeptons.at(i).lepMVA_TTH()     
                                            << std::endl;
            }

            std::cout << "==== Leptons ==== " << std::endl;

            for(int i=0; i<vLeptons.size() ; i++)
            {
                std::cout << "Lepton[" << i << "] : pt "        << vLeptons.at(i).pt()
                                            << "  eta  "        << vLeptons.at(i).eta()
                                            << "  phi  "        << vLeptons.at(i).phi()
                                            << "  isE  "        << vLeptons.at(i).isElectron()
                                            << "  isM  "        << vLeptons.at(i).isMuon()
                                            << "  istight  "    << vLeptons.at(i).isTightTTH()
                                            << "  lepMVA   "    << vLeptons.at(i).lepMVA_TTH()
                                            << std::endl;
            }

            //std::cout << "====  Taus   ==== " << std::endl;

            std::cout << "====  Jets   ==== " << std::endl;
            for(int i=0; i<vSelectedJets.size(); i++)
            {
                std::cout << "Jet[" << i << "] : pt "           << vSelectedJets.at(i).pt()
                                         << "  eta  "           << vSelectedJets.at(i).eta()
                                         << "  phi  "           << vSelectedJets.at(i).phi()
                                         << "  CSV  "           << vSelectedJets.at(i).CSVv2()
                                         << std::endl;
            }

        }

        if( test_stat )
        {
            std::cout << " ======================== "                               << std::endl;
            std::cout << " stat_2lss_SR_ee        = " << stat_2lss_SR_ee            << std::endl;
            std::cout << " stat_2lss_lepMVA_SB_ee = " << stat_2lss_lepMVA_SB_ee     << std::endl;
            std::cout << " stat_2lss_os_SB_ee     = " << stat_2lss_os_SB_ee         << std::endl;
            std::cout << " ======================== "                               << std::endl;
            std::cout << " stat_2lss_SR_em        = " << stat_2lss_SR_em            << std::endl;
            std::cout << " stat_2lss_lepMVA_SB_em = " << stat_2lss_lepMVA_SB_em     << std::endl;
            std::cout << " ======================== "                               << std::endl;
            std::cout << " stat_2lss_SR_mm        = " << stat_2lss_SR_mm            << std::endl;
            std::cout << " stat_2lss_lepMVA_SB_mm = " << stat_2lss_lepMVA_SB_mm     << std::endl;
            std::cout << " ======================== "                               << std::endl;
            std::cout << " stat_2lss_SR_tau       = " << stat_2lss_SR_tau           << std::endl;
            std::cout << " stat_2lss_fr_tau       = " << stat_2lss_fr_tau           << std::endl;
            std::cout << " ======================== "                               << std::endl;
            std::cout << " stat_3l_SR             = " << stat_3l_SR                 << std::endl;
            std::cout << " stat_3l_lepMVA_SB      = " << stat_3l_lepMVA_SB          << std::endl;
        }

        if(is_2lss_TTH_SR && cat_ee && produce_table_2lss_SR_ee)
        {
            std::cout << vEvent->at(0).id()     << " "
                      << 1.                     << " "
                      << max_Lep_eta            << " "
                      << nJet25_Recl            << " "
                      << mindr_lep1_jet         << " "
                      << mindr_lep2_jet         << " "
                      << met                    << " "
                      << avg_dr_jet             << " "
                      << MT_met_lep1            << " "
                      << LepGood_conePt0        << " "
                      << LepGood_conePt1        << " "
                      << signal_2lss_TT_MVA     << " "
                      << signal_2lss_TTV_MVA    << " "
                      << std::endl;
        }

        if(is_2lss_AppFakes_SR && cat_ee_2lss_FR && produce_table_2lss_lepMVA_SB_ee)
        {
            std::cout << vEvent->at(0).id()     << " "
                      << weight_FR_2lss         << " "
                      << max_Lep_eta            << " "
                      << nJet25_Recl            << " "
                      << mindr_lep1_jet         << " "
                      << mindr_lep2_jet         << " "
                      << met                    << " "
                      << avg_dr_jet             << " "
                      << MT_met_lep1            << " "
                      << LepGood_conePt0        << " "
                      << LepGood_conePt1        << " "
                      << signal_2lss_TT_MVA     << " "
                      << signal_2lss_TTV_MVA    << " "
                      << std::endl;
        }

        if(is_2lss_AppFlips_SR && cat_ee_2lss_QF && produce_table_2lss_os_SB_ee)
        {
            std::cout << vEvent->at(0).id()     << " "
                      << weight_QF_ee           << " "
                      << max_Lep_eta            << " "
                      << nJet25_Recl            << " "
                      << mindr_lep1_jet         << " "
                      << mindr_lep2_jet         << " "
                      << met                    << " "
                      << avg_dr_jet             << " "
                      << MT_met_lep1            << " "
                      << LepGood_conePt0        << " "
                      << LepGood_conePt1        << " "
                      << signal_2lss_TT_MVA     << " "
                      << signal_2lss_TTV_MVA    << " "
                      << std::endl;
        }

        if(is_2lss_TTH_SR && cat_em && produce_table_2lss_SR_em)
        {
            std::cout << vEvent->at(0).id()     << " "
                      << 1.                     << " "
                      << max_Lep_eta            << " "
                      << nJet25_Recl            << " "
                      << mindr_lep1_jet         << " "
                      << mindr_lep2_jet         << " "
                      << met                    << " "
                      << avg_dr_jet             << " "
                      << MT_met_lep1            << " "
                      << LepGood_conePt0        << " "
                      << LepGood_conePt1        << " "
                      << signal_2lss_TT_MVA     << " "
                      << signal_2lss_TTV_MVA    << " "
                      << std::endl;
        }

        if(is_2lss_AppFakes_SR && cat_em_2lss_FR && produce_table_2lss_lepMVA_SB_em)
        {
            std::cout << vEvent->at(0).id()     << " "
                      << weight_FR_2lss         << " "
                      << max_Lep_eta            << " "
                      << nJet25_Recl            << " "
                      << mindr_lep1_jet         << " "
                      << mindr_lep2_jet         << " "
                      << met                    << " "
                      << avg_dr_jet             << " "
                      << MT_met_lep1            << " "
                      << LepGood_conePt0        << " "
                      << LepGood_conePt1        << " "
                      << signal_2lss_TT_MVA     << " "
                      << signal_2lss_TTV_MVA    << " "
                      << std::endl;
        }

        if(is_2lss_TTH_SR && cat_mm && produce_table_2lss_SR_mm)
        {
            std::cout << vEvent->at(0).id()     << " "
                      << 1.                     << " "
                      << max_Lep_eta            << " "
                      << nJet25_Recl            << " "
                      << mindr_lep1_jet         << " "
                      << mindr_lep2_jet         << " "
                      << met                    << " "
                      << avg_dr_jet             << " "
                      << MT_met_lep1            << " "
                      << LepGood_conePt0        << " "
                      << LepGood_conePt1        << " "
                      << signal_2lss_TT_MVA     << " "
                      << signal_2lss_TTV_MVA    << " "
                      << std::endl;
        }

        if(is_2lss_AppFakes_SR && cat_mm_2lss_FR && produce_table_2lss_lepMVA_SB_mm)
        {
            std::cout << vEvent->at(0).id()     << " "
                      << weight_FR_2lss         << " "
                      << max_Lep_eta            << " "
                      << nJet25_Recl            << " "
                      << mindr_lep1_jet         << " "
                      << mindr_lep2_jet         << " "
                      << met                    << " "
                      << avg_dr_jet             << " "
                      << MT_met_lep1            << " "
                      << LepGood_conePt0        << " "
                      << LepGood_conePt1        << " "
                      << signal_2lss_TT_MVA     << " "
                      << signal_2lss_TTV_MVA    << " "
                      << std::endl;
        }

        if(is_3l_TTH_SR && produce_table_3l_SR)
        {
            std::cout << vEvent->at(0).id()     << " "
                      << 1.                     << " "
                      << max_Lep_eta            << " "
                      << nJet25_Recl            << " "
                      << mindr_lep1_jet         << " "
                      << mindr_lep2_jet         << " "
                      << met                    << " "
                      << avg_dr_jet             << " "
                      << MT_met_lep1            << " "
                      << LepGood_conePt0        << " "
                      << LepGood_conePt1        << " "
                      << signal_2lss_TT_MVA     << " "
                      << signal_2lss_TTV_MVA    << " "
                      << std::endl;
        }

        if(is_3l_AppFakes_SR && produce_table_3l_lepMVA_SB)
        {
            std::cout << vEvent->at(0).id()     << " "
                      << weight_FR_3l           << " "
                      << max_Lep_eta            << " "
                      << nJet25_Recl            << " "
                      << mindr_lep1_jet         << " "
                      << mindr_lep2_jet         << " "
                      << met                    << " "
                      << avg_dr_jet             << " "
                      << MT_met_lep1            << " "
                      << LepGood_conePt0        << " "
                      << LepGood_conePt1        << " "
                      << signal_2lss_TT_MVA     << " "
                      << signal_2lss_TTV_MVA    << " "
                      << std::endl;
        }

        // ########################################################################################################
        // #  __  __           _              _       _     _     _     _   _  ____ ___        _          __  __  #
        // # |  \/  | __ _  __| |_      _____(_) __ _| |__ | |_  | |   | | | |/ ___/ _ \   ___| |_ _   _ / _|/ _| #
        // # | |\/| |/ _` |/ _` \ \ /\ / / _ \ |/ _` | '_ \| __| | |   | |_| | |  | | | | / __| __| | | | |_| |_  #
        // # | |  | | (_| | (_| |\ V  V /  __/ | (_| | | | | |_  | |___|  _  | |__| |_| | \__ \ |_| |_| |  _|  _| #
        // # |_|  |_|\__,_|\__,_| \_/\_/ \___|_|\__, |_| |_|\__| |_____|_| |_|\____\___/  |___/\__|\__,_|_| |_|   #
        // #                                    |___/                                                             #
        // #                                                                                                      #
        // ########################################################################################################

        if ( !_isdata && _printLHCO_MC && ThreeLeptonSelection_TTH3l_MC()) PrintLHCOforMadweight_MC(jentry);

        // Common Selection:
        if ( !(vLeptons.size() >= 2
                    && vSelectedJets.size() >= 2) ) continue;

        float MET = vEvent->at(0).metpt();
        float METphi = vEvent->at(0).metphi();
        float METx = MET * TMath::Cos(METphi);
        float METy = MET * TMath::Sin(METphi);
        float METsum = vEvent->at(0).metsumet();
        int nlepsel = 0;
        float jet_px = 0, jet_py = 0, lep_px = 0, lep_py = 0, MHT = 0, met_ld = 0, jetht = 0;
        for (int i=0; i<vSelectedLeptons.size(); i++) {
            lep_px += vSelectedLeptons.at(i).p4().Px();
            lep_py += vSelectedLeptons.at(i).p4().Py();
            if ( vSelectedLeptons.at(i).pt() > 10. ) nlepsel++;
        }

        if ( nlepsel >= 2 && vSelectedLeptons.at(0).pt() < 20. ) nlepsel = -1.;
        TLorentzVector jetp4;
        for (int ijet=0; ijet < vSelectedJets.size(); ijet++) {
            jetp4.SetPtEtaPhiE(vSelectedJets.at(ijet).pt(), vSelectedJets.at(ijet).eta(), vSelectedJets.at(ijet).phi(), vSelectedJets.at(ijet).E());
            jet_px += jetp4.Px();
            jet_py += jetp4.Py();
            jetht += vSelectedJets.at(ijet).pt();
            theHistoManager->fillHisto("JetPt",  "Trig", "", "", vJet->at(ijet).pt(), weight);
        }


        MHT = sqrt( (jet_px+lep_px)*(jet_px+lep_px) + (jet_py+lep_py)*(jet_py+lep_py) );
        met_ld = 0.00397 * MET + 0.00265 * MHT;

        float Mllmin = 1000., Mllbest = 1000., Deltabest = 1000.;
        theHistoManager->fillHisto("nLep",   "PreSel", "", "", nlepsel, weight);
        theHistoManager->fillHisto("nLep loose", "PreSel", "", "", vLeptons.size(), weight);

        if ( nlepsel >= 2 ) {
            theHistoManager->fillHisto("lep1Pt", "PreSel", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep1Eta","PreSel", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Pt", "PreSel", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep2Eta","PreSel", "", "", vSelectedLeptons.at(1).eta(), weight);
            if ( nlepsel >= 3 ) theHistoManager->fillHisto("lep3Pt", "PreSel", "", "", vSelectedLeptons.at(2).pt(), weight);
            if ( nlepsel >= 3 ) theHistoManager->fillHisto("lep3Eta","PreSel", "", "", vSelectedLeptons.at(2).eta(), weight);
            float lepq = vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge();
            if ( nlepsel >= 3 ) lepq += vSelectedLeptons.at(2).charge();
            theHistoManager->fillHisto("lepQ", "PreSel", "", "", lepq, weight);
            for (int i=0; i<vSelectedLeptons.size()-1; i++) {
                for (int j=i+1; j<vSelectedLeptons.size(); j++) {
                    float mll = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M();
                    if ( mll < Mllmin ) Mllmin = mll;
                    if ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()) {
                        if ( fabs(mll - 91.188) < Deltabest ) {
                            Mllbest = mll;
                            Deltabest = fabs(mll - 91.188);
                        }
                        theHistoManager->fillHisto("Mll", "PreSel", "", "", mll, weight);
                        if ( nlepsel >= 3 ) theHistoManager->fillHisto("Mll", "PreSel 3l", "", "", mll, weight);
                        if ( nLooseBJets>=2 || nMediumBJets>=1 ) theHistoManager->fillHisto("Mll", "PreSel btag", "", "", mll, weight);
                        if ( (nLooseBJets>=2 || nMediumBJets>=1 ) && nlepsel >= 3 ) theHistoManager->fillHisto("Mll", "PreSel btag 3l", "", "", mll, weight);
                    }
                }
            }
        }

        theHistoManager->fillHisto("nJets",    "PreSel", "", "", vSelectedJets.size(), weight);
        theHistoManager->fillHisto("nLooseB",  "PreSel", "", "", nLooseBJets, weight);
        theHistoManager->fillHisto("nMediumB", "PreSel", "", "", nMediumBJets, weight);
        for (int ijet=0; ijet < vSelectedJets.size(); ijet++) {
            theHistoManager->fillHisto("JetPt",  "PreSel", "", "", vJet->at(ijet).pt(), weight);
            theHistoManager->fillHisto("JetEta", "PreSel", "", "", vJet->at(ijet).eta(), weight);
            theHistoManager->fillHisto("CSVv2",  "PreSel", "", "", vJet->at(ijet).CSVv2(), weight);
        }
        theHistoManager->fillHisto("METpx",  "PreSel", "", "", METx, weight);
        theHistoManager->fillHisto("METpy",  "PreSel", "", "", METy, weight);
        theHistoManager->fillHisto("MET"  ,  "PreSel", "", "", MET, weight);
        theHistoManager->fillHisto("METphi", "PreSel", "", "", METphi, weight);
        theHistoManager->fillHisto("METsum", "PreSel", "", "", METsum, weight);
        theHistoManager->fillHisto("MHT",    "PreSel", "", "", MHT, weight);
        theHistoManager->fillHisto("MET LD", "PreSel", "", "", met_ld, weight);

        int lepid2 = -1, lepid3 = -1;
        if ( nlepsel >= 2 ) {
            if ( abs(vSelectedLeptons.at(0).id()) == 11 
                    && abs(vSelectedLeptons.at(1).id()) == 11 ) lepid2 = 0; // ee
            if ( abs(vSelectedLeptons.at(0).id()) == 11 
                    && abs(vSelectedLeptons.at(1).id()) == 13 ) lepid2 = 1; // emu
            if ( abs(vSelectedLeptons.at(0).id()) == 13 
                    && abs(vSelectedLeptons.at(1).id()) == 11 ) lepid2 = 1;
            if ( abs(vSelectedLeptons.at(0).id()) == 13 
                    && abs(vSelectedLeptons.at(1).id()) == 13 ) lepid2 = 2; // mumu
        }
        if ( nlepsel >= 3 ) {
            if ( lepid2 == 0 && abs(vSelectedLeptons.at(2).id()) == 11 ) lepid3 = 0; // eee
            if ( lepid2 == 0 && abs(vSelectedLeptons.at(2).id()) == 13 ) lepid3 = 1; // eemu
            if ( lepid2 == 1 && abs(vSelectedLeptons.at(2).id()) == 11 ) lepid3 = 1;
            if ( lepid2 == 1 && abs(vSelectedLeptons.at(2).id()) == 13 ) lepid3 = 2; // emumu
            if ( lepid2 == 2 && abs(vSelectedLeptons.at(2).id()) == 11 ) lepid3 = 2;
            if ( lepid2 == 2 && abs(vSelectedLeptons.at(2).id()) == 13 ) lepid3 = 3; // mumumu
        }

        // ################################################################################
        // #  ____  ____    ____  ____ _____                   _       _     _            #
        // # |___ \|  _ \  | __ )|  _ \_   _| __   ____ _ _ __(_) __ _| |__ | | ___  ___  #
        // #   __) | | | | |  _ \| | | || |   \ \ / / _` | '__| |/ _` | '_ \| |/ _ \/ __| #
        // #  / __/| |_| | | |_) | |_| || |    \ V / (_| | |  | | (_| | |_) | |  __/\__ \ #
        // # |_____|____/  |____/|____/ |_|     \_/ \__,_|_|  |_|\__,_|_.__/|_|\___||___/ #
        // #                                                                              #
        // ################################################################################

        float dr;
        float lep1_mtw = -1., lep1_dr_min = 100., lep2_dr_min = 100., lep_eta_max = -1.;

        if (   vSelectedLeptons.size()     >= 2
                && vSelectedLeptons.at(0).pt() >  20 
                && vSelectedLeptons.at(1).pt() >  10 ) 
        {
            lep_eta_max = fabs(vSelectedLeptons.at(0).eta());

            if ( fabs(vSelectedLeptons.at(1).eta()) > lep_eta_max ) lep_eta_max = fabs(vSelectedLeptons.at(1).eta());

            for (int ijet=0; ijet < vJet->size() ; ijet++) 
            {
                dr = GetDeltaR( vSelectedLeptons.at(0).eta(), vSelectedLeptons.at(0).phi(), vJet->at(ijet).eta(), vJet->at(ijet).phi() );
                if ( dr < lep1_dr_min ) lep1_dr_min = dr;

                dr = GetDeltaR( vSelectedLeptons.at(1).eta(), vSelectedLeptons.at(1).phi(), vJet->at(ijet).eta(), vJet->at(ijet).phi() );
                if ( dr < lep2_dr_min ) lep2_dr_min = dr;
            }

            lep1_mtw = sqrt( 2 * vSelectedLeptons.at(0).p4().Pt() * MET * (1 - cos( vSelectedLeptons.at(0).phi() - METphi )));
        }	   

        int njj = 0; 
        float jet_dr_av = 0.;
        for (int ijet=0; ijet < vJet->size()-1 ; ijet++) 
        {
            for (int kjet=ijet+1; kjet < vJet->size() ; kjet++) 
            {
                jet_dr_av += GetDeltaR( vJet->at(ijet).eta(), vJet->at(ijet).phi(), vJet->at(kjet).eta(), vJet->at(kjet).phi() );
                njj++;
            }
        }
        if ( njj > 0 ) jet_dr_av = jet_dr_av / njj;

        // ########################################################
        // #  _     _     _                                       #
        // # | |__ (_)___| |_ ___   __ _ _ __ __ _ _ __ ___  ____ #
        // # | '_ \| / __| __/ _ \ / _` | '__/ _` | '_ ` _ \|_  / #
        // # | | | | \__ \ || (_) | (_| | | | (_| | | | | | |/ /  #
        // # |_| |_|_|___/\__\___/ \__, |_|  \__,_|_| |_| |_/___| #
        // #                       |___/                          #
        // #                                                      #
        // ########################################################

        // TTH3l
        if ( is_3l_TTH_SR ) {

            if ( _isdata ) {
                std::cout << "Run Event " << vEvent->at(0).run() <<  " " << vEvent->at(0).id() 
                    << " nlep/tight/fake " << vLeptons.size() <<  "/" << vSelectedLeptons.size() <<  "/" <<  vFakeLeptons.size() 
                    << " njet/loose/medium " << vSelectedJets.size() 
                    <<  "/" << nLooseBJets <<  "/" << nMediumBJets << std::endl;
                std::cout << "lep id " << lepid3 << " pT " << vSelectedLeptons.at(0).pt() <<  " " 
                    << vSelectedLeptons.at(1).pt() <<  " "
                    << vSelectedLeptons.at(2).pt() <<  " " << std::endl;
                std::cout << " " << std::endl;
            }

            theHistoManager->fillHisto("nLep",   "TTH3l", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "TTH3l", "", "", vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "TTH3l", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "TTH3l", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep3Pt", "TTH3l", "", "", vSelectedLeptons.at(2).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "TTH3l", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "TTH3l", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lep3Eta", "TTH3l", "", "", vSelectedLeptons.at(2).eta(), weight);
            if ( vSelectedLeptons.size() >= 4 ) {
                theHistoManager->fillHisto("lep4Pt",  "TTH3l", "", "", vSelectedLeptons.at(3).pt(), weight);
                theHistoManager->fillHisto("lep4Eta", "TTH3l", "", "", vSelectedLeptons.at(3).eta(), weight);
                theHistoManager->fillHisto("lepQ",    "TTH3l", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge()+vSelectedLeptons.at(2).charge()+vSelectedLeptons.at(3).charge(), weight);
            }
            else theHistoManager->fillHisto("lepQ", "TTH3l", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge()+vSelectedLeptons.at(2).charge(), weight);
            theHistoManager->fillHisto("nJets",     "TTH3l", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "TTH3l", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "TTH3l", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "TTH3l", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "TTH3l", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "TTH3l", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "TTH3l", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "TTH3l", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "TTH3l", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "TTH3l", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "TTH3l", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "TTH3l", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "TTH3l", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "TTH3l", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "TTH3l", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "TTH3l", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "TTH3l", "", "", jet_dr_av, weight);
            theHistoManager->fillHisto("LepId", "TTH3l", "", "", lepid3, weight);
            theHistoManager->fillHisto("HT",     "TTH3l", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","TTH3l", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","TTH3l", "", "", Mllbest, weight);
        }

        // WZ
        if ( is_3l_WZ_CR ) {
            theHistoManager->fillHisto("nLep",   "WZ", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "WZ", "", "", vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "WZ", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "WZ", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep3Pt", "WZ", "", "", vSelectedLeptons.at(2).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "WZ", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "WZ", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lep3Eta", "WZ", "", "", vSelectedLeptons.at(2).eta(), weight);
            theHistoManager->fillHisto("lepQ", "WZ", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge()+vSelectedLeptons.at(2).charge(), weight);
            int lepw = -1;
            float zpt = -1.;
            for (int i=0; i<vSelectedLeptons.size()-1; i++) {
                for (int j=i+1; j<vSelectedLeptons.size(); j++) {
                    if ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()) {
                        float mz = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M();
                        if ( fabs(mz - 91.188) < 10. ) {
                            if ( i == 0 && j == 1 ) lepw = 2;
                            if ( i == 0 && j == 2 ) lepw = 1;    
                            if ( i == 1 && j == 2 ) lepw = 0;
                            zpt = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).Pt();
                        }
                    }
                }
            }
            theHistoManager->fillHisto("nJets",     "WZ", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "WZ", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "WZ", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "WZ", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "WZ", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "WZ", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "WZ", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "WZ", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "WZ", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "WZ", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "WZ", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "WZ", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "WZ", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "WZ", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "WZ", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "WZ", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "WZ", "", "", jet_dr_av, weight);
            if ( lepw >= 0 ) {
                float mtw = sqrt( 2 * vSelectedLeptons.at(lepw).p4().Pt() * MET * (1 - cos( vSelectedLeptons.at(lepw).phi() - METphi )));
                theHistoManager->fillHisto("W Mt", "WZ", "", "", mtw, weight);
                theHistoManager->fillHisto("Z Pt", "WZ", "", "", zpt, weight);
            }
            TLorentzVector all_lep_invmass_p4;
            for (int i=0; i<vSelectedLeptons.size(); i++) {
                all_lep_invmass_p4 += vSelectedLeptons.at(i).p4();
            }
            float all_lep_invmass = all_lep_invmass_p4.M();
            float all_lep_sumofpt = sqrt( (lep_px*lep_px) + (lep_py*lep_py) );
            theHistoManager->fillHisto("inv.mass(l)", "WZ", "", "", all_lep_invmass, weight);
            theHistoManager->fillHisto("Pt Sum(l)",   "WZ", "", "", all_lep_sumofpt, weight);
            theHistoManager->fillHisto("LepId", "WZ", "", "", lepid3, weight);
            theHistoManager->fillHisto("HT",     "WZ", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","WZ", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","WZ", "", "", Mllbest, weight);
        }

        // WZrelaxed
        if ( is_3l_WZrel_CR ) {
            theHistoManager->fillHisto("nLep",   "WZrelaxed", "", "", vLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "WZrelaxed", "", "", vFakeLeptons.size()+vLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "WZrelaxed", "", "", vLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "WZrelaxed", "", "", vLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep3Pt", "WZrelaxed", "", "", vLeptons.at(2).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "WZrelaxed", "", "", vLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "WZrelaxed", "", "", vLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lep3Eta", "WZrelaxed", "", "", vLeptons.at(2).eta(), weight);
            theHistoManager->fillHisto("lepQ", "WZrelaxed", "", "", vLeptons.at(0).charge()+vLeptons.at(1).charge()+vLeptons.at(2).charge(), weight);
            int lepw = -1;
            float zpt = -1.;
            for (int i=0; i<vLeptons.size()-1; i++) {
                for (int j=i+1; j<vLeptons.size(); j++) {
                    if ( vLeptons.at(i).id() == -vLeptons.at(j).id()) {
                        float mz = ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M();
                        if ( fabs(mz - 91.188) < 10. ) {
                            if ( i == 0 && j == 1 ) lepw = 2;
                            if ( i == 0 && j == 2 ) lepw = 1;    
                            if ( i == 1 && j == 2 ) lepw = 0;
                            zpt = ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).Pt();
                        }
                    }
                }
            }
            theHistoManager->fillHisto("nJets",     "WZrelaxed", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "WZrelaxed", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "WZrelaxed", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "WZrelaxed", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "WZrelaxed", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "WZrelaxed", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "WZrelaxed", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "WZrelaxed", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "WZrelaxed", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "WZrelaxed", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "WZrelaxed", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "WZrelaxed", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "WZrelaxed", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "WZrelaxed", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "WZrelaxed", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "WZrelaxed", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "WZrelaxed", "", "", jet_dr_av, weight);
            if ( lepw >= 0 ) {
                float mtw = sqrt( 2 * vLeptons.at(lepw).p4().Pt() * MET * (1 - cos( vLeptons.at(lepw).phi() - METphi )));
                theHistoManager->fillHisto("W Mt", "WZrelaxed", "", "", mtw, weight);
                theHistoManager->fillHisto("Z Pt", "WZrelaxed", "", "", zpt, weight);
            }
            TLorentzVector all_lep_invmass_p4;
            for (int i=0; i<vLeptons.size(); i++) {
                all_lep_invmass_p4 += vLeptons.at(i).p4();
            }
            float all_lep_invmass = all_lep_invmass_p4.M();
            float all_lep_sumofpt = sqrt( (lep_px*lep_px) + (lep_py*lep_py) );
            theHistoManager->fillHisto("inv.mass(l)", "WZrelaxed", "", "", all_lep_invmass, weight);
            theHistoManager->fillHisto("Pt Sum(l)",   "WZrelaxed", "", "", all_lep_sumofpt, weight);
            theHistoManager->fillHisto("LepId", "WZrelaxed", "", "", lepid3, weight);
            theHistoManager->fillHisto("HT",     "WZrelaxed", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","WZrelaxed", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","WZrelaxed", "", "", Mllbest, weight);
        }

        // TTZ
        if ( is_3l_TTZ_CR ) {
            theHistoManager->fillHisto("nLep",    "TTZ", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",   "TTZ", "", "", vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt",  "TTZ", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt",  "TTZ", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep3Pt",  "TTZ", "", "", vSelectedLeptons.at(2).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "TTZ", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "TTZ", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lep3Eta", "TTZ", "", "", vSelectedLeptons.at(2).eta(), weight);
            if ( vSelectedLeptons.size() >= 4 ) {
                theHistoManager->fillHisto("lep4Pt",  "TTZ", "", "", vSelectedLeptons.at(3).pt(), weight);
                theHistoManager->fillHisto("lep4Eta", "TTZ", "", "", vSelectedLeptons.at(3).eta(), weight);
                theHistoManager->fillHisto("lepQ",    "TTZ", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge()+vSelectedLeptons.at(2).charge()+vSelectedLeptons.at(3).charge(), weight);
            }
            else theHistoManager->fillHisto("lepQ", "TTZ", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge()+vSelectedLeptons.at(2).charge(), weight);
            theHistoManager->fillHisto("nJets",     "TTZ", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "TTZ", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "TTZ", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "TTZ", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "TTZ", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "TTZ", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "TTZ", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "TTZ", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "TTZ", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "TTZ", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "TTZ", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "TTZ", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "TTZ", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "TTZ", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "TTZ", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "TTZ", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "TTZ", "", "", jet_dr_av, weight);
            theHistoManager->fillHisto("LepId", "TTZ", "", "", lepid3, weight);
            theHistoManager->fillHisto("HT",     "TTZ", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","TTZ", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","TTZ", "", "", Mllbest, weight);
        }

        // TTH2l
        if ( is_2lss_TTH_SR ) {
            theHistoManager->fillHisto("nLep",   "TTH2l", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "TTH2l","", "",  vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "TTH2l", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "TTH2l", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "TTH2l", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "TTH2l", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lepQ", "TTH2l", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge(), weight);
            theHistoManager->fillHisto("nJets",     "TTH2l", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "TTH2l", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "TTH2l", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "TTH2l", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "TTH2l", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "TTH2l", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "TTH2l", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "TTH2l", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "TTH2l", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "TTH2l", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "TTH2l", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "TTH2l", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "TTH2l", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "TTH2l", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "TTH2l", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "TTH2l", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "TTH2l", "", "", jet_dr_av, weight);
            theHistoManager->fillHisto("LepId", "TTH2l", "", "", lepid2, weight);
            theHistoManager->fillHisto("HT",     "TTH2l", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","TTH2l", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","TTH2l", "", "", Mllbest, weight);
        }

        // TTdilep
        if ( is_emu_TT_CR ) {
            theHistoManager->fillHisto("nLep",   "TTemu", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "TTemu","", "",  vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "TTemu", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "TTemu", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "TTemu", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "TTemu", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lepQ", "TTemu", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge(), weight);
            theHistoManager->fillHisto("nJets",     "TTemu", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "TTemu", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "TTemu", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "TTemu", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "TTemu", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "TTemu", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "TTemu", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "TTemu", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "TTemu", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "TTemu", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "TTemu", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "TTemu", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "TTemu", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "TTemu", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "TTemu", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "TTemu", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "TTemu", "", "", jet_dr_av, weight);
            theHistoManager->fillHisto("LepId", "TTemu", "", "", lepid2, weight);
            theHistoManager->fillHisto("HT",     "TTemu", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","TTemu", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","TTemu", "", "", Mllbest, weight);
        }

    }

}

// ##############################################################################
// #        _ _       _                   _                  _                  #
// #   __ _| | |  ___(_) __ _ _ __   __ _| |  _ __ ___  __ _(_) ___  _ __  ___  #
// #  / _` | | | / __| |/ _` | '_ \ / _` | | | '__/ _ \/ _` | |/ _ \| '_ \/ __| #
// # | (_| | | | \__ \ | (_| | | | | (_| | | | | |  __/ (_| | | (_) | | | \__ \ #
// #  \__,_|_|_| |___/_|\__, |_| |_|\__,_|_| |_|  \___|\__, |_|\___/|_| |_|___/ #
// #                    |___/                          |___/                    #
// #                                                                            #
// ##############################################################################

void TTbarHiggsMultileptonAnalysis::TwoLeptonsSameSignSelection_TTH2l(int evt)
{
    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "ttH_2lss",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "ttH_2lss",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    //std::cout << "weight: " << weight << "     taus: " << vSelectedTaus.size() << std::endl;

    if(weight==0)              return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    if(DEBUG)             std::cout << std::endl << " 2l ss SR ==================================================================" << std::endl;

    bool nTight             = ( n_tight == 2);
    if(!nTight)             return;

    if(DEBUG) std::cout << "nTight Ok... ";

    bool are_fullytight     =  ( vSelectedLeptons.at(0).isTightTTH()        && vSelectedLeptons.at(1).isTightTTH()
                              && vSelectedLeptons.at(0).passTightCharge()   && vSelectedLeptons.at(1).passTightCharge()
                              && vSelectedLeptons.at(0).cutEventSel()       && vSelectedLeptons.at(1).cutEventSel()
                              && vSelectedLeptons.at(0).noLostHits()        && vSelectedLeptons.at(1).noLostHits() );
    if(!are_fullytight)     return;
    if(DEBUG) std::cout << "are_tight Ok... ";

    bool leading_lep_pt     = (  ( (vSelectedLeptons.at(0).isElectron() ) && (vSelectedLeptons.at(0).pt() > 25) )
                              || ( (vSelectedLeptons.at(0).isMuon()     ) && (vSelectedLeptons.at(0).pt() > 20) ) );
    if(!leading_lep_pt)     return;

    if(DEBUG) std::cout << "leading_lep_pt Ok... ";

    bool following_lep_pt   = (  ( (vSelectedLeptons.at(1).isElectron() ) && (vSelectedLeptons.at(1).pt() > 15) )
                              || ( (vSelectedLeptons.at(1).isMuon()     ) && (vSelectedLeptons.at(1).pt() > 10) ) );
    if(!following_lep_pt)   return;

    if(DEBUG) std::cout << "following_lep_pt Ok... ";

    bool nJets              = ( vSelectedJets.size()        >= 4 );
    if(!nJets)              return;

    if(DEBUG) std::cout << "nJets Ok... ";

    bool nLooseBtag         = ( nLooseBJets                 >= 2 );
    bool nMediumBtag        = ( nMediumBJets                >= 1 );
    if(!nLooseBtag && !nMediumBtag)      return;

    if(DEBUG) std::cout << "Btag Ok... ";

    // Adding invariant mass cut on loose leptons pairs
    bool pass_invariantemasscut = true;
    for(int i=0; i<vLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLeptons.size(); j++)
        {
            if ( fabs( ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M() ) < 12 )
            { pass_invariantemasscut = false ;}
            //std::cout << std::endl << "lep[" << i << "] and lep[" << j << "]: " << fabs( ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M() ) << " " << std::endl; 
        }
    }
    if(!pass_invariantemasscut) return;

    if(DEBUG) std::cout << "invariantmasscut Ok... ";

    // ###############################
    // # Two leptons event selection #
    // ###############################

    bool same_sign      = ( vSelectedLeptons.at(0).charge() == vSelectedLeptons.at(1).charge() );
    if(!same_sign)      return;

    if(DEBUG) std::cout << "samesign Ok... ";

    // ##########
    // # Z veto #
    // ##########

    bool pass_Zveto = true;

    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).isElectron() &&  vSelectedLeptons.at(j).isElectron() )
               && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                "finalSel",   "ttH2lss",   "",   2, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "finalSel",   "ttH2lss", _sampleName.Data(),   3, weight);

    // ##################################################################################################################################

    // #################################
    // # b-tagging nominal reweighting #
    // #################################

    std::vector<double> jetPts;
    std::vector<double> jetEtas;
    std::vector<double> jetCSVs;
    std::vector<int>    jetFlavors;
    int iSys = 0;
    double wgt_csv, wgt_csv_def, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf, new_weight;

    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetPts.push_back(     vSelectedJets.at(i).pt()                );
        jetEtas.push_back(    vSelectedJets.at(i).eta()               );
        jetCSVs.push_back(    vSelectedJets.at(i).CSVv2()             );
        jetFlavors.push_back( vSelectedJets.at(i).jet_hadronFlavour() );
    }

    wgt_csv_def = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
    new_weight = weight * wgt_csv_def; // weight = weight * wgt_csv_def;

    // ##################################################################################################################################

    // ##################################
    // # b-tagging deriving systematics #
    // ##################################

    std::vector<double> weights_csv;
    double wgt_csv_def_sys = 0;

    for(int i=7; i<25; i++)
    {
        wgt_csv_def_sys = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, i, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf)/wgt_csv_def;
        weights_csv.push_back(wgt_csv_def_sys);
    }

    double min_weight_csv = *min_element(weights_csv.begin(),weights_csv.end());
    double max_weight_csv = *max_element(weights_csv.begin(),weights_csv.end());

    theHistoManager->fillHisto("WeightCSV_min", "finalSel", "ttH2l",    "",     min_weight_csv, 1);
    theHistoManager->fillHisto("WeightCSV_max", "finalSel", "ttH2l",    "",     max_weight_csv, 1);

    weight_csv_down = min_weight_csv;
    weight_csv_up   = max_weight_csv;    

    // ##################################################################################################################################

    if(  vSelectedLeptons.at(0).isElectron()
      && vSelectedLeptons.at(1).isElectron()
      && pass_Zveto
      && (met_ld > 0.2)                      )
    {
        if(DEBUG) std::cout << "ee + Zveto + metld Ok... SELECTED";
        
        if(vSelectedTaus.size() == 0)
        {
            cat_ee = vSelectedLeptons.at(0).charge();
            is_2lss_TTH_SR = true;
            stat_2lss_SR_ee = stat_2lss_SR_ee + 1;

            theHistoManager->fillHisto("CutFlow",                          "finalSel",   "ttH_2lss_ee",   "", 8                              , weight);

            theHistoManager->fillHisto("LeadingLeptonPt",                  "finalSel",   "ttH_2lss_ee",   "", vSelectedLeptons.at(0).pt()    , weight);
            theHistoManager->fillHisto("SubLeadingLeptonPt",               "finalSel",   "ttH_2lss_ee",   "", vSelectedLeptons.at(1).pt()    , weight);
            theHistoManager->fillHisto("MET",                              "finalSel",   "ttH_2lss_ee",   "", vEvent->at(0).metpt()          , weight);
            theHistoManager->fillHisto("MHT",                              "finalSel",   "ttH_2lss_ee",   "", MHT                            , weight);
            theHistoManager->fillHisto("MetLD",                            "finalSel",   "ttH_2lss_ee",   "", met_ld                         , weight);
            theHistoManager->fillHisto("TauMultiplicity",                  "finalSel",   "ttH_2lss_ee",   "", vSelectedTaus.size()           , weight);
            theHistoManager->fillHisto("JetMultiplicity",                  "finalSel",   "ttH_2lss_ee",   "", vSelectedJets.size()           , weight);
            theHistoManager->fillHisto("LooseBJetMultiplicity",            "finalSel",   "ttH_2lss_ee",   "", vSelectedBTagJets.size()       , weight);
            theHistoManager->fillHisto("MediumBJetMultiplicity",           "finalSel",   "ttH_2lss_ee",   "", vSelectedMediumBTagJets.size() , weight);
            theHistoManager->fillHisto("NumberOfSelectedLeptons",          "finalSel",   "ttH_2lss_ee",   "", vSelectedLeptons.size()        , weight);
        }

        if(vSelectedTaus.size() > 0)
        {
            cat_2ltau = vSelectedLeptons.at(0).charge();
            stat_2lss_SR_tau = stat_2lss_SR_tau + 1;
        }   
    }

    if( abs(vSelectedLeptons.at(0).id()) != abs(vSelectedLeptons.at(1).id()) )
    {
        if(DEBUG) std::cout << "em Ok... SELECTED";

        if(vSelectedTaus.size() == 0)
        {
            cat_em = vSelectedLeptons.at(0).charge();
            is_2lss_TTH_SR = true;
            stat_2lss_SR_em = stat_2lss_SR_em + 1;

            theHistoManager->fillHisto("CutFlow",                          "finalSel",   "ttH_2lss_em",   "", 8                              , weight);

            theHistoManager->fillHisto("LeadingLeptonPt",                  "finalSel",   "ttH_2lss_em",   "", vSelectedLeptons.at(0).pt()    , weight);
            theHistoManager->fillHisto("SubLeadingLeptonPt",               "finalSel",   "ttH_2lss_em",   "", vSelectedLeptons.at(1).pt()    , weight);
            theHistoManager->fillHisto("MET",                              "finalSel",   "ttH_2lss_em",   "", vEvent->at(0).metpt()          , weight);
            theHistoManager->fillHisto("MHT",                              "finalSel",   "ttH_2lss_em",   "", MHT                            , weight);
            theHistoManager->fillHisto("MetLD",                            "finalSel",   "ttH_2lss_em",   "", met_ld                         , weight);
            theHistoManager->fillHisto("TauMultiplicity",                  "finalSel",   "ttH_2lss_em",   "", vSelectedTaus.size()           , weight);
            theHistoManager->fillHisto("JetMultiplicity",                  "finalSel",   "ttH_2lss_em",   "", vSelectedJets.size()           , weight);
            theHistoManager->fillHisto("LooseBJetMultiplicity",            "finalSel",   "ttH_2lss_em",   "", vSelectedBTagJets.size()       , weight);
            theHistoManager->fillHisto("MediumBJetMultiplicity",           "finalSel",   "ttH_2lss_em",   "", vSelectedMediumBTagJets.size() , weight);
            theHistoManager->fillHisto("NumberOfSelectedLeptons",          "finalSel",   "ttH_2lss_em",   "", vSelectedLeptons.size()        , weight);
        }

        if(vSelectedTaus.size() > 0)
        {
            cat_2ltau = vSelectedLeptons.at(0).charge();
            stat_2lss_SR_tau = stat_2lss_SR_tau + 1;            
        }
    }

    if(  vSelectedLeptons.at(0).isMuon()
      && vSelectedLeptons.at(1).isMuon()  )
    {
        if(DEBUG) std::cout << "mm Ok... SELECTED";

        if(vSelectedTaus.size() == 0)
        {
            cat_mm = vSelectedLeptons.at(0).charge();
            is_2lss_TTH_SR = true;
            stat_2lss_SR_mm = stat_2lss_SR_mm + 1;

            theHistoManager->fillHisto("CutFlow",                          "finalSel",   "ttH_2lss_mm",   "", 8                              , weight);

            theHistoManager->fillHisto("LeadingLeptonPt",                  "finalSel",   "ttH_2lss_mm",   "", vSelectedLeptons.at(0).pt()    , weight);
            theHistoManager->fillHisto("SubLeadingLeptonPt",               "finalSel",   "ttH_2lss_mm",   "", vSelectedLeptons.at(1).pt()    , weight);
            theHistoManager->fillHisto("MET",                              "finalSel",   "ttH_2lss_mm",   "", vEvent->at(0).metpt()          , weight);
            theHistoManager->fillHisto("MHT",                              "finalSel",   "ttH_2lss_mm",   "", MHT                            , weight);
            theHistoManager->fillHisto("MetLD",                            "finalSel",   "ttH_2lss_mm",   "", met_ld                         , weight);
            theHistoManager->fillHisto("TauMultiplicity",                  "finalSel",   "ttH_2lss_mm",   "", vSelectedTaus.size()           , weight);
            theHistoManager->fillHisto("JetMultiplicity",                  "finalSel",   "ttH_2lss_mm",   "", vSelectedJets.size()           , weight);
            theHistoManager->fillHisto("LooseBJetMultiplicity",            "finalSel",   "ttH_2lss_mm",   "", vSelectedBTagJets.size()       , weight);
            theHistoManager->fillHisto("MediumBJetMultiplicity",           "finalSel",   "ttH_2lss_mm",   "", vSelectedMediumBTagJets.size() , weight);
            theHistoManager->fillHisto("NumberOfSelectedLeptons",          "finalSel",   "ttH_2lss_mm",   "", vSelectedLeptons.size()        , weight);
        }

        if(vSelectedTaus.size() > 0)
        {
            cat_2ltau = vSelectedLeptons.at(0).charge();
            is_2lss_TTH_SR = true;
            stat_2lss_SR_tau = stat_2lss_SR_tau + 1;
        }
    }

    // ####################################
    // #  ____  ____    ____  ____ _____  #
    // # |___ \|  _ \  | __ )|  _ \_   _| #
    // #   __) | | | | |  _ \| | | || |   #
    // #  / __/| |_| | | |_) | |_| || |   #
    // # |_____|____/  |____/|____/ |_|   #
    // #                                  #
    // ####################################

    // based on https://github.com/CERN-PH-CMG/cmgtools-lite/blob/0b47d4d1c50ea0e24ef0d9cf1c24c763e78c1bf0/TTHAnalysis/python/tools/kinMVA_2D_2lss_3l.py
    // and https://github.com/CERN-PH-CMG/cmgtools-lite/blob/0b47d4d1c50ea0e24ef0d9cf1c24c763e78c1bf0/TTHAnalysis/python/tools/eventVars_2lss.py

    // ======================================================================================================
    // variables against ttbar
    max_Lep_eta     = std::max( fabs(vSelectedLeptons.at(0).eta()), fabs(vSelectedLeptons.at(1).eta()) ) ;

    nJet25_Recl = vSelectedJets.size() ;     

    mindr_lep1_jet = 1000.;
    mindr_lep2_jet = 1000.;

    for (int i=0; i<vSelectedJets.size(); i++)
    {
        if( DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ) < mindr_lep1_jet )
        { mindr_lep1_jet = DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ); }
        if( DeltaRLeptonJet( vSelectedLeptons.at(1), vSelectedJets.at(i) ) < mindr_lep2_jet )
        { mindr_lep2_jet = DeltaRLeptonJet( vSelectedLeptons.at(1), vSelectedJets.at(i) ); }
    }

    float met_max   = 400;
    met             = std::min( vEvent->at(0).metpt(), met_max ) ;

    int njj     = 0 ;
    avg_dr_jet  = 0.;
    for (int ijet=0; ijet < vSelectedJets.size()-1 ; ijet++) 
    {
        for (int kjet=ijet+1; kjet < vSelectedJets.size() ; kjet++) 
        {
            avg_dr_jet += GetDeltaR( vSelectedJets.at(ijet).eta(), vSelectedJets.at(ijet).phi(), vSelectedJets.at(kjet).eta(), vSelectedJets.at(kjet).phi() );
            njj++;
        }
    }
    if ( njj > 0 ) avg_dr_jet = avg_dr_jet / njj;

    MT_met_lep1     = sqrt( 2. * vSelectedLeptons.at(0).pt() * vEvent->at(0).metpt() * (1 - cos( vSelectedLeptons.at(0).phi() - vEvent->at(0).metphi() )) );
    
    signal_2lss_TT_MVA  = mva_2lss_tt->EvaluateMVA("BDTG method");

    theHistoManager->fillHisto("Signal_2lss_TT_MVA",                       "finalSel", "ttH2lss",   "",  signal_2lss_TT_MVA,   weight);

    // ======================================================================================================
    // variables against ttV

    LepGood_conePt0 = vSelectedLeptons.at(0).pt() ;
    LepGood_conePt1 = vSelectedLeptons.at(1).pt() ;

    signal_2lss_TTV_MVA = mva_2lss_ttV->EvaluateMVA("BDTG method");

    theHistoManager->fillHisto("Signal_2lss_TTV_MVA",                      "finalSel", "ttH2lss",   "",  signal_2lss_TTV_MVA,  weight);

    //std::cout << " signal 2lss TT MVA: "  << signal_2lss_TT_MVA
    //    << " signal 2lss TTV MVA: " << signal_2lss_TTV_MVA << std::endl;

    // ======================================================================================================

    fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

void TTbarHiggsMultileptonAnalysis::TwoLeptonsSameSignSelection_ApplicationFakes(int evt)
{
    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "LepMVA_2l_SB",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "LepMVA_2l_SB",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    if(DEBUG)             std::cout << std::endl << " 2l ss FR ==================================================================" << std::endl;

    bool nLep               = ( vSelectedLeptons.size()                   >= 2 );
    if(!nLep)               return;
    if(DEBUG) std::cout << "nLep Ok... ";

    bool are_fullytight     =  ( vSelectedLeptons.at(0).isTightTTH()        && vSelectedLeptons.at(1).isTightTTH()
                              && vSelectedLeptons.at(0).passTightCharge()   && vSelectedLeptons.at(1).passTightCharge()
                              && vSelectedLeptons.at(0).cutEventSel()       && vSelectedLeptons.at(1).cutEventSel()
                              && vSelectedLeptons.at(0).noLostHits()        && vSelectedLeptons.at(1).noLostHits() );
    if(are_fullytight)     return;
    if(DEBUG) std::cout << "are_tight Ok... ";

    bool pass_nolosthits   =  ( vSelectedLeptons.at(0).noLostHits()        && vSelectedLeptons.at(1).noLostHits() );
    if(!pass_nolosthits)     return;
    if(DEBUG) std::cout << "pass_nolosthits Ok... ";

    bool leading_lep_pt     = (  ( (vSelectedLeptons.at(0).isElectron() ) && (vSelectedLeptons.at(0).pt() > 25) )
                              || ( (vSelectedLeptons.at(0).isMuon()     ) && (vSelectedLeptons.at(0).pt() > 20) ) );
    if(!leading_lep_pt)     return;
    if(DEBUG) std::cout << "leading_lep_pt Ok... ";

    bool following_lep_pt   = (  ( (vSelectedLeptons.at(1).isElectron() ) && (vSelectedLeptons.at(1).pt() > 15) )
                              || ( (vSelectedLeptons.at(1).isMuon()     ) && (vSelectedLeptons.at(1).pt() > 10) ) );
    if(!following_lep_pt)   return;
    if(DEBUG) std::cout << "following_lep_pt Ok... ";

    bool nJets              = ( vSelectedJets.size()        >= 4 );
    if(!nJets)              return;
    if(DEBUG) std::cout << "nJets Ok... ";

    bool nLooseBtag         = ( nLooseBJets                 >= 2 );
    bool nMediumBtag        = ( nMediumBJets                >= 1 );
    if(!nLooseBtag && !nMediumBtag)      return;
    if(DEBUG) std::cout << "Btag Ok... ";

    // Adding invariant mass cut on loose leptons pairs
    bool pass_invariantemasscut = true;
    for(int i=0; i<vLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLeptons.size(); j++)
        {
            if ( fabs( ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M() ) < 12 )
            { pass_invariantemasscut = false ;}
        }
    }
    if(!pass_invariantemasscut) return;
    if(DEBUG) std::cout << "invariantmasscut Ok... ";

    // ###############################
    // # Two leptons event selection #
    // ###############################

    bool same_sign      = ( vSelectedLeptons.at(0).charge() == vSelectedLeptons.at(1).charge() );
    if(!same_sign)      return;
    if(DEBUG) std::cout << "samesign Ok... ";

    // ##########
    // # Z veto #
    // ##########

    bool pass_Zveto = true;

    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).isElectron() &&  vSelectedLeptons.at(j).isElectron() )
               && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                "finalSel",   "LepMVA_2l_SB",   "",   2, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "finalSel",   "LepMVA_2l_SB", _sampleName.Data(),   3, weight);

    // ##################################################################################################################################

    // #########################
    // # fake rate reweighting #
    // #########################

    std::vector<double> leptonsPts;
    std::vector<double> leptonsEtas;
    std::vector<int>    leptonsIds;

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        if( !vSelectedLeptons.at(i).isTightTTH())
        {
            leptonsPts.push_back(     vSelectedLeptons.at(i).pt()                );
            leptonsEtas.push_back(    vSelectedLeptons.at(i).eta()               );
            leptonsIds.push_back(     vSelectedLeptons.at(i).id()                );
        }
    }

    weight_FR_2lss = get_FR_wgt_2l(leptonsPts, leptonsEtas, leptonsIds);

    // ##################################################################################################################################

    if(  vSelectedLeptons.at(0).isElectron()
      && vSelectedLeptons.at(1).isElectron()
      && pass_Zveto
      && (met_ld > 0.2)                      )
    {
        if(DEBUG) std::cout << "ee + Zveto + metld Ok... SELECTED";

        weightfake = get_FR_wgt_2l(leptonsPts, leptonsEtas, leptonsIds);
        //std::cout << "Fake 2lss weight: " << weightfake << std::endl; 

        if(vSelectedTaus.size() == 0)
        {
            stat_2lss_lepMVA_SB_ee = stat_2lss_lepMVA_SB_ee + 1;
            is_2lss_AppFakes_SR = true;
            cat_ee_2lss_FR = true;
            cat_ee_fake = vSelectedLeptons.at(0).charge();

            theHistoManager->fillHisto("CutFlow",                          "finalSel",   "LepMVA_2l_SB",   "", 8                              , weight);

            theHistoManager->fillHisto("LeadingLeptonPt",                  "finalSel",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
            theHistoManager->fillHisto("SubLeadingLeptonPt",               "finalSel",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
            theHistoManager->fillHisto("MET",                              "finalSel",   "LepMVA_2l_SB",   "", vEvent->at(0).metpt()          , weight);
            theHistoManager->fillHisto("MHT",                              "finalSel",   "LepMVA_2l_SB",   "", MHT                            , weight);
            theHistoManager->fillHisto("MetLD",                            "finalSel",   "LepMVA_2l_SB",   "", met_ld                         , weight);
            theHistoManager->fillHisto("TauMultiplicity",                  "finalSel",   "LepMVA_2l_SB",   "", vSelectedTaus.size()           , weight);
            theHistoManager->fillHisto("JetMultiplicity",                  "finalSel",   "LepMVA_2l_SB",   "", vSelectedJets.size()           , weight);
            theHistoManager->fillHisto("LooseBJetMultiplicity",            "finalSel",   "LepMVA_2l_SB",   "", vSelectedBTagJets.size()       , weight);
            theHistoManager->fillHisto("MediumBJetMultiplicity",           "finalSel",   "LepMVA_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
            theHistoManager->fillHisto("NumberOfSelectedLeptons",          "finalSel",   "LepMVA_2l_SB",   "", vSelectedLeptons.size()        , weight);
        }


        if(vSelectedTaus.size() > 0)
        {
            stat_2lss_fr_tau = stat_2lss_fr_tau + 1;
            cat_2ltau = true;
        }

    }

    if( abs(vSelectedLeptons.at(0).id()) != abs(vSelectedLeptons.at(1).id()) )
    {
        if(DEBUG) std::cout << "em Ok... SELECTED";

        weightfake = get_FR_wgt_2l(leptonsPts, leptonsEtas, leptonsIds);
        //std::cout << "Fake 2lss weight: " << weightfake << std::endl;

        if(vSelectedTaus.size() == 0)
        {
            stat_2lss_lepMVA_SB_em = stat_2lss_lepMVA_SB_em + 1;
            is_2lss_AppFakes_SR = true;
            cat_em_2lss_FR = true;
            cat_em_fake = vSelectedLeptons.at(0).charge();

            theHistoManager->fillHisto("CutFlow",                          "finalSel",   "LepMVA_2l_SB",   "", 8                              , weight);

            theHistoManager->fillHisto("LeadingLeptonPt",                  "finalSel",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
            theHistoManager->fillHisto("SubLeadingLeptonPt",               "finalSel",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
            theHistoManager->fillHisto("MET",                              "finalSel",   "LepMVA_2l_SB",   "", vEvent->at(0).metpt()          , weight);
            theHistoManager->fillHisto("MHT",                              "finalSel",   "LepMVA_2l_SB",   "", MHT                            , weight);
            theHistoManager->fillHisto("MetLD",                            "finalSel",   "LepMVA_2l_SB",   "", met_ld                         , weight);
            theHistoManager->fillHisto("TauMultiplicity",                  "finalSel",   "LepMVA_2l_SB",   "", vSelectedTaus.size()           , weight);
            theHistoManager->fillHisto("JetMultiplicity",                  "finalSel",   "LepMVA_2l_SB",   "", vSelectedJets.size()           , weight);
            theHistoManager->fillHisto("MediumBJetMultiplicity",           "finalSel",   "LepMVA_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
            theHistoManager->fillHisto("NumberOfSelectedLeptons",          "finalSel",   "LepMVA_2l_SB",   "", vSelectedLeptons.size()        , weight);
        }

        if(vSelectedTaus.size() > 0)
        {
            stat_2lss_fr_tau = stat_2lss_fr_tau + 1;
            cat_2ltau = true;
        }
    }

    if(  vSelectedLeptons.at(0).isMuon()
      && vSelectedLeptons.at(1).isMuon()  )
    {
        if(DEBUG) std::cout << "mm Ok... SELECTED";

        weightfake = get_FR_wgt_2l(leptonsPts, leptonsEtas, leptonsIds);
        //std::cout << "Fake 2lss weight: " << weightfake << std::endl;

        if(vSelectedTaus.size() == 0)
        {
            stat_2lss_lepMVA_SB_mm = stat_2lss_lepMVA_SB_mm + 1;
            is_2lss_AppFakes_SR = true;
            cat_mm_2lss_FR = true;
            cat_mm_fake = vSelectedLeptons.at(0).charge();

            theHistoManager->fillHisto("CutFlow",                          "finalSel",   "LepMVA_2l_SB",   "", 8                              , weight);

            theHistoManager->fillHisto("LeadingLeptonPt",                  "finalSel",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
            theHistoManager->fillHisto("SubLeadingLeptonPt",               "finalSel",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
            theHistoManager->fillHisto("MET",                              "finalSel",   "LepMVA_2l_SB",   "", vEvent->at(0).metpt()          , weight);
            theHistoManager->fillHisto("MHT",                              "finalSel",   "LepMVA_2l_SB",   "", MHT                            , weight);
            theHistoManager->fillHisto("MetLD",                            "finalSel",   "LepMVA_2l_SB",   "", met_ld                         , weight);
            theHistoManager->fillHisto("TauMultiplicity",                  "finalSel",   "LepMVA_2l_SB",   "", vSelectedTaus.size()           , weight);
            theHistoManager->fillHisto("JetMultiplicity",                  "finalSel",   "LepMVA_2l_SB",   "", vSelectedJets.size()           , weight);
            theHistoManager->fillHisto("LooseBJetMultiplicity",            "finalSel",   "LepMVA_2l_SB",   "", vSelectedBTagJets.size()       , weight);
            theHistoManager->fillHisto("MediumBJetMultiplicity",           "finalSel",   "LepMVA_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
            theHistoManager->fillHisto("NumberOfSelectedLeptons",          "finalSel",   "LepMVA_2l_SB",   "", vSelectedLeptons.size()        , weight);
        }

        if(vSelectedTaus.size() > 0)
        {
            stat_2lss_fr_tau = stat_2lss_fr_tau + 1;
            cat_2ltau = true;
        }
    }

    // ####################################
    // #  ____  ____    ____  ____ _____  #
    // # |___ \|  _ \  | __ )|  _ \_   _| #
    // #   __) | | | | |  _ \| | | || |   #
    // #  / __/| |_| | | |_) | |_| || |   #
    // # |_____|____/  |____/|____/ |_|   #
    // #                                  #
    // ####################################

    // based on https://github.com/CERN-PH-CMG/cmgtools-lite/blob/0b47d4d1c50ea0e24ef0d9cf1c24c763e78c1bf0/TTHAnalysis/python/tools/kinMVA_2D_2lss_3l.py
    // and https://github.com/CERN-PH-CMG/cmgtools-lite/blob/0b47d4d1c50ea0e24ef0d9cf1c24c763e78c1bf0/TTHAnalysis/python/tools/eventVars_2lss.py

    // ======================================================================================================
    // variables against ttbar
    max_Lep_eta     = std::max( fabs(vSelectedLeptons.at(0).eta()), fabs(vSelectedLeptons.at(1).eta()) ) ;

    nJet25_Recl = vSelectedJets.size() ;

    mindr_lep1_jet = 1000.;
    mindr_lep2_jet = 1000.;

    for (int i=0; i<vSelectedJets.size(); i++)
    {
        if( DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ) < mindr_lep1_jet )
        { mindr_lep1_jet = DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ); }
        if( DeltaRLeptonJet( vSelectedLeptons.at(1), vSelectedJets.at(i) ) < mindr_lep2_jet )
        { mindr_lep2_jet = DeltaRLeptonJet( vSelectedLeptons.at(1), vSelectedJets.at(i) ); }
    }

    float met_max   = 400;
    met             = std::min( vEvent->at(0).metpt(), met_max ) ;

    int njj = 0;
    avg_dr_jet = 0.;
    for (int ijet=0; ijet < vSelectedJets.size()-1 ; ijet++)
    {
        for (int kjet=ijet+1; kjet < vSelectedJets.size() ; kjet++)
        {
            avg_dr_jet += GetDeltaR( vSelectedJets.at(ijet).eta(), vSelectedJets.at(ijet).phi(), vSelectedJets.at(kjet).eta(), vSelectedJets.at(kjet).phi() );
            njj++;
        }
    }
    if ( njj > 0 ) avg_dr_jet = avg_dr_jet / njj;

    MT_met_lep1     = sqrt( 2. * vSelectedLeptons.at(0).pt() * vEvent->at(0).metpt() * (1 - cos( vSelectedLeptons.at(0).phi() - vEvent->at(0).metphi() )));

    signal_2lss_TT_MVA  = mva_2lss_tt->EvaluateMVA("BDTG method");

    theHistoManager->fillHisto("Signal_2lss_TT_MVA",                       "finalSel", "ttH2lss",   "",  signal_2lss_TT_MVA,   weight);

    // ======================================================================================================
    // variables against ttV

    LepGood_conePt0 = vSelectedLeptons.at(0).pt() ;
    LepGood_conePt1 = vSelectedLeptons.at(1).pt() ;

    signal_2lss_TTV_MVA = mva_2lss_ttV->EvaluateMVA("BDTG method");

    theHistoManager->fillHisto("Signal_2lss_TTV_MVA",                      "finalSel", "ttH2lss",   "",  signal_2lss_TTV_MVA,  weight);

    //std::cout << " signal 2lss TT MVA: "  << signal_2lss_TT_MVA
    //    << " signal 2lss TTV MVA: " << signal_2lss_TTV_MVA << std::endl;

    // ======================================================================================================

    fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

void TTbarHiggsMultileptonAnalysis::TwoLeptonsSameSignSelection_ApplicationFlips(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "ttH2lss",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "ttH2lss",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    if(DEBUG) std::cout << std::endl << " 2l ss OS ==================================================================" << std::endl;

    // ####################
    // # Common selection #
    // ####################

    bool nTight             = ( n_tight == 2);
    if(!nTight)             return;
    if(DEBUG) std::cout << "nTight Ok... ";

    bool are_fullytight     =  ( vSelectedLeptons.at(0).isTightTTH()        && vSelectedLeptons.at(1).isTightTTH()
                              && vSelectedLeptons.at(0).passTightCharge()   && vSelectedLeptons.at(1).passTightCharge()
                              && vSelectedLeptons.at(0).cutEventSel()       && vSelectedLeptons.at(1).cutEventSel()
                              && vSelectedLeptons.at(0).noLostHits()        && vSelectedLeptons.at(1).noLostHits() );
    if(!are_fullytight)     return;
    if(DEBUG) std::cout << "are_tight Ok... ";

    bool are_electrons      = ( vSelectedLeptons.at(0).isElectron() && vSelectedLeptons.at(1).isElectron() );
    if(!are_electrons)      return;

    bool leading_lep_pt     = (vSelectedLeptons.at(0).pt() > 25);
    if(!leading_lep_pt)     return;
    if(DEBUG) std::cout << "leading_lep_pt Ok... ";

    bool following_lep_pt   = (vSelectedLeptons.at(1).pt() > 15);
    if(!following_lep_pt)   return;
    if(DEBUG) std::cout << "following_lep_pt Ok... ";

    bool nJets              = ( vSelectedJets.size()        >= 4 );
    if(!nJets)              return;
    if(DEBUG) std::cout << "nJets Ok... ";

    bool nLooseBtag         = ( nLooseBJets                 >= 2 );
    bool nMediumBtag        = ( nMediumBJets                >= 1 );
    if(!nLooseBtag && !nMediumBtag)      return;
    if(DEBUG) std::cout << "Btag Ok... ";

    // Adding invariant mass cut on loose leptons pairs
    bool pass_invariantemasscut = true;
    for(int i=0; i<vLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLeptons.size(); j++)
        {
            if ( fabs( ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M() ) < 12 )
            { pass_invariantemasscut = false ;}
        }
    }
    if(!pass_invariantemasscut) return;
    if(DEBUG) std::cout << "invariantmasscut Ok... ";

    // ###############################
    // # Two leptons event selection #
    // ###############################

    bool same_sign      = ( vSelectedLeptons.at(0).charge() == vSelectedLeptons.at(1).charge() );
    if(same_sign)      return;
    if(DEBUG) std::cout << "oppositesign Ok... ";

    // ##########
    // # Z veto # here for leptons of same charge only !
    // ##########

    bool pass_Zveto = true;

    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).isElectron() &&  vSelectedLeptons.at(j).isElectron() )
               && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                "finalSel",   "ttH2lss",   "",   2, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "finalSel",   "ttH2lss", _sampleName.Data(),   3, weight);

    // #################################
    // # b-tagging nominal reweighting #
    // #################################

    std::vector<double> jetPts;
    std::vector<double> jetEtas;
    std::vector<double> jetCSVs;
    std::vector<int>    jetFlavors;
    int iSys = 0;
    double wgt_csv, wgt_csv_def, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf, new_weight;

    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetPts.push_back(     vSelectedJets.at(i).pt()                );
        jetEtas.push_back(    vSelectedJets.at(i).eta()               );
        jetCSVs.push_back(    vSelectedJets.at(i).CSVv2()             );
        jetFlavors.push_back( vSelectedJets.at(i).jet_hadronFlavour() );
    }

    wgt_csv_def = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
    new_weight = weight * wgt_csv_def; // weight = weight * wgt_csv_def;

    // ##################################################################################################################################

    // ##################################
    // # b-tagging deriving systematics #
    // ##################################

    std::vector<double> weights_csv;
    double wgt_csv_def_sys = 0;

    for(int i=7; i<25; i++)
    {
        wgt_csv_def_sys = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, i, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf)/wgt_csv_def;
        weights_csv.push_back(wgt_csv_def_sys);
    }

    double min_weight_csv = *min_element(weights_csv.begin(),weights_csv.end());
    double max_weight_csv = *max_element(weights_csv.begin(),weights_csv.end());
    theHistoManager->fillHisto("WeightCSV_min",  "",   "ttH2l",   "", min_weight_csv, 1);
    theHistoManager->fillHisto("WeightCSV_max",  "",   "ttH2l",   "", max_weight_csv, 1);

    weight_csv_down = min_weight_csv;
    weight_csv_up   = max_weight_csv;    

    // ##################################################################################################################################

    // ####################
    // # flip reweighting #
    // ####################

    std::vector<double> leptonsPts;
    std::vector<double> leptonsEtas;

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        if( !vSelectedLeptons.at(i).isTightTTH())
        {
            leptonsPts.push_back(     vSelectedLeptons.at(i).pt()                );
            leptonsEtas.push_back(    vSelectedLeptons.at(i).eta()               );
        }
    }

    // ##################################################################################################################################

    if(  pass_Zveto
      && (met_ld > 0.2) )
    {
        if(DEBUG) std::cout << "ee + Zveto + metld Ok... SELECTED";

        weightflip = get_QF_wgt_2l(leptonsPts, leptonsEtas);

        //std::cout << " Flip 2lss weight: " << weightflip << std::endl;

        if(vSelectedTaus.size() == 0)
        {
            stat_2lss_os_SB_ee = stat_2lss_os_SB_ee + 1;
            is_2lss_AppFlips_SR = true;
            cat_ee_2lss_QF = true;
            cat_ee_flip = vSelectedLeptons.at(0).charge();

            theHistoManager->fillHisto("CutFlow",                          "finalSel",   "ttH_2lss_ee",   "", 8                              , weight);

            theHistoManager->fillHisto("LeadingLeptonPt",                  "finalSel",   "ttH_2lss_ee",   "", vSelectedLeptons.at(0).pt()    , weight);
            theHistoManager->fillHisto("SubLeadingLeptonPt",               "finalSel",   "ttH_2lss_ee",   "", vSelectedLeptons.at(1).pt()    , weight);
            theHistoManager->fillHisto("MET",                              "finalSel",   "ttH_2lss_ee",   "", vEvent->at(0).metpt()          , weight);
            theHistoManager->fillHisto("MHT",                              "finalSel",   "ttH_2lss_ee",   "", MHT                            , weight);
            theHistoManager->fillHisto("MetLD",                            "finalSel",   "ttH_2lss_ee",   "", met_ld                         , weight);
            theHistoManager->fillHisto("TauMultiplicity",                  "finalSel",   "ttH_2lss_ee",   "", vSelectedTaus.size()           , weight);
            theHistoManager->fillHisto("JetMultiplicity",                  "finalSel",   "ttH_2lss_ee",   "", vSelectedJets.size()           , weight);
            theHistoManager->fillHisto("LooseBJetMultiplicity",            "finalSel",   "ttH_2lss_ee",   "", vSelectedBTagJets.size()       , weight);
            theHistoManager->fillHisto("MediumBJetMultiplicity",           "finalSel",   "ttH_2lss_ee",   "", vSelectedMediumBTagJets.size() , weight);
            theHistoManager->fillHisto("NumberOfSelectedLeptons",          "finalSel",   "ttH_2lss_ee",   "", vSelectedLeptons.size()        , weight);
        }

        if(vSelectedTaus.size() > 0)
        {
            cat_2ltau = true;
        }
    }

    // ####################################
    // #  ____  ____    ____  ____ _____  #
    // # |___ \|  _ \  | __ )|  _ \_   _| #
    // #   __) | | | | |  _ \| | | || |   #
    // #  / __/| |_| | | |_) | |_| || |   #
    // # |_____|____/  |____/|____/ |_|   #
    // #                                  #
    // ####################################

    // based on https://github.com/CERN-PH-CMG/cmgtools-lite/blob/0b47d4d1c50ea0e24ef0d9cf1c24c763e78c1bf0/TTHAnalysis/python/tools/kinMVA_2D_2lss_3l.py
    // and https://github.com/CERN-PH-CMG/cmgtools-lite/blob/0b47d4d1c50ea0e24ef0d9cf1c24c763e78c1bf0/TTHAnalysis/python/tools/eventVars_2lss.py

    // ======================================================================================================
    // variables against ttbar
    max_Lep_eta     = std::max( fabs(vSelectedLeptons.at(0).eta()), fabs(vSelectedLeptons.at(1).eta()) ) ;

    nJet25_Recl = vSelectedJets.size() ;

    mindr_lep1_jet = 1000.;
    mindr_lep2_jet = 1000.;

    for (int i=0; i<vSelectedJets.size(); i++)
    {
        if( DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ) < mindr_lep1_jet )
        { mindr_lep1_jet = DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ); }
        if( DeltaRLeptonJet( vSelectedLeptons.at(1), vSelectedJets.at(i) ) < mindr_lep2_jet )
        { mindr_lep2_jet = DeltaRLeptonJet( vSelectedLeptons.at(1), vSelectedJets.at(i) ); }
    }

    float met_max = 400;
    met             = std::min( vEvent->at(0).metpt(), met_max ) ;

    int njj = 0;
    avg_dr_jet = 0.;
    for (int ijet=0; ijet < vSelectedJets.size()-1 ; ijet++)
    {
        for (int kjet=ijet+1; kjet < vSelectedJets.size() ; kjet++)
        {
            avg_dr_jet += GetDeltaR( vSelectedJets.at(ijet).eta(), vSelectedJets.at(ijet).phi(), vSelectedJets.at(kjet).eta(), vSelectedJets.at(kjet).phi() );
            njj++;
        }
    }
    if ( njj > 0 ) avg_dr_jet = avg_dr_jet / njj;

    MT_met_lep1     = sqrt( 2 * vSelectedLeptons.at(0).pt() * vEvent->at(0).metpt() * (1 - cos( vSelectedLeptons.at(0).phi() - vEvent->at(0).metphi() )));

    //std::cout << "MT Calculated in 2lss flip" << std::endl;

    signal_2lss_TT_MVA  = mva_2lss_tt->EvaluateMVA("BDTG method");

    theHistoManager->fillHisto("Signal_2lss_TT_MVA",                       "finalSel", "ttH2lss",   "",  signal_2lss_TT_MVA,   weight);

    // ======================================================================================================
    // variables against ttV

    LepGood_conePt0 = vSelectedLeptons.at(0).pt() ;
    LepGood_conePt1 = vSelectedLeptons.at(1).pt() ;

    signal_2lss_TTV_MVA = mva_2lss_ttV->EvaluateMVA("BDTG method");

    theHistoManager->fillHisto("Signal_2lss_TTV_MVA",                      "finalSel", "ttH2lss",   "",  signal_2lss_TTV_MVA,  weight);

    //std::cout << " signal 2lss TT MVA: "  << signal_2lss_TT_MVA
    //    << " signal 2lss TTV MVA: " << signal_2lss_TTV_MVA << std::endl;

    // ======================================================================================================

    fillOutputTree();
}

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_TTH3l(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "ttH_3l",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "ttH_3l",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)
    //if(vSelectedTaus.size()>0) return;

    // ####################
    // # Common selection #
    // ####################

    if(DEBUG)             std::cout << std::endl << " 3l ss SR ==================================================================" << std::endl;

    bool nLep               = ( vSelectedLeptons.size()     >= 3 );
    if(!nLep)               return;

    if(DEBUG) std::cout << "nLep Ok... ";

    bool nTight             = ( n_tight                     >= 3 );
    if(!nTight)             return;

    bool are_fullytight     =  ( vSelectedLeptons.at(0).isTightTTH()        && vSelectedLeptons.at(1).isTightTTH()      && vSelectedLeptons.at(2).isTightTTH()
                              //&& vSelectedLeptons.at(0).passTightCharge()   && vSelectedLeptons.at(1).passTightCharge() && vSelectedLeptons.at(2).passTightCharge()
                              && vSelectedLeptons.at(0).cutEventSel()       && vSelectedLeptons.at(1).cutEventSel()     && vSelectedLeptons.at(2).cutEventSel()
                              && vSelectedLeptons.at(0).noLostHits()        && vSelectedLeptons.at(1).noLostHits()      && vSelectedLeptons.at(2).noLostHits()      );
    if(!are_fullytight)     return;
    if(DEBUG) std::cout << "nTight Ok... ";

    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt() > 20 );
    if(!leading_lep_pt)     return;

    if(DEBUG) std::cout << "leading_lep_pt Ok... ";

    bool following_lep_pt   = ( (vSelectedLeptons.at(1).pt() > 10) && (vSelectedLeptons.at(2).pt() > 10) );
    if(!following_lep_pt)   return;

    if(DEBUG) std::cout << "following_lep_pt Ok... ";

    bool nJets              = ( vSelectedJets.size()        >= 2 );
    if(!nJets)              return;

    if(DEBUG) std::cout << "nJets Ok... ";

    bool nLooseBtag         = ( nLooseBJets                 >= 2 );
    bool nMediumBtag        = ( nMediumBJets                >= 1 );
    if(!nLooseBtag && !nMediumBtag)      return;

    if(DEBUG) std::cout << "Btag Ok... ";

    // Adding invariant mass cut on loose leptons pairs
    bool pass_invariantemasscut = true;
    for(int i=0; i<vLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLeptons.size(); j++)
        {
            if ( fabs( ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M() ) < 12 )
            { pass_invariantemasscut = false ;}
        }
    }
    if(!pass_invariantemasscut) return;    

    if(DEBUG) std::cout << "invariantmasscut Ok... ";

    // ##########
    // # Z veto #
    // ##########

    bool pass_Zveto = true;
    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id()     == -vSelectedLeptons.at(j).id()                             )
               && ( vSelectedLeptons.at(i).charge() == -vSelectedLeptons.at(j).charge()                         )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
    }
    if(!pass_Zveto)       return;

    theHistoManager->fillHisto("CutFlow",                "finalSel",   "ttH_3l",   "",   2, weight);

    if(DEBUG) std::cout << "Zveto Ok... ";

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, tau_px = 0, tau_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    TLorentzVector taup4;
    for(int i=0; i<vSelectedTaus.size(); i++)
    {
        taup4.SetPtEtaPhiE(vSelectedTaus.at(i).pt(), vSelectedTaus.at(i).eta(), vSelectedTaus.at(i).phi(), vSelectedTaus.at(i).E());
        tau_px = tau_px + taup4.Px();
        tau_py = tau_py + taup4.Py();
    }

    MHT = sqrt( (jet_px + lepton_px + tau_px) * (jet_px + lepton_px + tau_px) + (jet_py + lepton_py + tau_py) * (jet_py + lepton_py + tau_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "finalSel",   "ttH_3l", _sampleName.Data(),   3, weight);
    if( vSelectedJets.size() >= 4  )  theHistoManager->fillHisto("CutFlow",                "finalSel",   "ttH_3l", _sampleName.Data(),   4, weight);

    if(DEBUG) std::cout << " MHT =  " << MHT << "MET = " << vEvent->at(0).metpt() << " met_ld = " << met_ld ;

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        for(int j=0; j<vSelectedLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    if( met_ld               > 0.3 ) theHistoManager->fillHisto("CutFlow",                "finalSel",   "ttH_3l",   "",   6, weight);
    if( isSFOS                     ) theHistoManager->fillHisto("CutFlow",                "finalSel",   "ttH_3l",   "",   7, weight);

    if(DEBUG) std::cout << std::endl << "nJets: " << vSelectedJets.size() << " met_ld: " << met_ld << " isSFOS " << isSFOS << std::endl;

    if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return;

    if(DEBUG) std::cout << "nJets and met_ld Ok... ";

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vSelectedLeptons.at(i).charge();
    }
    if( fabs(sum_charges_3l) != 1 ) return;

    if(DEBUG) std::cout << "sumOfCharges Ok... ";

    // #################################
    // # b-tagging nominal reweighting #
    // #################################

    std::vector<double> jetPts;
    std::vector<double> jetEtas;
    std::vector<double> jetCSVs;
    std::vector<int>    jetFlavors;
    int iSys = 0;
    double wgt_csv, wgt_csv_def, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf, new_weight;

    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetPts.push_back(     vSelectedJets.at(i).pt()                );
        jetEtas.push_back(    vSelectedJets.at(i).eta()               );
        jetCSVs.push_back(    vSelectedJets.at(i).CSVv2()             );
        jetFlavors.push_back( vSelectedJets.at(i).jet_hadronFlavour() );
    }

    wgt_csv_def = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
    new_weight = weight * wgt_csv_def; // weight = weight * wgt_csv_def;

    // ##################################################################################################################################

    // ##################################
    // # b-tagging deriving systematics #
    // ##################################

    std::vector<double> weights_csv;
    double wgt_csv_def_sys = 0;

    for(int i=7; i<25; i++)
    {
        wgt_csv_def_sys = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, i, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf)/wgt_csv_def;
        weights_csv.push_back(wgt_csv_def_sys);
    }

    double min_weight_csv = *min_element(weights_csv.begin(),weights_csv.end()); 
    double max_weight_csv = *max_element(weights_csv.begin(),weights_csv.end()); 
    theHistoManager->fillHisto("WeightCSV_min",  "",   "ttH2l",   "", min_weight_csv, 1);
    theHistoManager->fillHisto("WeightCSV_max",  "",   "ttH2l",   "", max_weight_csv, 1);

    weight_csv_down = min_weight_csv;
    weight_csv_up   = max_weight_csv;

    // ##################################################################################################################################

    stat_3l_SR  = stat_3l_SR + 1;
    cat_3l      = sum_charges_3l;
    is_3l_TTH_SR = true;

    theHistoManager->fillHisto("CutFlow",                          "finalSel",   "ttH_3l",   "", 8                              , weight);

    theHistoManager->fillHisto("LeadingLeptonPt",                  "finalSel",   "ttH_3l",   "", vSelectedLeptons.at(0).pt()    , weight);
    theHistoManager->fillHisto("SubLeadingLeptonPt",               "finalSel",   "ttH_3l",   "", vSelectedLeptons.at(1).pt()    , weight);
    theHistoManager->fillHisto("ThirdLeptonPt",                    "finalSel",   "ttH_3l",   "", vSelectedLeptons.at(2).pt()    , weight);
    theHistoManager->fillHisto("MET",                              "finalSel",   "ttH_3l",   "", vEvent->at(0).metpt()          , weight);
    theHistoManager->fillHisto("MHT",                              "finalSel",   "ttH_3l",   "", MHT                            , weight);
    theHistoManager->fillHisto("MetLD",                            "finalSel",   "ttH_3l",   "", met_ld                         , weight);
    theHistoManager->fillHisto("TauMultiplicity",                  "finalSel",   "ttH_3l",   "", vSelectedTaus.size()           , weight);
    theHistoManager->fillHisto("JetMultiplicity",                  "finalSel",   "ttH_3l",   "", vSelectedJets.size()           , weight);
    theHistoManager->fillHisto("LooseBJetMultiplicity",            "finalSel",   "ttH_3l",   "", vSelectedBTagJets.size()       , weight);
    theHistoManager->fillHisto("MediumBJetMultiplicity",           "finalSel",   "ttH_3l",   "", vSelectedMediumBTagJets.size() , weight);
    theHistoManager->fillHisto("SumOfThreeLeptonsCharges",         "finalSel",   "ttH_3l",   "", sum_charges_3l                 , weight);
    theHistoManager->fillHisto("NumberOfSelectedLeptons",          "finalSel",   "ttH_3l",   "", vSelectedLeptons.size()        , weight);

    if(DEBUG) std::cout << "Ok... SELECTED";

    // ####################################
    // #  ____  ____    ____  ____ _____  #
    // # |___ \|  _ \  | __ )|  _ \_   _| #
    // #   __) | | | | |  _ \| | | || |   #
    // #  / __/| |_| | | |_) | |_| || |   #
    // # |_____|____/  |____/|____/ |_|   #
    // #                                  #
    // ####################################

    // based on https://github.com/CERN-PH-CMG/cmgtools-lite/blob/0b47d4d1c50ea0e24ef0d9cf1c24c763e78c1bf0/TTHAnalysis/python/tools/kinMVA_2D_2lss_3l.py

    // ======================================================================================================
    // variables against ttbar

    max_Lep_eta     = std::max( fabs(vSelectedLeptons.at(0).eta()), fabs(vSelectedLeptons.at(1).eta()) ) ;

    MT_met_lep1     = sqrt( 2 * vSelectedLeptons.at(0).pt() * vEvent->at(0).metpt() * (1 - cos( vSelectedLeptons.at(0).phi() - vEvent->at(0).metphi() )));

    nJet25_Recl     = vSelectedJets.size() ;
    mhtJet25_Recl   = sqrt( (jet_px*jet_px) + (jet_py*jet_py) );

    int njj = 0;
    avg_dr_jet = 0.;
    for (int ijet=0; ijet < vSelectedJets.size()-1 ; ijet++)
    {
        for (int kjet=ijet+1; kjet < vSelectedJets.size() ; kjet++)
        {
            avg_dr_jet += GetDeltaR( vSelectedJets.at(ijet).eta(), vSelectedJets.at(ijet).phi(), vSelectedJets.at(kjet).eta(), vSelectedJets.at(kjet).phi() );
            njj++;
        }
    }
    if ( njj > 0 ) avg_dr_jet = avg_dr_jet / njj;

    mindr_lep1_jet = 1000.;
    mindr_lep2_jet = 1000.;

    for (int i=0; i<vSelectedJets.size(); i++)
    {
        if( DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ) < mindr_lep1_jet )
        { mindr_lep1_jet = DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ); }
        if( DeltaRLeptonJet( vSelectedLeptons.at(2), vSelectedJets.at(i) ) < mindr_lep2_jet )
        { mindr_lep2_jet = DeltaRLeptonJet( vSelectedLeptons.at(2), vSelectedJets.at(i) ); }
    }

    signal_3l_TT_MVA    = mva_3l_tt->EvaluateMVA("BDTG method");

    theHistoManager->fillHisto("Signal_3l_TT_MVA",                         "finalSel",   "ttH_3l",   "",  signal_3l_TT_MVA,   weight);

    // ======================================================================================================
    // variables against ttV

    met             = vEvent->at(0).metpt() ;

    LepGood_conePt0 = vSelectedLeptons.at(0).pt() ;
    LepGood_conePt1 = vSelectedLeptons.at(2).pt() ;

    signal_3l_TTV_MVA   = mva_3l_ttV->EvaluateMVA("BDTG method");

    theHistoManager->fillHisto("Signal_3l_TTV_MVA",                        "finalSel",   "ttH_3l",   "",  signal_3l_TTV_MVA,  weight);

    //std::cout << " signal 3l   TT MVA: "  << signal_3l_TT_MVA
    //          << " signal 3l   TTV MVA: " << signal_3l_TTV_MVA << std::endl;

    // ======================================================================================================

    fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_ApplicationFakes(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "ttH_3l",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "ttH_3l",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)
    //if(vSelectedTaus.size()>0) return;

    // ####################
    // # Common selection #
    // ####################

    if(DEBUG)             std::cout << std::endl << " 3l ss FR ==================================================================" << std::endl;

    bool nLep               = ( vSelectedLeptons.size()     >= 3 );
    if(!nLep)               return;

    if(DEBUG) std::cout << "nLep Ok... ";

    bool are_fullytight     =  ( vSelectedLeptons.at(0).isTightTTH()        && vSelectedLeptons.at(1).isTightTTH()      && vSelectedLeptons.at(2).isTightTTH()
                              //&& vSelectedLeptons.at(0).passTightCharge()   && vSelectedLeptons.at(1).passTightCharge() && vSelectedLeptons.at(2).passTightCharge()
                              && vSelectedLeptons.at(0).cutEventSel()       && vSelectedLeptons.at(1).cutEventSel()     && vSelectedLeptons.at(2).cutEventSel()
                              && vSelectedLeptons.at(0).noLostHits()        && vSelectedLeptons.at(1).noLostHits()      && vSelectedLeptons.at(2).noLostHits()      );
    if(are_fullytight)     return;
    if(DEBUG) std::cout << "nNonTight Ok... ";

    bool pass_nolosthits    =  ( vSelectedLeptons.at(0).noLostHits()        && vSelectedLeptons.at(1).noLostHits()       && vSelectedLeptons.at(2).noLostHits()     );
    if(!pass_nolosthits)     return;
    if(DEBUG) std::cout << "pass_nolosthits Ok... ";

    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt() > 20 );
    if(!leading_lep_pt)     return;

    if(DEBUG) std::cout << "leading_lep_pt Ok... ";

    bool following_lep_pt   = ( (vSelectedLeptons.at(1).pt() > 10) && (vSelectedLeptons.at(2).pt() > 10) );
    if(!following_lep_pt)   return;

    if(DEBUG) std::cout << "following_lep_pt Ok... ";

    bool nJets              = ( vSelectedJets.size()        >= 2 );
    if(!nJets)              return;

    if(DEBUG) std::cout << "nJets Ok... ";

    bool nLooseBtag         = ( nLooseBJets                 >= 2 );
    bool nMediumBtag        = ( nMediumBJets                >= 1 );
    if(!nLooseBtag && !nMediumBtag)      return;

    if(DEBUG) std::cout << "Btag Ok...";

    // Adding invariant mass cut on loose leptons pairs
    bool pass_invariantemasscut = true;
    for(int i=0; i<vLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLeptons.size(); j++)
        {
            if ( fabs( ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M() ) < 12 )
            { pass_invariantemasscut = false ;}
        }
    }
    if(!pass_invariantemasscut) return;

    if(DEBUG) std::cout << "invariantmasscut Ok... ";

    // ##########
    // # Z veto #
    // ##########

    bool pass_Zveto = true;
    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id()     == -vSelectedLeptons.at(j).id()                             )
               && ( vSelectedLeptons.at(i).charge() == -vSelectedLeptons.at(j).charge()                         )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
    }
    if(!pass_Zveto)       return;

    if(DEBUG) std::cout << "Zveto Ok... ";

    theHistoManager->fillHisto("CutFlow",                "finalSel",   "ttH_3l",   "",   2, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "finalSel",   "ttH_3l", _sampleName.Data(),   3, weight);
    if( vSelectedJets.size() >= 4  )  theHistoManager->fillHisto("CutFlow",                "finalSel",   "ttH_3l", _sampleName.Data(),   4, weight);

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vLeptons.size(); i++)
    {
        for(int j=0; j<vLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vLeptons.at(i).id() == -vLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    if( met_ld               > 0.3 ) theHistoManager->fillHisto("CutFlow",                "finalSel",   "ttH_3l",   "",   6, weight);
    if( isSFOS                     ) theHistoManager->fillHisto("CutFlow",                "finalSel",   "ttH_3l",   "",   7, weight);

    if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return;

    if(DEBUG) std::cout << "nJets and met_ld Ok... ";

    int sum_charges = 0;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        sum_charges = sum_charges + vSelectedLeptons.at(i).charge();
    }

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vSelectedLeptons.at(i).charge();
    }
    if( fabs(sum_charges_3l) != 1 ) return;

    if(DEBUG) std::cout << "sumOfCharges Ok... ";

    // ##################################################################################################################################

    // #################################
    // # b-tagging nominal reweighting #
    // #################################

    std::vector<double> jetPts;
    std::vector<double> jetEtas;
    std::vector<double> jetCSVs;
    std::vector<int>    jetFlavors;
    int iSys = 0;
    double wgt_csv, wgt_csv_def, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf, new_weight;

    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetPts.push_back(     vSelectedJets.at(i).pt()                );
        jetEtas.push_back(    vSelectedJets.at(i).eta()               );
        jetCSVs.push_back(    vSelectedJets.at(i).CSVv2()             );
        jetFlavors.push_back( vSelectedJets.at(i).jet_hadronFlavour() );
    }

    wgt_csv_def = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
    new_weight = weight * wgt_csv_def; // weight = weight * wgt_csv_def;

    // ##################################################################################################################################

    // ##################################
    // # b-tagging deriving systematics #
    // ##################################

    std::vector<double> weights_csv;
    double wgt_csv_def_sys = 0;

    for(int i=7; i<25; i++)
    {
        wgt_csv_def_sys = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, i, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf)/wgt_csv_def;
        weights_csv.push_back(wgt_csv_def_sys);
    }

    double min_weight_csv = *min_element(weights_csv.begin(),weights_csv.end()); 
    double max_weight_csv = *max_element(weights_csv.begin(),weights_csv.end()); 
    theHistoManager->fillHisto("WeightCSV_min",  "",   "ttH2l",   "", min_weight_csv, 1);
    theHistoManager->fillHisto("WeightCSV_max",  "",   "ttH2l",   "", max_weight_csv, 1);

    weight_csv_down = min_weight_csv;
    weight_csv_up   = max_weight_csv;

    // ##################################################################################################################################

    // #########################
    // # fake rate reweighting #
    // #########################

    std::vector<double> leptonsPts;
    std::vector<double> leptonsEtas;
    std::vector<int>    leptonsIds;

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        if( !vSelectedLeptons.at(i).isTightTTH())
        {
            leptonsPts.push_back(     vSelectedLeptons.at(i).pt()                );
            leptonsEtas.push_back(    vSelectedLeptons.at(i).eta()               );
            leptonsIds.push_back(     vSelectedLeptons.at(i).id()                );
        }
    }

    weightfake = get_FR_wgt_3l(leptonsPts, leptonsEtas, leptonsIds);

    // ##################################################################################################################################

    stat_3l_lepMVA_SB   = stat_3l_lepMVA_SB + 1;
    cat_3l_fake         = sum_charges_3l;
    is_3l_AppFakes_SR = true;

    //std::cout << "Fake 3l weight: " << weightfake << std::endl;

    theHistoManager->fillHisto("CutFlow",                          "finalSel",   "ttH_3l",   "", 8                              , weight);

    theHistoManager->fillHisto("LeadingLeptonPt",                  "finalSel",   "ttH_3l",   "", vSelectedLeptons.at(0).pt()    , weight);
    theHistoManager->fillHisto("SubLeadingLeptonPt",               "finalSel",   "ttH_3l",   "", vSelectedLeptons.at(1).pt()    , weight);
    theHistoManager->fillHisto("ThirdLeptonPt",                    "finalSel",   "ttH_3l",   "", vSelectedLeptons.at(2).pt()    , weight);
    theHistoManager->fillHisto("MET",                              "finalSel",   "ttH_3l",   "", vEvent->at(0).metpt()          , weight);
    theHistoManager->fillHisto("MHT",                              "finalSel",   "ttH_3l",   "", MHT                            , weight);
    theHistoManager->fillHisto("MetLD",                            "finalSel",   "ttH_3l",   "", met_ld                         , weight);
    theHistoManager->fillHisto("TauMultiplicity",                  "finalSel",   "ttH_3l",   "", vSelectedTaus.size()           , weight);
    theHistoManager->fillHisto("JetMultiplicity",                  "finalSel",   "ttH_3l",   "", vSelectedJets.size()           , weight);
    theHistoManager->fillHisto("LooseBJetMultiplicity",            "finalSel",   "ttH_3l",   "", vSelectedBTagJets.size()       , weight);
    theHistoManager->fillHisto("MediumBJetMultiplicity",           "finalSel",   "ttH_3l",   "", vSelectedMediumBTagJets.size() , weight);
    theHistoManager->fillHisto("SumOfLeptonsCharges",              "finalSel",   "ttH_3l",   "", sum_charges                    , weight);
    theHistoManager->fillHisto("SumOfThreeLeptonsCharges",         "finalSel",   "ttH_3l",   "", sum_charges_3l                 , weight);
    theHistoManager->fillHisto("NumberOfSelectedLeptons",          "finalSel",   "ttH_3l",   "", vSelectedLeptons.size()        , weight);

    // ####################################
    // #  ____  ____    ____  ____ _____  #
    // # |___ \|  _ \  | __ )|  _ \_   _| #
    // #   __) | | | | |  _ \| | | || |   #
    // #  / __/| |_| | | |_) | |_| || |   #
    // # |_____|____/  |____/|____/ |_|   #
    // #                                  #
    // ####################################

    // based on https://github.com/CERN-PH-CMG/cmgtools-lite/blob/0b47d4d1c50ea0e24ef0d9cf1c24c763e78c1bf0/TTHAnalysis/python/tools/kinMVA_2D_2lss_3l.py

    // ======================================================================================================
    // variables against ttbar

    max_Lep_eta     = std::max( fabs(vSelectedLeptons.at(0).eta()), fabs(vSelectedLeptons.at(1).eta()) ) ;

    MT_met_lep1     = sqrt( 2 * vSelectedLeptons.at(0).p4().Pt() * vEvent->at(0).metpt() * (1 - cos( vSelectedLeptons.at(0).phi() - vEvent->at(0).metphi() )));

    nJet25_Recl     = vSelectedJets.size() ;
    mhtJet25_Recl   = sqrt( (jet_px*jet_px) + (jet_py*jet_py) );

    int njj = 0;
    avg_dr_jet = 0.;
    for (int ijet=0; ijet < vSelectedJets.size()-1 ; ijet++)
    {
        for (int kjet=ijet+1; kjet < vSelectedJets.size() ; kjet++)
        {
            avg_dr_jet += GetDeltaR( vSelectedJets.at(ijet).eta(), vSelectedJets.at(ijet).phi(), vSelectedJets.at(kjet).eta(), vSelectedJets.at(kjet).phi() );
            njj++;
        }
    }
    if ( njj > 0 ) avg_dr_jet = avg_dr_jet / njj;

    mindr_lep1_jet = 1000.;
    mindr_lep2_jet = 1000.;

    for (int i=0; i<vSelectedJets.size(); i++)
    {
        if( DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ) < mindr_lep1_jet )
        { mindr_lep1_jet = DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ); }
        if( DeltaRLeptonJet( vSelectedLeptons.at(2), vSelectedJets.at(i) ) < mindr_lep2_jet )
        { mindr_lep2_jet = DeltaRLeptonJet( vSelectedLeptons.at(2), vSelectedJets.at(i) ); }
    }

    signal_3l_TT_MVA    = mva_3l_tt->EvaluateMVA("BDTG method");

    theHistoManager->fillHisto("Signal_3l_TT_MVA",                         "finalSel",   "ttH_3l",   "",  signal_3l_TT_MVA,   weight);

    // ======================================================================================================
    // variables against ttV

    met             = vEvent->at(0).metpt() ;

    LepGood_conePt0 = vSelectedLeptons.at(0).pt() ;
    LepGood_conePt1 = vSelectedLeptons.at(2).pt() ;

    signal_3l_TTV_MVA   = mva_3l_ttV->EvaluateMVA("BDTG method");

    theHistoManager->fillHisto("Signal_3l_TTV_MVA",                        "finalSel",   "ttH_3l",   "",  signal_3l_TTV_MVA,  weight);

    //std::cout << " signal 3l   TT MVA: "  << signal_3l_TT_MVA
    //          << " signal 3l   TTV MVA: " << signal_3l_TTV_MVA << std::endl;

    // ======================================================================================================

    fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

// # ######################################################################
// #                  _             _                  _                  #
// #   ___ ___  _ __ | |_ _ __ ___ | |  _ __ ___  __ _(_) ___  _ __  ___  #
// #  / __/ _ \| '_ \| __| '__/ _ \| | | '__/ _ \/ _` | |/ _ \| '_ \/ __| #
// # | (_| (_) | | | | |_| | | (_) | | | | |  __/ (_| | | (_) | | | \__ \ #
// #  \___\___/|_| |_|\__|_|  \___/|_| |_|  \___|\__, |_|\___/|_| |_|___/ #
// #                                             |___/                    #
// #                                                                      #
// ########################################################################

void TTbarHiggsMultileptonAnalysis::TwoLeptonsSameSignSelection_LepMVA_sideband(int evt) // TO BE UPDATED
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "LepMVA_2l_SB",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "LepMVA_2l_SB",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep               = ( vInclusiveFakeLeptons.size()              == 2 );
    bool nLepTight          = ( vSelectedLeptons.size()                   == 1 );
    if(!nLep)               return;
    if(!nLepTight)          return;

    bool leading_lep_pt     = ( vInclusiveFakeLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)     return;

    //bool following_lep_pt   = ( vSelectedLeptons.at(1).pt()               > 10 );
    bool following_lep_pt   = (  ( (abs(vInclusiveFakeLeptons.at(1).id()) == 11) && (vInclusiveFakeLeptons.at(1).pt() > 15) )
            || ( (abs(vInclusiveFakeLeptons.at(1).id()) == 13) && (vInclusiveFakeLeptons.at(1).pt() > 10) ) );
    if(!following_lep_pt)   return;

    bool passMll12Gt12      = ( ( vInclusiveFakeLeptons.at(0).p4() + vInclusiveFakeLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)      return;

    bool nJets              = ( vSelectedJets.size()                      >= 4 );
    if(!nJets)              return;

    bool nLooseBtag         = ( nLooseBJets                               >= 2 );
    bool nMediumBtag        = ( nMediumBJets                              >= 1 );
    if(!nLooseBtag && !nMediumBtag)      return;


    // ###############################
    // # Two leptons event selection #
    // ###############################

    bool same_sign      = ( vInclusiveFakeLeptons.at(0).charge() == vInclusiveFakeLeptons.at(1).charge() );
    if(!same_sign)      return;

    // ##########
    // # Z veto #
    // ##########

    bool pass_Zveto = true;
    for(int i=0; i<vInclusiveFakeLeptons.size(); i++)
    {
        for(int j=0; j<vInclusiveFakeLeptons.size(); j++)
        {
            if (  ( i                           != j                                                       )
                    && ( vInclusiveFakeLeptons.at(i).id() == -vInclusiveFakeLeptons.at(j).id()                            )
                    && ( fabs( ( vInclusiveFakeLeptons.at(i).p4() + vInclusiveFakeLeptons.at(j).p4() ).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                "finalSel",   "LepMVA_2l_SB",   "",   2, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vInclusiveFakeLeptons.size(); i++)
    {
        lepton_px = lepton_px + vInclusiveFakeLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vInclusiveFakeLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "finalSel",   "LepMVA_2l_SB", _sampleName.Data(),   3, weight);

    if(  (abs(vInclusiveFakeLeptons.at(0).id()) == 11)
            && (abs(vInclusiveFakeLeptons.at(1).id()) == 11)
            && (pass_Zveto                            )
            && (met_ld                           > 0.2) // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "finalSel",   "LepMVA_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "finalSel",   "LepMVA_2l_SB",   "", vInclusiveFakeLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "finalSel",   "LepMVA_2l_SB",   "", vInclusiveFakeLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "finalSel",   "LepMVA_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "finalSel",   "LepMVA_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "finalSel",   "LepMVA_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "finalSel",   "LepMVA_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "finalSel",   "LepMVA_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "finalSel",   "LepMVA_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "finalSel",   "LepMVA_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "finalSel",   "LepMVA_2l_SB",   "", vInclusiveFakeLeptons.size()        , weight);
    }

    if(  ( (abs(vInclusiveFakeLeptons.at(0).id()) == 11) && (abs(vInclusiveFakeLeptons.at(1).id()) == 13) )
            || ( (abs(vInclusiveFakeLeptons.at(0).id()) == 13) && (abs(vInclusiveFakeLeptons.at(1).id()) == 11) ) // smarter way to do is probably exists
            && true                                      // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "finalSel",   "LepMVA_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "finalSel",   "LepMVA_2l_SB",   "", vInclusiveFakeLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "finalSel",   "LepMVA_2l_SB",   "", vInclusiveFakeLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "finalSel",   "LepMVA_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "finalSel",   "LepMVA_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "finalSel",   "LepMVA_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "finalSel",   "LepMVA_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "finalSel",   "LepMVA_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "finalSel",   "LepMVA_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "finalSel",   "LepMVA_2l_SB",   "", vInclusiveFakeLeptons.size()        , weight);
    }

    if(  (abs(vInclusiveFakeLeptons.at(0).id()) == 13)
            && (abs(vInclusiveFakeLeptons.at(1).id()) == 13)
            && true                                      // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "finalSel",   "LepMVA_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "finalSel",   "LepMVA_2l_SB",   "", vInclusiveFakeLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "finalSel",   "LepMVA_2l_SB",   "", vInclusiveFakeLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "finalSel",   "LepMVA_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "finalSel",   "LepMVA_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "finalSel",   "LepMVA_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "finalSel",   "LepMVA_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "finalSel",   "LepMVA_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "finalSel",   "LepMVA_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "finalSel",   "LepMVA_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "finalSel",   "LepMVA_2l_SB",   "", vInclusiveFakeLeptons.size()        , weight);
    }

    is_2lss_LepMVA_SB = true;   

    //fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

void TTbarHiggsMultileptonAnalysis::TwoLeptonsSameSignSelection_JetMultiplicity_sideband(int evt) // TO BE UPDATED
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "JM_2l_SB",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "JM_2l_SB",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep               = ( vSelectedLeptons.size()                   == 2 );
    if(!nLep)               return;

    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)     return;

    bool following_lep_pt   = (  ( (abs(vSelectedLeptons.at(1).id()) == 11) && (vSelectedLeptons.at(1).pt() > 15) )
            || ( (abs(vSelectedLeptons.at(1).id()) == 13) && (vSelectedLeptons.at(1).pt() > 10) ) );
    if(!following_lep_pt)   return;

    bool passMll12Gt12      = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)      return;

    bool nJets              = ( vSelectedJets.size()                      == 3 );
    if(!nJets)              return;

    bool nLooseBtag         = ( nLooseBJets                               >= 2 );
    bool nMediumBtag        = ( nMediumBJets                              >= 1 );
    if(!nLooseBtag && !nMediumBtag)      return;

    // ###############################
    // # Two leptons event selection #
    // ###############################

    bool same_sign      = ( vSelectedLeptons.at(0).charge() == vSelectedLeptons.at(1).charge() );
    if(!same_sign)      return;

    // ##########
    // # Z veto #
    // ##########

    bool pass_Zveto = true;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        for(int j=0; j<vSelectedLeptons.size(); j++)
        {
            if (  ( i                           != j                                                       )
                    && ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()                            )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                "finalSel",   "JM_2l_SB",   "",   2, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "finalSel",   "JM_2l_SB", _sampleName.Data(),   3, weight);

    if(  (abs(vSelectedLeptons.at(0).id()) == 11)
            && (abs(vSelectedLeptons.at(1).id()) == 11)
            && (pass_Zveto                            )
            && (met_ld                           > 0.2) // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "finalSel",   "JM_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "finalSel",   "JM_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "finalSel",   "JM_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "finalSel",   "JM_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "finalSel",   "JM_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "finalSel",   "JM_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "finalSel",   "JM_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "finalSel",   "JM_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "finalSel",   "JM_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "finalSel",   "JM_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "finalSel",   "JM_2l_SB",   "", vSelectedLeptons.size()        , weight);
    }

    if(  ( (abs(vSelectedLeptons.at(0).id()) == 11) && (abs(vSelectedLeptons.at(1).id()) == 13) )
            || ( (abs(vSelectedLeptons.at(0).id()) == 13) && (abs(vSelectedLeptons.at(1).id()) == 11) ) // smarter way to do is probably exists
            && true                                      // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "finalSel",   "JM_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "finalSel",   "JM_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "finalSel",   "JM_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "finalSel",   "JM_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "finalSel",   "JM_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "finalSel",   "JM_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "finalSel",   "JM_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "finalSel",   "JM_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "finalSel",   "JM_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "finalSel",   "JM_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "finalSel",   "JM_2l_SB",   "", vSelectedLeptons.size()        , weight);
    }

    if(  (abs(vSelectedLeptons.at(0).id()) == 13)
            && (abs(vSelectedLeptons.at(1).id()) == 13)
            && true                                      // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "finalSel",   "JM_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "finalSel",   "JM_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "finalSel",   "JM_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "finalSel",   "JM_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "finalSel",   "JM_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "finalSel",   "JM_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "finalSel",   "JM_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "finalSel",   "JM_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "finalSel",   "JM_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "finalSel",   "JM_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "finalSel",   "JM_2l_SB",   "", vSelectedLeptons.size()        , weight);
    }

    is_2lss_JM_SB = true;   

    //fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

void TTbarHiggsMultileptonAnalysis::TwoLeptonsSameSignSelection_TTbar(int evt) // TO BE UPDATED
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "TT_2l_CR",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "TT_2l_CR",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep               = ( vSelectedLeptons.size()                   == 2 );
    if(!nLep)               return;

    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)     return;

    bool following_lep_pt   = (  ( (abs(vSelectedLeptons.at(1).id()) == 11) && (vSelectedLeptons.at(1).pt() > 15) )
            || ( (abs(vSelectedLeptons.at(1).id()) == 13) && (vSelectedLeptons.at(1).pt() > 10) ) );
    if(!following_lep_pt)   return;

    bool passMll12Gt12      = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)      return;

    bool nJets              = ( vSelectedJets.size()                      >= 2 );
    if(!nJets)              return;

    bool nMediumBtag        = ( nMediumBJets                              >= 1 );
    if(!nMediumBtag)      return;

    // ###############################
    // #        e+mu- selection      #
    // ###############################

    if ( vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge() != 0 ) return;
    if ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(1).id() ) return;

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if(  ( (abs(vSelectedLeptons.at(0).id()) == 11) && (abs(vSelectedLeptons.at(1).id()) == 13) )
            || ( (abs(vSelectedLeptons.at(0).id()) == 13) && (abs(vSelectedLeptons.at(1).id()) == 11) ) // smarter way to do this probably exists
            && (  vSelectedLeptons.at(0).charge()         == -vSelectedLeptons.at(1).charge()         )
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "finalSel",   "TT_2l_CR",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "finalSel",   "TT_2l_CR",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "finalSel",   "TT_2l_CR",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "finalSel",   "TT_2l_CR",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MetLD",                            "finalSel",   "TT_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "finalSel",   "TT_2l_CR",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "finalSel",   "TT_2l_CR",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "finalSel",   "TT_2l_CR",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "finalSel",   "TT_2l_CR",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "finalSel",   "TT_2l_CR",   "", vSelectedLeptons.size()        , weight);
    }

    is_emu_TT_CR = true;   

    //fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_CR_WZ(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel", "WZ_CR",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel", "WZ_CR",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    if(DEBUG)             std::cout << std::endl << " 3l ss SR ==================================================================" << std::endl;

    bool nLep               = ( vSelectedLeptons.size()     >= 3 );
    if(!nLep)               return;

    if(DEBUG) std::cout << "nLep Ok... ";

    bool nTight             = ( n_tight                     >= 3 );
    if(!nTight)             return;

    bool are_fullytight     =  ( vSelectedLeptons.at(0).isTightTTH()        && vSelectedLeptons.at(1).isTightTTH()      && vSelectedLeptons.at(2).isTightTTH()
                              && vSelectedLeptons.at(0).cutEventSel()       && vSelectedLeptons.at(1).cutEventSel()     && vSelectedLeptons.at(2).cutEventSel()
                              && vSelectedLeptons.at(0).noLostHits()        && vSelectedLeptons.at(1).noLostHits()      && vSelectedLeptons.at(2).noLostHits()      );
    if(!are_fullytight)     return;
    if(DEBUG) std::cout << "nTight Ok... ";

    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt() > 20 );
    if(!leading_lep_pt)     return;

    if(DEBUG) std::cout << "leading_lep_pt Ok... ";

    bool following_lep_pt   = ( (vSelectedLeptons.at(1).pt() > 10) && (vSelectedLeptons.at(2).pt() > 10) );
    if(!following_lep_pt)   return;

    if(DEBUG) std::cout << "following_lep_pt Ok... ";

    bool nJets              = ( vSelectedJets.size()        >= 2 );
    if(!nJets)              return;

    if(DEBUG) std::cout << "nJets Ok... ";

    // Adding invariant mass cut on loose leptons pairs
    bool pass_invariantemasscut = true;
    for(int i=0; i<vLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLeptons.size(); j++)
        {
            if ( fabs( ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M() ) < 12 )
            { pass_invariantemasscut = false ;}
        }
    }
    if(!pass_invariantemasscut) return;    

    if(DEBUG) std::cout << "invariantmasscut Ok... ";

    // ##########
    // # Z peak #
    // ##########

    bool pass_Zpeak = false;
    float Delta = 10, ZM = -1., Zpt = -1., Lep1Z = -1, Lep2Z = -1, LepW = -1;

    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()                               )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < Delta ) )
            {
                Lep1Z      = i;
                Lep2Z      = j;
                ZM         = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M();
                Zpt        = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).Pt();
                Delta      = fabs(ZM - 91.188);
                pass_Zpeak = true;
            }
        }
    }
    if(!pass_Zpeak)       return;

    // 
    if( ( (Lep1Z == 0) && (Lep2Z == 1) ) || ( (Lep1Z == 1) && (Lep2Z == 0) ) ) LepW = 2;
    if( ( (Lep1Z == 0) && (Lep2Z == 2) ) || ( (Lep1Z == 2) && (Lep2Z == 0) ) ) LepW = 1;    
    if( ( (Lep1Z == 1) && (Lep2Z == 2) ) || ( (Lep1Z == 2) && (Lep2Z == 1) ) ) LepW = 0;
    //AC considering cases with 4 leptons, taking lepW w/ highest pT
    if( ( (Lep1Z == 0) && (Lep2Z == 3) ) || ( (Lep1Z == 3) && (Lep2Z == 0) ) ) LepW = 1;
    if( ( (Lep1Z == 1) && (Lep2Z == 3) ) || ( (Lep1Z == 3) && (Lep2Z == 1) ) ) LepW = 0;
    if( ( (Lep1Z == 2) && (Lep2Z == 3) ) || ( (Lep1Z == 3) && (Lep2Z == 2) ) ) LepW = 0;


    float MTW = 0. ;
    if (LepW >=0) MTW = sqrt( 2 * vSelectedLeptons.at(LepW).p4().Pt() * vEvent->at(0).metpt() * (1 - cos( vSelectedLeptons.at(LepW).phi() - vEvent->at(0).metphi() )));

    theHistoManager->fillHisto("CutFlow",                           "noSel", "WZ_CR",   "",   3, weight);
    theHistoManager->fillHisto("CutFlow",                        "PassingZ", "WZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass",        "PassingZ", "WZ_CR",   "",  ZM, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, tau_px = 0, tau_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    TLorentzVector taup4;
    for(int i=0; i<vSelectedTaus.size(); i++)
    {
        taup4.SetPtEtaPhiE(vSelectedTaus.at(i).pt(), vSelectedTaus.at(i).eta(), vSelectedTaus.at(i).phi(), vSelectedTaus.at(i).E());
        tau_px = tau_px + taup4.Px();
        tau_py = tau_py + taup4.Py();
    }

    MHT = sqrt( (jet_px + lepton_px + tau_px) * (jet_px + lepton_px + tau_px) + (jet_py + lepton_py + tau_py) * (jet_py + lepton_py + tau_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "WZ_CR",   "", _sampleName.Data(),   3, weight);
    if( vSelectedJets.size() >= 4  )  theHistoManager->fillHisto("CutFlow",                "WZ_CR",   "", _sampleName.Data(),   4, weight);

    if(DEBUG) std::cout << " MHT =  " << MHT << "MET = " << vEvent->at(0).metpt() << " met_ld = " << met_ld ;

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        for(int j=0; j<vSelectedLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    if( met_ld               > 0.3 ) theHistoManager->fillHisto("CutFlow",                "WZ_CR",   "",   "",   6, weight);
    if( isSFOS                     ) theHistoManager->fillHisto("CutFlow",                "WZ_CR",   "",   "",   7, weight);

    if(DEBUG) std::cout << std::endl << "nJets: " << vSelectedJets.size() << " met_ld: " << met_ld << " isSFOS " << isSFOS << std::endl;

    if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return;

    if(DEBUG) std::cout << "nJets and met_ld Ok... ";

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vSelectedLeptons.at(i).charge();
    }
    if( fabs(sum_charges_3l) != 1 ) return;

    if(DEBUG) std::cout << "sumOfCharges Ok... ";
    
    // ##################################################################################################################################

    // ###########################
    // # b-tagging Scale Factors #
    // ###########################

    std::vector<double> jetPts;
    std::vector<double> jetEtas;
    std::vector<double> jetCSVs;
    std::vector<int>    jetFlavors;
    int iSys = 0;
    double wgt_csv, wgt_csv_def, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf, new_weight;
    
    zero_btagSF    = 1.;        
    zero_btagSF_up = 1.;        
    zero_btagSF_do = 1.;        

    if(!_isdata)
    {
        for(int i=0; i<vSelectedJets.size(); i++)
        {
            //std::cout   << "pt: "       << vSelectedJets.at(i).pt() 
            //            << " eta: "     << vSelectedJets.at(i).eta() 
            //            << " flavor: "  << vSelectedJets.at(i).jet_hadronFlavour() << std::endl;

            float  pt_test              = vSelectedJets.at(i).pt();
            float  eta_test             = vSelectedJets.at(i).eta();
            float  disc_test            = vSelectedJets.at(i).CSVv2();

            float btag_eff = 1.;
            float jet_scalefactor = 1., jet_scalefactor_up = 1., jet_scalefactor_do = 1.;

            if ( vSelectedJets.at(i).jet_hadronFlavour() == 5 )
            {
                btag_eff = get_eff_btagging( pt_test, eta_test, 5);
                jet_scalefactor      = reader.eval_auto_bounds(  "central",  BTagEntryX::FLAV_B,     eta_test,   pt_test );
                jet_scalefactor_up   = reader.eval_auto_bounds(  "up",       BTagEntryX::FLAV_B,     eta_test,   pt_test );
                jet_scalefactor_do   = reader.eval_auto_bounds(  "down",     BTagEntryX::FLAV_B,     eta_test,   pt_test );  
            }

            if ( vSelectedJets.at(i).jet_hadronFlavour() == 4 )
            {
                btag_eff = get_eff_btagging( pt_test, eta_test, 4);
                jet_scalefactor      = reader.eval_auto_bounds(  "central",  BTagEntryX::FLAV_C,     eta_test,   pt_test );
                jet_scalefactor_up   = reader.eval_auto_bounds(  "up",       BTagEntryX::FLAV_C,     eta_test,   pt_test );
                jet_scalefactor_do   = reader.eval_auto_bounds(  "down",     BTagEntryX::FLAV_C,     eta_test,   pt_test );  
            }

            if ( vSelectedJets.at(i).jet_hadronFlavour() == 0 )
            {
                btag_eff = get_eff_btagging( pt_test, eta_test, 0);
                jet_scalefactor      = reader.eval_auto_bounds(  "central",  BTagEntryX::FLAV_UDSG,     eta_test,   pt_test );
                jet_scalefactor_up   = reader.eval_auto_bounds(  "up",       BTagEntryX::FLAV_UDSG,     eta_test,   pt_test );
                jet_scalefactor_do   = reader.eval_auto_bounds(  "down",     BTagEntryX::FLAV_UDSG,     eta_test,   pt_test );  
            }

            zero_btagSF    *= ( 1. - (btag_eff * jet_scalefactor) );
            zero_btagSF_up *= ( 1. - (btag_eff * jet_scalefactor_up) );
            zero_btagSF_do *= ( 1. - (btag_eff * jet_scalefactor_do) );

            //std::cout << "Jet[" << i << "]: pt << " << pt_test << " eta " << eta_test <<  " flavor " << vSelectedJets.at(i).jet_hadronFlavour() << std::endl;
            //std::cout << "btag_eff: " << btag_eff << std::endl;
            //std::cout << "central: " << jet_scalefactor << " up: " << jet_scalefactor_up << " down: " << jet_scalefactor_do << std::endl;
            //std::cout << "prob: " << prob << " btag_SF: " << zero_btagSF << std::endl;
        }
    }

    // weight *= zero_btagSF;
     
    //wgt_csv_def = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
    //new_weight = weight * wgt_csv_def; // weight = weight * wgt_csv_def;

    // ##############
    // # b-jet veto #
    // ##############

    //bool nLooseBtag       = ( nLooseBJets                               == 0 );
    //bool nMediumBtag      = ( nMediumBJets                              == 0 );
    //if(!nLooseBtag || !nMediumBtag)      return;

    theHistoManager->fillHisto("CutFlow",                 "PassingbJetsVeto", "WZ_CR",   "",   1, weight*zero_btagSF);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingbJetsVeto", "WZ_CR",   "",  ZM, weight*zero_btagSF);

    // Building arbitrary variables with all leptons...

    TLorentzVector all_lep_invmass_p4;

    int   all_lep_sumofcharges = 0;
    float all_lep_invmass      = 0;
    float all_lep_sumofpt      = 0;

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        all_lep_sumofcharges = all_lep_sumofcharges + vSelectedLeptons.at(i).charge();
        all_lep_invmass_p4   = all_lep_invmass_p4   + vSelectedLeptons.at(i).p4();
    }

    all_lep_invmass = all_lep_invmass_p4.M();
    all_lep_sumofpt = sqrt( (lepton_px*lepton_px) + (lepton_py*lepton_py) );

    // Building arbitrary variables with remaining (after Z peak determination) leptons...

    TLorentzVector rem_lep_invmass_p4;

    int   rem_lep_sumofcharges = 0;
    float rem_lep_invmass      = 0;
    float rem_lep_sumofpt      = 0;
    float rem_lep_px           = 0;
    float rem_lep_py           = 0;


    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        if( i!=Lep1Z && i!=Lep2Z)
        {
            rem_lep_sumofcharges = rem_lep_sumofcharges + vSelectedLeptons.at(i).charge();
            rem_lep_px           = rem_lep_px           + vSelectedLeptons.at(i).p4().Px();
            rem_lep_py           = rem_lep_py           + vSelectedLeptons.at(i).p4().Py();
            rem_lep_invmass_p4   = rem_lep_invmass_p4   + vSelectedLeptons.at(i).p4();
        }
    }

    rem_lep_invmass = rem_lep_invmass_p4.M();
    rem_lep_sumofpt = sqrt( (rem_lep_px*rem_lep_px) + (rem_lep_py*rem_lep_py) );

   
    // ##################################################################################################################################

    theHistoManager->fillHisto("ZCandidateInvariantMass",        "finalSel", "WZ_CR",   "",   ZM                    , weight*zero_btagSF);
    theHistoManager->fillHisto("ZCandidateTransverseMomentum",   "finalSel", "WZ_CR",   "",   Zpt                   , weight);//ero_btagSF);
    theHistoManager->fillHisto("MET",                            "finalSel", "WZ_CR",   "",   vEvent->at(0).metpt() , weight*zero_btagSF);
    theHistoManager->fillHisto("MHT",                            "finalSel", "WZ_CR",   "",   MHT                   , weight*zero_btagSF);
    theHistoManager->fillHisto("MetLD",                          "finalSel", "WZ_CR",   "",   met_ld                , weight*zero_btagSF);
    theHistoManager->fillHisto("TauMultiplicity",                "finalSel", "WZ_CR",   "",   vSelectedTaus.size()  , weight*zero_btagSF);
    theHistoManager->fillHisto("JetMultiplicity",                "finalSel", "WZ_CR",   "",   vSelectedJets.size()  , weight*zero_btagSF);

    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR",   "",   MTW                   , weight*zero_btagSF);
    if(weights_pdf.size()> 110)
    {
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_1002",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(1));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_1003",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(2));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_1004",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(3));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_1005",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(4));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_1006",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(5));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_1007",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(6));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_1008",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(7));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_1009",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(8));
    
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2001",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(9));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2002",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(10));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2003",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(11));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2004",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(12));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2005",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(13));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2006",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(14));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2007",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(15));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2008",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(16));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2009",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(17));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2010",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(18));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2011",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(19));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2012",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(20));
   
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2013",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(21));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2014",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(22));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2015",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(23));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2016",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(24));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2017",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(25));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2018",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(26));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2019",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(27));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2020",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(28));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2021",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(29));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2022",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(30));
  
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2023",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(31));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2024",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(32));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2025",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(33));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2026",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(34));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2027",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(35));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2028",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(36));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2029",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(37));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2030",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(38));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2031",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(39));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2032",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(40));
    
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2033",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(41));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2034",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(42));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2035",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(43));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2036",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(44));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2037",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(45));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2038",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(46));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2039",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(47));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2040",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(48));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2041",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(49));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2042",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(50));
    
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2043",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(51));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2044",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(52));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2045",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(53));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2046",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(54));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2047",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(55));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2048",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(56));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2049",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(57));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2050",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(58));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2051",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(59));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2052",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(60));
    
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2053",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(61));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2054",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(62));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2055",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(63));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2056",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(64));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2057",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(65));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2058",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(66));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2059",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(67));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2060",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(68));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2061",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(69));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2062",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(70));
   
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2063",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(71));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2064",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(72));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2065",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(73));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2066",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(74));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2067",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(75));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2068",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(76));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2069",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(77));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2070",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(78));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2071",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(79));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2072",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(80));
   
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2073",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(81));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2074",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(82));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2075",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(83));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2076",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(84));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2077",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(85));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2078",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(86));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2079",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(87));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2080",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(88));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2081",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(89));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2082",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(90));
    
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2083",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(91));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2084",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(92));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2085",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(93));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2086",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(94));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2087",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(95));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2088",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(96));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2089",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(97));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2090",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(98));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2091",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(99));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2092",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(100));
  
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2093",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(101));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2094",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(102));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2095",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(103));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2096",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(104));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2097",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(105));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2098",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(106));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2099",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(107));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2100",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(108));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2101",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(109));
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_2102",   "",   MTW              , weight*zero_btagSF*weights_pdf.at(110));
    }
    
    theHistoManager->fillHisto("InvariantMassOfSelectedLeptons", "finalSel", "WZ_CR",   "",   all_lep_invmass       , weight*zero_btagSF);
    theHistoManager->fillHisto("SumOfSelectedLeptonsCharges",    "finalSel", "WZ_CR",   "",   all_lep_sumofcharges  , weight*zero_btagSF);
    theHistoManager->fillHisto("SumVecPtSelectedLeptons",        "finalSel", "WZ_CR",   "",   all_lep_sumofpt       , weight*zero_btagSF);

    theHistoManager->fillHisto("InvMassRemainingLepton",         "finalSel", "WZ_CR",   "",   rem_lep_invmass       , weight*zero_btagSF);
    theHistoManager->fillHisto("SumOfRemainingLeptonsCharges",   "finalSel", "WZ_CR",   "",   rem_lep_sumofcharges  , weight*zero_btagSF);
    theHistoManager->fillHisto("SumVecPtRemainingLeptons",       "finalSel", "WZ_CR",   "",   rem_lep_sumofpt       , weight*zero_btagSF);

    theHistoManager->fillHisto2D("InvMassLastLeptonVSZMass",     "finalSel", "WZ_CR",   "",   rem_lep_invmass,    ZM, weight*zero_btagSF);
    theHistoManager->fillHisto2D("SumPtLepVSZMass",              "finalSel", "WZ_CR",   "",   rem_lep_sumofpt,    ZM, weight*zero_btagSF);
    theHistoManager->fillHisto2D("METLDVSZMass",                 "finalSel", "WZ_CR",   "",   met_ld,             ZM, weight*zero_btagSF);
    
    theHistoManager->fillHisto("JetMultiplicity",                "finalSel", "WZ_CR_BSFUp",   "",   vSelectedJets.size()  , weight*zero_btagSF_up);
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_BSFUp",   "",   MTW                   , weight*zero_btagSF_up);
    theHistoManager->fillHisto("SumOfSelectedLeptonsCharges",    "finalSel", "WZ_CR_BSFUp",   "",   all_lep_sumofcharges  , weight*zero_btagSF_up);
    theHistoManager->fillHisto("ZCandidateInvariantMass",        "finalSel", "WZ_CR_BSFUp",   "",   ZM                    , weight*zero_btagSF_up);
    theHistoManager->fillHisto("ZCandidateTransverseMomentum",   "finalSel", "WZ_CR_BSFUp",   "",   Zpt                   , weight*zero_btagSF_up);
    
    theHistoManager->fillHisto("JetMultiplicity",                "finalSel", "WZ_CR_BSFDo",   "",   vSelectedJets.size()  , weight*zero_btagSF_do);
    theHistoManager->fillHisto("MTW",                            "finalSel", "WZ_CR_BSFDo",   "",   MTW                   , weight*zero_btagSF_do);
    theHistoManager->fillHisto("SumOfSelectedLeptonsCharges",    "finalSel", "WZ_CR_BSFDo",   "",   all_lep_sumofcharges  , weight*zero_btagSF_do);
    theHistoManager->fillHisto("ZCandidateInvariantMass",        "finalSel", "WZ_CR_BSFDo",   "",   ZM                    , weight*zero_btagSF_do);
    theHistoManager->fillHisto("ZCandidateTransverseMomentum",   "finalSel", "WZ_CR_BSFDo",   "",   Zpt                   , weight*zero_btagSF_do);
 
   //std::cout << "LepW "<< LepW <<std::endl;
   //std::cout << "vSelectedLeptons.size() " <<vSelectedLeptons.size() << std::endl;
   //std::cout << "phi " << vSelectedLeptons.at(LepW).phi() <<" "<< vEvent->at(0).metphi()<< std::endl;

   //std::cout <<"evt " << vEvent->at(0).id() <<  std::endl;

    is_3l_WZ_CR = true;
}

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_CR_WZrelaxed(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel", "WZrel_CR",   "",    vLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel", "WZrel_CR",   "",    vLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    if(DEBUG)             std::cout << std::endl << " 3l ss SR ==================================================================" << std::endl;

    bool nLep               = ( vSelectedLeptons.size()     >= 3 );
    if(!nLep)               return;

    if(DEBUG) std::cout << "nLep Ok... ";

    //bool nTight             = ( n_tight                     >= 3 );
    //if(!nTight)             return;

    //bool are_fullytight     =  ( vSelectedLeptons.at(0).isTightTTH()        && vSelectedLeptons.at(1).isTightTTH()      && vSelectedLeptons.at(2).isTightTTH()
    //                          && vSelectedLeptons.at(0).cutEventSel()       && vSelectedLeptons.at(1).cutEventSel()     && vSelectedLeptons.at(2).cutEventSel()
    //                          && vSelectedLeptons.at(0).noLostHits()        && vSelectedLeptons.at(1).noLostHits()      && vSelectedLeptons.at(2).noLostHits()      );
    //if(!are_fullytight)     return;
    if(DEBUG) std::cout << "nTight Ok... ";

    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt() > 20 );
    if(!leading_lep_pt)     return;

    if(DEBUG) std::cout << "leading_lep_pt Ok... ";

    bool following_lep_pt   = ( (vSelectedLeptons.at(1).pt() > 10) && (vSelectedLeptons.at(2).pt() > 10) );
    if(!following_lep_pt)   return;

    if(DEBUG) std::cout << "following_lep_pt Ok... ";

    bool nJets              = ( vSelectedJets.size()        >= 2 );
    if(!nJets)              return;

    if(DEBUG) std::cout << "nJets Ok... ";

    // Adding invariant mass cut on loose leptons pairs
    bool pass_invariantemasscut = true;
    for(int i=0; i<vLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLeptons.size(); j++)
        {
            if ( fabs( ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M() ) < 12 )
            { pass_invariantemasscut = false ;}
        }
    }
    if(!pass_invariantemasscut) return;    

    if(DEBUG) std::cout << "invariantmasscut Ok... ";

    // ##########
    // # Z peak #
    // ##########

    bool pass_Zpeak = false;
    float Delta = 10, ZM = -1., Zpt = -1., Lep1Z = -1, Lep2Z = -1, LepW = -1;

    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()                               )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < Delta ) )
            {
                Lep1Z      = i;
                Lep2Z      = j;
                ZM         = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M();
                Zpt        = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).Pt();
                Delta      = fabs(ZM - 91.188);
                pass_Zpeak = true;
            }
        }
    }
    if(!pass_Zpeak)       return;

    // 
    if( ( (Lep1Z == 0) && (Lep2Z == 1) ) || ( (Lep1Z == 1) && (Lep2Z == 0) ) ) LepW = 2;
    if( ( (Lep1Z == 0) && (Lep2Z == 2) ) || ( (Lep1Z == 2) && (Lep2Z == 0) ) ) LepW = 1;    
    if( ( (Lep1Z == 1) && (Lep2Z == 2) ) || ( (Lep1Z == 2) && (Lep2Z == 1) ) ) LepW = 0;
    //AC considering cases with 4 leptons, taking lepW w/ highest pT
    if( ( (Lep1Z == 0) && (Lep2Z == 3) ) || ( (Lep1Z == 3) && (Lep2Z == 0) ) ) LepW = 1;
    if( ( (Lep1Z == 1) && (Lep2Z == 3) ) || ( (Lep1Z == 3) && (Lep2Z == 1) ) ) LepW = 0;
    if( ( (Lep1Z == 2) && (Lep2Z == 3) ) || ( (Lep1Z == 3) && (Lep2Z == 2) ) ) LepW = 0;


    float MTW = 0. ;
    if (LepW >=0) MTW = sqrt( 2 * vSelectedLeptons.at(LepW).p4().Pt() * vEvent->at(0).metpt() * (1 - cos( vSelectedLeptons.at(LepW).phi() - vEvent->at(0).metphi() )));

    theHistoManager->fillHisto("CutFlow",                           "noSel", "WZ_CR",   "",   3, weight);
    theHistoManager->fillHisto("CutFlow",                        "PassingZ", "WZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass",        "PassingZ", "WZ_CR",   "",  ZM, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, tau_px = 0, tau_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    TLorentzVector taup4;
    for(int i=0; i<vSelectedTaus.size(); i++)
    {
        taup4.SetPtEtaPhiE(vSelectedTaus.at(i).pt(), vSelectedTaus.at(i).eta(), vSelectedTaus.at(i).phi(), vSelectedTaus.at(i).E());
        tau_px = tau_px + taup4.Px();
        tau_py = tau_py + taup4.Py();
    }

    MHT = sqrt( (jet_px + lepton_px + tau_px) * (jet_px + lepton_px + tau_px) + (jet_py + lepton_py + tau_py) * (jet_py + lepton_py + tau_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "WZ_CR",   "", _sampleName.Data(),   3, weight);
    if( vSelectedJets.size() >= 4  )  theHistoManager->fillHisto("CutFlow",                "WZ_CR",   "", _sampleName.Data(),   4, weight);

    if(DEBUG) std::cout << " MHT =  " << MHT << "MET = " << vEvent->at(0).metpt() << " met_ld = " << met_ld ;

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        for(int j=0; j<vSelectedLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    if( met_ld               > 0.3 ) theHistoManager->fillHisto("CutFlow",                "WZ_CR",   "",   "",   6, weight);
    if( isSFOS                     ) theHistoManager->fillHisto("CutFlow",                "WZ_CR",   "",   "",   7, weight);

    if(DEBUG) std::cout << std::endl << "nJets: " << vSelectedJets.size() << " met_ld: " << met_ld << " isSFOS " << isSFOS << std::endl;

    if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return;

    if(DEBUG) std::cout << "nJets and met_ld Ok... ";

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vSelectedLeptons.at(i).charge();
    }
    if( fabs(sum_charges_3l) != 1 ) return;

    if(DEBUG) std::cout << "sumOfCharges Ok... ";

    // ##############
    // # b-jet veto #
    // ##############

    bool nLooseBtag       = ( nLooseBJets                               == 0 );
    bool nMediumBtag      = ( nMediumBJets                              == 0 );
    if(!nLooseBtag || !nMediumBtag)      return;

    theHistoManager->fillHisto("CutFlow",                 "PassingbJetsVeto", "WZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingbJetsVeto", "WZ_CR",   "",  ZM, weight);

    // Building arbitrary variables with all leptons...

    TLorentzVector all_lep_invmass_p4;

    int   all_lep_sumofcharges = 0;
    float all_lep_invmass      = 0;
    float all_lep_sumofpt      = 0;

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        all_lep_sumofcharges = all_lep_sumofcharges + vSelectedLeptons.at(i).charge();
        all_lep_invmass_p4   = all_lep_invmass_p4   + vSelectedLeptons.at(i).p4();
    }

    all_lep_invmass = all_lep_invmass_p4.M();
    all_lep_sumofpt = sqrt( (lepton_px*lepton_px) + (lepton_py*lepton_py) );

    // Building arbitrary variables with remaining (after Z peak determination) leptons...

    TLorentzVector rem_lep_invmass_p4;

    int   rem_lep_sumofcharges = 0;
    float rem_lep_invmass      = 0;
    float rem_lep_sumofpt      = 0;
    float rem_lep_px           = 0;
    float rem_lep_py           = 0;


    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        if( i!=Lep1Z && i!=Lep2Z)
        {
            rem_lep_sumofcharges = rem_lep_sumofcharges + vSelectedLeptons.at(i).charge();
            rem_lep_px           = rem_lep_px           + vSelectedLeptons.at(i).p4().Px();
            rem_lep_py           = rem_lep_py           + vSelectedLeptons.at(i).p4().Py();
            rem_lep_invmass_p4   = rem_lep_invmass_p4   + vSelectedLeptons.at(i).p4();
        }
    }

    rem_lep_invmass = rem_lep_invmass_p4.M();
    rem_lep_sumofpt = sqrt( (rem_lep_px*rem_lep_px) + (rem_lep_py*rem_lep_py) );

    theHistoManager->fillHisto("ZCandidateInvariantMass",        "finalSel", "WZrel_CR",   "",   ZM                    , weight);
    theHistoManager->fillHisto("ZCandidateTransverseMomentum",   "finalSel", "WZrel_CR",   "",   Zpt                   , weight);
    theHistoManager->fillHisto("MET",                            "finalSel", "WZrel_CR",   "",   vEvent->at(0).metpt() , weight);
    theHistoManager->fillHisto("MHT",                            "finalSel", "WZrel_CR",   "",   MHT                   , weight);
    theHistoManager->fillHisto("MetLD",                          "finalSel", "WZrel_CR",   "",   met_ld                , weight);
    theHistoManager->fillHisto("TauMultiplicity",                "finalSel", "WZrel_CR",   "",   vSelectedTaus.size()  , weight);
    theHistoManager->fillHisto("JetMultiplicity",                "finalSel", "WZrel_CR",   "",   vSelectedJets.size()  , weight);

    theHistoManager->fillHisto("InvariantMassOfSelectedLeptons", "finalSel", "WZrel_CR",   "",   all_lep_invmass       , weight);
    theHistoManager->fillHisto("SumOfSelectedLeptonsCharges",    "finalSel", "WZrel_CR",   "",   all_lep_sumofcharges  , weight);
    theHistoManager->fillHisto("SumVecPtSelectedLeptons",        "finalSel", "WZrel_CR",   "",   all_lep_sumofpt       , weight);

    theHistoManager->fillHisto("InvMassRemainingLepton",         "finalSel", "WZrel_CR",   "",   rem_lep_invmass       , weight);
    theHistoManager->fillHisto("SumOfRemainingLeptonsCharges",   "finalSel", "WZrel_CR",   "",   rem_lep_sumofcharges  , weight);
    theHistoManager->fillHisto("SumVecPtRemainingLeptons",       "finalSel", "WZrel_CR",   "",   rem_lep_sumofpt       , weight);

    theHistoManager->fillHisto2D("InvMassLastLeptonVSZMass",     "finalSel", "WZrel_CR",   "",   rem_lep_invmass,    ZM, weight);
    theHistoManager->fillHisto2D("SumPtLepVSZMass",              "finalSel", "WZrel_CR",   "",   rem_lep_sumofpt,    ZM, weight);
    theHistoManager->fillHisto2D("METLDVSZMass",                 "finalSel", "WZrel_CR",   "",   met_ld,             ZM, weight);

    is_3l_WZrel_CR = true;
}

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_TTZ(int evt)
{

    theHistoManager->fillHisto2D("LeptonsVsJets",           "noSel", "TTZ_CR",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("LeptonsVsBJets",          "noSel", "TTZ_CR",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    if(DEBUG)             std::cout << std::endl << " 3l ss SR ==================================================================" << std::endl;

    bool nLep               = ( vSelectedLeptons.size()     >= 3 );
    if(!nLep)               return;

    if(DEBUG) std::cout << "nLep Ok... ";

    bool nTight             = ( n_tight                     >= 3 );
    if(!nTight)             return;

    bool are_fullytight     =  ( vSelectedLeptons.at(0).isTightTTH()        && vSelectedLeptons.at(1).isTightTTH()      && vSelectedLeptons.at(2).isTightTTH()
                              && vSelectedLeptons.at(0).cutEventSel()       && vSelectedLeptons.at(1).cutEventSel()     && vSelectedLeptons.at(2).cutEventSel()
                              && vSelectedLeptons.at(0).noLostHits()        && vSelectedLeptons.at(1).noLostHits()      && vSelectedLeptons.at(2).noLostHits()      );
    if(!are_fullytight)     return;
    if(DEBUG) std::cout << "nTight Ok... ";

    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt() > 20 );
    if(!leading_lep_pt)     return;

    if(DEBUG) std::cout << "leading_lep_pt Ok... ";

    bool following_lep_pt   = ( (vSelectedLeptons.at(1).pt() > 10) && (vSelectedLeptons.at(2).pt() > 10) );
    if(!following_lep_pt)   return;

    if(DEBUG) std::cout << "following_lep_pt Ok... ";

    bool nJets              = ( vSelectedJets.size()        >= 2 );
    if(!nJets)              return;

    if(DEBUG) std::cout << "nJets Ok... ";

    // Adding invariant mass cut on loose leptons pairs
    bool pass_invariantemasscut = true;
    for(int i=0; i<vLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLeptons.size(); j++)
        {
            if ( fabs( ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M() ) < 12 )
            { pass_invariantemasscut = false ;}
        }
    }
    if(!pass_invariantemasscut) return;    

    if(DEBUG) std::cout << "invariantmasscut Ok... ";

    // ##########
    // # Z peak #
    // ##########

    bool pass_Zpeak = false;
    float Delta = 10, ZM = -1., Zpt = -1., Lep1Z = -1, Lep2Z = -1, LepW = -1;

    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()                               )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < Delta ) )
            {
                Lep1Z      = i;
                Lep2Z      = j;
                ZM         = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M();
                Zpt        = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).Pt();
                Delta      = fabs(ZM - 91.188);
                pass_Zpeak = true;
            }
        }
    }
    if(!pass_Zpeak)       return;

    // 
    if( ( (Lep1Z == 0) && (Lep2Z == 1) ) || ( (Lep1Z == 1) && (Lep2Z == 0) ) ) LepW = 2;
    if( ( (Lep1Z == 0) && (Lep2Z == 2) ) || ( (Lep1Z == 2) && (Lep2Z == 0) ) ) LepW = 1;    
    if( ( (Lep1Z == 1) && (Lep2Z == 2) ) || ( (Lep1Z == 2) && (Lep2Z == 1) ) ) LepW = 0;
    //AC considering cases with 4 leptons, taking lepW w/ highest pT
    if( ( (Lep1Z == 0) && (Lep2Z == 3) ) || ( (Lep1Z == 3) && (Lep2Z == 0) ) ) LepW = 1;
    if( ( (Lep1Z == 1) && (Lep2Z == 3) ) || ( (Lep1Z == 3) && (Lep2Z == 1) ) ) LepW = 0;
    if( ( (Lep1Z == 2) && (Lep2Z == 3) ) || ( (Lep1Z == 3) && (Lep2Z == 2) ) ) LepW = 0;


    float MTW = 0. ;
    if (LepW >=0) MTW = sqrt( 2 * vSelectedLeptons.at(LepW).p4().Pt() * vEvent->at(0).metpt() * (1 - cos( vSelectedLeptons.at(LepW).phi() - vEvent->at(0).metphi() )));

    theHistoManager->fillHisto("CutFlow",                           "noSel", "WZ_CR",   "",   3, weight);
    theHistoManager->fillHisto("CutFlow",                        "PassingZ", "WZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass",        "PassingZ", "WZ_CR",   "",  ZM, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, tau_px = 0, tau_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    TLorentzVector taup4;
    for(int i=0; i<vSelectedTaus.size(); i++)
    {
        taup4.SetPtEtaPhiE(vSelectedTaus.at(i).pt(), vSelectedTaus.at(i).eta(), vSelectedTaus.at(i).phi(), vSelectedTaus.at(i).E());
        tau_px = tau_px + taup4.Px();
        tau_py = tau_py + taup4.Py();
    }

    MHT = sqrt( (jet_px + lepton_px + tau_px) * (jet_px + lepton_px + tau_px) + (jet_py + lepton_py + tau_py) * (jet_py + lepton_py + tau_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "WZ_CR",   "", _sampleName.Data(),   3, weight);
    if( vSelectedJets.size() >= 4  )  theHistoManager->fillHisto("CutFlow",                "WZ_CR",   "", _sampleName.Data(),   4, weight);

    if(DEBUG) std::cout << " MHT =  " << MHT << "MET = " << vEvent->at(0).metpt() << " met_ld = " << met_ld ;

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        for(int j=0; j<vSelectedLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    if( met_ld               > 0.3 ) theHistoManager->fillHisto("CutFlow",                "WZ_CR",   "",   "",   6, weight);
    if( isSFOS                     ) theHistoManager->fillHisto("CutFlow",                "WZ_CR",   "",   "",   7, weight);

    if(DEBUG) std::cout << std::endl << "nJets: " << vSelectedJets.size() << " met_ld: " << met_ld << " isSFOS " << isSFOS << std::endl;

    if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return;

    if(DEBUG) std::cout << "nJets and met_ld Ok... ";

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vSelectedLeptons.at(i).charge();
    }
    if( fabs(sum_charges_3l) != 1 ) return;

    if(DEBUG) std::cout << "sumOfCharges Ok... ";

    // #########
    // # b-jet #
    // #########

    bool nLooseBtag       = ( nLooseBJets                               >= 2 );
    bool nMediumBtag      = ( nMediumBJets                              >= 1 );
    if(!nLooseBtag || !nMediumBtag)      return;

    theHistoManager->fillHisto("CutFlow",                 "PassingbJetsVeto", "WZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingbJetsVeto", "WZ_CR",   "",  ZM, weight);

    // Building arbitrary variables with all leptons...

    TLorentzVector all_lep_invmass_p4;

    int   all_lep_sumofcharges = 0;
    float all_lep_invmass      = 0;
    float all_lep_sumofpt      = 0;

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        all_lep_sumofcharges = all_lep_sumofcharges + vSelectedLeptons.at(i).charge();
        all_lep_invmass_p4   = all_lep_invmass_p4   + vSelectedLeptons.at(i).p4();
    }

    all_lep_invmass = all_lep_invmass_p4.M();
    all_lep_sumofpt = sqrt( (lepton_px*lepton_px) + (lepton_py*lepton_py) );

    // Building arbitrary variables with remaining (after Z peak determination) leptons...

    TLorentzVector rem_lep_invmass_p4;

    int   rem_lep_sumofcharges = 0;
    float rem_lep_invmass      = 0;
    float rem_lep_sumofpt      = 0;
    float rem_lep_px           = 0;
    float rem_lep_py           = 0;


    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        if( i!=Lep1Z && i!=Lep2Z)
        {
            rem_lep_sumofcharges = rem_lep_sumofcharges + vSelectedLeptons.at(i).charge();
            rem_lep_px           = rem_lep_px           + vSelectedLeptons.at(i).p4().Px();
            rem_lep_py           = rem_lep_py           + vSelectedLeptons.at(i).p4().Py();
            rem_lep_invmass_p4   = rem_lep_invmass_p4   + vSelectedLeptons.at(i).p4();
        }
    }

    rem_lep_invmass = rem_lep_invmass_p4.M();
    rem_lep_sumofpt = sqrt( (rem_lep_px*rem_lep_px) + (rem_lep_py*rem_lep_py) );

    theHistoManager->fillHisto("LeadingLeptonPt",                "finalSel", "TTZ_CR",   "",   vSelectedLeptons.at(0).pt()   , weight);
    theHistoManager->fillHisto("SubLeadingLeptonPt",             "finalSel", "TTZ_CR",   "",   vSelectedLeptons.at(1).pt()   , weight);    

    theHistoManager->fillHisto("ZCandidateInvariantMass",        "finalSel", "TTZ_CR",   "",   ZM                    , weight);
    theHistoManager->fillHisto("ZCandidateTransverseMomentum",   "finalSel", "TTZ_CR",   "",   Zpt                   , weight);
    theHistoManager->fillHisto("MET",                            "finalSel", "TTZ_CR",   "",   vEvent->at(0).metpt() , weight);
    theHistoManager->fillHisto("MHT",                            "finalSel", "TTZ_CR",   "",   MHT                   , weight);
    theHistoManager->fillHisto("MetLD",                          "finalSel", "TTZ_CR",   "",   met_ld                , weight);
    theHistoManager->fillHisto("TauMultiplicity",                "finalSel", "TTZ_CR",   "",   vSelectedTaus.size()  , weight);
    theHistoManager->fillHisto("JetMultiplicity",                "finalSel", "TTZ_CR",   "",   vSelectedJets.size()  , weight);

    theHistoManager->fillHisto("InvariantMassOfSelectedLeptons", "finalSel", "TTZ_CR",   "",   all_lep_invmass       , weight);
    theHistoManager->fillHisto("SumOfSelectedLeptonsCharges",    "finalSel", "TTZ_CR",   "",   all_lep_sumofcharges  , weight);
    theHistoManager->fillHisto("SumVecPtSelectedLeptons",        "finalSel", "TTZ_CR",   "",   all_lep_sumofpt       , weight);

    theHistoManager->fillHisto("InvMassRemainingLepton",         "finalSel", "TTZ_CR",   "",   rem_lep_invmass       , weight);
    theHistoManager->fillHisto("SumOfRemainingLeptonsCharges",   "finalSel", "TTZ_CR",   "",   rem_lep_sumofcharges  , weight);
    theHistoManager->fillHisto("SumVecPtRemainingLeptons",       "finalSel", "TTZ_CR",   "",   rem_lep_sumofpt       , weight);

    theHistoManager->fillHisto2D("InvMassLastLeptonVSZMass",     "finalSel", "TTZ_CR",   "",   rem_lep_invmass,    ZM, weight);
    theHistoManager->fillHisto2D("SumPtLepVSZMass",              "finalSel", "TTZ_CR",   "",   rem_lep_sumofpt,    ZM, weight);
    theHistoManager->fillHisto2D("METLDVSZMass",                 "finalSel", "TTZ_CR",   "",   met_ld,             ZM, weight);

    if(vSelectedJets.size() >= 4)
    {
        theHistoManager->fillHisto("LeadingLeptonPt",                "finalSel4j", "TTZ_CR",   "",   vSelectedLeptons.at(0).pt()  , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",             "finalSel4j", "TTZ_CR",   "",   vSelectedLeptons.at(1).pt()  , weight);

        theHistoManager->fillHisto("MET",                            "finalSel4j", "TTZ_CR",   "",   vEvent->at(0).metpt() , weight);
        theHistoManager->fillHisto("JetMultiplicity",                "finalSel4j", "TTZ_CR",   "",   vSelectedJets.size()  , weight);
        theHistoManager->fillHisto("ZCandidateInvariantMass",        "finalSel4j", "TTZ_CR",   "",   ZM                    , weight);
    }

    is_3l_TTZ_CR = true;
}

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_CR_Zl(int evt)
{

}

bool TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_TTH3l_MC() 
{ 
    bool sel_MC = true;

    //Check decays 
    if (!((vTruth->at(0).ttbar_decay()==1 && vTruth->at(0).boson_decay()==0) || //ttH
                (vTruth->at(0).ttbar_decay()==2 && vTruth->at(0).boson_decay()==1) || //ttH
                (vTruth->at(0).ttbar_decay()==1 && vTruth->at(0).boson_decay()==2) || //tt semi-lep, ttZ
                (vTruth->at(0).ttbar_decay()==2 && vTruth->at(0).boson_decay()==3)    //tt di-lep, ttW
         )) 
    { 
        sel_MC = false; 
        return sel_MC;}

        // 
        if (vTruth->at(0).JetsHighestPt_phi().size()<2 || vTruth->at(0).JetsClosestMw_phi().size()<2 || vTruth->at(0).JetsLowestMjj_phi().size()<2) 
        { 
            sel_MC = false; 
            return sel_MC;}

            //pt, eta of leptons	
            if (!(vTruth->at(0).Leptons_pt().at(0)> 10 && 
                        vTruth->at(0).Leptons_pt().at(1)> 10 &&
                        vTruth->at(0).Leptons_pt().at(2)> 10   )) sel_MC = false; 


            if (!(fabs(vTruth->at(0).Leptons_eta().at(0)) <2.5 && 
                        fabs(vTruth->at(0).Leptons_eta().at(1)) <2.5 &&
                        fabs(vTruth->at(0).Leptons_eta().at(2)) <2.5   )) sel_MC = false; 


            //lead. lepton
            if (!(vTruth->at(0).Leptons_pt().at(0) > 20 || 
                        vTruth->at(0).Leptons_pt().at(1) > 20 ||
                        vTruth->at(0).Leptons_pt().at(2) > 20  )) sel_MC = false; 


            //SFOS && M(ll) not in 81-101 ??? 
            int SFOSpair = -1;
            if ((vTruth->at(0).Leptons_id().at(0)==-vTruth->at(0).Leptons_id().at(1))) SFOSpair = 0;
            if ((vTruth->at(0).Leptons_id().at(1)==-vTruth->at(0).Leptons_id().at(2))) SFOSpair = 1;
            if ((vTruth->at(0).Leptons_id().at(0)==-vTruth->at(0).Leptons_id().at(2))) SFOSpair = 2;


            TLorentzVector Lep1;
            Lep1.SetPtEtaPhiE(vTruth->at(0).Leptons_pt().at(0),  vTruth->at(0).Leptons_eta().at(0), vTruth->at(0).Leptons_phi().at(0), vTruth->at(0).Leptons_E().at(0));
            TLorentzVector Lep2;
            Lep2.SetPtEtaPhiE(vTruth->at(0).Leptons_pt().at(1),  vTruth->at(0).Leptons_eta().at(1), vTruth->at(0).Leptons_phi().at(1), vTruth->at(0).Leptons_E().at(1));
            TLorentzVector Lep3;
            Lep3.SetPtEtaPhiE(vTruth->at(0).Leptons_pt().at(2),  vTruth->at(0).Leptons_eta().at(2), vTruth->at(0).Leptons_phi().at(2), vTruth->at(0).Leptons_E().at(2));


            if ( !(( Lep1+Lep2 ).M()  > 12 && ( Lep1+Lep3 ).M()  > 12 && ( Lep2+Lep3 ).M()  > 12 )) sel_MC = false;

            if ( (SFOSpair == 0 && fabs( (Lep1+Lep2 ).M()-91.188 ) < 10. ) || 
                    (SFOSpair == 1 && fabs( (Lep2+Lep3 ).M()-91.188 ) < 10. ) ||  
                    (SFOSpair == 2 && fabs( (Lep1+Lep3 ).M()-91.188 ) < 10. )    ) sel_MC = false;

            return sel_MC;

}

void TTbarHiggsMultileptonAnalysis::initializeOutputTree()
{

    outputfile->cd();
    tOutput = new TTree("Tree", "Tree");

    tOutput->Branch("mc_event",&mc_event,"mc_event/I");
    tOutput->Branch("weight",&weight,"weight/F");
    tOutput->Branch("weightfake",&weightfake,"weightfake/F");
    tOutput->Branch("weightflip",&weightflip,"weightflip/F");
    tOutput->Branch("mc_weight",&mc_weight,"mc_weight/F"); 
    tOutput->Branch("weight_scale_muF0p5",&weight_scale_muF0p5,"weight_scale_muF0p5/F");
    tOutput->Branch("weight_scale_muF2",&weight_scale_muF2,"weight_scale_muF2/F");
    tOutput->Branch("weight_scale_muR0p5",&weight_scale_muR0p5,"weight_scale_muR0p5/F");
    tOutput->Branch("weight_scale_muR2",&weight_scale_muR2,"weight_scale_muR2/F");
    tOutput->Branch("weight_csv_down",&weight_csv_down,"weight_csv_down/F");
    tOutput->Branch("weight_csv_up",&weight_csv_up,"weight_csv_up/F");
    tOutput->Branch("weights_pdf","std::vector<float>",&weights_pdf);
    tOutput->Branch("ids_pdf","std::vector<std::string>",&ids_pdf);

    tOutput->Branch("PV_weight",&weight_PV,"PV_weight/F");
    tOutput->Branch("mc_3l_category",&mc_3l_category,"mc_3l_category/I");
    tOutput->Branch("mc_ttbar_decay",&mc_ttbar_decay,"mc_ttbar_decay/I");
    tOutput->Branch("mc_boson_decay",&mc_boson_decay,"mc_boson_decay/I");
    tOutput->Branch("mc_ttZhypAllowed",&mc_ttZhypAllowed,"mc_ttZhypAllowed/I");
    tOutput->Branch("mc_nJets25",&mc_nJets25,"mc_nJets25/I");
    tOutput->Branch("mc_nBtagJets25",&mc_nBtagJets25,"mc_nBtagJets25/I");
    tOutput->Branch("mc_nMediumBtagJets25",&mc_nMediumBtagJets25,"mc_nMediumBtagJets25/I");
    tOutput->Branch("mc_nNonBtagJets25",&mc_nNonBtagJets25,"mc_nNonBtagJets25/I");

    tOutput->Branch("catJets",&catJets,"catJets/I");

    tOutput->Branch("is_2lss_TTH_SR",&is_2lss_TTH_SR,"is_2lss_TTH_SR/B");
    tOutput->Branch("is_3l_TTH_SR",&is_3l_TTH_SR,"is_3l_TTH_SR/B");

    tOutput->Branch("is_emu_TT_CR",&is_emu_TT_CR,"is_emu_TT_CR/B");
    tOutput->Branch("is_3l_WZrel_CR",&is_3l_WZrel_CR,"is_3l_WZrel_CR/B");
    tOutput->Branch("is_3l_TTZ_CR",&is_3l_TTZ_CR,"is_3l_TTZ_CR/B");
     
    tOutput->Branch("is_2bTight",&is_2bTight,"is_2bTight/I");

    tOutput->Branch("cat_ee",&cat_ee,"cat_ee/I");
    tOutput->Branch("cat_ee_fake",&cat_ee_fake,"cat_ee_fake/I");
    tOutput->Branch("cat_ee_flip",&cat_ee_flip,"cat_ee_flip/I");
    tOutput->Branch("cat_em",&cat_em,"cat_em/I");
    tOutput->Branch("cat_em_fake",&cat_em_fake,"cat_em_fake/I");
    tOutput->Branch("cat_em_flip",&cat_em_flip,"cat_em_flip/I");
    tOutput->Branch("cat_mm",&cat_mm,"cat_mm/I");
    tOutput->Branch("cat_mm_fake",&cat_mm_fake,"cat_mm_fake/I");
    tOutput->Branch("cat_2ltau",&cat_2ltau,"cat_2ltau/I");
    tOutput->Branch("cat_3l",&cat_3l,"cat_3l/I");
    tOutput->Branch("cat_3l_fake",&cat_3l_fake,"cat_3l_fake/I");

    tOutput->Branch("cat_HtoWW",&cat_HtoWW,"cat_HtoWW/B");
    tOutput->Branch("cat_HtoZZ",&cat_HtoZZ,"cat_HtoZZ/B");
    tOutput->Branch("cat_Htott",&cat_Htott,"cat_Htott/B");

    tOutput->Branch("is_trigger",&is_trigger,"is_trigger/B");

    tOutput->Branch("max_Lep_eta", &max_Lep_eta, "max_Lep_eta/F");
    tOutput->Branch("MT_met_lep1",&MT_met_lep1,"MT_met_lep1/F");
    tOutput->Branch("nJet25_Recl",&nJet25_Recl,"nJet25_Recl/F");
    tOutput->Branch("mindr_lep1_jet",&mindr_lep1_jet,"mindr_lep1_jet/F");
    tOutput->Branch("mindr_lep2_jet",&mindr_lep2_jet,"mindr_lep2_jet/F");
    tOutput->Branch("LepGood_conePt0",&LepGood_conePt0,"LepGood_conePt0/F");
    tOutput->Branch("LepGood_conePt1",&LepGood_conePt1,"LepGood_conePt1/F");
    tOutput->Branch("met",&met,"met/F");
    tOutput->Branch("avg_dr_jet",&avg_dr_jet,"avg_dr_jet/F");

    tOutput->Branch("signal_2lss_TT_MVA",&signal_2lss_TT_MVA,"signal_2lss_TT_MVA/F");
    tOutput->Branch("signal_2lss_TTV_MVA",&signal_2lss_TTV_MVA,"signal_2lss_TTV_MVA/F");

    tOutput->Branch("signal_3l_TT_MVA",&signal_3l_TT_MVA,"signal_3l_TT_MVA/F");
    tOutput->Branch("signal_3l_TTV_MVA",&signal_3l_TTV_MVA,"signal_3l_TTV_MVA/F");

    tOutput->Branch("multilepton_Lepton1_Id",               &multilepton_Lepton1_Id,                "multilepton_Lepton1_Id/I");
    tOutput->Branch("multilepton_Lepton1_P4",               "TLorentzVector",                       &multilepton_Lepton1_P4);
    tOutput->Branch("multilepton_Lepton1_DeltaR_Matched",   &multilepton_Lepton1_DeltaR_Matched,    "multilepton_Lepton1_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Lepton1_Label_Matched",    &multilepton_Lepton1_Label_Matched,     "multilepton_Lepton1_Label_Matched/I");
    tOutput->Branch("multilepton_Lepton1_Id_Matched",       &multilepton_Lepton1_Id_Matched,        "multilepton_Lepton1_Id_Matched/I");
    tOutput->Branch("multilepton_Lepton1_P4_Matched",       "TLorentzVector",                       &multilepton_Lepton1_P4_Matched);
    tOutput->Branch("multilepton_Lepton2_Id",               &multilepton_Lepton2_Id,                "multilepton_Lepton2_Id/I");
    tOutput->Branch("multilepton_Lepton2_P4",               "TLorentzVector",                       &multilepton_Lepton2_P4);
    tOutput->Branch("multilepton_Lepton2_DeltaR_Matched",   &multilepton_Lepton2_DeltaR_Matched,    "multilepton_Lepton2_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Lepton2_Label_Matched",    &multilepton_Lepton2_Label_Matched,     "multilepton_Lepton2_Label_Matched/I");
    tOutput->Branch("multilepton_Lepton2_Id_Matched",       &multilepton_Lepton2_Id_Matched,        "multilepton_Lepton2_Id_Matched/I");
    tOutput->Branch("multilepton_Lepton2_P4_Matched",       "TLorentzVector",                       &multilepton_Lepton2_P4_Matched);
    tOutput->Branch("multilepton_Lepton3_Id",               &multilepton_Lepton3_Id,                "multilepton_Lepton3_Id/I");
    tOutput->Branch("multilepton_Lepton3_P4",               "TLorentzVector",                       &multilepton_Lepton3_P4);
    tOutput->Branch("multilepton_Lepton3_DeltaR_Matched",   &multilepton_Lepton3_DeltaR_Matched,    "multilepton_Lepton3_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Lepton3_Label_Matched",    &multilepton_Lepton3_Label_Matched,     "multilepton_Lepton3_Label_Matched/I");
    tOutput->Branch("multilepton_Lepton3_Id_Matched",       &multilepton_Lepton3_Id_Matched,        "multilepton_Lepton3_Id_Matched/I");
    tOutput->Branch("multilepton_Lepton3_P4_Matched",       "TLorentzVector",                       &multilepton_Lepton3_P4_Matched);
    tOutput->Branch("multilepton_Lepton4_Id",               &multilepton_Lepton4_Id,                "multilepton_Lepton4_Id/I");
    tOutput->Branch("multilepton_Lepton4_P4",               "TLorentzVector",                       &multilepton_Lepton4_P4);
    tOutput->Branch("multilepton_Lepton4_DeltaR_Matched",   &multilepton_Lepton4_DeltaR_Matched,    "multilepton_Lepton4_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Lepton4_Label_Matched",    &multilepton_Lepton4_Label_Matched,     "multilepton_Lepton4_Label_Matched/I");
    tOutput->Branch("multilepton_Lepton4_Id_Matched",       &multilepton_Lepton4_Id_Matched,        "multilepton_Lepton4_Id_Matched/I");
    tOutput->Branch("multilepton_Lepton4_P4_Matched",       "TLorentzVector",                       &multilepton_Lepton4_P4_Matched);

    tOutput->Branch("multilepton_Bjet1_Id",                 &multilepton_Bjet1_Id,                  "multilepton_Bjet1_Id/I");
    tOutput->Branch("multilepton_Bjet1_P4",                 "TLorentzVector",                       &multilepton_Bjet1_P4);
    tOutput->Branch("multilepton_Bjet1_CSV",                &multilepton_Bjet1_CSV,                 "multilepton_Bjet1_CSV/F");
    tOutput->Branch("multilepton_Bjet1_JEC_Up",             &multilepton_Bjet1_JEC_Up,              "multilepton_Bjet1_JEC_Up/F");
    tOutput->Branch("multilepton_Bjet1_JEC_Down",           &multilepton_Bjet1_JEC_Down,            "multilepton_Bjet1_JEC_Down/F");
    tOutput->Branch("multilepton_Bjet1_JER_Up",             &multilepton_Bjet1_JER_Up,              "multilepton_Bjet1_JER_Up/F");
    tOutput->Branch("multilepton_Bjet1_JER_Down",           &multilepton_Bjet1_JER_Down,            "multilepton_Bjet1_JER_Down/F");
    tOutput->Branch("multilepton_Bjet1_DeltaR_Matched",     &multilepton_Bjet1_DeltaR_Matched,      "multilepton_Bjet1_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Bjet1_Label_Matched",      &multilepton_Bjet1_Label_Matched,       "multilepton_Bjet1_Label_Matched/I");
    tOutput->Branch("multilepton_Bjet1_Id_Matched",         &multilepton_Bjet1_Id_Matched,          "multilepton_Bjet1_Id_Matched/I");
    tOutput->Branch("multilepton_Bjet1_P4_Matched",         "TLorentzVector",                       &multilepton_Bjet1_P4_Matched);

    tOutput->Branch("multilepton_Bjet2_Id",                 &multilepton_Bjet2_Id,                  "multilepton_Bjet2_Id/I");
    tOutput->Branch("multilepton_Bjet2_P4",                 "TLorentzVector",                       &multilepton_Bjet2_P4);
    tOutput->Branch("multilepton_Bjet2_CSV",                &multilepton_Bjet2_CSV,                 "multilepton_Bjet2_CSV/F");
    tOutput->Branch("multilepton_Bjet2_JEC_Up",             &multilepton_Bjet2_JEC_Up,              "multilepton_Bjet2_JEC_Up/F");
    tOutput->Branch("multilepton_Bjet2_JEC_Down",           &multilepton_Bjet2_JEC_Down,            "multilepton_Bjet2_JEC_Down/F");
    tOutput->Branch("multilepton_Bjet2_JER_Up",             &multilepton_Bjet2_JER_Up,              "multilepton_Bjet2_JER_Up/F");
    tOutput->Branch("multilepton_Bjet2_JER_Down",           &multilepton_Bjet2_JER_Down,            "multilepton_Bjet2_JER_Down/F");
    tOutput->Branch("multilepton_Bjet2_DeltaR_Matched",     &multilepton_Bjet2_DeltaR_Matched,      "multilepton_Bjet2_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Bjet2_Label_Matched",      &multilepton_Bjet2_Label_Matched,       "multilepton_Bjet2_Label_Matched/I");
    tOutput->Branch("multilepton_Bjet2_Id_Matched",         &multilepton_Bjet2_Id_Matched,          "multilepton_Bjet2_Id_Matched/I");
    tOutput->Branch("multilepton_Bjet2_P4_Matched",         "TLorentzVector",                       &multilepton_Bjet2_P4_Matched);

    tOutput->Branch("multilepton_JetHighestPt1_Id",&multilepton_JetHighestPt1_Id,"multilepton_JetHighestPt1_Id/I");
    tOutput->Branch("multilepton_JetHighestPt1_P4","TLorentzVector",&multilepton_JetHighestPt1_P4);
    tOutput->Branch("multilepton_JetHighestPt1_CSV",&multilepton_JetHighestPt1_CSV,"multilepton_JetHighestPt1_CSV/F");
    tOutput->Branch("multilepton_JetHighestPt1_JEC_Up",&multilepton_JetHighestPt1_JEC_Up,"multilepton_JetHighestPt1_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestPt1_JEC_Down",&multilepton_JetHighestPt1_JEC_Down,"multilepton_JetHighestPt1_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestPt1_JER_Up",&multilepton_JetHighestPt1_JER_Up,"multilepton_JetHighestPt1_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestPt1_JER_Down",&multilepton_JetHighestPt1_JER_Down,"multilepton_JetHighestPt1_JER_Down/F");

    tOutput->Branch("multilepton_JetHighestPt2_Id",&multilepton_JetHighestPt2_Id,"multilepton_JetHighestPt2_Id/I");
    tOutput->Branch("multilepton_JetHighestPt2_P4","TLorentzVector",&multilepton_JetHighestPt2_P4);
    tOutput->Branch("multilepton_JetHighestPt2_CSV",&multilepton_JetHighestPt2_CSV,"multilepton_JetHighestPt2_CSV/F");
    tOutput->Branch("multilepton_JetHighestPt2_JEC_Up",&multilepton_JetHighestPt2_JEC_Up,"multilepton_JetHighestPt2_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestPt2_JEC_Down",&multilepton_JetHighestPt2_JEC_Down,"multilepton_JetHighestPt2_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestPt2_JER_Up",&multilepton_JetHighestPt2_JER_Up,"multilepton_JetHighestPt2_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestPt2_JER_Down",&multilepton_JetHighestPt2_JER_Down,"multilepton_JetHighestPt2_JER_Down/F");

    tOutput->Branch("multilepton_JetClosestMw1_Id",&multilepton_JetClosestMw1_Id,"multilepton_JetClosestMw1_Id/I");
    tOutput->Branch("multilepton_JetClosestMw1_P4","TLorentzVector",&multilepton_JetClosestMw1_P4);
    tOutput->Branch("multilepton_JetClosestMw1_CSV",&multilepton_JetClosestMw1_CSV,"multilepton_JetClosestMw1_CSV/F");
    tOutput->Branch("multilepton_JetClosestMw1_JEC_Up",&multilepton_JetClosestMw1_JEC_Up,"multilepton_JetClosestMw1_JEC_Up/F");
    tOutput->Branch("multilepton_JetClosestMw1_JEC_Down",&multilepton_JetClosestMw1_JEC_Down,"multilepton_JetClosestMw1_JEC_Down/F");
    tOutput->Branch("multilepton_JetClosestMw1_JER_Up",&multilepton_JetClosestMw1_JER_Up,"multilepton_JetClosestMw1_JER_Up/F");
    tOutput->Branch("multilepton_JetClosestMw1_JER_Down",&multilepton_JetClosestMw1_JER_Down,"multilepton_JetClosestMw1_JER_Down/F");

    tOutput->Branch("multilepton_JetClosestMw2_Id",&multilepton_JetClosestMw2_Id,"multilepton_JetClosestMw2_Id/I");
    tOutput->Branch("multilepton_JetClosestMw2_P4","TLorentzVector",&multilepton_JetClosestMw2_P4);
    tOutput->Branch("multilepton_JetClosestMw2_CSV",&multilepton_JetClosestMw2_CSV,"multilepton_JetClosestMw2_CSV/F");
    tOutput->Branch("multilepton_JetClosestMw2_JEC_Up",&multilepton_JetClosestMw2_JEC_Up,"multilepton_JetClosestMw2_JEC_Up/F");
    tOutput->Branch("multilepton_JetClosestMw2_JEC_Down",&multilepton_JetClosestMw2_JEC_Down,"multilepton_JetClosestMw2_JEC_Down/F");
    tOutput->Branch("multilepton_JetClosestMw2_JER_Up",&multilepton_JetClosestMw2_JER_Up,"multilepton_JetClosestMw2_JER_Up/F");
    tOutput->Branch("multilepton_JetClosestMw2_JER_Down",&multilepton_JetClosestMw2_JER_Down,"multilepton_JetClosestMw2_JER_Down/F");

    tOutput->Branch("multilepton_JetLowestMjj1_Id",&multilepton_JetLowestMjj1_Id,"multilepton_JetLowestMjj1_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj1_P4","TLorentzVector",&multilepton_JetLowestMjj1_P4);
    tOutput->Branch("multilepton_JetLowestMjj1_CSV",&multilepton_JetLowestMjj1_CSV,"multilepton_JetLowestMjj1_CSV/F");
    tOutput->Branch("multilepton_JetLowestMjj1_JEC_Up",&multilepton_JetLowestMjj1_JEC_Up,"multilepton_JetLowestMjj1_JEC_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj1_JEC_Down",&multilepton_JetLowestMjj1_JEC_Down,"multilepton_JetLowestMjj1_JEC_Down/F");
    tOutput->Branch("multilepton_JetLowestMjj1_JER_Up",&multilepton_JetLowestMjj1_JER_Up,"multilepton_JetLowestMjj1_JER_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj1_JER_Down",&multilepton_JetLowestMjj1_JER_Down,"multilepton_JetLowestMjj1_JER_Down/F");

    tOutput->Branch("multilepton_JetLowestMjj2_Id",&multilepton_JetLowestMjj2_Id,"multilepton_JetLowestMjj2_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj2_P4","TLorentzVector",&multilepton_JetLowestMjj2_P4);
    tOutput->Branch("multilepton_JetLowestMjj2_CSV",&multilepton_JetLowestMjj2_CSV,"multilepton_JetLowestMjj2_CSV/F");
    tOutput->Branch("multilepton_JetLowestMjj2_JEC_Up",&multilepton_JetLowestMjj2_JEC_Up,"multilepton_JetLowestMjj2_JEC_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj2_JEC_Down",&multilepton_JetLowestMjj2_JEC_Down,"multilepton_JetLowestMjj2_JEC_Down/F");
    tOutput->Branch("multilepton_JetLowestMjj2_JER_Up",&multilepton_JetLowestMjj2_JER_Up,"multilepton_JetLowestMjj2_JER_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj2_JER_Down",&multilepton_JetLowestMjj2_JER_Down,"multilepton_JetLowestMjj2_JER_Down/F");

    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_Id",&multilepton_JetHighestPt1_2ndPair_Id,"multilepton_JetHighestPt1_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_P4","TLorentzVector",&multilepton_JetHighestPt1_2ndPair_P4);
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_CSV",&multilepton_JetHighestPt1_2ndPair_CSV,"multilepton_JetHighestPt1_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_JEC_Up",&multilepton_JetHighestPt1_2ndPair_JEC_Up,"multilepton_JetHighestPt1_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_JEC_Down",&multilepton_JetHighestPt1_2ndPair_JEC_Down,"multilepton_JetHighestPt1_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_JER_Up",&multilepton_JetHighestPt1_2ndPair_JER_Up,"multilepton_JetHighestPt1_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_JER_Down",&multilepton_JetHighestPt1_2ndPair_JER_Down,"multilepton_JetHighestPt1_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_Id",&multilepton_JetHighestPt2_2ndPair_Id,"multilepton_JetHighestPt2_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_P4","TLorentzVector",&multilepton_JetHighestPt2_2ndPair_P4);
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_CSV",&multilepton_JetHighestPt2_2ndPair_CSV,"multilepton_JetHighestPt2_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_JEC_Up",&multilepton_JetHighestPt2_2ndPair_JEC_Up,"multilepton_JetHighestPt2_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_JEC_Down",&multilepton_JetHighestPt2_2ndPair_JEC_Down,"multilepton_JetHighestPt2_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_JER_Up",&multilepton_JetHighestPt2_2ndPair_JER_Up,"multilepton_JetHighestPt2_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_JER_Down",&multilepton_JetHighestPt2_2ndPair_JER_Down,"multilepton_JetHighestPt2_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_Id",&multilepton_JetClosestMw1_2ndPair_Id,"multilepton_JetClosestMw1_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_P4","TLorentzVector",&multilepton_JetClosestMw1_2ndPair_P4);
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_CSV",&multilepton_JetClosestMw1_2ndPair_CSV,"multilepton_JetClosestMw1_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_JEC_Up",&multilepton_JetClosestMw1_2ndPair_JEC_Up,"multilepton_JetClosestMw1_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_JEC_Down",&multilepton_JetClosestMw1_2ndPair_JEC_Down,"multilepton_JetClosestMw1_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_JER_Up",&multilepton_JetClosestMw1_2ndPair_JER_Up,"multilepton_JetClosestMw1_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_JER_Down",&multilepton_JetClosestMw1_2ndPair_JER_Down,"multilepton_JetClosestMw1_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_Id",&multilepton_JetClosestMw2_2ndPair_Id,"multilepton_JetClosestMw2_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_P4","TLorentzVector",&multilepton_JetClosestMw2_2ndPair_P4);
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_CSV",&multilepton_JetClosestMw2_2ndPair_CSV,"multilepton_JetClosestMw2_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_JEC_Up",&multilepton_JetClosestMw2_2ndPair_JEC_Up,"multilepton_JetClosestMw2_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_JEC_Down",&multilepton_JetClosestMw2_2ndPair_JEC_Down,"multilepton_JetClosestMw2_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_JER_Up",&multilepton_JetClosestMw2_2ndPair_JER_Up,"multilepton_JetClosestMw2_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_JER_Down",&multilepton_JetClosestMw2_2ndPair_JER_Down,"multilepton_JetClosestMw2_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_Id",&multilepton_JetLowestMjj1_2ndPair_Id,"multilepton_JetLowestMjj1_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_P4","TLorentzVector",&multilepton_JetLowestMjj1_2ndPair_P4);
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_CSV",&multilepton_JetLowestMjj1_2ndPair_CSV,"multilepton_JetLowestMjj1_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_JEC_Up",&multilepton_JetLowestMjj1_2ndPair_JEC_Up,"multilepton_JetLowestMjj1_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_JEC_Down",&multilepton_JetLowestMjj1_2ndPair_JEC_Down,"multilepton_JetLowestMjj1_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_JER_Up",&multilepton_JetLowestMjj1_2ndPair_JER_Up,"multilepton_JetLowestMjj1_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_JER_Down",&multilepton_JetLowestMjj1_2ndPair_JER_Down,"multilepton_JetLowestMjj1_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_Id",&multilepton_JetLowestMjj2_2ndPair_Id,"multilepton_JetLowestMjj2_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_P4","TLorentzVector",&multilepton_JetLowestMjj2_2ndPair_P4);
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_CSV",&multilepton_JetLowestMjj2_2ndPair_CSV,"multilepton_JetLowestMjj2_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_JEC_Up",&multilepton_JetLowestMjj2_2ndPair_JEC_Up,"multilepton_JetLowestMjj2_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_JEC_Down",&multilepton_JetLowestMjj2_2ndPair_JEC_Down,"multilepton_JetLowestMjj2_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_JER_Up",&multilepton_JetLowestMjj2_2ndPair_JER_Up,"multilepton_JetLowestMjj2_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_JER_Down",&multilepton_JetLowestMjj2_2ndPair_JER_Down,"multilepton_JetLowestMjj2_2ndPair_JER_Down/F");

    // Test adding truth information
    
    tOutput->Branch("multilepton_h0_Id",                    &multilepton_h0_Id,                 "multilepton_h0_Id/I");
    tOutput->Branch("multilepton_h0_P4",                    "TLorentzVector",                   &multilepton_h0_P4);
    tOutput->Branch("multilepton_t1_Id",                    &multilepton_t1_Id,                 "multilepton_t1_Id/I");
    tOutput->Branch("multilepton_t1_P4",                    "TLorentzVector",                   &multilepton_t1_P4);
    tOutput->Branch("multilepton_t2_Id",                    &multilepton_t2_Id,                 "multilepton_t2_Id/I");
    tOutput->Branch("multilepton_t2_P4",                    "TLorentzVector",                   &multilepton_t2_P4);

    // End test adding truth information

    tOutput->Branch("multilepton_mET","TLorentzVector",&multilepton_mET);
    tOutput->Branch("multilepton_mETcov00",&multilepton_mETcov00,"multilepton_mETcov00/D");
    tOutput->Branch("multilepton_mETcov01",&multilepton_mETcov01,"multilepton_mETcov01/D");
    tOutput->Branch("multilepton_mETcov10",&multilepton_mETcov10,"multilepton_mETcov10/D");
    tOutput->Branch("multilepton_mETcov11",&multilepton_mETcov11,"multilepton_mETcov11/D");
    tOutput->Branch("multilepton_mHT",&multilepton_mHT,"multilepton_mHT/F");
    tOutput->Branch("multilepton_Ptot","TLorentzVector",&multilepton_Ptot);

    return;
}

void TTbarHiggsMultileptonAnalysis::fillOutputTree(){

    bool is2lss=false, is3l=false, is4l=false;

    int tot_charge = 0;
    int tot_id = 0;
    if (vSelectedLeptons.size()>=4) {
        for (unsigned int i=0; i<4; i++) {
            tot_charge += vSelectedLeptons.at(i).charge();
            tot_id += vSelectedLeptons.at(i).id();
        }
    }
    if ((is_3l_TTH_SR || is_3l_TTZ_CR || is_3l_WZrel_CR) && vSelectedLeptons.size()>=4 && tot_charge==0 && (tot_id==0 || abs(tot_id)==2)) is4l = true;
    else if ( ((is_3l_TTH_SR || is_3l_TTZ_CR || is_3l_WZrel_CR) && vSelectedLeptons.size()==3)
            || ((is_3l_TTH_SR || is_3l_TTZ_CR || is_3l_WZrel_CR) && vSelectedLeptons.size()>=4) ) 
        is3l = true;
    else if ( is_2lss_TTH_SR && vSelectedLeptons.size()==2 && vSelectedLeptons.at(0).charge()==vSelectedLeptons.at(1).charge()) is2lss = true;
    if (!is2lss && !is3l && !is4l) return;

    if (vSelectedJets.size()<2) return;
    if (!(vSelectedBTagJets.size()>=2 || (vSelectedMediumBTagJets.size()==1))) return; 

    //if (vSelectedLeptons.size()<2) return; // 2lss only at the moment

    //std::cout << "lept="<<vSelectedLeptons.size()<<" fake="<<vFakeLeptons.size()<<std::endl;
    //std::cout << "btag="<<vSelectedBTagJets.size()<<" nonbtag="<<vSelectedNonBTagJets.size()<<std::endl;

    //Choosing 2 b-jets
    bool doSelectOnlyBjets = false;

    TLorentzVector Bjet1, Bjet2;
    int ib1=-1, ib2=-1;
    selectBjets("HighestBtagDiscrim", &ib1, &ib2, doSelectOnlyBjets);
    if (ib1!=-1) Bjet1.SetPtEtaPhiE(vSelectedJets.at(ib1).pt(), vSelectedJets.at(ib1).eta(), vSelectedJets.at(ib1).phi(), vSelectedJets.at(ib1).E());
    if (ib2!=-1) Bjet2.SetPtEtaPhiE(vSelectedJets.at(ib2).pt(), vSelectedJets.at(ib2).eta(), vSelectedJets.at(ib2).phi(), vSelectedJets.at(ib2).E());

    //2lss
    if (is2lss && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2>=4) catJets = kCat_2lss_2b_4j;
    else if (is2lss && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1>=4) catJets = kCat_2lss_1b_4j;
    else if (is2lss && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==3) catJets = kCat_2lss_2b_3j;
    else if (is2lss && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1==3) catJets = kCat_2lss_1b_3j;
    else if (is2lss && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==2) catJets = kCat_2lss_2b_2j;
    //4l 
    else if (is4l && ib1!=-1 && ib2!=-1) catJets = kCat_4l_2b;
    else if (is4l && ib1!=-1 && ib2==-1) catJets = kCat_4l_1b;
    //3l
    else if (is3l && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2>=2) catJets = kCat_3l_2b_2j;
    else if (is3l && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1>=2) catJets = kCat_3l_1b_2j;
    else if (is3l && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==1) catJets = kCat_3l_2b_1j;
    else if (is3l && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1==1) catJets = kCat_3l_1b_1j;
    else if (is3l && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==0) catJets = kCat_3l_2b_0j;
    else catJets = -1;

    //std::cout << "catJets="<<catJets<<std::endl;

    multilepton_Lepton1_Id = -999;
    multilepton_Lepton2_Id = -999;
    multilepton_Lepton3_Id = -999;
    multilepton_Lepton4_Id = -999;

    if (vSelectedLeptons.size()>=2){
        multilepton_Lepton1_P4 = vSelectedLeptons.at(0).p4();
        multilepton_Lepton1_Id = vSelectedLeptons.at(0).id();
        multilepton_Lepton2_P4 = vSelectedLeptons.at(1).p4();
        multilepton_Lepton2_Id = vSelectedLeptons.at(1).id();
    }

    if (vSelectedLeptons.size()>=3)
    {
        multilepton_Lepton3_P4 = vSelectedLeptons.at(2).p4();
        multilepton_Lepton3_Id = vSelectedLeptons.at(2).id();
    }

    if (vSelectedLeptons.size()>=4 && is4l)
    {
        multilepton_Lepton4_P4 = vSelectedLeptons.at(3).p4();
        multilepton_Lepton4_Id = vSelectedLeptons.at(3).id();
    }

    multilepton_Bjet1_Id = -999;
    if (ib1!=-1){
        FillJetInfoOutputTree(&multilepton_Bjet1_Id, 5, &multilepton_Bjet1_P4, Bjet1, &multilepton_Bjet1_CSV, vSelectedJets.at(ib1).CSVv2(), &multilepton_Bjet1_JEC_Up, &multilepton_Bjet1_JEC_Down, vSelectedJets.at(ib1).JES_uncert(), &multilepton_Bjet1_JER_Up, &multilepton_Bjet1_JER_Down, vSelectedJets.at(ib1).pt_JER(), vSelectedJets.at(ib1).pt_JER_up(), vSelectedJets.at(ib1).pt_JER_down());
        //multilepton_Bjet1_P4 = Bjet1;
        //multilepton_Bjet1_Id = 5;
    }
    multilepton_Bjet2_Id = -999;
    if (ib2!=-1){
        FillJetInfoOutputTree(&multilepton_Bjet2_Id, 5, &multilepton_Bjet2_P4, Bjet2, &multilepton_Bjet2_CSV, vSelectedJets.at(ib2).CSVv2(), &multilepton_Bjet2_JEC_Up, &multilepton_Bjet2_JEC_Down, vSelectedJets.at(ib2).JES_uncert(), &multilepton_Bjet2_JER_Up, &multilepton_Bjet2_JER_Down, vSelectedJets.at(ib2).pt_JER(), vSelectedJets.at(ib2).pt_JER_up(), vSelectedJets.at(ib2).pt_JER_down());
        //multilepton_Bjet2_P4 = Bjet2;
        //multilepton_Bjet2_Id = 5;
    }

    // ###############################################################################
    // #                  _       _     _               _                            #
    // #  _ __ ___   __ _| |_ ___| |__ (_)_ __   __ _  | |_ ___     __ _  ___ _ __   #
    // # | '_ ` _ \ / _` | __/ __| '_ \| | '_ \ / _` | | __/ _ \   / _` |/ _ \ '_ \  #
    // # | | | | | | (_| | || (__| | | | | | | | (_| | | || (_) | | (_| |  __/ | | | #
    // # |_| |_| |_|\__,_|\__\___|_| |_|_|_| |_|\__, |  \__\___/   \__, |\___|_| |_| #
    // #                                        |___/              |___/             #
    // #                                                                             #
    // ###############################################################################

    float lep1_dr_gen       = 100.,     lep2_dr_gen     = 100.,     lep3_dr_gen     = 100.,     lep4_dr_gen     = 100. ;
    float jet1_dr_gen       = 100.,     jet2_dr_gen     = 100.;
    float lep1_dr_gen_min   = 100.,     lep2_dr_gen_min = 100.,     lep3_dr_gen_min = 100.,     lep4_dr_gen_min = 100. ;
    float jet1_dr_gen_min   = 100.,     jet2_dr_gen_min = 100.;
    int   lep1_matched      = -1,       lep2_matched    = -1,       lep3_matched   = -1,       lep4_matched    = -1;
    int   jet1_matched      = -1,       jet2_matched    = -1;  

    TLorentzVector LeptonX;

    for(unsigned int itruth = 0; itruth < vTruth->at(0).mc_truth_label().size() ; itruth++)
    {

        if( abs(vTruth->at(0).mc_truth_id().at(itruth)) < 18 )
        {
            lep1_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta().at(itruth),  vTruth->at(0).mc_truth_phi().at(itruth), vSelectedLeptons.at(0).eta(), vSelectedLeptons.at(0).phi() );
            if( lep1_dr_gen < lep1_dr_gen_min)
            {   lep1_dr_gen_min = lep1_dr_gen;  lep1_matched = itruth;  }

            lep2_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta().at(itruth),  vTruth->at(0).mc_truth_phi().at(itruth), vSelectedLeptons.at(1).eta(), vSelectedLeptons.at(1).phi() );
            if( lep2_dr_gen < lep2_dr_gen_min)
            {   lep2_dr_gen_min = lep2_dr_gen;  lep2_matched = itruth;  }

            if(vSelectedLeptons.size()>=3)
            {    
                lep3_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta().at(itruth),  vTruth->at(0).mc_truth_phi().at(itruth), vSelectedLeptons.at(2).eta(), vSelectedLeptons.at(2).phi() );
                if( lep3_dr_gen < lep3_dr_gen_min)
                {   lep3_dr_gen_min = lep3_dr_gen;  lep3_matched = itruth;  }
            }

            if(vSelectedLeptons.size()>=4)
            {
                lep4_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta().at(itruth),  vTruth->at(0).mc_truth_phi().at(itruth), vSelectedLeptons.at(3).eta(), vSelectedLeptons.at(3).phi() );
                if( lep4_dr_gen < lep4_dr_gen_min)
                {   lep4_dr_gen_min = lep4_dr_gen;  lep4_matched = itruth;  }
            }

            jet1_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta().at(itruth),  vTruth->at(0).mc_truth_phi().at(itruth), vSelectedJets.at(ib1).eta(), vSelectedJets.at(ib1).phi() );
            if( jet1_dr_gen < jet1_dr_gen_min)
            {   jet1_dr_gen_min = jet1_dr_gen;  jet1_matched = itruth;  }        

            jet2_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta().at(itruth),  vTruth->at(0).mc_truth_phi().at(itruth), vSelectedJets.at(ib2).eta(), vSelectedJets.at(ib2).phi() );
            if( jet2_dr_gen < jet2_dr_gen_min)
            {   jet2_dr_gen_min = jet2_dr_gen;  jet2_matched = itruth;  }
        }
    }

    if(false)
    {
        std::cout << "lep1_matched: "   << lep1_matched                                     << std::endl;
        std::cout << "pt: "             << vTruth->at(0).mc_truth_pt().at(lep1_matched)     << std::endl;
        std::cout << "eta: "            << vTruth->at(0).mc_truth_eta().at(lep1_matched)    << std::endl;
        std::cout << "phi: "            << vTruth->at(0).mc_truth_phi().at(lep1_matched)    << std::endl;
        std::cout << "E: "              << vTruth->at(0).mc_truth_E().at(lep1_matched)      << std::endl;
        std::cout << "Id: "             << vTruth->at(0).mc_truth_id().at(lep1_matched)     << std::endl;
        std::cout << "Label: "          << vTruth->at(0).mc_truth_label().at(lep1_matched)  << std::endl;
    }

    if(lep1_matched >= 0)
    {
        LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(lep1_matched),       vTruth->at(0).mc_truth_eta().at(lep1_matched), 
                                vTruth->at(0).mc_truth_phi().at(lep1_matched),      vTruth->at(0).mc_truth_E().at(lep1_matched)     );
        multilepton_Lepton1_P4_Matched      = LeptonX;
        multilepton_Lepton1_Id_Matched      = vTruth->at(0).mc_truth_id().at(lep1_matched);
        multilepton_Lepton1_Label_Matched   = vTruth->at(0).mc_truth_label().at(lep1_matched);
        multilepton_Lepton1_DeltaR_Matched  = lep1_dr_gen_min;
    }

    if(lep2_matched >= 0)
    {    
        LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(lep2_matched),       vTruth->at(0).mc_truth_eta().at(lep2_matched),
                                vTruth->at(0).mc_truth_phi().at(lep2_matched),      vTruth->at(0).mc_truth_E().at(lep2_matched)     );
        multilepton_Lepton2_P4_Matched      = LeptonX;
        multilepton_Lepton2_Id_Matched      = vTruth->at(0).mc_truth_id().at(lep2_matched);
        multilepton_Lepton2_Label_Matched   = vTruth->at(0).mc_truth_label().at(lep2_matched);
        multilepton_Lepton2_DeltaR_Matched  = lep2_dr_gen_min;
    }

    if(lep3_matched >= 0)
    {
        if(vSelectedLeptons.size()>=3)
        {
            LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(lep3_matched),       vTruth->at(0).mc_truth_eta().at(lep3_matched),
                                    vTruth->at(0).mc_truth_phi().at(lep3_matched),      vTruth->at(0).mc_truth_E().at(lep3_matched)     );
            multilepton_Lepton3_P4_Matched      = LeptonX;
            multilepton_Lepton3_Id_Matched      = vTruth->at(0).mc_truth_id().at(lep3_matched);
            multilepton_Lepton3_Label_Matched   = vTruth->at(0).mc_truth_label().at(lep3_matched);
            multilepton_Lepton3_DeltaR_Matched  = lep3_dr_gen_min;
        }
    }

    if(lep4_matched >= 0)
    {
        if(vSelectedLeptons.size()>=4)
        {
            LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(lep4_matched),       vTruth->at(0).mc_truth_eta().at(lep4_matched),
                                    vTruth->at(0).mc_truth_phi().at(lep4_matched),      vTruth->at(0).mc_truth_E().at(lep4_matched)     );
            multilepton_Lepton4_P4_Matched      = LeptonX;
            multilepton_Lepton4_Id_Matched      = vTruth->at(0).mc_truth_id().at(lep4_matched);
            multilepton_Lepton4_Label_Matched   = vTruth->at(0).mc_truth_label().at(lep4_matched);
            multilepton_Lepton4_DeltaR_Matched  = lep4_dr_gen_min;
        }
    }

    if(false)
    {
        std::cout << " ============ "   << std::endl;
        std::cout << "jet1_matched: "   << jet1_matched                                     << std::endl;
        std::cout << "pt: "             << vTruth->at(0).mc_truth_pt().at(jet1_matched)     << std::endl;
        std::cout << "eta: "            << vTruth->at(0).mc_truth_eta().at(jet1_matched)    << std::endl;
        std::cout << "phi: "            << vTruth->at(0).mc_truth_phi().at(jet1_matched)    << std::endl;
        std::cout << "E: "              << vTruth->at(0).mc_truth_E().at(jet1_matched)      << std::endl;
        std::cout << "Id: "             << vTruth->at(0).mc_truth_id().at(jet1_matched)     << std::endl;
        std::cout << "Label: "          << vTruth->at(0).mc_truth_label().at(jet1_matched)  << std::endl;    
    }

    if(jet1_matched >= 0)
    {
        LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(jet1_matched),       vTruth->at(0).mc_truth_eta().at(jet1_matched),
                                vTruth->at(0).mc_truth_phi().at(jet1_matched),      vTruth->at(0).mc_truth_E().at(jet1_matched)     );
        multilepton_Bjet1_P4_Matched        = LeptonX;
        multilepton_Bjet1_Id_Matched        = vTruth->at(0).mc_truth_id().at(jet1_matched);
        multilepton_Bjet1_Label_Matched     = vTruth->at(0).mc_truth_label().at(jet1_matched);
        multilepton_Bjet1_DeltaR_Matched    = jet1_dr_gen_min;
    }

    if(jet2_matched >= 0)
    {
        LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(jet2_matched),       vTruth->at(0).mc_truth_eta().at(jet2_matched),
                                vTruth->at(0).mc_truth_phi().at(jet2_matched),      vTruth->at(0).mc_truth_E().at(jet2_matched)     );
        multilepton_Bjet2_P4_Matched        = LeptonX;
        multilepton_Bjet2_Id_Matched        = vTruth->at(0).mc_truth_id().at(jet2_matched);
        multilepton_Bjet2_Label_Matched     = vTruth->at(0).mc_truth_label().at(jet2_matched);
        multilepton_Bjet2_DeltaR_Matched    = jet2_dr_gen_min;
    }

    // ========================

    //Choose 2 jets
    TLorentzVector Pjet1, Pjet2;
    float pt_max=0, pt_max2=0; int ij1=-1, ij2=-1;
    float diffmass_min = 10000, mass_min = 10000; int ik1=-1, ik2=-1, il1=-1, il2=-1;
    for (unsigned int ij=0; ij<vSelectedJets.size(); ij++){
        if (ij==ib1 || ij==ib2) continue;
        if (vSelectedJets.at(ij).pt() > pt_max ) {
            pt_max2 = pt_max;
            ij2 = ij1;
            pt_max = vSelectedJets.at(ij).pt();
            ij1 = ij;
        } 
        if (vSelectedJets.at(ij).pt() < pt_max && vSelectedJets.at(ij).pt() > pt_max2){
            pt_max2 = vSelectedJets.at(ij).pt(); 
            ij2 = ij; 
        } 
        for (unsigned int ik=0; ik<vSelectedJets.size(); ik++){
            if (ik==ij) continue;
            if (ik==ib1 || ik==ib2) continue;
            Pjet1.SetPtEtaPhiE(vSelectedJets.at(ij).pt(), vSelectedJets.at(ij).eta(), vSelectedJets.at(ij).phi(), vSelectedJets.at(ij).E());
            Pjet2.SetPtEtaPhiE(vSelectedJets.at(ik).pt(), vSelectedJets.at(ik).eta(), vSelectedJets.at(ik).phi(), vSelectedJets.at(ik).E()); 
            if (TMath::Abs((Pjet1+Pjet2).M()-80.419)<diffmass_min){
                ik1=ij;
                ik2=ik;
                diffmass_min = TMath::Abs((Pjet1+Pjet2).M()-80.419);
            } 
            if ((Pjet1+Pjet2).M()<mass_min){
                il1=ij;
                il2=ik;
                mass_min = (Pjet1+Pjet2).M();
            } 
        } 
    }  

    //Choose 2 more jets
    int io1=-1, io2=-1, ip1=-1, ip2=-1, im1=-1, im2=-1;
    diffmass_min = 10000, mass_min = 10000, pt_max2 = 0, pt_max = 0;
    for (unsigned int im=0; im<vSelectedJets.size(); im++){
        if (im==ib1 || im==ib2 || im==ik1 || im==ik2) continue;
        if (vSelectedJets.at(im).pt() > pt_max ) {
            pt_max2 = pt_max;
            im2 = im1;
            pt_max = vSelectedJets.at(im).pt();
            im1 = im;
        }
        if (vSelectedJets.at(im).pt() < pt_max && vSelectedJets.at(im).pt() > pt_max2){
            pt_max2 = vSelectedJets.at(im).pt();
            im2 = im;
        }
        for (unsigned int in=0; in<vSelectedJets.size(); in++){
            if (in==ib1 || in==ib2 || in==ik1 || in==ik2 || in==im) continue;
            Pjet1.SetPtEtaPhiE(vSelectedJets.at(im).pt(), vSelectedJets.at(im).eta(), vSelectedJets.at(im).phi(), vSelectedJets.at(im).E());
            Pjet2.SetPtEtaPhiE(vSelectedJets.at(in).pt(), vSelectedJets.at(in).eta(), vSelectedJets.at(in).phi(), vSelectedJets.at(in).E());
            if (TMath::Abs((Pjet1+Pjet2).M()-80.419)<diffmass_min){
                io1=im;
                io2=in;
                diffmass_min = TMath::Abs((Pjet1+Pjet2).M()-80.419);
            }
            if ((Pjet1+Pjet2).M()<mass_min){
                ip1=im;
                ip2=in;
                mass_min = (Pjet1+Pjet2).M();
            }
        }
    }

    multilepton_JetHighestPt1_Id = -999;
    multilepton_JetHighestPt2_Id = -999;
    multilepton_JetClosestMw1_Id = -999;
    multilepton_JetClosestMw2_Id = -999;
    multilepton_JetLowestMjj1_Id = -999;
    multilepton_JetLowestMjj2_Id = -999;

    multilepton_JetHighestPt1_2ndPair_Id = -999;
    multilepton_JetHighestPt2_2ndPair_Id = -999;
    multilepton_JetClosestMw1_2ndPair_Id = -999;
    multilepton_JetClosestMw2_2ndPair_Id = -999;
    multilepton_JetLowestMjj1_2ndPair_Id = -999;
    multilepton_JetLowestMjj2_2ndPair_Id = -999;

    TLorentzVector Jet1, Jet2;

    if (ij1!=-1 && ij2==-1){
        Jet1.SetPtEtaPhiE(vSelectedJets.at(ij1).pt(), vSelectedJets.at(ij1).eta(), vSelectedJets.at(ij1).phi(), vSelectedJets.at(ij1).E());
        FillJetInfoOutputTree(&multilepton_JetHighestPt1_Id, 1, &multilepton_JetHighestPt1_P4, Jet1, &multilepton_JetHighestPt1_CSV, vSelectedJets.at(ij1).CSVv2(), &multilepton_JetHighestPt1_JEC_Up, &multilepton_JetHighestPt1_JEC_Down, vSelectedJets.at(ij1).JES_uncert(), &multilepton_JetHighestPt1_JER_Up, &multilepton_JetHighestPt1_JER_Down, vSelectedJets.at(ij1).pt_JER(), vSelectedJets.at(ij1).pt_JER_up(), vSelectedJets.at(ij1).pt_JER_down());
        //multilepton_JetHighestPt1_Id = 1;
        //multilepton_JetHighestPt1_P4.SetPtEtaPhiE(vSelectedJets.at(ij1).pt(), vSelectedJets.at(ij1).eta(), vSelectedJets.at(ij1).phi(), vSelectedJets.at(ij1).E());
    }   
    if (ij1!=-1 && ij2!=-1) {
        Jet1.SetPtEtaPhiE(vSelectedJets.at(ij1).pt(), vSelectedJets.at(ij1).eta(), vSelectedJets.at(ij1).phi(), vSelectedJets.at(ij1).E());
        Jet2.SetPtEtaPhiE(vSelectedJets.at(ij2).pt(), vSelectedJets.at(ij2).eta(), vSelectedJets.at(ij2).phi(), vSelectedJets.at(ij2).E());
        FillJetInfoOutputTree(&multilepton_JetHighestPt1_Id, 1, &multilepton_JetHighestPt1_P4, Jet1, &multilepton_JetHighestPt1_CSV, vSelectedJets.at(ij1).CSVv2(), &multilepton_JetHighestPt1_JEC_Up, &multilepton_JetHighestPt1_JEC_Down, vSelectedJets.at(ij1).JES_uncert(), &multilepton_JetHighestPt1_JER_Up, &multilepton_JetHighestPt1_JER_Down, vSelectedJets.at(ij1).pt_JER(), vSelectedJets.at(ij1).pt_JER_up(), vSelectedJets.at(ij1).pt_JER_down());
        FillJetInfoOutputTree(&multilepton_JetHighestPt2_Id, 1, &multilepton_JetHighestPt2_P4, Jet2, &multilepton_JetHighestPt2_CSV, vSelectedJets.at(ij2).CSVv2(), &multilepton_JetHighestPt2_JEC_Up, &multilepton_JetHighestPt2_JEC_Down, vSelectedJets.at(ij2).JES_uncert(), &multilepton_JetHighestPt2_JER_Up, &multilepton_JetHighestPt2_JER_Down, vSelectedJets.at(ij2).pt_JER(), vSelectedJets.at(ij2).pt_JER_up(), vSelectedJets.at(ij2).pt_JER_down());	
        //multilepton_JetHighestPt1_Id = 1;
        //multilepton_JetHighestPt2_Id = 1;
        //multilepton_JetHighestPt1_P4.SetPtEtaPhiE(vSelectedJets.at(ij1).pt(), vSelectedJets.at(ij1).eta(), vSelectedJets.at(ij1).phi(), vSelectedJets.at(ij1).E());
        //multilepton_JetHighestPt2_P4.SetPtEtaPhiE(vSelectedJets.at(ij2).pt(), vSelectedJets.at(ij2).eta(), vSelectedJets.at(ij2).phi(), vSelectedJets.at(ij2).E());
    }
    if (ik1!=-1 && ik2!=-1){
        Jet1.SetPtEtaPhiE(vSelectedJets.at(ik1).pt(), vSelectedJets.at(ik1).eta(), vSelectedJets.at(ik1).phi(), vSelectedJets.at(ik1).E());
        Jet2.SetPtEtaPhiE(vSelectedJets.at(ik2).pt(), vSelectedJets.at(ik2).eta(), vSelectedJets.at(ik2).phi(), vSelectedJets.at(ik2).E());
        FillJetInfoOutputTree(&multilepton_JetClosestMw1_Id, 2, &multilepton_JetClosestMw1_P4, Jet1, &multilepton_JetClosestMw1_CSV, vSelectedJets.at(ik1).CSVv2(), &multilepton_JetClosestMw1_JEC_Up, &multilepton_JetClosestMw1_JEC_Down, vSelectedJets.at(ik1).JES_uncert(), &multilepton_JetClosestMw1_JER_Up, &multilepton_JetClosestMw1_JER_Down, vSelectedJets.at(ik1).pt_JER(), vSelectedJets.at(ik1).pt_JER_up(), vSelectedJets.at(ik1).pt_JER_down());
        FillJetInfoOutputTree(&multilepton_JetClosestMw2_Id, 2, &multilepton_JetClosestMw2_P4, Jet2, &multilepton_JetClosestMw2_CSV, vSelectedJets.at(ik2).CSVv2(), &multilepton_JetClosestMw2_JEC_Up, &multilepton_JetClosestMw2_JEC_Down, vSelectedJets.at(ik2).JES_uncert(), &multilepton_JetClosestMw2_JER_Up, &multilepton_JetClosestMw2_JER_Down, vSelectedJets.at(ik2).pt_JER(), vSelectedJets.at(ik2).pt_JER_up(), vSelectedJets.at(ik2).pt_JER_down());
        //multilepton_JetClosestMw1_Id = 2;
        //multilepton_JetClosestMw2_Id = 2;
        //multilepton_JetClosestMw1_P4.SetPtEtaPhiE(vSelectedJets.at(ik1).pt(), vSelectedJets.at(ik1).eta(), vSelectedJets.at(ik1).phi(), vSelectedJets.at(ik1).E());
        //multilepton_JetClosestMw2_P4.SetPtEtaPhiE(vSelectedJets.at(ik2).pt(), vSelectedJets.at(ik2).eta(), vSelectedJets.at(ik2).phi(), vSelectedJets.at(ik2).E());
    }
    if (il1!=-1 && il2!=-1){
        Jet1.SetPtEtaPhiE(vSelectedJets.at(il1).pt(), vSelectedJets.at(il1).eta(), vSelectedJets.at(il1).phi(), vSelectedJets.at(il1).E());
        Jet2.SetPtEtaPhiE(vSelectedJets.at(il2).pt(), vSelectedJets.at(il2).eta(), vSelectedJets.at(il2).phi(), vSelectedJets.at(il2).E());
        FillJetInfoOutputTree(&multilepton_JetLowestMjj1_Id, 3, &multilepton_JetLowestMjj1_P4, Jet1, &multilepton_JetLowestMjj1_CSV, vSelectedJets.at(il1).CSVv2(), &multilepton_JetLowestMjj1_JEC_Up, &multilepton_JetLowestMjj1_JEC_Down, vSelectedJets.at(il1).JES_uncert(), &multilepton_JetLowestMjj1_JER_Up, &multilepton_JetLowestMjj1_JER_Down, vSelectedJets.at(il1).pt_JER(), vSelectedJets.at(il1).pt_JER_up(), vSelectedJets.at(il1).pt_JER_down());
        FillJetInfoOutputTree(&multilepton_JetLowestMjj2_Id, 3, &multilepton_JetLowestMjj2_P4, Jet2, &multilepton_JetLowestMjj2_CSV, vSelectedJets.at(il2).CSVv2(), &multilepton_JetLowestMjj2_JEC_Up, &multilepton_JetLowestMjj2_JEC_Down, vSelectedJets.at(il2).JES_uncert(), &multilepton_JetLowestMjj2_JER_Up, &multilepton_JetLowestMjj2_JER_Down, vSelectedJets.at(il2).pt_JER(), vSelectedJets.at(il2).pt_JER_up(), vSelectedJets.at(il2).pt_JER_down());
        //multilepton_JetLowestMjj1_Id = 3;
        //multilepton_JetLowestMjj2_Id = 3;
        //multilepton_JetLowestMjj1_P4.SetPtEtaPhiE(vSelectedJets.at(il1).pt(), vSelectedJets.at(il1).eta(), vSelectedJets.at(il1).phi(), vSelectedJets.at(il1).E());
        //multilepton_JetLowestMjj2_P4.SetPtEtaPhiE(vSelectedJets.at(il2).pt(), vSelectedJets.at(il2).eta(), vSelectedJets.at(il2).phi(), vSelectedJets.at(il2).E());
    }

    //2nd pair (first one: closest to Mw)
    if (is2lss && ij1!=-1 && ij2!=-1){
        if (im1!=-1 && im2==-1){
            Jet1.SetPtEtaPhiE(vSelectedJets.at(im1).pt(), vSelectedJets.at(im1).eta(), vSelectedJets.at(im1).phi(), vSelectedJets.at(im1).E());
            FillJetInfoOutputTree(&multilepton_JetHighestPt1_2ndPair_Id, 1, &multilepton_JetHighestPt1_2ndPair_P4, Jet1, &multilepton_JetHighestPt1_2ndPair_CSV, vSelectedJets.at(im1).CSVv2(), &multilepton_JetHighestPt1_2ndPair_JEC_Up, &multilepton_JetHighestPt1_2ndPair_JEC_Down, vSelectedJets.at(im1).JES_uncert(), &multilepton_JetHighestPt1_2ndPair_JER_Up, &multilepton_JetHighestPt1_2ndPair_JER_Down, vSelectedJets.at(im1).pt_JER(), vSelectedJets.at(im1).pt_JER_up(), vSelectedJets.at(im1).pt_JER_down());
            //multilepton_JetHighestPt1_2ndPair_Id = 1;
            //multilepton_JetHighestPt1_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(im1).pt(), vSelectedJets.at(im1).eta(), vSelectedJets.at(im1).phi(), vSelectedJets.at(im1).E());
        }
        if (im1!=-1 && im2!=-1){
            Jet1.SetPtEtaPhiE(vSelectedJets.at(im1).pt(), vSelectedJets.at(im1).eta(), vSelectedJets.at(im1).phi(), vSelectedJets.at(im1).E());
            Jet2.SetPtEtaPhiE(vSelectedJets.at(im2).pt(), vSelectedJets.at(im2).eta(), vSelectedJets.at(im2).phi(), vSelectedJets.at(im2).E());
            FillJetInfoOutputTree(&multilepton_JetHighestPt1_2ndPair_Id, 1, &multilepton_JetHighestPt1_2ndPair_P4, Jet1, &multilepton_JetHighestPt1_2ndPair_CSV, vSelectedJets.at(im1).CSVv2(), &multilepton_JetHighestPt1_2ndPair_JEC_Up, &multilepton_JetHighestPt1_2ndPair_JEC_Down, vSelectedJets.at(im1).JES_uncert(), &multilepton_JetHighestPt1_2ndPair_JER_Up, &multilepton_JetHighestPt1_2ndPair_JER_Down, vSelectedJets.at(im1).pt_JER(), vSelectedJets.at(im1).pt_JER_up(), vSelectedJets.at(im1).pt_JER_down());
            FillJetInfoOutputTree(&multilepton_JetHighestPt2_2ndPair_Id, 1, &multilepton_JetHighestPt2_2ndPair_P4, Jet2, &multilepton_JetHighestPt2_2ndPair_CSV, vSelectedJets.at(im2).CSVv2(), &multilepton_JetHighestPt2_2ndPair_JEC_Up, &multilepton_JetHighestPt2_2ndPair_JEC_Down, vSelectedJets.at(im2).JES_uncert(), &multilepton_JetHighestPt2_2ndPair_JER_Up, &multilepton_JetHighestPt2_2ndPair_JER_Down, vSelectedJets.at(im2).pt_JER(), vSelectedJets.at(im2).pt_JER_up(), vSelectedJets.at(im2).pt_JER_down());
            //multilepton_JetHighestPt1_2ndPair_Id = 1;
            //multilepton_JetHighestPt1_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(im1).pt(), vSelectedJets.at(im1).eta(), vSelectedJets.at(im1).phi(), vSelectedJets.at(im1).E());
            //multilepton_JetHighestPt2_2ndPair_Id = 1;
            //multilepton_JetHighestPt2_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(im2).pt(), vSelectedJets.at(im2).eta(), vSelectedJets.at(im2).phi(), vSelectedJets.at(im2).E());   
        }
        if (io1!=-1 && io2!=-1){
            Jet1.SetPtEtaPhiE(vSelectedJets.at(ip1).pt(), vSelectedJets.at(ip1).eta(), vSelectedJets.at(ip1).phi(), vSelectedJets.at(ip1).E());
            Jet2.SetPtEtaPhiE(vSelectedJets.at(io2).pt(), vSelectedJets.at(io2).eta(), vSelectedJets.at(io2).phi(), vSelectedJets.at(io2).E());
            FillJetInfoOutputTree(&multilepton_JetClosestMw1_2ndPair_Id, 2, &multilepton_JetClosestMw1_2ndPair_P4, Jet1, &multilepton_JetClosestMw1_2ndPair_CSV, vSelectedJets.at(io1).CSVv2(), &multilepton_JetClosestMw1_2ndPair_JEC_Up, &multilepton_JetClosestMw1_2ndPair_JEC_Down, vSelectedJets.at(io1).JES_uncert(), &multilepton_JetClosestMw1_2ndPair_JER_Up, &multilepton_JetClosestMw1_2ndPair_JER_Down, vSelectedJets.at(io1).pt_JER(), vSelectedJets.at(io1).pt_JER_up(), vSelectedJets.at(io1).pt_JER_down());
            FillJetInfoOutputTree(&multilepton_JetClosestMw2_2ndPair_Id, 2, &multilepton_JetClosestMw2_2ndPair_P4, Jet2, &multilepton_JetClosestMw2_2ndPair_CSV, vSelectedJets.at(io2).CSVv2(), &multilepton_JetClosestMw2_2ndPair_JEC_Up, &multilepton_JetClosestMw2_2ndPair_JEC_Down, vSelectedJets.at(io2).JES_uncert(), &multilepton_JetClosestMw2_2ndPair_JER_Up, &multilepton_JetClosestMw2_2ndPair_JER_Down, vSelectedJets.at(io2).pt_JER(), vSelectedJets.at(io2).pt_JER_up(), vSelectedJets.at(io2).pt_JER_down());
            //multilepton_JetClosestMw1_2ndPair_Id = 2;
            //multilepton_JetClosestMw2_2ndPair_Id = 2;
            //multilepton_JetClosestMw1_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(io1).pt(), vSelectedJets.at(io1).eta(), vSelectedJets.at(io1).phi(), vSelectedJets.at(io1).E());
            //multilepton_JetClosestMw2_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(io2).pt(), vSelectedJets.at(io2).eta(), vSelectedJets.at(io2).phi(), vSelectedJets.at(io2).E());
        }
        if (ip1!=-1 && ip2!=-1){
            Jet1.SetPtEtaPhiE(vSelectedJets.at(ip1).pt(), vSelectedJets.at(ip1).eta(), vSelectedJets.at(ip1).phi(), vSelectedJets.at(ip1).E());
            Jet2.SetPtEtaPhiE(vSelectedJets.at(ip2).pt(), vSelectedJets.at(ip2).eta(), vSelectedJets.at(ip2).phi(), vSelectedJets.at(ip2).E());
            FillJetInfoOutputTree(&multilepton_JetLowestMjj1_2ndPair_Id, 3, &multilepton_JetLowestMjj1_2ndPair_P4, Jet1, &multilepton_JetLowestMjj1_2ndPair_CSV, vSelectedJets.at(ip1).CSVv2(), &multilepton_JetLowestMjj1_2ndPair_JEC_Up, &multilepton_JetLowestMjj1_2ndPair_JEC_Down, vSelectedJets.at(ip1).JES_uncert(), &multilepton_JetLowestMjj1_2ndPair_JER_Up, &multilepton_JetLowestMjj1_2ndPair_JER_Down, vSelectedJets.at(ip1).pt_JER(), vSelectedJets.at(ip1).pt_JER_up(), vSelectedJets.at(ip1).pt_JER_down());
            FillJetInfoOutputTree(&multilepton_JetLowestMjj2_2ndPair_Id, 3, &multilepton_JetLowestMjj2_2ndPair_P4, Jet2, &multilepton_JetLowestMjj2_2ndPair_CSV, vSelectedJets.at(ip2).CSVv2(), &multilepton_JetLowestMjj2_2ndPair_JEC_Up, &multilepton_JetLowestMjj2_2ndPair_JEC_Down, vSelectedJets.at(ip2).JES_uncert(), &multilepton_JetLowestMjj2_2ndPair_JER_Up, &multilepton_JetLowestMjj2_2ndPair_JER_Down, vSelectedJets.at(ip2).pt_JER(), vSelectedJets.at(ip2).pt_JER_up(), vSelectedJets.at(ip2).pt_JER_down());
            //multilepton_JetLowestMjj1_2ndPair_Id = 3;
            //multilepton_JetLowestMjj2_2ndPair_Id = 3;
            //multilepton_JetLowestMjj1_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(ip1).pt(), vSelectedJets.at(ip1).eta(), vSelectedJets.at(ip1).phi(), vSelectedJets.at(ip1).E());
            //multilepton_JetLowestMjj2_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(ip2).pt(), vSelectedJets.at(ip2).eta(), vSelectedJets.at(ip2).phi(), vSelectedJets.at(ip2).E());
        }
    }

    // ##########################################################################
    // #      _                  _       _                                      #
    // #  ___| |_ __ _ _ __   __| | __ _| | ___  _ __   ___    __ _  ___ _ __   #
    // # / __| __/ _` | '_ \ / _` |/ _` | |/ _ \| '_ \ / _ \  / _` |/ _ \ '_ \  #
    // # \__ \ || (_| | | | | (_| | (_| | | (_) | | | |  __/ | (_| |  __/ | | | #
    // # |___/\__\__,_|_| |_|\__,_|\__,_|_|\___/|_| |_|\___|  \__, |\___|_| |_| #
    // #                                                      |___/             #
    // #                                                                        #
    // ##########################################################################

    for(unsigned int itruth = 0; itruth < vTruth->at(0).mc_truth_label().size() ; itruth++)
    {
        TLorentzVector LeptonX;

        if( vTruth->at(0).mc_truth_label().at(itruth) == 1 )
        { 
            LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(itruth), 
                                    vTruth->at(0).mc_truth_eta().at(itruth), 
                                    vTruth->at(0).mc_truth_phi().at(itruth), 
                                    vTruth->at(0).mc_truth_E().at(itruth) );
           
            multilepton_h0_P4 = LeptonX;
            multilepton_h0_Id = vTruth->at(0).mc_truth_id().at(itruth);
        }

        if( vTruth->at(0).mc_truth_label().at(itruth) == 2 )
        {
            LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(itruth),
                                    vTruth->at(0).mc_truth_eta().at(itruth), 
                                    vTruth->at(0).mc_truth_phi().at(itruth), 
                                    vTruth->at(0).mc_truth_E().at(itruth) );

            multilepton_t1_P4 = LeptonX;
            multilepton_t1_Id = vTruth->at(0).mc_truth_id().at(itruth);
        }        

        if( vTruth->at(0).mc_truth_label().at(itruth) == 3 )
        {
            LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(itruth),
                                    vTruth->at(0).mc_truth_eta().at(itruth),
                                    vTruth->at(0).mc_truth_phi().at(itruth),
                                    vTruth->at(0).mc_truth_E().at(itruth) );

            multilepton_t2_P4 = LeptonX;
            multilepton_t2_Id = vTruth->at(0).mc_truth_id().at(itruth);
        }
    }

    multilepton_mET.SetPtEtaPhiE(vEvent->at(0).metpt(), 0, vEvent->at(0).metphi(), vEvent->at(0).metpt());
    multilepton_mETcov00 = vEvent->at(0).metcov00();
    multilepton_mETcov01 = vEvent->at(0).metcov01();
    multilepton_mETcov10 = vEvent->at(0).metcov10();
    multilepton_mETcov11 = vEvent->at(0).metcov11();
    multilepton_mHT = vEvent->at(0).metsumet();

    mc_ttZhypAllowed = 0;
    /*
        if(vSelectedLeptons.size()==3) 
        {
            if ( vSelectedLeptons.at(0).charge()==vSelectedLeptons.at(1).charge() && vSelectedLeptons.at(1).charge()==vSelectedLeptons.at(2).charge() ) mc_ttZhypAllowed =-1;
            else if (  ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(1).id() ) 
            || ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(2).id() ) 
            || ( vSelectedLeptons.at(1).id() == -vSelectedLeptons.at(2).id() ))
            mc_ttZhypAllowed = 1; 
        }
    */

    if (multilepton_Lepton1_Id!=-999 && multilepton_Lepton2_Id!=-999 && multilepton_Lepton3_Id!=-999){
        if (multilepton_Lepton1_Id*multilepton_Lepton2_Id>0 && multilepton_Lepton2_Id*multilepton_Lepton3_Id>0) mc_ttZhypAllowed =-1;
        else if ( (multilepton_Lepton1_Id==-multilepton_Lepton2_Id)
                || (multilepton_Lepton1_Id==-multilepton_Lepton3_Id)
                || (multilepton_Lepton2_Id==-multilepton_Lepton3_Id))
            mc_ttZhypAllowed = 1; }

    mc_nJets25              = vSelectedJets.size();
    mc_nBtagJets25          = vSelectedBTagJets.size();
    mc_nMediumBtagJets25    = vSelectedMediumBTagJets.size();
    mc_nNonBtagJets25       = vSelectedNonBTagJets.size();
    
    is_2bTight = (vSelectedMediumBTagJets.size()>=2)?1:0;

    tOutput->Fill();

    //if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

void TTbarHiggsMultileptonAnalysis::FillJetInfoOutputTree(int* tree_Id, int Id, TLorentzVector* tree_P4, TLorentzVector P4, float* tree_CSV, float CSV, float* tree_JEC_Up, float* tree_JEC_Down, float JEC_value, float* tree_JER_Up, float* tree_JER_Down, float JER, float JER_Up, float JER_Down){

    *tree_Id = Id;
    *tree_P4 = P4;

    *tree_CSV = CSV;

    *tree_JEC_Up = P4.E()*(1.+JEC_value);
    *tree_JEC_Down = P4.E()*(1.-JEC_value);

    *tree_JER_Up = P4.E()*JER_Up/JER;
    *tree_JER_Down = P4.E()*JER_Down/JER;

    return;
}

void TTbarHiggsMultileptonAnalysis::PrintLHCOforMadweight_MC(int evt)
{  

    //std::cout <<" _MC 1"<< std::endl;

    /*if( proc<-1 || proc > 6 )f
      {
      std::cout << "proc can only take following values: -1,1,2,3,4,5,6" << std::endl;
      std::cout << "3l final state specific" << std::endl;    
      std::cout << "1,2,3,4: specific to the ttH final state, cf patches provided to madweight" << std::endl;
      std::cout << "5: specific to the ttZ with l+l-l+ final state" << std::endl;
      std::cout << "6: specific to the ttZ with l+l-l- final state" << std::endl;
      std::cout << "-1: no selection on the final state applied" << std::endl;

      return;
      }*/

    fline0 = "0      " + std::string(Form("%d",evt)) + "     " + trig;

    int nobj = 1;	      

    int multilepton_Lepton1_Id_LHCO = -666;
    int multilepton_Lepton2_Id_LHCO = -666;
    int multilepton_Lepton3_Id_LHCO = -666;

    //
    // LHCO lepton ID convention
    //  
    if (abs(vTruth->at(0).Leptons_id().at(0))==11) multilepton_Lepton1_Id_LHCO = 1 ;
    else if (abs(vTruth->at(0).Leptons_id().at(0))==13) multilepton_Lepton1_Id_LHCO = 2 ;
    if (abs(vTruth->at(0).Leptons_id().at(1))==11) multilepton_Lepton2_Id_LHCO = 1 ;
    else if (abs(vTruth->at(0).Leptons_id().at(1))==13) multilepton_Lepton2_Id_LHCO = 2 ;
    if (abs(vTruth->at(0).Leptons_id().at(2))==11) multilepton_Lepton3_Id_LHCO = 1 ;
    else if (abs(vTruth->at(0).Leptons_id().at(2))==13) multilepton_Lepton3_Id_LHCO = 2 ;

    //
    // LHCO phi convention
    //
    float multilepton_Lepton1_phi       = Phi_0_2Pi(vTruth->at(0).Leptons_phi().at(0));
    float multilepton_Lepton2_phi       = Phi_0_2Pi(vTruth->at(0).Leptons_phi().at(1));
    float multilepton_Lepton3_phi       = Phi_0_2Pi(vTruth->at(0).Leptons_phi().at(2));
    float multilepton_Bjet1_phi         = Phi_0_2Pi(vTruth->at(0).Bjets_phi().at(0));
    float multilepton_Bjet2_phi         = Phi_0_2Pi(vTruth->at(0).Bjets_phi().at(1));	  
    float multilepton_JetHighestPt1_phi = Phi_0_2Pi(vTruth->at(0).JetsHighestPt_phi().at(0));
    float multilepton_JetHighestPt2_phi = Phi_0_2Pi(vTruth->at(0).JetsHighestPt_phi().at(1));
    float multilepton_JetClosestMw1_phi = Phi_0_2Pi(vTruth->at(0).JetsClosestMw_phi().at(0));
    float multilepton_JetClosestMw2_phi = Phi_0_2Pi(vTruth->at(0).JetsClosestMw_phi().at(1));
    float multilepton_JetLowestMjj1_phi = Phi_0_2Pi(vTruth->at(0).JetsLowestMjj_phi().at(0));
    float multilepton_JetLowestMjj2_phi = Phi_0_2Pi(vTruth->at(0).JetsLowestMjj_phi().at(1));

    //std::cout <<" _MC 22"<< std::endl;

    // l1
    std::string l1_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d", 
                nobj,multilepton_Lepton1_Id_LHCO,vTruth->at(0).Leptons_eta().at(0),multilepton_Lepton1_phi,vTruth->at(0).Leptons_pt().at(0),0.0,vTruth->at(0).Leptons_id().at(0)/abs(vTruth->at(0).Leptons_id().at(0)),0,0,0,0));
    nobj++; 
    //std::cout <<" _MC 23"<< std::endl;
    // l2
    std::string l2_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
                nobj,multilepton_Lepton2_Id_LHCO,vTruth->at(0).Leptons_eta().at(1),multilepton_Lepton2_phi,vTruth->at(1).Leptons_pt().at(1),0.0,vTruth->at(0).Leptons_id().at(1)/abs(vTruth->at(0).Leptons_id().at(1)),0,0,0,0));
    nobj++;	
    //std::cout <<" _MC 24"<< std::endl;
    // l3
    std::string l3_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
                nobj,multilepton_Lepton3_Id_LHCO,vTruth->at(0).Leptons_eta().at(2),multilepton_Lepton3_phi,vTruth->at(2).Leptons_pt().at(2),0.0,vTruth->at(0).Leptons_id().at(2)/abs(vTruth->at(0).Leptons_id().at(2)),0,0,0,0));
    nobj++;
    //std::cout <<" _MC 3"<< std::endl;

    //										    
    std::string j1_fline;
    std::string j2_fline;

    if ( _processLHCO_MC == 5 || _processLHCO_MC == 6 || _processLHCO_MC == 4 || _processLHCO_MC == 3 )
    {			     
        // j1
        j1_fline = std::string(Form("%d	  %d	 %.2f	    %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
                    nobj,4,vTruth->at(0).JetsClosestMw_eta().at(0),multilepton_JetClosestMw1_phi,vTruth->at(0).JetsClosestMw_pt().at(0),0.0,1,0,0,0,0));
        nobj++;

        // j2
        j2_fline = std::string(Form("%d	  %d	  %.2f      %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
                    nobj,4,vTruth->at(0).JetsClosestMw_eta().at(1),multilepton_JetClosestMw2_phi,vTruth->at(0).JetsClosestMw_pt().at(1),0.0,1,0,0,0,0));
        nobj++;
    }
    else if ( _processLHCO_MC == 1 || _processLHCO_MC == 2)
    {
        // j1
        j1_fline = std::string(Form("%d	  %d	  %.2f      %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
                    nobj,4,vTruth->at(0).JetsLowestMjj_eta().at(0),multilepton_JetLowestMjj1_phi,vTruth->at(0).JetsLowestMjj_pt().at(0),0.0,1,0,0,0,0));
        nobj++;

        // j2
        j2_fline = std::string(Form("%d	  %d	   %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
                    nobj,4,vTruth->at(0).JetsLowestMjj_eta().at(1),multilepton_JetLowestMjj2_phi,vTruth->at(0).JetsLowestMjj_pt().at(1),0.0,1,0,0,0,0));
        nobj++;
    }

    //
    TLorentzVector BJet1;
    BJet1.SetPtEtaPhiE(vTruth->at(0).Bjets_pt().at(0), vTruth->at(0).Bjets_eta().at(0), vTruth->at(0).Bjets_phi().at(0), vTruth->at(0).Bjets_E().at(0));

    TLorentzVector BJet2;
    BJet2.SetPtEtaPhiE(vTruth->at(0).Bjets_pt().at(1), vTruth->at(0).Bjets_eta().at(1), vTruth->at(0).Bjets_phi().at(1), vTruth->at(0).Bjets_E().at(1));


    // bj1
    std::string b1_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
                nobj,4,vTruth->at(0).Bjets_eta().at(0),multilepton_Bjet1_phi,vTruth->at(0).Bjets_pt().at(0),BJet1.M(),1,2,0,0,0));
    nobj++;


    // bj2 
    std::string b2_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
                nobj,4,vTruth->at(0).Bjets_eta().at(1),multilepton_Bjet2_phi,vTruth->at(0).Bjets_pt().at(1),BJet2.M(),1,2,0,0,0));
    nobj++;

    // met
    std::string met_fline = std::string(Form("%d     %d	  %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d	  %d",
                nobj,6,0.,Phi_0_2Pi(vTruth->at(0).metGen_phi()),vTruth->at(0).metGen_pt(),0.,0,0,0,0,0));
    nobj++;


    //    
    fout_MC << fline00   << std::endl;
    fout_MC << fline0    << std::endl;
    fout_MC << l1_fline  << std::endl;
    fout_MC << l2_fline  << std::endl;
    fout_MC << l3_fline  << std::endl;	  
    if (!(_processLHCO_MC == 7 || _processLHCO_MC == 8)) fout_MC << j1_fline  << std::endl;// don't print jets for ttW hypothesis
    if (!(_processLHCO_MC == 7 || _processLHCO_MC == 8)) fout_MC << j2_fline  << std::endl;// don't print jets for ttW hypothesis
    fout_MC << b1_fline  << std::endl;
    fout_MC << b2_fline  << std::endl;
    fout_MC << met_fline << std::endl;      

}

void TTbarHiggsMultileptonAnalysis::PrintLHCOforMadweight_RECO(int evt)
{

    if (vSelectedLeptons.size()!=3 || vSelectedBTagJets.size()<2 || vSelectedJets.size()<4 ) return; 

    fline0 = "0      " + std::string(Form("%d",evt)) + "     " + trig;

    int nobj = 1;	      

    int multilepton_Lepton1_Id_LHCO = -666;
    int multilepton_Lepton2_Id_LHCO = -666;
    int multilepton_Lepton3_Id_LHCO = -666;

    //
    // LHCO lepton ID convention
    //  
    if(abs(multilepton_Lepton1_Id)==11) multilepton_Lepton1_Id_LHCO = 1 ;
    else if(abs(multilepton_Lepton1_Id)==13) multilepton_Lepton1_Id_LHCO = 2 ;
    if(abs(multilepton_Lepton2_Id)==11) multilepton_Lepton2_Id_LHCO = 1 ;
    else if(abs(multilepton_Lepton2_Id)==13) multilepton_Lepton2_Id_LHCO = 2 ;
    if(abs(multilepton_Lepton3_Id)==11) multilepton_Lepton3_Id_LHCO = 1 ;
    else if(abs(multilepton_Lepton3_Id)==13) multilepton_Lepton3_Id_LHCO = 2 ;

    //
    // LHCO phi convention
    //

    float multilepton_Lepton1_phi       = Phi_0_2Pi(multilepton_Lepton1_P4.Phi());
    float multilepton_Lepton2_phi       = Phi_0_2Pi(multilepton_Lepton1_P4.Phi());
    float multilepton_Lepton3_phi       = Phi_0_2Pi(multilepton_Lepton1_P4.Phi());
    float multilepton_Bjet1_phi         = Phi_0_2Pi(multilepton_Bjet1_P4.Phi());
    float multilepton_Bjet2_phi         = Phi_0_2Pi(multilepton_Bjet2_P4.Phi());	 
    float multilepton_JetHighestPt1_phi = Phi_0_2Pi(multilepton_JetHighestPt1_P4.Phi());
    float multilepton_JetHighestPt2_phi = Phi_0_2Pi(multilepton_JetHighestPt1_P4.Phi());
    float multilepton_JetClosestMw1_phi = Phi_0_2Pi(multilepton_JetClosestMw1_P4.Phi());
    float multilepton_JetClosestMw2_phi = Phi_0_2Pi(multilepton_JetClosestMw2_P4.Phi());
    float multilepton_JetLowestMjj1_phi = Phi_0_2Pi(multilepton_JetLowestMjj1_P4.Phi());
    float multilepton_JetLowestMjj2_phi = Phi_0_2Pi(multilepton_JetLowestMjj2_P4.Phi());


    // l1
    std::string l1_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d", nobj,multilepton_Lepton1_Id_LHCO,multilepton_Lepton1_P4.Eta(),multilepton_Lepton1_phi,multilepton_Lepton1_P4.Pt(),0.0,multilepton_Lepton1_Id/abs(multilepton_Lepton1_Id),0,0,0,0));
    nobj++; 

    // l2
    std::string l2_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d", nobj,multilepton_Lepton2_Id_LHCO,multilepton_Lepton2_P4.Eta(),multilepton_Lepton2_phi,multilepton_Lepton2_P4.Pt(),0.0,multilepton_Lepton2_Id/abs(multilepton_Lepton2_Id),0,0,0,0));
    nobj++;	

    // l3
    std::string l3_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d", nobj,multilepton_Lepton3_Id_LHCO,multilepton_Lepton3_P4.Eta(),multilepton_Lepton3_phi,multilepton_Lepton3_P4.Pt(),0.0,multilepton_Lepton3_Id/abs(multilepton_Lepton3_Id),0,0,0,0));
    nobj++;

    //										    
    std::string j1_fline;
    std::string j2_fline;

    if ( _processLHCO_RECO == 5 || _processLHCO_RECO == 6 || _processLHCO_RECO == 4 || _processLHCO_RECO == 3 )
    {			     
        // j1
        j1_fline = std::string(Form("%d	  %d	 %.2f	    %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
                    nobj,4,multilepton_JetClosestMw1_P4.Eta(),multilepton_JetClosestMw1_phi,multilepton_JetClosestMw1_P4.Pt(),0.0,1,0,0,0,0));
        nobj++;

        // j2
        j2_fline = std::string(Form("%d	  %d	  %.2f      %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
                    nobj,4,multilepton_JetClosestMw2_P4.Eta(),multilepton_JetClosestMw2_phi,multilepton_JetClosestMw2_P4.Pt(),0.0,1,0,0,0,0));
        nobj++;
    }
    else if ( _processLHCO_RECO == 1 || _processLHCO_RECO == 2 ) 
    {
        // j1
        j1_fline = std::string(Form("%d	  %d	  %.2f      %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
                    nobj,4,multilepton_JetLowestMjj1_P4.Eta(),multilepton_JetLowestMjj1_phi,multilepton_JetLowestMjj1_P4.Pt(),0.0,1,0,0,0,0));
        nobj++;

        // j2
        j2_fline = std::string(Form("%d	  %d	   %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
                    nobj,4,multilepton_JetLowestMjj2_P4.Eta(),multilepton_JetLowestMjj2_phi,multilepton_JetLowestMjj2_P4.Pt(),0.0,1,0,0,0,0));
        nobj++;
    }

    // for B-jet mass
    TLorentzVector BJet1;
    BJet1.SetPtEtaPhiE(multilepton_Bjet1_P4.Pt(), multilepton_Bjet1_P4.Eta(), multilepton_Bjet1_P4.Phi(), multilepton_Bjet1_P4.E());

    TLorentzVector BJet2;
    BJet2.SetPtEtaPhiE(multilepton_Bjet2_P4.Pt(), multilepton_Bjet2_P4.Eta(), multilepton_Bjet2_P4.Phi(), multilepton_Bjet2_P4.E());


    // bj1
    std::string b1_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
                nobj,4,multilepton_Bjet1_P4.Eta(),multilepton_Bjet1_phi,multilepton_Bjet1_P4.Pt(),BJet1.M(),1,2,0,0,0));
    nobj++;


    // bj2 
    std::string b2_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
                nobj,4,multilepton_Bjet2_P4.Eta(),multilepton_Bjet2_phi,multilepton_Bjet2_P4.Pt(),BJet2.M(),1,2,0,0,0));
    nobj++;

    // met
    std::string met_fline = std::string(Form("%d     %d	  %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d	  %d",
                nobj,6,0.,Phi_0_2Pi(multilepton_mET.Phi()),multilepton_mET.Pt(),0.,0,0,0,0,0));
    nobj++;


    //    
    fout_RECO << fline00	<< std::endl;
    fout_RECO << fline0	<< std::endl;
    fout_RECO << l1_fline  << std::endl;
    fout_RECO << l2_fline  << std::endl;
    fout_RECO << l3_fline  << std::endl;	
    if (!(_processLHCO_RECO == 7 || _processLHCO_RECO == 8)) fout_RECO << j1_fline  << std::endl;// don't print jets for ttW hypothesis
    if (!(_processLHCO_RECO == 7 || _processLHCO_RECO == 8)) fout_RECO << j2_fline  << std::endl;// don't print jets for ttW hypothesis
    fout_RECO << b1_fline  << std::endl;
    fout_RECO << b2_fline  << std::endl;
    fout_RECO << met_fline << std::endl;      

}

/*void TTbarHiggsMultileptonAnalysis::ProducePUweight()
  {

// TTbarHiggsMultileptonAnalysis::TTbarHiggsMultileptonAnalysis(
// TString inputFileName, TChain *tree, TString the_sampleName, TString treeName, TString outputFile, bool isdata, float xsec, float lumi, int nowe, int nmax)

std::cout << "Initializing PU from MC..." << std::endl;

const char *fname_str1  = "input_WZJets_MC.txt";
const char *stream_str1 = "Nt";
TString inputFileName1 = *fname_str1;
TChain *tree1 = 0;
TString treeName1 = *stream_str1;
int nmax1 = -1;

std::ifstream infile1;
infile1.open(inputFileName1);
std::string ifile1 = "";
while( getline(infile1, ifile1) )
{
std::string fnameStr1 = std::string(ifile1);
tree1->Add(fnameStr1.c_str());
std::cout << "file: " << fnameStr1 << std::endl;
}
infile1.close();

Init(tree1);

TH1F * PU_MC = new TH1F("PU_MC","PU_MC",100,0,99);
//TH1F PU_MC("PU_MC","PU_MC",100,0,99); 

Long64_t nentries1 = fChain->GetEntries();
int nentries_max1 = nentries1;

for (Long64_t jentry=0; jentry<nentries_max1;jentry++)
{
std::cout << "n[" << jentry << "] / " << nentries_max1 << std::endl;
int PU_M = 0;
PU_M = vEvent->at(0).pv_n();
std::cout << "Number of PV in MC: " << PU_M << std::endl;
PU_MC->Fill(PU_M,1); 
}

std::cout << "Initializing PU from Data..." << std::endl;

const char *fname_str2  = "input_DoubleEG_DATA.txt";
const char *stream_str2 = "Nt";
TString inputFileName2 = *fname_str2;
TChain *tree2 = 0;
TString treeName2 = *stream_str2;
int nmax2 = -1;

std::ifstream infile2;
infile2.open(inputFileName2);
std::string ifile2 = "";
while( getline(infile2, ifile2) )
{
std::string fnameStr2 = std::string(ifile2);
tree2->Add(fnameStr2.c_str());
std::cout << "file: " << fnameStr2 << std::endl;
}
infile2.close();

Init(tree2);

TH1F * PU_DATA = new TH1F("PU_DATA","PU_DATA",100,0,99);
//TH1F PU_DATA("PU_DATA","PU_DATA",100,0,99); 

Long64_t nentries2 = fChain->GetEntries();
int nentries_max2 = nentries2;

for (Long64_t jentry=0; jentry<nentries_max2;jentry++)
{
    if (fChain == 0) break;
    std::cout << "n[" << jentry << "] / " << nentries_max2 << std::endl;
    int PU_D = 0;
    //PU_D = vEvent->at(0).pv_n();
    std::cout << "Number of PV in DATA: " << PU_D << std::endl;
    PU_DATA->Fill(PU_D,1);
}

TH1F * PU_weight = new TH1F("PU_weight","PU_weight",100,0,99);
PU_DATA->Divide(PU_MC);
PU_weight = PU_DATA;

std::cout << "PU weight produced..." << std::endl;
}


float TTbarHiggsMultileptonAnalysis::PUweight()
{
    float weight = 1;
    return weight;
}*/

void TTbarHiggsMultileptonAnalysis::selectBjets(std::string BjetSel, int* ibsel1, int* ibsel2, bool doSelectOnlyBjets){

    //Selects the two highest b-tag jets. If only one b-tag select just this one.
    int ib1=-1, ib2=-1;

    if (BjetSel=="HighestBtagDiscrim"){
        float btag_max=-999, btag_max2=-999;
        for (unsigned int ib=0; ib<vSelectedJets.size(); ib++){
            if (doSelectOnlyBjets && (vSelectedJets.at(ib).CSVv2()<0.46)) continue;
            if (vSelectedJets.at(ib).CSVv2()>btag_max){
                btag_max2 = btag_max;
                ib2 = ib1;
                btag_max = vSelectedJets.at(ib).CSVv2();
                ib1 = ib;
            }
            if (vSelectedJets.at(ib).CSVv2()<btag_max && vSelectedJets.at(ib).CSVv2()>btag_max2){
                btag_max2 = vSelectedJets.at(ib).CSVv2();
                ib2 = ib;
            }
        }
    }
    if (BjetSel=="BtagHighestPt"){
        float pt_max=0, pt_max2=0;
        for (unsigned int ib=0; ib<vSelectedJets.size(); ib++){
            //            if (vSelectedJets.at(ib).CSVv2()<0.423) continue;
            if (vSelectedJets.at(ib).pt()>pt_max){
                pt_max2 = pt_max;
                ib2 = ib1;
                pt_max = vSelectedJets.at(ib).pt();
                ib1 = ib;
            }
            if (vSelectedJets.at(ib).pt()<pt_max && vSelectedJets.at(ib).pt()>pt_max2){
                pt_max2 = vSelectedJets.at(ib).pt();
                ib2 = ib;
            }
        }
    }

    *ibsel1 = ib1;
    *ibsel2 = ib2;

}

float TTbarHiggsMultileptonAnalysis::Phi_0_2Pi(float phi)
{
    float phi_0_2pi = phi;
    if (phi>= TMath::Pi()) phi_0_2pi  -= 2.*TMath::Pi();
    if (phi<0.)            phi_0_2pi  += 2.*TMath::Pi();
    return phi_0_2pi;
}

float TTbarHiggsMultileptonAnalysis::GetDeltaR(float eta1,float phi1,float eta2,float phi2)
{
    float DeltaPhi = TMath::Abs(phi2 - phi1);
    if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
    return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}
