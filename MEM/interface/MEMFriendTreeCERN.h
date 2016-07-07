#ifndef _MEMFriendTreeCERN_h
#define _MEMFriendTreeCERN_h

#include <ctime>
#include <map>
#include <string>
#include <vector>
#include <TLorentzVector.h>

struct MultiLepton;
struct ConfigParser;
struct HypIntegrator;
struct Permutations;


class  MEMFriendTreeCERN {
    public:
	MEMFriendTreeCERN() : multiLepton(nullptr), cfgParser(nullptr), hypIntegrator(nullptr) {}
        ~MEMFriendTreeCERN() ;
        MEMFriendTreeCERN(const MEMFriendTreeCERN &other) = delete;
        MEMFriendTreeCERN(MEMFriendTreeCERN &&other) = delete;
        MEMFriendTreeCERN & operator=(const MEMFriendTreeCERN &other) = delete;
        MEMFriendTreeCERN & operator=(MEMFriendTreeCERN &&other) = delete;

        void init(const std::string &config);

        void clear();
        void addLepton(const TLorentzVector &p4, int pdgId);
        void addJet(const std::string &what,  const TLorentzVector &p4, float CSV);
        void addBJet(const TLorentzVector &p4, float CSV);
        void setMET(const TLorentzVector &p4, double cov00, double cov01, double cov10, double cov11, double mHT) ;
        bool setCategory(const std::string &cat);
        std::map<std::string,float> compute() ;
    protected:
        MultiLepton  *multiLepton;
        ConfigParser *cfgParser;
        HypIntegrator *hypIntegrator;
        std::vector<Permutations*> MEMpermutations;

        // after config parsing
        int nhyp;
        std::string* shyp;
        int* hyp;
        int* nPointsHyp;
        int* index_hyp;

        int index[7];

        int doOptim;
        int doOptimTopHad, doOptimTopLep, doOptimHiggs, doOptimW;
        int doMinimization;
        std::string JetChoice;
        int nPermutationJetSyst;
        double xsTTH;
        double xsTTLL;
        double xsTTW;
        double xsTTbar;

        std::string index_CatJets[12];
        Double_t index_XS[7];
};


#endif
