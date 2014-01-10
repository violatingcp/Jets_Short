// Template class: MetPlotsExample
// Description:  Example of simple EDAnalyzer for jets.
// Author: K. Kousouris
// Date:  25 - August - 2008
#ifndef JetPlotsExample_h
#define JetPlotsExample_h
#include <TH1.h>
#include <TFile.h>
#include "TNamed.h"
#include <vector>
#include <map>
#include "FWCore/Framework/interface/EDAnalyzer.h"

template<class MET>
class METPlotsExample : public edm::EDAnalyzer 
   {
     public:
       METPlotsExample(edm::ParameterSet const& cfg);
     private:
       typedef std::vector<MET> METCollection;
       void FillHist1D(const TString& histName, const Double_t& x);
       void beginJob();
       void analyze(edm::Event const& e, edm::EventSetup const& iSetup);
       void endJob();
       std::map<TString, TH1*> m_HistNames1D;  
       TFile* m_file;
       /////// Configurable parameters /////////////////////////////////////
       /////// MET algorithm: it can be any  PF Gen algorithm //////
       std::string METAlgorithm;
       /////// Leptons to determine the Recoil                //////
       std::string recoilLeptons;
       /////// use leptons to determint the recoil            //////
       bool useRecoil;
       /////// Histogram where the plots are stored //////////////////////// 
       std::string HistoFileName;
   };
#endif
