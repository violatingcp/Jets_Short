// Implementation of template class: JetPlotsExample
// Description:  Example of simple EDAnalyzer for jets.
// Author: K. Kousouris, K. Mishra, S. Rappoccio
// Date:  25 - August - 2008
#include "RecoJets/JetAnalyzers/interface/JetPlotsExample.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/PileupJetIdentifier.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <TFile.h>
#include <cmath>
using namespace edm;
using namespace reco;
using namespace std;
////////////////////////////////////////////////////////////////////////////////////////
template<class Jet>
JetPlotsExample<Jet>::JetPlotsExample(edm::ParameterSet const& cfg)
{
  JetAlgorithm  = cfg.getParameter<std::string> ("JetAlgorithm"); 
  HistoFileName = cfg.getParameter<std::string> ("HistoFileName");
  NJets         = cfg.getParameter<int> ("NJets");
  JetPtMin      = cfg.getParameter<double> ("JetPtMin");
  useJecLevels  = cfg.exists("jecLevels");
  PUJetIdDisc   = cfg.getParameter<edm::InputTag> ("PUJetDiscriminant"); 
  PUJetId       = cfg.getParameter<edm::InputTag> ("PUJetId"); 
  if ( useJecLevels )
    jecLevels     = cfg.getParameter<std::string> ("jecLevels");
}
////////////////////////////////////////////////////////////////////////////////////////
template<class Jet>
void JetPlotsExample<Jet>::beginJob() 
{
  TString hname;
  m_file = new TFile(HistoFileName.c_str(),"RECREATE"); 
  /////////// Booking histograms //////////////////////////
  hname = "JetPt";
  m_HistNames1D[hname] = new TH1F(hname,hname,100,0,1000);
  hname = "JetEta";
  m_HistNames1D[hname] = new TH1F(hname,hname,120,-6,6);
  hname = "JetPhi";
  m_HistNames1D[hname] = new TH1F(hname,hname,100,-M_PI,M_PI);
  hname = "JetArea";
  m_HistNames1D[hname] = new TH1F(hname,hname,100,0.0,5.0);
  hname = "NumberOfJets";
  m_HistNames1D[hname] = new TH1F(hname,hname,100,0,100);
  hname = "ChargedHadronEnergyFraction";
  m_HistNames1D[hname] = new TH1F(hname,hname,100,0,1); 
  hname = "NeutralHadronEnergyFraction";
  m_HistNames1D[hname] = new TH1F(hname,hname,100,0,1); 
  hname = "PhotonEnergyFraction";
  m_HistNames1D[hname] = new TH1F(hname,hname,100,0,1); 
  hname = "ElectronEnergyFraction";
  m_HistNames1D[hname] = new TH1F(hname,hname,100,0,1); 
  hname = "MuonEnergyFraction";
  m_HistNames1D[hname] = new TH1F(hname,hname,100,0,1); 
  hname = "PUDiscriminant";
  m_HistNames1D[hname] = new TH1F(hname,hname,100,-1,1); 
  hname = "PUId";
  m_HistNames1D[hname] = new TH1F(hname,hname,10,0,10); 
}
////////////////////////////////////////////////////////////////////////////////////////
template<class Jet>
void JetPlotsExample<Jet>::analyze(edm::Event const& evt, edm::EventSetup const& iSetup) 
{
  /////////// Get the jet collection //////////////////////
  //Handle<JetCollection> jets;
   edm::Handle<edm::View<reco::Jet> > jets;
  evt.getByLabel(JetAlgorithm,jets);

  /////////// Get pu jet id map //////////////////////
  //Handle<JetCollection> jets;
  edm::Handle<edm::ValueMap<float> > PUJetIdMVA;
  edm::Handle<edm::ValueMap<int> >   PUJetIdFlag;
  if(PUJetIdDisc.label().size() != 0 && PUJetId.label().size() != 0) { 
    evt.getByLabel(PUJetIdDisc,PUJetIdMVA);
    evt.getByLabel(PUJetId    ,PUJetIdFlag);
  }
  /////////// Get Muons for Cleaning ////////////
  edm::Handle<CandidateView> leptons;
  evt.getByLabel("muons", leptons);

  // typename JetCollection::const_iterator i_jet;
  int index = 0;
  TString hname; 
  /////////// Count the jets in the event /////////////////
  hname = "NumberOfJets";
  FillHist1D(hname,jets->size()); 
  /////////// Fill Histograms for the leading NJet jets ///
  edm::View<reco::Jet>::const_iterator i_jet, endpjets = jets->end(); 
  for (i_jet = jets->begin();  i_jet != endpjets && index < NJets;  ++i_jet) {
    if(i_jet->pt() < JetPtMin) continue;
    bool pIsClean = true;
    //Remove all jets near a muon
    for ( CandidateView::const_iterator lepton = leptons->begin();
	   lepton != leptons->end(); ++lepton ) {
      if(lepton->pt() < 5.) continue;
      if(deltaR(i_jet->eta(),i_jet->phi(),lepton->eta(),lepton->phi()) < 0.5) pIsClean = false;
    }
    if(!pIsClean) continue;
    index++;
    double jecFactor = 1.0;
    double chf = 0.0;
    double nhf = 0.0;
    double pef = 0.0;
    double eef = 0.0;
    double mef = 0;
    
    edm::Ptr<reco::Jet> ptrToJet = jets->ptrAt( i_jet - jets->begin() );
    if ( ptrToJet.isNonnull() && ptrToJet.isAvailable() ) {
      reco::PFJet const * pfJet = dynamic_cast<reco::PFJet const *>( ptrToJet.get() );
      if ( pfJet != 0) { 
	chf = pfJet->chargedHadronEnergyFraction();
	nhf = pfJet->neutralHadronEnergyFraction();
	pef = pfJet->photonEnergyFraction();
	eef = pfJet->electronEnergy() / pfJet->energy();
	mef = pfJet->muonEnergyFraction();
	//if ( useJecLevels ) jecFactor = pfJet->jecFactor( jecLevels );
      }
      }
    bool passPU = true;
    float mva     = 0;
    int   idflag  = 0;
    if(PUJetIdDisc.label().size() != 0 && PUJetId.label().size() != 0) { 
      mva     = (*PUJetIdMVA) [ptrToJet];
      idflag = (*PUJetIdFlag)[ptrToJet];
      if(! PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose )) passPU  = false;
      }
    if(!passPU) continue;
    hname = "JetPt";
    FillHist1D(hname,(*i_jet).pt() * jecFactor );   
    hname = "JetEta";
    FillHist1D(hname,(*i_jet).eta());
    hname = "JetPhi";
    FillHist1D(hname,(*i_jet).phi());
    hname = "JetArea";
    FillHist1D(hname,(*i_jet).jetArea());
    hname = "ChargedHadronEnergyFraction";
    FillHist1D(hname, chf) ;
    hname = "NeutralHadronEnergyFraction";
    FillHist1D(hname, nhf);
    hname = "PhotonEnergyFraction";
    FillHist1D(hname, pef);
    hname = "ElectronEnergyFraction";
    FillHist1D(hname, eef);  
    hname = "MuonEnergyFraction";
    FillHist1D(hname, mef);  
    //PU Jet Id
    if(PUJetIdDisc.label().size() != 0 && PUJetId.label().size() != 0) { 
	hname = "PUDiscriminant";
	FillHist1D(hname, mva);  
	hname = "PUId";
	FillHist1D(hname, idflag);  
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////
template<class Jet>
void JetPlotsExample<Jet>::endJob() 
{
  /////////// Write Histograms in output ROOT file ////////
  if (m_file !=0) 
    {
      m_file->cd();
      for (std::map<TString, TH1*>::iterator hid = m_HistNames1D.begin(); hid != m_HistNames1D.end(); hid++)
        hid->second->Write();
      delete m_file;
      m_file = 0;      
    }
}
////////////////////////////////////////////////////////////////////////////////////////
template<class Jet>
void JetPlotsExample<Jet>::FillHist1D(const TString& histName,const Double_t& value) 
{
  std::map<TString, TH1*>::iterator hid=m_HistNames1D.find(histName);
  if (hid==m_HistNames1D.end())
    std::cout << "%fillHist -- Could not find histogram with name: " << histName << std::endl;
  else
    hid->second->Fill(value);
}
/////////// Register Modules ////////
#include "FWCore/Framework/interface/MakerMacros.h"
/////////// Calo Jet Instance ////////
typedef JetPlotsExample<CaloJet> CaloJetPlotsExample;
DEFINE_FWK_MODULE(CaloJetPlotsExample);
/////////// Cen Jet Instance ////////
typedef JetPlotsExample<GenJet> GenJetPlotsExample;
DEFINE_FWK_MODULE(GenJetPlotsExample);
/////////// PF Jet Instance ////////
typedef JetPlotsExample<PFJet> PFJetPlotsExample;
DEFINE_FWK_MODULE(PFJetPlotsExample);
/////////// JPT Jet Instance ////////
typedef JetPlotsExample<JPTJet> JPTJetPlotsExample;
DEFINE_FWK_MODULE(JPTJetPlotsExample);
