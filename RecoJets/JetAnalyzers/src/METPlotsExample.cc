#include "RecoJets/JetAnalyzers/interface/METPlotsExample.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <TFile.h>
#include <TLorentzVector.h>
#include <cmath>
using namespace edm;
using namespace reco;
using namespace std;
////////////////////////////////////////////////////////////////////////////////////////
template<class MET>
METPlotsExample<MET>::METPlotsExample(edm::ParameterSet const& cfg)
{
  METAlgorithm  = cfg.getParameter<std::string> ("MetAlgorithm"); 
  HistoFileName = cfg.getParameter<std::string> ("HistoFileName");
  useRecoil     = cfg.getParameter<bool>        ("useRecoil");  
  recoilLeptons = cfg.getParameter<std::string> ("recoilLeptons"); 
}
////////////////////////////////////////////////////////////////////////////////////////
template<class MET>
void METPlotsExample<MET>::beginJob() 
{
  TString hname;
  m_file = new TFile(HistoFileName.c_str(),"RECREATE"); 
  /////////// Booking histograms //////////////////////////
  hname = "MET";
  m_HistNames1D[hname] = new TH1F(hname,hname,50,0,200);
  hname = "METPhi";
  m_HistNames1D[hname] = new TH1F(hname,hname,40,-M_PI,M_PI);
  hname = "METX";
  m_HistNames1D[hname] = new TH1F(hname,hname,50,-100,100);
  hname = "METY";
  m_HistNames1D[hname] = new TH1F(hname,hname,50,-100,100);
  hname = "sumEt";
  m_HistNames1D[hname] = new TH1F(hname,hname,50,0,3000);
  hname = "UPara";
  m_HistNames1D[hname] = new TH1F(hname,hname,20,-100,200);
  hname = "UPerp";
  m_HistNames1D[hname] = new TH1F(hname,hname,20,-100,100);
  hname = "UParaPlusPt";
  m_HistNames1D[hname] = new TH1F(hname,hname,20,-100,100);
  hname = "Response";
  m_HistNames1D[hname] = new TH1F(hname,hname,10,-0.5,2);
}
////////////////////////////////////////////////////////////////////////////////////////
template<class MET>
void METPlotsExample<MET>::analyze(edm::Event const& evt, edm::EventSetup const& iSetup) 
{
  /////////// Get the MET collection //////////////////////
  edm::Handle<edm::View<reco::MET> > METs;
  evt.getByLabel(METAlgorithm,METs);

  /////////// Get the Leptonic Recoil //////////////////////
  TLorentzVector recoil;
  if(useRecoil) { 
    edm::Handle<CandidateView> leptons;
    evt.getByLabel(recoilLeptons, leptons);
    std::vector<TLorentzVector> leptonVecs;
    int massId1  = -1; 
    int massId2  = -1; 
    int leptonId = 0;
    double lMass = 0.;
    for ( CandidateView::const_iterator lepton = leptons->begin();
          lepton != leptons->end(); ++lepton ) {
      if(lepton->pt() < 10.) continue;
      TLorentzVector pVec; pVec.SetPtEtaPhiM(lepton->pt(),lepton->eta(),lepton->phi(),lepton->mass());
      for(unsigned int i0 = 0; i0 < leptonVecs.size(); i0++) {
	double pMass = (pVec+leptonVecs[i0]).M();
	if(pMass > lMass) {
	  massId1 = leptonId; 
	  massId2 = i0; 
	  lMass = pMass;
	}
      }	
      leptonVecs.push_back(pVec);
      leptonId++;
    }
    if(massId1 > -1) recoil+=leptonVecs[massId1];
    if(massId2 > -1) recoil+=leptonVecs[massId2];
  }
  int index = 0;
  TString hname; 
  /////////// Count the Met in the event /////////////////
  edm::View<reco::MET>::const_iterator i_met, endpmets = METs->end(); 
  for (i_met = METs->begin();  i_met != endpmets;  ++i_met) {
    hname = "MET";
    FillHist1D(hname,(*i_met).pt()  );   
    hname = "METPhi";
    FillHist1D(hname,(*i_met).phi());
    hname = "METX";
    FillHist1D(hname,(*i_met).pt()*cos((*i_met).phi()));
    hname = "METY";
    FillHist1D(hname,(*i_met).pt()*sin((*i_met).phi()));
    hname = "sumEt";
    FillHist1D(hname,(*i_met).sumEt());
    TLorentzVector pU;
    pU.SetPtEtaPhiM(i_met->pt(),0.,i_met->phi(),0.);
    pU+=recoil;
    double deltaPhi = (recoil.Phi()-pU.Phi()); 
    if(deltaPhi >  M_PI) deltaPhi -= 2.*M_PI;
    if(deltaPhi < -M_PI) deltaPhi += 2.*M_PI;
    hname = "UPara";
    FillHist1D(hname,pU.Pt()*cos(deltaPhi));
    hname = "UPerp";
    FillHist1D(hname,pU.Pt()*sin(deltaPhi));
    hname = "UParaPlusPt";
    FillHist1D(hname,pU.Pt()*cos(deltaPhi)-recoil.Pt());
    hname = "Response";
    if(recoil.Pt() > 40) FillHist1D(hname,(pU.Pt()*cos(deltaPhi)/recoil.Pt()));
    index++;
  }
}
////////////////////////////////////////////////////////////////////////////////////////
template<class MET>
void METPlotsExample<MET>::endJob() 
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
template<class MET>
void METPlotsExample<MET>::FillHist1D(const TString& histName,const Double_t& value) 
{
  std::map<TString, TH1*>::iterator hid=m_HistNames1D.find(histName);
  if (hid==m_HistNames1D.end())
    std::cout << "%fillHist -- Could not find histogram with name: " << histName << std::endl;
  else
    hid->second->Fill(value);
}
/////////// Register Modules ////////
#include "FWCore/Framework/interface/MakerMacros.h"
/////////// PF MET Instance ////////
typedef METPlotsExample<PFMET> PFMETPlotsExample;
DEFINE_FWK_MODULE(PFMETPlotsExample);
typedef METPlotsExample<CaloMET> CaloMETPlotsExample;
DEFINE_FWK_MODULE(CaloMETPlotsExample);
