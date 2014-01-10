#include "PhysicsTools/UtilAlgos/interface/EDFilterObjectWrapper.h"
#include "RecoJets/JetAnalyzers/interface/PFJetIDSelectionFunctorBasic.h"
typedef edm::FilterObjectWrapper<PFJetIDSelectionFunctorBasic, std::vector<reco::PFJet> > PFJetIDSelectionFunctorBasicFilter;
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFJetIDSelectionFunctorBasicFilter);
