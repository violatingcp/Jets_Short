#ifndef _interface_PFJetIDSelectionFunctorBasic_h
#define _interface_PFJetIDSelectionFunctorBasic_h

#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "PhysicsTools/SelectorUtils/interface/Selector.h"
 
#include <TMath.h>
class PFJetIDSelectionFunctorBasic : public Selector<reco::PFJet> {

 public: // interface
  
  enum Version_t { FIRSTDATA, N_VERSIONS };
  enum Quality_t { LOOSE, TIGHT, N_QUALITY};
  
  PFJetIDSelectionFunctorBasic() {}
  
  PFJetIDSelectionFunctorBasic( edm::ParameterSet const & params ) 
    {
      std::string versionStr = params.getParameter<std::string>("version");
      std::string qualityStr = params.getParameter<std::string>("quality");
      
      if ( versionStr == "FIRSTDATA" ) 
	version_ = FIRSTDATA;
      else
	version_ = FIRSTDATA;  /// will have other options eventually, most likely
      
      if      ( qualityStr == "LOOSE") quality_ = LOOSE;
      else if ( qualityStr == "TIGHT") quality_ = TIGHT;
      else quality_ = LOOSE;
      
      push_back("CHF" );
      push_back("NHF" );
      push_back("CEF" );
      push_back("NEF" );
      push_back("NCH" );
      push_back("nConstituents");
      
      // Set some default cuts for LOOSE, TIGHT
      if ( quality_ == LOOSE ) {
	set("CHF", 0.0);
	set("NHF", 0.99);
	set("CEF", 0.99);
	set("NEF", 0.99);
	set("NCH", 0);
	set("nConstituents", 1);
      } else if ( quality_ == TIGHT ) {
	set("CHF", 0.0);
	set("NHF", 0.9);
	set("CEF", 0.99);
	set("NEF", 0.9);
	set("NCH", 0);
	set("nConstituents", 1);      
      }
      
      
    // Now check the configuration to see if the user changed anything
      if ( params.exists("CHF") ) set("CHF", params.getParameter<double>("CHF") );
      if ( params.exists("NHF") ) set("NHF", params.getParameter<double>("NHF") );
      if ( params.exists("CEF") ) set("CEF", params.getParameter<double>("CEF") );
      if ( params.exists("NEF") ) set("NEF", params.getParameter<double>("NEF") );
      if ( params.exists("NCH") ) set("NCH", params.getParameter<int>   ("NCH") );
      if ( params.exists("nConstuents") ) set("nConstituents", params.getParameter<int> ("nConstituents") );
      
      if ( params.exists("cutsToIgnore") )
	setIgnoredCuts( params.getParameter<std::vector<std::string> >("cutsToIgnore") );
      
      indexNConstituents_ = index_type (&bits_, "nConstituents");
      indexNEF_ = index_type (&bits_, "NEF");
      indexNHF_ = index_type (&bits_, "NHF");
      indexCEF_ = index_type (&bits_, "CEF");
      indexCHF_ = index_type (&bits_, "CHF");
      indexNCH_ = index_type (&bits_, "NCH");
      
      retInternal_ = getBitTemplate();
     
    }
  
  PFJetIDSelectionFunctorBasic( Version_t version,
				Quality_t quality ) :
    version_(version), quality_(quality)
    {
      push_back("CHF" );
      push_back("NHF" );
      push_back("CEF" );
      push_back("NEF" );
      push_back("NCH" );
      push_back("nConstituents");
      
      
      // Set some default cuts for LOOSE, TIGHT
      if ( quality_ == LOOSE ) {
	set("CHF", 0.0);
	set("NHF", 0.99);
	set("CEF", 0.99);
	set("NEF", 0.99);
	set("NCH", 0);
	set("nConstituents", 1);
      } else if ( quality_ == TIGHT ) {
	set("CHF", 0.0);
	set("NHF", 0.9);
	set("CEF", 0.99);
	set("NEF", 0.9);
	set("NCH", 0);
	set("nConstituents", 1);      
      }
      
      indexNConstituents_ = index_type (&bits_, "nConstituents");
      indexNEF_ = index_type (&bits_, "NEF");
      indexNHF_ = index_type (&bits_, "NHF");
      indexCEF_ = index_type (&bits_, "CEF");
      indexCHF_ = index_type (&bits_, "CHF");
      indexNCH_ = index_type (&bits_, "NCH");
      
      retInternal_ = getBitTemplate();   
    }
                            
 
    using Selector<reco::PFJet>::operator();
    
    // 
    // Accessor from *CORRECTED* 4-vector, EMF, and Jet ID. 
    // This can be used with reco quantities. 
    // 
    bool operator()( reco::PFJet const & jet, 
		     pat::strbitset & ret )  
    {
      if ( version_ == FIRSTDATA ) return firstDataCuts( jet, ret );
      else {
	return false;
      }
    }
    
    // 
    // cuts based on craft 08 analysis. 
    // 
    bool firstDataCuts( reco::Jet const & jet,
			pat::strbitset & ret) 
    {    
      
      ret.set(false);
      
      // cache some variables
      double chf = 0.0;
      double nhf = 0.0;
      double cef = 0.0;
      double nef = 0.0;
      int    nch = 0;
      int    nconstituents = 0;
      
      reco::PFJet const    * pfJet    = dynamic_cast<reco::PFJet const *>(&jet);
      reco::BasicJet const * basicJet = dynamic_cast<reco::BasicJet const *>(&jet);
      
      if ( pfJet != 0 ) {
	chf = pfJet->chargedHadronEnergyFraction();
	nhf = ( pfJet->neutralHadronEnergy() + pfJet->HFHadronEnergy() ) / pfJet->energy();
	cef = pfJet->chargedEmEnergyFraction();
	nef = pfJet->neutralEmEnergyFraction();
	nch = pfJet->chargedMultiplicity();
	nconstituents = pfJet->numberOfDaughters();
      } // end if PF jet
      // Handle the special case where this is a composed jet for
      // subjet analyses
      else if ( basicJet != 0 ) {
	double e_chf = 0.0;
	double e_nhf = 0.0;
	double e_cef = 0.0;
	double e_nef = 0.0;
	nch = 0;
	nconstituents = 0;
	
	for ( reco::Jet::const_iterator ibegin = basicJet->begin(),
		iend = pfJet->end(), isub = ibegin; isub != iend; ++isub ) {
	  reco::PFJet const * pfsub = dynamic_cast<reco::PFJet const *>( &*isub );
	  e_cef += pfsub->chargedEmEnergy();
	  e_nef += pfsub->neutralEmEnergy();
	  nch += pfsub->chargedMultiplicity();
	  nconstituents += pfsub->numberOfDaughters();
	}
	double e = basicJet->energy();
	if ( e > 0.000001 ) {
	  chf = e_chf / e;
	  nhf = e_nhf / e;
	  cef = e_cef / e;
	  nef = e_nef / e;
	}
      }// end if basic jet
      
      // Cuts for all |eta|:
      if ( ignoreCut(indexNConstituents_) || nconstituents > cut(indexNConstituents_, int() ) ) passCut( ret, indexNConstituents_);
      if ( ignoreCut(indexNEF_)           || ( nef < cut(indexNEF_, double()) ) ) passCut( ret, indexNEF_);
      if ( ignoreCut(indexNHF_)           || ( nhf < cut(indexNHF_, double()) ) ) passCut( ret, indexNHF_);    
      // Cuts for |eta| < 2.4:
      if ( ignoreCut(indexCEF_)           || ( cef < cut(indexCEF_, double()) || std::abs(jet.eta()) > 2.4 ) ) passCut( ret, indexCEF_);
      if ( ignoreCut(indexCHF_)           || ( chf > cut(indexCHF_, double()) || std::abs(jet.eta()) > 2.4 ) ) passCut( ret, indexCHF_);
      if ( ignoreCut(indexNCH_)           || ( nch > cut(indexNCH_, int())    || std::abs(jet.eta()) > 2.4 ) ) passCut( ret, indexNCH_);    
      setIgnored( ret );
      return (bool)ret;
    }
      
    private: // member variables
	
	Version_t version_;
      Quality_t quality_;
      
      index_type indexNConstituents_;
      index_type indexNEF_;
      index_type indexNHF_;
      index_type indexCEF_;
      index_type indexCHF_;
      index_type indexNCH_;
      
    };
 
#endif
