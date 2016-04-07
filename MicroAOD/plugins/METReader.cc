// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
//
// class declaration
//

class METReader : public edm::EDAnalyzer {
   public:
      explicit METReader(const edm::ParameterSet&);
      ~METReader();

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::METCollection> metToken_;
};

METReader::METReader(const edm::ParameterSet& iConfig):
    metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets")))
{
}

METReader::~METReader()
{
}


void
METReader::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
    edm::Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken_, mets);
    const pat::MET &met = mets->front();
    std::cout<<"************************************"<<std::endl;
    printf(" MET: pt %5.1f, pt raw %5.1f, pt t1 %5.1f, pt t1p2 %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f.\n MET with JES up/down: %.1f/%.1f\n MET with JER up/down: %.1f/%.1f\n MET with MUON up/down: %.1f/%.1f\n MET with ELECTRON up/down: %.1f/%.1f\n MET with TAU up/down: %.1f/%.1f\n MET with UNCL up/down: %.1f/%.1f\n",
	   met.pt(), met.corPt(pat::MET::Raw), met.corPt(pat::MET::Type1),met.corPt(pat::MET::Type1XY),  met.phi(), met.sumEt(), 
	   met.genMET()->pt(),
	   met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown),
	   met.shiftedPt(pat::MET::JetResUp), met.shiftedPt(pat::MET::JetResDown),
	   met.shiftedPt(pat::MET::MuonEnUp), met.shiftedPt(pat::MET::MuonEnDown),
	   met.shiftedPt(pat::MET::ElectronEnUp), met.shiftedPt(pat::MET::ElectronEnDown),
	   met.shiftedPt(pat::MET::TauEnUp), met.shiftedPt(pat::MET::TauEnDown),
	   met.shiftedPt(pat::MET::UnclusteredEnUp), met.shiftedPt(pat::MET::UnclusteredEnDown));


    


    printf("\n");
}

//define this as a plug-in
DEFINE_FWK_MODULE(METReader);
