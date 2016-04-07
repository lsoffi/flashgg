#include "flashgg/Taggers/interface/LeptonSelection.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "TLorentzVector.h"

using namespace std;
using namespace edm;


namespace flashgg {

    std::vector<edm::Ptr<flashgg::Muon> > selectMuons( const std::vector<edm::Ptr<flashgg::Muon> > &muonPointers, Ptr<flashgg::DiPhotonCandidate> dipho,
            const std::vector<edm::Ptr<reco::Vertex> > &vertexPointers, double muonEtaThreshold, double muonPtThreshold, double muPFIsoSumRelThreshold,
                                                       double dRPhoLeadMuonThreshold, double dRPhoSubLeadMuonThreshold )
    {

        std::vector<edm::Ptr<flashgg::Muon> > goodMuons;

        for( unsigned int muonIndex = 0; muonIndex < muonPointers.size(); muonIndex++ ) {
            Ptr<flashgg::Muon> muon = muonPointers[muonIndex];

            /*
            std::cout << " Muon index " << muonIndex << " has pt eta weight: "
                      << muon->pt() << " " << muon->eta() << " "
                      << muon->centralWeight() << std::endl;
            auto weightList = muon->weightList();
            for( unsigned int i = 0 ; i < weightList.size() ; i++ ) {
                std::cout << "    " << weightList[i] << " " << muon->weight( weightList[i] );
            }
            std::cout << std::endl;
            */

            if( fabs( muon->eta() ) > muonEtaThreshold ) { continue; }
            if( muon->pt() < muonPtThreshold ) { continue; }

            int vtxInd = 0;
            double dzmin = 9999;

            for( size_t ivtx = 0 ; ivtx < vertexPointers.size(); ivtx++ ) {

                Ptr<reco::Vertex> vtx = vertexPointers[ivtx];

                if( !muon->innerTrack() ) { continue; }

                if( fabs( muon->innerTrack()->vz() - vtx->position().z() ) < dzmin ) {

                    dzmin = fabs( muon->innerTrack()->vz() - vtx->position().z() );

                    vtxInd = ivtx;
                }

            }

            Ptr<reco::Vertex> best_vtx = vertexPointers[vtxInd];

            //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Tight_Muon and https://cmssdt.cern.ch/SDT/lxr/source/RecoBTag/SoftLepton/plugins/SoftPFMuonTagInfoProducer.cc#0135

            if( !muon::isTightMuon( *muon, *best_vtx ) ) continue; 

            //I = [sumChargedHadronPt+ max(0.,sumNeutralHadronPt+sumPhotonPt-0.5sumPUPt]/pt
            //https://cmssdt.cern.ch/SDT/doxygen/CMSSW_5_3_14/doc/html/df/d33/structreco_1_1MuonPFIsolation.html
            // deltaBeta correction

            double muPFIsoSumRel = ( muon->pfIsolationR04().sumChargedHadronPt + max( 0.,
                                     muon->pfIsolationR04().sumNeutralHadronEt + muon->pfIsolationR04().sumPhotonEt - 0.5 * muon->pfIsolationR04().sumPUPt ) ) / ( muon->pt() );
            if( muPFIsoSumRel > muPFIsoSumRelThreshold ) { continue; }

            float dRPhoLeadMuon = deltaR( muon->eta(), muon->phi(), dipho->leadingPhoton()->superCluster()->eta(), dipho->leadingPhoton()->superCluster()->phi() ) ;
            float dRPhoSubLeadMuon = deltaR( muon->eta(), muon->phi(), dipho->subLeadingPhoton()->superCluster()->eta(), dipho->subLeadingPhoton()->superCluster()->phi() );

            //https://github.com/h2gglobe/h2gglobe/blob/master/PhotonAnalysis/src/PhotonAnalysis.cc#L5369
            if( dRPhoLeadMuon < dRPhoLeadMuonThreshold || dRPhoSubLeadMuon < dRPhoSubLeadMuonThreshold ) { continue; }

            goodMuons.push_back( muon );
        }
        return goodMuons;
    }

    std::vector<edm::Ptr<Electron> > selectElectrons( const std::vector<edm::Ptr<flashgg::Electron> > &ElectronPointers,
            const std::vector<edm::Ptr<reco::Vertex> > &vertexPointers , double ElectronPtThreshold, double DeltaRTrkElec, double TransverseImpactParam,
            double LongitudinalImapctParam, double NonTrigMVAThreshold, double IsoThreshold, double NumOfMissingHitsThreshold, vector<double> EtaCuts )
    {

        std::vector<edm::Ptr<flashgg::Electron> > goodElectrons;

        // std::cout << " LC DEBUG (LeptonSlection.cc) nElectrons " << ElectronPointers.size() << std::endl;
        for( unsigned int ElectronIndex = 0; ElectronIndex < ElectronPointers.size(); ElectronIndex++ ) {

            Ptr<flashgg::Electron> Electron = ElectronPointers[ElectronIndex];

            /*
            std::cout << " LeptonSelection Electron index " << ElectronIndex << " has pt eta weight: "
                      << Electron->pt() << " " << Electron->eta() << " "
                      << Electron->centralWeight() << std::endl;
            auto weightList = Electron->weightList();
            for( unsigned int i = 0 ; i < weightList.size() ; i++ ) {
                std::cout << "    " << weightList[i] << " " << Electron->weight( weightList[i] );
            }
            std::cout << std::endl;
            */

            float Electron_eta = fabs( Electron->superCluster()->eta() );
            Ptr<reco::Vertex> Electron_vtx = chooseElectronVertex( Electron, vertexPointers );
            float dxy = Electron->gsfTrack()->dxy( Electron_vtx->position() );
            float dz = Electron->gsfTrack()->dz( Electron_vtx->position() );
            //std::cout << "[DEBUG] LeptonSelection.cc Electron eta " << Electron_eta <<  "  -- cuts <"<< EtaCuts[0] << ", <" << EtaCuts[1] << ", >" << EtaCuts[2] << std::endl;
            //std::cout << "[DEBUG] LeptonSelection.cc Electron pt " << Electron->pt() <<", pt cut  >" << ElectronPtThreshold<<  std::endl;
            //std::cout << "[DEBUG] LeptonSelection.cc non trig MVA " << Electron->nonTrigMVA() <<", mva cut  >" << NonTrigMVAThreshold <<  std::endl;
            //std::cout << "[DEBUG] LeptonSelection.cc iso" << Electron->standardHggIso()  <<", iso cut  <" << IsoThreshold <<  std::endl;
            //std::cout << "[DEBUG] LeptonSelection.cc missing Hits" << Electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS )  <<", cut  <" << NumOfMissingHitsThreshold <<  std::endl;
            //std::cout << "[DEBUG] LeptonSelection.cc hasConversion " << Electron->hasMatchedConversion() << " needed 0  to pass " << std::endl;
            //std::cout << "[DEBUG] LeptonSelection.cc dxy " << dxy << " -- cut  < "<<TransverseImpactParam << std::endl;
            //std::cout << "[DEBUG] LeptonSelection.cc dz " << dz << " -- cut  < "<< LongitudinalImapctParam << std::endl;

            if( Electron_eta > EtaCuts[2] && ( Electron_eta > EtaCuts[0] && Electron_eta < EtaCuts[1] ) ) { continue; }
            //fixed the ||  to become and && in above, otherwise all electrons above EtaCuts[2] (which is ~1.44) are thrown out.
            //ie it would have ignored all electrons on the endcaps.. 
            //std::cout << "[DEBUG] passed eta cuts. " << std::endl;
            if( Electron->pt() < ElectronPtThreshold ) { continue; }
            //std::cout << "[DEBUG] passed pt cuts." << std::endl;



            if( Electron->nonTrigMVA() < NonTrigMVAThreshold ) { continue; }
            //std::cout << "[DEBUG] passed non trig mva." << std::endl;
         
            //if( Electron->standardHggIso() > IsoThreshold ) { continue; } //FIXME
            //std::cout << "[DEBUG] passed iso." << std::endl; //FIXME
            //The default standardHggIso is delivering values which are an order or magnitude larger than the cut. Maybe need to add rho correction? 
            std::cout << "[WARNING] Isolation cut for flashggElectrons temporarily disabled pending optimisation." << std::endl; //FIXME


            if( Electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS ) > NumOfMissingHitsThreshold ) { continue; }
            //std::cout << "[DEBUG] passed  missing hits cut." << std::endl;

            if( Electron->hasMatchedConversion() ) { continue; }
            //std::cout << "[DEBUG] passed hasConversion." << std::endl;

            if( fabs(dxy) > TransverseImpactParam ) { continue; }
            //std::cout << "[DEBUG] passed dxy cut " << std::endl;
            if( fabs(dz) > LongitudinalImapctParam ) { continue; }
            //std::cout << "[DEBUG] passed dz cut." << std::endl;

            //std::cout << "[DEBUG] LeptonSelection.c   ... pushing back Electron index " << ElectronIndex << std::endl;
            goodElectrons.push_back( Electron );
        }

        return goodElectrons;
    }

    std::vector<edm::Ptr<Electron> > selectMediumElectrons( const std::vector<edm::Ptr<flashgg::Electron> > &ElectronPointers,
                                                            const std::vector<edm::Ptr<reco::Vertex> > &vertexPointers,
                                                            Ptr<flashgg::DiPhotonCandidate> dipho,  
                                                            float rho,
                                                            double ElectronPtThreshold,
                                                            double dRPhoLeadEleThreshold, double dRPhoSubLeadEleThreshold )
    {

        std::vector<edm::Ptr<flashgg::Electron> > goodElectrons;

        for( unsigned int ElectronIndex = 0; ElectronIndex < ElectronPointers.size(); ElectronIndex++ ) {
            
            Ptr<flashgg::Electron> Electron = ElectronPointers[ElectronIndex];

            float Electron_eta = fabs( Electron->superCluster()->eta() );
            if( Electron_eta > 2.5) continue; 
            if( Electron_eta > 1.442 && Electron_eta< 1.566) continue;   

            if( Electron->pt() < ElectronPtThreshold ) continue;   

            // hardcoded cuts
            // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns
            float sieie_cut = 0.0101;
            float detacut = 0.0103;
            float dphicut = 0.0336;
            float hoecut = 0.0876;
            float relisocut = 0.0766;
            float ooemoopcut = 0.0174;
            float dxycut = 0.0118;
            float dzcut = 0.373;
            int misshcut = 2;  
            // + conversion veto 
            if (Electron_eta>1.5) {
                sieie_cut = 0.0283;
                detacut = 0.00733;
                dphicut = 0.114;
                hoecut = 0.0678;
                relisocut = 0.0678;
                ooemoopcut = 0.0898;
                dxycut = 0.0739;
                dzcut = 0.602;
                misshcut = 1;  
            }

            // ID selection
            float HoE          = Electron->hcalOverEcal();   
            float DeltaPhiIn   = fabs(Electron->deltaPhiSuperClusterTrackAtVtx());    
            float DeltaEtaIn   = fabs(Electron->deltaEtaSuperClusterTrackAtVtx());          
            float Full5x5Sieie = fabs(Electron->full5x5_sigmaIetaIeta());   
            float ecalEne      = Electron->ecalEnergy();  
            float OneOverEoP;    
            if (ecalEne==0) {   
                OneOverEoP = 1000000.; 
            } else {
                OneOverEoP = fabs(1.0/ecalEne - (Electron->eSuperClusterOverP())/ecalEne);     
            }
            bool passHoE     = HoE<hoecut;
            bool passDPhi    = DeltaPhiIn<dphicut;
            bool passDEta    = DeltaEtaIn<detacut;
            bool passSieie   = Full5x5Sieie<sieie_cut;
            bool passOoemoop = OneOverEoP<ooemoopcut;
            bool passId = passHoE && passDPhi && passDEta && passSieie && passOoemoop;
            
            // Isolation 
            float scEta = fabs(Electron->superCluster()->eta());
            reco::GsfElectron::PflowIsolationVariables pfIso = Electron->pfIsolationVariables();    
            float corrHadPlusPho = pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho*effectiveAreaEle03(scEta);   
            if (corrHadPlusPho<=0) corrHadPlusPho = 0.; 
            float absIsoWeffArea = pfIso.sumChargedHadronPt + corrHadPlusPho; 
            float relIso = absIsoWeffArea/(Electron->pt());         
            bool passIso = relIso<relisocut;

            // IP
            Ptr<reco::Vertex> Electron_vtx = chooseElectronVertex( Electron, vertexPointers );
            float dxy = fabs(Electron->gsfTrack()->dxy(Electron_vtx->position()));
            float dz  = fabs(Electron->gsfTrack()->dz(Electron_vtx->position()));
            bool passdz  = dz<dzcut;
            bool passdxy = dxy<dxycut;
            bool passIp  = passdz && passdxy;

            // conversion rejection
            int missHits   = Electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS );
            bool matchConv = Electron->hasMatchedConversion();
            bool passmhits = missHits<=misshcut; 
            bool passConv  = !matchConv && passmhits;

            // all
            bool selected = passId && passIp && passConv && passIso;
            if (!selected) continue;

            // should be far from the photons
            float scPhi = fabs(Electron->superCluster()->phi());
            float dRPhoLeadEle = deltaR( scEta, scPhi, dipho->leadingPhoton()->superCluster()->eta(), dipho->leadingPhoton()->superCluster()->phi() ) ;
            float dRPhoSubLeadEle = deltaR( scEta, scPhi, dipho->subLeadingPhoton()->superCluster()->eta(), dipho->subLeadingPhoton()->superCluster()->phi() );
            if( dRPhoLeadEle < dRPhoLeadEleThreshold || dRPhoSubLeadEle < dRPhoSubLeadEleThreshold ) continue;

            goodElectrons.push_back( Electron );
        }
        
        return goodElectrons;
    }


    Ptr<reco::Vertex>  chooseElectronVertex( Ptr<flashgg::Electron> &elec, const std::vector<edm::Ptr<reco::Vertex> > &vertices )
    {
        double vtx_dz = 1000000;
        unsigned int min_dz_vtx = -1;
        for( unsigned int vtxi = 0; vtxi < vertices.size(); vtxi++ ) {

            Ptr<reco::Vertex> vtx = vertices[vtxi];

            if( fabs( vtx_dz ) > fabs( elec->gsfTrack()->dz( vtx->position() ) ) ) {
                
                vtx_dz = elec->gsfTrack()->dz( vtx->position() );
                min_dz_vtx = vtxi;
            }
        }
        return vertices[min_dz_vtx];
    }

    float effectiveAreaEle03(float theEta) {
        
        float theEA = -999;
        if(fabs(theEta) < 1) theEA = 0.1752;
        else if(fabs(theEta) < 1.479) theEA = 0.1862;
        else if(fabs(theEta) < 2.0) theEA = 0.1411;
        else if(fabs(theEta) < 2.2) theEA = 0.1534;
        else if(fabs(theEta) < 2.3) theEA = 0.1903;
        else if(fabs(theEta) < 2.4) theEA = 0.2243;
        else if(fabs(theEta) < 2.5) theEA = 0.2687;

        return theEA;
    }
}
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

