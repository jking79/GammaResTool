// -*- C++ -*-
//
// Package:    GammaResTool
// Class:      GammaResTool
// gammaResTool_cc.wip
/**\class GammaResTool GammaResTool.cc LLPgammaAnalyzer/plugins/GammaResTool.cc

 		Description: [one line class summary]

 		Implementation: [Notes on implementation]

*/
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//

//----------------------------------------  cc file   --------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------

#include "GammaResTool/gammaResTool/plugins/gammaResTool.hh"
using namespace std;

//#define DEBUG true
#define DEBUG false

//
// constructors and destructor
//
GammaResTool::GammaResTool(const edm::ParameterSet& iConfig) :

// -- declare tags ----------------------------------------------------------

	// flags
	doTwoTier(iConfig.existsAs<bool>("doTwoTier")  ? iConfig.getParameter<bool>("doTwoTier")  : false),

	// tracks
	tracksTag(iConfig.getParameter<edm::InputTag>("tracks")),

    // vertices
    verticesTag(iConfig.getParameter<edm::InputTag>("vertices")),
	
	// electrons
	electronsTag(iConfig.getParameter<edm::InputTag>("electrons")),  

	// recHits
	recHitsEBTag(iConfig.getParameter<edm::InputTag>("recHitsEB")),  
	recHitsEETag(iConfig.getParameter<edm::InputTag>("recHitsEE")),

    kuCCStcRecHitsEBTag(iConfig.getParameter<edm::InputTag>("kuCCStcRecHitsEB")),
    kuCCStcRecHitsEETag(iConfig.getParameter<edm::InputTag>("kuCCStcRecHitsEE")),

	// gedphotons
	gedPhotonsTag(iConfig.getParameter<edm::InputTag>("gedPhotons")),

	// ootPhotons
	ootPhotonsTag(iConfig.getParameter<edm::InputTag>("ootPhotons")),

	// ECAL RECORDS 
    caloGeometryToken_(esConsumes()),
    ecalLaserDbServiceToken_(esConsumes()),
    ecalIntercalibConstantsToken_(esConsumes()),
    ecalADCToGeVConstantToken_(esConsumes()),
    EcalPedestalsToken_(esConsumes())


// -- end of tag declarations ---------------------------------------
{ //<<<< GammaResTool::GammaResTool(const edm::ParameterSet& iConfig) :

	usesResource();
	usesResource("TFileService");

// -- consume tags ------------------------------------------------------------
	if( DEBUG ) std::cout << "In constructor for GammaResTool - tag and tokens" << std::endl;

	// tracks 
	tracksToken_ = consumes<std::vector<reco::Track>>(tracksTag);

	// vertices
	verticesToken_ = consumes<std::vector<reco::Vertex>>(verticesTag);

	// leptons
	electronsToken_	= consumes<std::vector<pat::Electron>>(electronsTag);

	// rechits
	recHitsEBToken_	= consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(recHitsEBTag);
	recHitsEEToken_	= consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(recHitsEETag);

	if( doTwoTier ){
    kuCCStcRecHitsEBToken_ = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(kuCCStcRecHitsEBTag);
    kuCCStcRecHitsEEToken_ = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(kuCCStcRecHitsEETag);
	}//if( doTwoTier )

	// photons
	gedPhotonsToken_ = consumes<std::vector<pat::Photon>>(gedPhotonsTag);
	ootPhotonsToken_ = consumes<std::vector<pat::Photon>>(ootPhotonsTag);

// ---------------------------------------------------------------------------------

}//>>>>GammaResTool::GammaResTool(const edm::ParameterSet& iConfig)


GammaResTool::~GammaResTool(){
	///////////////////////////////////////////////////////////////////
	// do anything here that needs to be done at desctruction time   //
	// (e.g. close files, deallocate resources etc.)                 //
	///////////////////////////////////////////////////////////////////
}//>>>>GammaResTool::~GammaResTool()


//
// member functions
//

float GammaResTool::getPhotonSeedTime( pat::Photon photon ){

	const auto & phosc = photon.superCluster().isNonnull() ? photon.superCluster() : photon.parentSuperCluster();
	const auto & seedDetId = phosc->seed()->seed(); // get seed detid
	const auto recHits = ((seedDetId.subdetId() == EcalSubdetector::EcalBarrel) ? recHitsEB_ : recHitsEE_); // which recHits to use
	const auto seedHit = recHits->find(seedDetId); // get the underlying rechit
	const auto seedTime = ((seedHit != recHits->end()) ? seedHit->time() : -9999.f);
	if( DEBUG && seedTime == -9999.f ) std::cout << "Bad Photon seed time !!!! " << std::endl;
	return seedTime;

}//<<>>float GammaResTool::getPhotonSeedTime( pat::Photon )

int GammaResTool::getRhIdx( uInt rhDetID, std::vector<uInt> rhID ){

    //b_rhID->GetEntry(entry);
	int nRhIds = rhID.size();
    for( int idx = 0; idx < nRhIds; idx++ ){ if( rhDetID == rhID[idx] ) return idx; }
    //std::cout << " -- !! no rhDetID to (*rhID)[idx] match !! ---------------------- " << std::endl;
    return -1;

}//<<>>int makehists::getRhIdx( int rhDetID )


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ------------ method called for each event	------------
void GammaResTool::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

	using namespace edm;

// -- Consume Tokens --------------------------------------------
	if( DEBUG ) std::cout << "Consume Tokens -------------------------------------------- " << std::endl;

	// TRACKS
	iEvent.getByToken(tracksToken_, tracks_);

	// VERTICES
	iEvent.getByToken(verticesToken_, vertices_);

	// LEPTONS & PHOTONS
	iEvent.getByToken(electronsToken_, electrons_);

	// PHOTONS
	iEvent.getByToken(gedPhotonsToken_, gedPhotons_);
	iEvent.getByToken(ootPhotonsToken_, ootPhotons_);

	// ECAL RECHITS
	iEvent.getByToken(recHitsEBToken_, recHitsEB_);
	iEvent.getByToken(recHitsEEToken_, recHitsEE_);

	if( doTwoTier ){
    iEvent.getByToken(kuCCStcRecHitsEBToken_, kuCCStcRecHitsEB_);
    iEvent.getByToken(kuCCStcRecHitsEEToken_, kuCCStcRecHitsEE_);
	}//if( doTwoTier )

	// GEOMETRY : https://gitlab.cern.ch/shervin/ECALELF
	caloGeo_ = iSetup.getHandle(caloGeometryToken_); 
	barrelGeometry = caloGeo_->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalBarrel);
	endcapGeometry = caloGeo_->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalEndcap); 

    // Laser constants : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h
    laser_ = iSetup.getHandle(ecalLaserDbServiceToken_);
    evTime = iEvent.time();

	// Intercalibration constants : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h
  	interCalib_ = iSetup.getHandle(ecalIntercalibConstantsToken_);
	interCalibMap = &interCalib_->getMap();

	// ADCToGeV : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h
  	adcToGeV_ = iSetup.getHandle(ecalADCToGeVConstantToken_);
  	adcToGeVEB = adcToGeV_->getEBValue();
  	adcToGeVEE = adcToGeV_->getEEValue();

	// Pedestals : https://github.com/ferriff/usercode/blob/master/DBDump/plugins/DBDump.cc
  	pedestals_ = iSetup.getHandle(EcalPedestalsToken_);


// -- Process Objects ------------------------------------------

// -- event information 

   	run   = iEvent.id().run();
   	lumi  = iEvent.luminosityBlock();
   	event = iEvent.id().event();
    if( DEBUG ) std::cout << "******************************************************************************************************" << std::endl;
	if( DEBUG ) std::cout << "Processing event: " << event << " in run: " << run << " and lumiblock: " << lumi << std::endl;

// -- Process Prime Vertix
	const auto & primevtx = vertices_->front();
	
	auto vtxX = primevtx.position().x();
	auto vtxY = primevtx.position().y();
	auto vtxZ = primevtx.position().z();

    std::vector<EcalRecHit>         frtrechits;
    std::vector<EcalRecHit>         fccrechits;
    std::vector<pat::Photon>        fphotons;
    std::vector<pat::Electron>  	felectrons;

	float minRecHitEnergy = 0.0;	
	if( DEBUG ) std::cout << "Processing RecHits" << std::endl;
	for (const auto recHit : *recHitsEB_ ){ if( recHit.energy() > minRecHitEnergy ) frtrechits.push_back(recHit); }
    for (const auto recHit : *recHitsEE_ ){ if( recHit.energy() > minRecHitEnergy ) frtrechits.push_back(recHit); }
	if( doTwoTier ){
    for (const auto recHit : *kuCCStcRecHitsEB_ ){ if( recHit.energy() > minRecHitEnergy ) fccrechits.push_back(recHit); }
    for (const auto recHit : *kuCCStcRecHitsEE_ ){ if( recHit.energy() > minRecHitEnergy ) fccrechits.push_back(recHit); }
	}//if( doTwoTier )

    if( DEBUG ) std::cout << "Processing " << gedPhotons_->size() << " gedPhotons" << std::endl;
    string phoMvaWp80("mvaPhoID-RunIIFall17-v1-wp80'");
    for( const auto photon : *gedPhotons_ ){

		//auto passIdCut = photon.photonID(phoMvaWp80);
		auto passIdCut = true;
		auto timecut = getPhotonSeedTime(photon) > -25.0;
        if( passIdCut && timecut ){ fphotons.push_back(photon); }
		//if( classcut && timecut ){ fphotons.push_back(photon); phoOOT.push_back(false); }

    }//<<>>for( const auto photon : *gedPhotons_ )
    if( DEBUG ) std::cout << "Selected " << fphotons.size() << " photons." << std::endl;

/*
    if( DEBUG ) std::cout << "Processing ootPhotons" << std::endl;
    for( const auto ootphoton : *ootPhotons_ ){
        if( getPhotonSeedTime(ootphoton) > -25.0 ){ fphotons.push_back(photon); phoOOT.push_back(true); }
    }//<<>>for( const auto photon : *gedPhotons_ )

	const auto nPhotons = fphotons.size();
    vector<bool> phoExcluded(nPhotons+1,false);
    for( int io = 0; io < nPhotons; io++ ){
        double minDr(0.1);
		double lowestDr(10.0);
        int match(-1);
		auto ofPhoton = (*fphotons)[io];
        auto oEta = ofPhoton.eta();
        auto oPhi = ofPhoton.phi();
        for( int ip = io; ip < nPhotons; ip++ ){
			auto pfPhoton = (*fphotons)[ip];
            auto pEta = pfPhoton.eta();
            auto pPhi = pfPhoton.phi();
            auto dRmatch = deltaR( pEta, oEta, pPhi, oPhi );
            if( dRmatch < lowestDr ){ lowestDr = dRmatch; match = ip; }
        }//<<>>for( int ip; ip < nPhotons; ip++ )
        if( lowestDr < minDr ){
            auto pPt = ((*fphotons)[match]).pt();
            auto oPt = ((*fphotons)[io]).pt();
            if( oPt > pPt ) phoExcluded[match] = true;
            else phoExcluded[io] = true;
        }//<<>>if( dRmatch < 0.3 )
    }//<<>>for( int io = 0; io < nOotPhotons; io++ )
*/

    if( DEBUG ) std::cout << "Processing Electrons" << std::endl;
	for( const auto electron : *electrons_ ){
		felectrons.push_back(electron);
	}//<<>>for( const auto electron : *electrons_ )

    //------------------------------------------------------------------------------------
    if( DEBUG ) std::cout << "Processing RecHits" << std::endl;

    rhCaliID.clear();
	rhCaliRtTime.clear();
	rhCaliCCTime.clear();

    std::vector<float> rhEnergy, rhAmp;
    std::vector<float> rhRtTime, rhCCTime, rhTOF;
    std::vector<uInt>  rhID;

    //if( DEBUG ) std::cout << " - enetering RecHit loop" << std::endl;
    for (const auto recHit : frtrechits ){

        //if( DEBUG ) std::cout << " -- proccesing ID info" << std::endl;
        // something in this section is seg faluting after several rechits for crab jobs
        const auto recHitID = getRawID(recHit);
        const auto isEB = getIsEB(recHit);
        //if( DEBUG ) std::cout << " -- proccesing EBEE info" << std::endl;
        //if( DEBUG ) std::cout << " -- proccesing GEO info" << std::endl;
        const auto geometry( isEB ? barrelGeometry : endcapGeometry );
        //if( DEBUG ) std::cout << " -- proccesing POSITION info" << std::endl;
		const auto recHitPos = geometry->getGeometry(recHitID)->getPosition();	
        const auto rhX = recHitPos.x();
        const auto rhY = recHitPos.y();
        const auto rhZ = recHitPos.z();
        //if( DEBUG ) std::cout << " -- proccesing TOF info" << std::endl;
        const auto d_rh = hypo(rhX,rhY,rhZ);
        const auto d_pv = hypo(rhX-vtxX,rhY-vtxY,rhZ-vtxZ);
        const auto tof = (d_rh-d_pv)/SOL;
        //if( DEBUG ) std::cout << " -- proccesing SWISSCROSS info" << std::endl;
        //float swisscross(0.0);
        //if( isEB ) swisscross = EcalTools::swissCross(recHitID, *recHitsEB_, 0.0, true);
        //else swisscross = EcalTools::swissCross(recHitID, *recHitsEE_, 0.0, true);

        //if( DEBUG ) std::cout << " -- proccesing LASER info" << std::endl;
        // adcToGeVInfo : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/src/EcalClusterLazyTools.cc#0204
        const auto laser = laser_->getLaserCorrection(recHitID,evTime);
        const auto interCalibIter = interCalibMap->find(recHitID);
        const auto interCalib = ((interCalibIter != interCalibMap->end()) ? (*interCalibIter) : - 1.f);
        //if( DEBUG ) std::cout << " -- proccesing ADC info" << std::endl;
        //if ((laser > 0.f) && (interCalib > 0.f) && (adcToGeV > 0.f)) rhadcToGeV[pos] = (laser*interCalib*adcToGeV);
        const float adcToGeV( isEB ? adcToGeVEB : adcToGeVEE );
        //if( DEBUG ) std::cout << " -- proccesing PED info" << std::endl;
        // pedestal info
        const auto & pediter = pedestals_->find(recHitID);
		const auto pedrms12 = (pediter != pedestals_->end()) ? pediter->rms(1) : 0.0;
		const auto rhadcToGeV = laser*interCalib*adcToGeV;
		const auto amplitude = ( pedrms12 != 0 && rhadcToGeV != 0 ) ? (recHit.energy()/rhadcToGeV)/pedrms12 : 0;

		auto cctime = -99.0;
		if( doTwoTier ){
		for (const auto ccRecHit : fccrechits ){ if( getRawID(ccRecHit) == recHitID ){ cctime = ccRecHit.time(); break; } }
		}//if( doTwoTier )

        //if( DEBUG ) std::cout << " -- storing values BASE" << std::endl;
        rhID.push_back(recHitID);
        //rhPosX.push_back(rhX);
        //rhPosY.push_back(rhY);
        //rhPosZ.push_back(rhZ);
        rhTOF.push_back(tof);
        //rhPosEta.push_back(recHitPos.eta());
        //rhPosPhi.push_back(recHitPos.phi());
        rhRtTime.push_back(recHit.time());
        rhCCTime.push_back(cctime);
        //if( DEBUG ) std::cout << " -- storing values FLAGS" << std::endl;
        //rhisOOT.push_back(recHit.checkFlag(EcalRecHit::kOutOfTime));
        rhEnergy.push_back(recHit.energy());
		rhAmp.push_back(amplitude);
        //energyError()
        //rhSwCross.push_back(swisscross);
        //rhisWeird.push_back(recHit.checkFlag(EcalRecHit::kWeird));
        //rhisDiWeird.push_back(recHit.checkFlag(EcalRecHit::kDiWeird));
        //rhisGS6.push_back(recHit.checkFlag(EcalRecHit::kHasSwitchToGain6));
        //rhisGS1.push_back(recHit.checkFlag(EcalRecHit::kHasSwitchToGain1));
        //if ((laser > 0.f) && (interCalib > 0.f) && (adcToGeV > 0.f)) 
        //if( DEBUG ) std::cout << " -- storing values PED" << std::endl;
        //rhadcToGeV.push_back(laser*interCalib*adcToGeV);
        //if( DEBUG ) std::cout << " -- next rechit" << std::endl;

		if( recHit.energy() > 5.0 ){
			rhCaliID.push_back(recHitID);
			rhCaliRtTime.push_back(recHit.time());
			rhCaliCCTime.push_back(cctime);
		}//<<>>if( recHit.energy() > 5.0 )

    }//<<>>for (const auto recHit : *recHitsEB_ )   

    if( DEBUG ) std::cout << " - enetering Photon loop" << std::endl;

    std::vector<uInt> locRHCands;
	std::vector<pat::Photon> gloPhotons;
    std::vector<uInt> locSeedRHs{0,0};
    std::vector<uInt> gloSeedRHs{0,0};
	vector<vector<int>> offsets{{1,0},{-1,0},{0,1},{0,-1},{1,1},{1,-1},{-1,1},{-1,-1}};
    for (const auto photon : fphotons ){

		//auto passIdCut = photon.photonID(phomMvaWp80);
		//if( not passIdCut ) continue;
		const auto &phosc = photon.superCluster().isNonnull() ? photon.superCluster() : photon.parentSuperCluster();
		const auto scptr = phosc.get();
		scGroup phoSCGroup{*scptr};
		const auto seedDetId = scptr->seed()->seed(); // seed detid
		const auto isEB = (seedDetId.subdetId() == EcalSubdetector::EcalBarrel); // which subdet
		const auto phoRecHits = (isEB ? *recHitsEB_ : *recHitsEE_ );

		// select global photons
		if ( photon.hasPixelSeed() ){ 
			float bestmatch(0.1);
			float elTrackZ(1000.0);
    		for (const auto electron : felectrons){
				auto match = reco::deltaR(electron,photon);
          		if ( match < bestmatch ){ elTrackZ = electron.trackPositionAtVtx().Z(); bestmatch = match; }
    		}//<<>>for (const auto electron : felectrons){
			auto eleMatch = elTrackZ < 1000.0;
			auto trackMatch = abs( elTrackZ - vtxZ ) < 1.0;
			if( eleMatch && trackMatch ) gloPhotons.push_back(photon); 
		}//<<>>if (not inpho.hasPixSeed)

		// select local photons rechit+neighbor crystals 

      	const auto & ph2ndMoments = noZS::EcalClusterTools::cluster2ndMoments( *scptr, phoRecHits);
    	const auto smaj  = ph2ndMoments.sMaj;
    	const auto smin  = ph2ndMoments.sMin;
        if ( smin < 0.3 && smaj < 0.5){
			auto phoRhGroup = getRHGroup( phoSCGroup, 1.0 );
			if( DEBUG ) std::cout << " Examining Photon with " << phoRhGroup.size()	<< " rechits." << std::endl;
			for( auto rechit : phoRhGroup ){
				const auto rhDetId = rechit.detid();
				const auto rhEnergy = recHitE( rhDetId, phoRecHits );
				for( auto offset : offsets ){ 
					const auto nbDetId = ( isEB ) ? EBDetId::offsetBy( rhDetId, offset[0], offset[1] ) : EEDetId::offsetBy( rhDetId, offset[0], offset[1] );
					auto neighborEnergy = recHitE( nbDetId, phoRecHits );
					auto ordered = rhEnergy > neighborEnergy;
					auto close = rhEnergy < 1.20*neighborEnergy;
					//if( DEBUG ) std::cout << " Examining rechit pair with " << rhEnergy << " & " << neighborEnergy << " energies/" << std::endl;
					if( ordered && close ){  // need to be within 20% of energy
						if( DEBUG ) std::cout << " Matching loc rechit pair with " << rhEnergy << " & " << neighborEnergy << " energies/" << std::endl;
						locRHCands.push_back( rhDetId.rawId() ); 
						locRHCands.push_back( nbDetId.rawId() );
					}//<<>>if( high < 1.20*low )
				}//<<>>for( auto offset : offsets )
			}//<<>>for( rechit : phoRhGroup )
		}//<<>>if ( smin < 0.3 && smaj < 0.5)

	}//<<>>for (const auto photon : fphotons )

	if( DEBUG ) std::cout << " Selected " << gloPhotons.size() << " global photons and " << locRHCands.size() << " local rechits. " << std::endl;

	// select rhs for global
	int nGloPhos = gloPhotons.size();
	if( nGloPhos > 1 ){
		float zMassMatch(35.00);
		float zMass(91.1876);
		vector<int> phoIndx{0,0};
		for( int first(0); first < nGloPhos; first++ ){
            auto pho1Eta = gloPhotons[first].eta();
            auto pho1Phi = gloPhotons[first].phi();
            auto pho1Pt = gloPhotons[first].pt();
            auto pho1E = gloPhotons[first].energy();
			TLorentzVector pho1vec;
			pho1vec.SetPtEtaPhiE(pho1Pt, pho1Eta, pho1Phi, pho1E);
			for( int second(first+1); second < nGloPhos; second++ ){
				auto pho2Eta = gloPhotons[second].eta();
				auto pho2Phi = gloPhotons[second].phi();
                auto pho2Pt = gloPhotons[second].pt();
                auto pho2E = gloPhotons[second].energy();
				TLorentzVector pho2vec; 
				pho2vec.SetPtEtaPhiE(pho2Pt, pho2Eta, pho2Phi, pho2E);
				//const auto dr12 = pho1vec.DeltaR(pho2vec);
				pho1vec += pho2vec;
				auto pairMass = pho1vec.M();
				if( pairMass > 60.0 && pairMass < 120.0 ){ 
					auto zMassDiff = abs(pairMass-zMass); 
					if( zMassDiff < zMassMatch ){ phoIndx[0] = first; phoIndx[1] = second; zMassMatch = zMassDiff; }
				}//<<>>if( pairMass > 60.0 && pairMass < 120.0 )
			}//<<>>for( int second(first+1); second < nGloPhos; second++ )
		}//<<>>for( int first(0); first < nGloPhos; first++ )
		if( zMassMatch < 35.00 ){
			auto pho0 = gloPhotons[phoIndx[0]];
			auto pho1 = gloPhotons[phoIndx[1]];
        	const auto &phosc0 = pho0.superCluster().isNonnull() ? pho0.superCluster() : pho0.parentSuperCluster();
            const auto &phosc1 = pho1.superCluster().isNonnull() ? pho1.superCluster() : pho1.parentSuperCluster();
			gloSeedRHs[0] = ((phosc0.get())->seed()->seed()).rawId();
			gloSeedRHs[1] = ((phosc1.get())->seed()->seed()).rawId();
			if( DEBUG ) std::cout << " Matching glo photon pair with : " << gloSeedRHs[0] << " & " << gloSeedRHs[1] << " with dZmass : " << zMassMatch << std::endl;
		}//<<>>if( zMassMatch < 35.00 )
	}//<<>>if( gloPhotons.size() > 1 )i

	int nLocRHCands = locRHCands.size();
	if( nLocRHCands > 1 ){
		int lead(0);
		for( int next(2); next+1 < nLocRHCands; next+=2 ){
			auto leadE = rhEnergy[getRhIdx(locRHCands[lead],rhID)];
			auto nextRhE = rhEnergy[getRhIdx(locRHCands[next],rhID)];
			if( nextRhE > leadE ) lead = next;
		}//<<>>for( int it(0); it+1 < nLocRHCands; it += 2; )
		locSeedRHs[0] = locRHCands[lead];
		locSeedRHs[1] = locRHCands[lead+1];
	}//<<>>if( nLocRHCands > 1 )

	if( DEBUG ) std::cout << " Storing ids Global : " << gloSeedRHs[0] << ", " << gloSeedRHs[1] << " and Local: " << locSeedRHs[0] << ", " << locSeedRHs[1] << std::endl;

	resRhID.clear();
	resAmp.clear();
	resE.clear();
	resRtTime.clear();
	resCCTime.clear();
	resTOF.clear();

	if( locSeedRHs[0] ) locMatches++;
    if( gloSeedRHs[0] ) gloMatches++;
	nEvents++;
	std::vector<uInt> lgRhIds{locSeedRHs[0],locSeedRHs[1],gloSeedRHs[0],gloSeedRHs[1]};
	for( auto lgRhId : lgRhIds ){
			
		auto idx = ( lgRhId > 0 ) ? getRhIdx(lgRhId,rhID) : 0; 
		resRhID.push_back( ( lgRhId > 0 ) ? lgRhId : 0 );
		resAmp.push_back( ( lgRhId > 0 ) ? rhAmp[idx] : 0 );
        resE.push_back( ( lgRhId > 0 ) ? rhEnergy[idx] : 0 );
        resRtTime.push_back( ( lgRhId > 0 ) ? rhRtTime[idx] : 0 );
        resCCTime.push_back( ( lgRhId > 0 ) ? rhCCTime[idx] : 0 );
		resTOF.push_back( ( lgRhId > 0 ) ? rhTOF[idx] : 0 );

	}//<<>>for( auto lgRhId : lgRhIds )
    if( DEBUG ) std::cout << " Stored times Global : " << resRtTime[2] << ", " << resRtTime[3] << " and Local: " << resRtTime[0] << ", " << resRtTime[1] << std::endl;

	// -- Fill output trees ------------------------------------------
	if( DEBUG ) std::cout << "---------- Next Event -----" << std::endl;
	outTree->Fill();

	// -- EOFun ------------------------------------------------------
	//#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	//	 ESHandle<SetupData> pSetup;
	//	 iSetup.get<SetupRecord>().get(pSetup);
	//#endif
}//>>>>void GammaResTool::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)


// ------------ method called once each job just before starting event loop	------------
void GammaResTool::beginJob(){

    locMatches = 0;
    gloMatches = 0;
	nEvents = 0;

	// Book output files and trees
	edm::Service<TFileService> fs;
	outTree = fs->make<TTree>("llpgtree","llpgtree");

	// Book histograms
	
    //------ 1D Hists --------------------------------------------------------------------------

	//------ 2D Hists --------------------------------------------------------------------------

	std::cout << "Histograms Booked" << std::endl;

	// Create output Tree branches -----------------------------

	// Run, Lumi, Event info
	outTree->Branch("run", &run);
	outTree->Branch("lumi", &lumi);
	//outTree->Branch("event", &event, "event/l");

    outTree->Branch("rhCaliID", &rhCaliID);
    outTree->Branch("rhCaliRtTime", &rhCaliRtTime);
    outTree->Branch("rhCaliCCTime", &rhCaliCCTime);

    outTree->Branch("resRhID", &resRhID);
    outTree->Branch("resAmp", &resAmp);
    outTree->Branch("resE", &resE);
    outTree->Branch("resRtTime", &resRtTime);
    outTree->Branch("resCCTime", &resCCTime);
    outTree->Branch("resTOF", &resTOF);	

}//>>>>void GammaResTool::beginJob()


// ------------ method called once each job just after ending the event loop	------------
void GammaResTool::endJob(){

	if( nEvents == 0 ) nEvents = 1;
	//if( DEBUG ) 
	std::cout << " Found Global % " << gloMatches/nEvents << " and Local % " << locMatches/nEvents << " of " << nEvents << std::endl;

}//>>>>void GammaResTool::endJob()


// ------------ method fills 'descriptions' with the allowed parameters for the module	------------
void GammaResTool::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);

	//Specify that only 'tracks' is allowed
	//To use, remove the default given above and uncomment below
	//ParameterSetDescription desc;
	//desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
	//descriptions.addDefault(desc);
}//>>>>void GammaResTool::fillDescriptions(edm::ConfigurationDescriptions& descriptions)


//define this as a plug-in
DEFINE_FWK_MODULE(GammaResTool);
