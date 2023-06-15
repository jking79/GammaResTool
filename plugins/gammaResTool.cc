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

//#define DIAG true   replaced with tag below
//#define DIAG false

//
// constructors and destructor
//
GammaResTool::GammaResTool(const edm::ParameterSet& iConfig) :

// -- declare tags ----------------------------------------------------------

	// flags
	doTwoTier(iConfig.existsAs<bool>("doTwoTier")  ? iConfig.getParameter<bool>("doTwoTier")  : false),
    doDiag(iConfig.existsAs<bool>("doDiag")  ? iConfig.getParameter<bool>("doDiag")  : false),

	// tracks
	tracksTag(iConfig.getParameter<edm::InputTag>("tracks")),

    // vertices
    verticesTag(iConfig.getParameter<edm::InputTag>("vertices")),
	
	// electrons
	electronsTag(iConfig.getParameter<edm::InputTag>("electrons")),  

	// recHits
	recHitsEBTag(iConfig.getParameter<edm::InputTag>("recHitsEB")),  
	recHitsEETag(iConfig.getParameter<edm::InputTag>("recHitsEE")),

    kuRtStcRecHitsEBTag(iConfig.getParameter<edm::InputTag>("kuRtStcRecHitsEB")),
    kuRtStcRecHitsEETag(iConfig.getParameter<edm::InputTag>("kuRtStcRecHitsEE")),

    kuCCStcRecHitsEBTag(iConfig.getParameter<edm::InputTag>("kuCCStcRecHitsEB")),
    kuCCStcRecHitsEETag(iConfig.getParameter<edm::InputTag>("kuCCStcRecHitsEE")),

    unCCRecHitsEBTag(iConfig.getParameter<edm::InputTag>("ku_uncalibratedRecHitsEB")),
    unCCRecHitsEETag(iConfig.getParameter<edm::InputTag>("ku_uncalibratedRecHitsEE")),

	// gedphotons
	gedPhotonsTag(iConfig.getParameter<edm::InputTag>("gedPhotons")),

	// ootPhotons
	ootPhotonsTag(iConfig.getParameter<edm::InputTag>("ootPhotons")),

	// ECAL RECORDS 
    caloGeometryToken_(esConsumes()),
	caloTopologyToken_{esConsumes()},
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

    if( doDiag ){
    	kuRtStcRecHitsEBToken_ = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(kuRtStcRecHitsEBTag);
    	kuRtStcRecHitsEEToken_ = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(kuRtStcRecHitsEETag);
    }//if( doTwoTier )

	if( doTwoTier ){
    	kuCCStcRecHitsEBToken_ = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(kuCCStcRecHitsEBTag);
    	kuCCStcRecHitsEEToken_ = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(kuCCStcRecHitsEETag);
		unCCRecHitsEBToken_ = 
			consumes<edm::SortedCollection<EcalUncalibratedRecHit,edm::StrictWeakOrdering<EcalUncalibratedRecHit>>>(unCCRecHitsEBTag);
    	unCCRecHitsEEToken_ = 
			consumes<edm::SortedCollection<EcalUncalibratedRecHit,edm::StrictWeakOrdering<EcalUncalibratedRecHit>>>(unCCRecHitsEETag);
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
	//if( DEBUG && seedTime == -9999.f ) std::cout << "Bad Photon seed time !!!! " << std::endl;
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

	bool storeEvent(true);

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

    if( doDiag ){
    	iEvent.getByToken(kuRtStcRecHitsEBToken_, kuRtStcRecHitsEB_);
    	iEvent.getByToken(kuRtStcRecHitsEEToken_, kuRtStcRecHitsEE_);
    }//if( doTwoTier )

	if( doTwoTier ){
    	iEvent.getByToken(kuCCStcRecHitsEBToken_, kuCCStcRecHitsEB_);
    	iEvent.getByToken(kuCCStcRecHitsEEToken_, kuCCStcRecHitsEE_);
    	iEvent.getByToken(unCCRecHitsEBToken_, unCCRecHitsEB_);
    	iEvent.getByToken(unCCRecHitsEEToken_, unCCRecHitsEE_);
	}//if( doTwoTier )

	// GEOMETRY : https://gitlab.cern.ch/shervin/ECALELF
	caloGeo_ = iSetup.getHandle(caloGeometryToken_); 
	barrelGeometry = caloGeo_->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalBarrel);
	endcapGeometry = caloGeo_->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalEndcap); 

	// CaloTopology
	caloTopo_ = iSetup.getHandle(caloTopologyToken_);
	topology = &iSetup.getData(caloTopologyToken_);

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

// -- Process ECAL RecHits
    std::vector<EcalRecHit>         frechits;
    std::vector<EcalRecHit>         frtrechits;
    std::vector<EcalRecHit>         fccrechits;
    std::vector<EcalUncalibratedRecHit> funccrechits;
    std::vector<pat::Photon>        fphotons;
    std::vector<pat::Electron>  	felectrons;

	float minRecHitEnergy = 0.3;
    float minRecHitAmp = 5;	
	if( DEBUG ) std::cout << "Processing RecHits" << std::endl;
	for (const auto recHit : *recHitsEB_ ){ if( recHit.energy() > minRecHitEnergy ) frechits.push_back(recHit); }
    for (const auto recHit : *recHitsEE_ ){ if( recHit.energy() > minRecHitEnergy ) frechits.push_back(recHit); }
    if( doDiag ){
        for (const auto recHit : *kuRtStcRecHitsEB_ ){ if( recHit.energy() > minRecHitEnergy ) frtrechits.push_back(recHit); }
        for (const auto recHit : *kuRtStcRecHitsEE_ ){ if( recHit.energy() > minRecHitEnergy ) frtrechits.push_back(recHit); }
    }//<<>>if( doDiag )
	if( doTwoTier ){
    	for (const auto recHit : *kuCCStcRecHitsEB_ ){ if( recHit.energy() > minRecHitEnergy ) fccrechits.push_back(recHit); }
    	for (const auto recHit : *kuCCStcRecHitsEE_ ){ if( recHit.energy() > minRecHitEnergy ) fccrechits.push_back(recHit); }
        for (const auto recHit : *unCCRecHitsEB_ ){ if( recHit.amplitude() > minRecHitAmp ) funccrechits.push_back(recHit); }
        for (const auto recHit : *unCCRecHitsEE_ ){ if( recHit.amplitude() > minRecHitAmp ) funccrechits.push_back(recHit); }
	}//<<>>if( doTwoTier )

// -- Process gedPhotons
    if( DEBUG ) std::cout << "Processing " << gedPhotons_->size() << " gedPhotons" << std::endl;
	string phoMvaWp80("mvaPhoID-RunIIFall17-v1-wp80");//2018
    //string phoMvaWp80("mvaPhoID-RunIIFall17-v2-wp80");//r3
    string phoCutTight("cutBasedPhotonID-Fall17-94X-V2-tight");//2022
    string phoCutLoose("cutBasedPhotonID-Fall17-94X-V2-loose");//2022
    //string phoCutLoose("cutBasedPhotonID-Fall17-94X-V1-loose");//2018
	float phoMinPt = 5.0;
	float phoMinSeedTime = -25.0;
    for( const auto photon : *gedPhotons_ ){

		auto passIdCut = photon.photonID(phoCutLoose);
        //auto passIdCut = true;
		auto minPhoPt = photon.pt() > phoMinPt;
		auto phoSeedTime = getPhotonSeedTime(photon);
		auto timecut = phoSeedTime > phoMinSeedTime;
		
		if( DEBUG ) std::cout << " - pho time cut : " << phoSeedTime << " passID: " << passIdCut << " pt: " << photon.pt() << std::endl;
        if( passIdCut && timecut && minPhoPt ){ fphotons.push_back(photon); if( DEBUG ) std::cout << " - >> Photon Passed " << std::endl; }
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


// -- Process electrons

    //std::vector<pat::Electron> gloPhotons;
    if( DEBUG ) std::cout << "Processing Electrons" << std::endl;
    string eleMvaWp80("mvaEleID-Fall17-noIso-V2-wp80");
    string eleMvaWpLoose("mvaEleID-Fall17-noIso-V2-wpLoose");//r3
    string eleCutTight("cutBasedElectronID-Fall17-94X-V2-tight");
    string eleCutLoose("cutBasedElectronID-Fall17-94X-V2-loose");//2022
    //string eleCutLoose("cutBasedElectronID-Fall17-94X-V1-loose");//2018
	for( const auto electron : *electrons_ ){

        auto passIdCut = electron.electronID(eleCutLoose);
        //auto passIdCut = true;
		if( passIdCut ) felectrons.push_back(electron);

	}//<<>>for( const auto electron : *electrons_ )

    //------------------------------------------------------------------------------------
    if( DEBUG ) std::cout << "Processing Mini & Cali RecHits" << std::endl;

    rhCaliID.clear();
	rhCaliEnergy.clear();
	rhCaliRtTime.clear();
	rhCaliCCTime.clear();

	std::vector<uInt> miniRhId;
    std::vector<float> miniRhAmp;
    std::vector<float> miniRhEnergy;
    std::vector<float> miniRhRtTime;
    std::vector<float> miniRhCCTime;
    std::vector<float> miniRhTOF;

    //if( DEBUG ) std::cout << " - enetering RecHit loop" << std::endl;
    for (const auto recHit : frechits ){

		const auto recHitID = getRawID(recHit);
        auto cctime = -999.0;
        if( doTwoTier ){
            for (const auto ccRecHit : fccrechits ){ if( getRawID(ccRecHit) == recHitID ){ cctime = ccRecHit.time(); break; } }
        }//if( doTwoTier )

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
        //if( DEBUG ) std::cout << " -- proccesing LASER info" << std::endl;
        // adcToGeVInfo : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/src/EcalClusterLazyTools.cc#0204
        const auto laser = laser_->getLaserCorrection(recHitID,evTime);
        const auto interCalibIter = interCalibMap->find(recHitID);
        const auto interCalib = ((interCalibIter != interCalibMap->end()) ? (*interCalibIter) : - 1.f);
        //if( DEBUG ) std::cout << " -- proccesing ADC info" << std::endl;
        //if ((laser > 0.f) && (interCalib > 0.f) && (adcToGeV > 0.f)) rhadcToGeV[pos] = (laser*interCalib*adcToGeV);
        const float adcToGeV0( isEB ? adcToGeVEB : adcToGeVEE );
        //if( DEBUG ) std::cout << " -- proccesing PED info" << std::endl;
        // pedestal info
        const auto & pediter = pedestals_->find(recHitID);
        const auto pedrms12 = (pediter != pedestals_->end()) ? pediter->rms(1) : 0.0;
        const auto adcToGeV = laser*interCalib*adcToGeV0;
        const auto amplitude = ( pedrms12 != 0 && adcToGeV != 0 ) ? (recHit.energy()/adcToGeV)/pedrms12 : 0;

		miniRhId.push_back(recHitID);
		miniRhAmp.push_back(amplitude);
		miniRhEnergy.push_back(recHit.energy());
		miniRhRtTime.push_back(recHit.time());
		miniRhCCTime.push_back(cctime);
		miniRhTOF.push_back(tof);

        float rhECut(5.0);
        if( recHit.energy() > rhECut ){

            rhCaliID.push_back(recHitID);
			rhCaliEnergy.push_back(recHit.energy());
            rhCaliRtTime.push_back(recHit.time());
            rhCaliCCTime.push_back(cctime);

        }//<<>>if( recHit.energy() > 5.0 )

    }//<<>>for (const auto recHit : *recHitsEB_ ) 

	if( DEBUG ) std::cout << "Processing Diag RecHits" << std::endl;

	rhID.clear();
	rhRtTime.clear();
	rhCCTime.clear();
	//rhTimeErr, 
	rhTOF.clear();
    rhEnergy.clear();
	rhAmp.clear();
	rhRtisOOT.clear();
    rhCCisOOT.clear();
	rhisGS6.clear(); 
	rhisGS1.clear();
	rhisWeird.clear();
	rhisDiWeird.clear();
	rhadcToGeV.clear();
	rhSwCross.clear();
	rhpedrms12.clear();
    //rhped12, rhped6, rhped1;
    //rhpedrms12, rhpedrms6, rhpedrms1;

    //if( DEBUG ) std::cout << " - enetering Diag RecHit loop" << std::endl;
    if( doDiag ){
	    for (const auto recHit : frtrechits ){

		if( recHit.energy() > 2.0 ){
	
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
	        float swisscross(10.0);
	        if( isEB ) swisscross = EcalTools::swissCross(recHitID, *kuRtStcRecHitsEB_, 1.0, true);
	        else swisscross = EcalTools::swissCross(recHitID, *kuRtStcRecHitsEE_, 1.0, true);
	
	        //if( DEBUG ) std::cout << " -- proccesing LASER info" << std::endl;
	        // adcToGeVInfo : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/src/EcalClusterLazyTools.cc#0204
	        const auto laser = laser_->getLaserCorrection(recHitID,evTime);
	        const auto interCalibIter = interCalibMap->find(recHitID);
	        const auto interCalib = ((interCalibIter != interCalibMap->end()) ? (*interCalibIter) : - 1.f);
	        //if( DEBUG ) std::cout << " -- proccesing ADC info" << std::endl;
	        //if ((laser > 0.f) && (interCalib > 0.f) && (adcToGeV > 0.f)) rhadcToGeV[pos] = (laser*interCalib*adcToGeV);
	        const float adcToGeV0( isEB ? adcToGeVEB : adcToGeVEE );
	        //if( DEBUG ) std::cout << " -- proccesing PED info" << std::endl;
	        // pedestal info
	        const auto & pediter = pedestals_->find(recHitID);
			const auto pedrms12 = (pediter != pedestals_->end()) ? pediter->rms(1) : 0.0;
			const auto adcToGeV = laser*interCalib*adcToGeV0;
			const auto amplitude = ( pedrms12 != 0 && adcToGeV != 0 ) ? (recHit.energy()/adcToGeV)/pedrms12 : 0;

			auto cctime = -999.0;
			auto ccisoot = false;
			if( doTwoTier ){
				for (const auto ccRecHit : fccrechits ){ 
					if( getRawID(ccRecHit) == recHitID ){ 
						cctime = ccRecHit.time(); 
						ccisoot = ccRecHit.checkFlag(EcalRecHit::kOutOfTime);
						break; 
					}//<<>>for (const auto ccRecHit : fccrechits )
				}//<<>>for (const auto ccRecHit : fccrechits )
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
	        rhEnergy.push_back(recHit.energy());
			rhAmp.push_back(amplitude);
	        //energyError()
	        rhRtisOOT.push_back(recHit.checkFlag(EcalRecHit::kOutOfTime));
            rhCCisOOT.push_back(ccisoot);
			//rhisGood.push_back(recHit.checkFlag(EcalRecHit::kGood));
	        rhSwCross.push_back(swisscross);
	        rhisWeird.push_back(recHit.checkFlag(EcalRecHit::kWeird));
	        rhisDiWeird.push_back(recHit.checkFlag(EcalRecHit::kDiWeird));
	        rhisGS6.push_back(recHit.checkFlag(EcalRecHit::kHasSwitchToGain6));
	        rhisGS1.push_back(recHit.checkFlag(EcalRecHit::kHasSwitchToGain1));
	        //if ((laser > 0.f) && (interCalib > 0.f) && (adcToGeV > 0.f)) 
	        //if( DEBUG ) std::cout << " -- storing values PED" << std::endl;
	        rhadcToGeV.push_back(adcToGeV);
			rhpedrms12.push_back(pedrms12);
	        //if( DEBUG ) std::cout << " -- next rechit" << std::endl;
	
		}//<<>> if( recHit.energy() > 2.0 )

	    }//<<>>for (const auto recHit : *recHitsEB_ )   
    }//<<>>if( doDiag )

	if( DEBUG ) std::cout << "Processing Uncalibrated RecHits" << std::endl;

    unrhJitter.clear();
    unrhNonJitter.clear();
    unrhEncNonJitter.clear();
    unrhEnergy.clear();

    if( doTwoTier ){
        for (const auto recHit : funccrechits ){

        	const auto recHitID = getRawID(recHit);
        	const auto isEB = getIsEB(recHit);

			if( isEB ){

				float energy = -9.0;
				float rhNonTime = -99.0;
        		for (const auto ccRecHit : fccrechits ){ 
					if( getRawID(ccRecHit) == recHitID ){ 
						energy = ccRecHit.energy();
						rhNonTime = ccRecHit.nonCorrectedTime(); 
						break; 
					}//<<>>if( getRawID(ccRecHit) == recHitID ) 
				}//<<>>for (const auto ccRecHit : fccrechits )

				auto jitter = recHit.jitter();
				auto nonjitter = recHit.chi2();
				auto encNonTime = recHit.nonCorrectedTime();

				std::cout << "NonCorrectedTime : rh = " << rhNonTime << " unrh = " << encNonTime; 
				std::cout << " diff = " << rhNonTime - encNonTime << std::endl;

				//std::cout << " Storing cc encoding info " << std::endl;
				unrhJitter.push_back(jitter);
				unrhNonJitter.push_back(nonjitter);
				unrhEncNonJitter.push_back(encNonTime);
				unrhEnergy.push_back(energy);

/*
				auto amp = recHit.amplitude();
				auto time = 25*recHit.jitter();
				auto untime = 25*recHit.chi2();
				auto enctime = 25*recHit.nonCorrectedTime();
				hist2d[0]->Fill(amp,energy);
				auto offset = 0.0;
				auto slope = 1.3;
				auto actdif = time - untime;
				auto adjtime = slope*time + offset;
				//auto adjtime = (1.2 + 0.4*(1.0-std::exp(-1*amp/20)))*time + offset; 
				//auto fac = (std::abs(time) - 12.5)*(std::abs(time) - 12.5) - 156.25;
				//auto adjtime = (1.0 + 0.8*(1.0-std::exp(-1*fac*fac)))*time + offset;
				auto encdif = adjtime - untime;
				auto reso = enctime - untime;

				hist2d[7]->Fill( encdif, energy );
                totrhs0++;//hist1d[0]->Fill(16);
                if( energy > 0.5 ) totrhs05++;//hist1d[1]->Fill(16); //if( inwindow ) hist1d[0]->Fill(i); }
                if( energy > 1.0 ) totrhs1++; //hist1d[2]->Fill(16); //if( inwindow ) hist1d[1]->Fill(i); }
                if( energy > 2.0 ) totrhs2++;//hist1d[3]->Fill(16); //if( inwindow ) hist1d[2]->Fill(i); }
                if( energy > 5.0 ) totrhs5++;//hist1d[4]->Fill(16); //if( inwindow ) hist1d[3]->Fill(i); }
                if( energy > 10.0 ) totrhs10++;//hist1d[5]->Fill(16); //if( inwindow ) hist1d[4]->Fill(i); }
				for( int i = 0; i < 15; i++ ){
					bool inwindow = std::abs(adjtime) < i; 
					if( inwindow ) hist1d[0]->Fill(i);
					if( energy > 0.5 ){ if( inwindow ) hist1d[1]->Fill(i); }
					if( energy > 1.0 ){ if( inwindow ) hist1d[2]->Fill(i); }
					if( energy > 2.0 ){ if( inwindow ) hist1d[3]->Fill(i); }
                    if( energy > 5.0 ){ if( inwindow ) hist1d[4]->Fill(i); }
                    if( energy > 10.0 ){ if( inwindow ) hist1d[5]->Fill(i); }
				}//<<>>for( int i = 0; i < 13; i++ )

                float ethres = 10.0;
                float maxrange = 8.0;
                float minrange = -8.0;
                if( energy > ethres ){

					totrhs++;
					if( adjtime < maxrange && adjtime > minrange ) encrhs++;
					hist2d[2]->Fill( untime, actdif );
                    hist2d[3]->Fill( untime, encdif );
                    hist2d[4]->Fill( untime, reso );
					hist2d[5]->Fill( encdif, amp );
                    hist2d[6]->Fill( actdif, amp );

				}//<<>>if( energy > ethres )
*/

			}//<<>>if( isEB )

        }//<<>>for (const auto recHit : *funccrechits )   
    }//<<>>if( doTwoTier )


    if( DEBUG ) std::cout << " - enetering Photon loop" << std::endl;

    std::vector<uInt> locRHCands;
	std::vector<pat::Photon> gloPhotons;
    std::vector<pat::Photon> locPhotons;
    std::vector<pat::Photon> selPhotons;
    std::vector<int> selPhoType;
    std::vector<uInt> locSeedRHs{0,0};
    std::vector<uInt> gloSeedRHs{0,0};
	float gloDiMass(-1), gloDiAngle(-1), gloDiDr(-1), gloDiPhi(-1), gloDiEta(-1);
	//std::vector<int> selPhoIndx{-1,-1,-1};
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
    		}//<<>>for (const auto electron : felectrons)
			auto eleMatch = elTrackZ < 1000.0;
			auto dz = abs( elTrackZ - vtxZ );
			auto trackMatch = dz < 1.0;
			if( eleMatch && trackMatch ){ 
				gloPhotons.push_back(photon); 
				if( DEBUG ) std::cout << " Electron track match for Glo photon with " << dz << " dz" << std::endl;
			}//<<<>>if( eleMatch && trackMatch )
		}//<<>>if (not inpho.hasPixSeed)

		// select local photons rechit+neighbor crystals 

      	const auto & ph2ndMoments = noZS::EcalClusterTools::cluster2ndMoments( *scptr, phoRecHits);
    	const auto smaj  = ph2ndMoments.sMaj;
    	const auto smin  = ph2ndMoments.sMin;
        if ( smin < 0.3 && smaj < 0.5){
			auto phoRhGroup = getRHGroup( phoSCGroup, 0.0 );
			if( DEBUG ) std::cout << " Examining Photon with " << phoRhGroup.size()	<< " rechits." << std::endl;
			for( auto rechit : phoRhGroup ){
				const auto rhDetId = rechit.detid();
				const auto lrhEnergy = recHitE( rhDetId, phoRecHits );
				for( auto offset : offsets ){ 
					const auto nbDetId = ( isEB ) ? EBDetId::offsetBy( rhDetId, offset[0], offset[1] ) : EEDetId::offsetBy( rhDetId, offset[0], offset[1] );
					auto neighborEnergy = recHitE( nbDetId, phoRecHits );
					auto ordered = lrhEnergy > neighborEnergy;
					auto close = lrhEnergy < 1.20*neighborEnergy;
					//if( DEBUG ) std::cout << " Examining rechit pair with " << lrhEnergy << " & " << neighborEnergy << " energies" << std::endl;
					if( ordered && close ){  // need to be within 20% of energy
						if( DEBUG ) std::cout << " Matching loc rechit pair with e: " << lrhEnergy << " : " << neighborEnergy;
                        if( DEBUG ) std::cout << " (" << lrhEnergy/neighborEnergy << ")"<< " for lead rh id: " << rhDetId.rawId() <<std::endl;
						locPhotons.push_back( photon );
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
				TVector3 pho2vec3(pho2Pt, pho2Eta, pho2Phi); 
				pho2vec.SetPtEtaPhiE(pho2Pt, pho2Eta, pho2Phi, pho2E);
				const auto dr12 = pho1vec.DeltaR(pho2vec);
				const auto ang12 = pho1vec.Angle(pho2vec3);
				const auto dphi12 = pho1vec.DeltaPhi(pho2vec);
                const auto deta12 = std::abs(pho1Eta - pho2Eta);
				pho1vec += pho2vec;
				auto pairMass = pho1vec.M();
				if( pairMass > 60.0 && pairMass < 120.0 ){ 
					auto zMassDiff = std::abs(pairMass-zMass); 
					if( zMassDiff < zMassMatch ){ 
						phoIndx[0] = first; phoIndx[1] = second; zMassMatch = zMassDiff; 
						gloDiMass = pairMass; gloDiAngle = ang12; gloDiDr = dr12; gloDiPhi = dphi12; gloDiEta = deta12;
						if( DEBUG ) std::cout << " Matching glo pho pair with " << zMassDiff << " mass diff" << std::endl;	
					}//<<>>if( zMassDiff < zMassMatch )}
				}//<<>>if( pairMass > 60.0 && pairMass < 120.0 )
			}//<<>>for( int second(first+1); second < nGloPhos; second++ )
		}//<<>>for( int first(0); first < nGloPhos; first++ )
		if( zMassMatch < 35.00 ){
			
			//selPhoIndx[1] = phoIndx[0];
			auto pho0 = gloPhotons[phoIndx[0]];
			selPhotons.push_back(pho0);
			selPhoType.push_back(1);
			//selPhoIndx[2] = phoIndx[1];
			auto pho1 = gloPhotons[phoIndx[1]];
            selPhotons.push_back(pho1);
            selPhoType.push_back(2);
        	const auto &phosc0 = pho0.superCluster().isNonnull() ? pho0.superCluster() : pho0.parentSuperCluster();
            const auto &phosc1 = pho1.superCluster().isNonnull() ? pho1.superCluster() : pho1.parentSuperCluster();
			gloSeedRHs[0] = ((phosc0.get())->seed()->seed()).rawId();
			gloSeedRHs[1] = ((phosc1.get())->seed()->seed()).rawId();
			if( DEBUG ) std::cout << " Selecting matching glo photon pair with : " << gloSeedRHs[0] << " & " << gloSeedRHs[1];
            if( DEBUG ) std::cout << " and dZmass : " << zMassMatch << std::endl;
			

		}//<<>>if( zMassMatch < 35.00 )
	}//<<>>if( gloPhotons.size() > 1 )i

	if( DEBUG ) std::cout << "Sorting Local RecHits " << std::endl;
	int nLocRHCands = locRHCands.size();
	if( nLocRHCands > 0 ){
		if( nLocRHCands > 2 ){
			int lead(0);
			for( int next(2); next < nLocRHCands; next+=2 ){
				auto leadE = miniRhEnergy[getRhIdx(locRHCands[lead],miniRhId)];
				auto nextRhE = miniRhEnergy[getRhIdx(locRHCands[next],miniRhId)];
				if( nextRhE > leadE ) lead = next;
			}//<<>>for( int it(0); it+1 < nLocRHCands; it += 2; )
			locSeedRHs[0] = locRHCands[lead];
			locSeedRHs[1] = locRHCands[lead+1];
			//selPhoIndx[0] = lead/2;
            selPhotons.push_back(locPhotons[lead/2]);
            selPhoType.push_back(0);
		} else { //<<>>if( nLocRHCands > 2 )
            locSeedRHs[0] = locRHCands[0];
            locSeedRHs[1] = locRHCands[1];
            //selPhoIndx[0] = 0;
            selPhotons.push_back(locPhotons[0]);
            selPhoType.push_back(0);
		}//<<>>if( nLocRHCands > 2 )
    }//<<>>if( nLocRHCands > 0 )

	if( DEBUG ) std::cout << " Storing ids Global : " << gloSeedRHs[0] << ", " << gloSeedRHs[1]; 
    if( DEBUG ) std::cout << " and Local: " << locSeedRHs[0] << ", " << locSeedRHs[1] << std::endl;

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
	bool hasResRHs( false );
	for( auto lgRhId : lgRhIds ){
					
		auto goodrh = lgRhId > 0;
		if( goodrh ) hasResRHs = true;
        auto idx = goodrh ? getRhIdx(lgRhId,miniRhId) : 0;
        if( DEBUG ) std::cout << " - Storing : " << " id: " << lgRhId << " idx: " << idx << std::endl;
		resRhID.push_back( goodrh ? lgRhId : 0 );
		resAmp.push_back( goodrh ? miniRhAmp[idx] : -999 );
        resE.push_back( goodrh ? miniRhEnergy[idx] : -999 );
        resRtTime.push_back( goodrh ? miniRhRtTime[idx] : -999 );
        resCCTime.push_back( goodrh ? miniRhCCTime[idx] : -999 );
		resTOF.push_back( goodrh ? miniRhTOF[idx] : -999 );

	}//<<>>for( auto lgRhId : lgRhIds )
	storeEvent = storeEvent || hasResRHs;
    if( DEBUG ) std::cout << " Stored times Global : " << resRtTime[2] << ", " << resRtTime[3];
    if( DEBUG ) std::cout << " and Local: " << resRtTime[0] << ", " << resRtTime[1] << std::endl;

	//gloDiMass = pairMass; gloDiAngle = ang12; gloDiDr = dr12;

	if( DEBUG ) std::cout << "Skimming Selected Photons " << std::endl;

	phoDiMass = gloDiMass;
	phoDiAngle = gloDiAngle;
	phoDiDr = gloDiDr;
    phoDiPhi = gloDiPhi;
    phoDiEta = gloDiEta;

	phoEnergy.clear();
	phoRhIds.clear();
    phoPt.clear();
    phoEta.clear();
    phoPhi.clear();
	phoHadOverEM.clear();
	phoSigmaIEtaIEta.clear();
	phoCov2IEtaIEta.clear();
	phoCov2IEtaIPhi.clear();
	phoCov2IPhiIPhi.clear();
	phoEcalRHSumEtConeDR04.clear();
	phoHcalTwrSumEtConeDR04.clear();
	phoTrkSumPtSolidConeDR04.clear();
	phoTrkSumPtHollowConeDR04.clear();
	phoR9.clear();
	phoSelType.clear();

	int idx(0);
	for( auto photon : selPhotons ){
		
		if( DEBUG ) std::cout << " - skimming photon with idx: " << idx << std::endl;
        const auto &phosc = photon.superCluster().isNonnull() ? photon.superCluster() : photon.parentSuperCluster();
        const auto scptr = phosc.get();
        scGroup phoSCGroup{*scptr};
		auto phoRhGroup = getRHGroup( phoSCGroup, 0.0 );
		auto phoRhIdsGroup = getRhGrpIDs( phoRhGroup );
        const auto seedDetId = scptr->seed()->seed(); // seed detid
        const auto isEB = (seedDetId.subdetId() == EcalSubdetector::EcalBarrel); // which subdet
        const auto phoRecHits = (isEB ? *recHitsEB_ : *recHitsEE_ );
		auto lCov = EcalClusterTools::localCovariances( *(phosc->seed()), &phoRecHits, &(*topology));
		phoCov2IEtaIEta.push_back(std::sqrt(lCov[0]));
		phoCov2IEtaIPhi.push_back(std::sqrt(lCov[1]));
		phoCov2IPhiIPhi.push_back(std::sqrt(lCov[2]));
		phoRhIds.push_back(phoRhIdsGroup);
		phoEnergy.push_back(photon.energy());  
        phoPt.push_back(photon.pt());
        phoEta.push_back(photon.eta());
        phoPhi.push_back(photon.phi());
		phoHadOverEM.push_back(photon.hadronicOverEm());
		phoSigmaIEtaIEta.push_back(photon.sigmaIetaIeta() );
		phoEcalRHSumEtConeDR04.push_back(photon.ecalRecHitSumEtConeDR04());
		phoHcalTwrSumEtConeDR04.push_back(photon.hcalTowerSumEtConeDR04());
		phoTrkSumPtSolidConeDR04.push_back(photon.trkSumPtSolidConeDR04());
		phoTrkSumPtHollowConeDR04.push_back(photon.trkSumPtHollowConeDR04());
		phoR9.push_back(photon.r9());
		phoSelType.push_back(selPhoType[idx]);
		idx++;

	}//<<>>for( auto idx : selPhoIndx )

	// -- Fill output trees ------------------------------------------
	if( DEBUG ) std::cout << "---------- Next Event -----" << std::endl;
	outTree->Fill();
	//if( storeEvent ) outTree->Fill();

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

	//totrhs = 0;
    //totrhs0 = 0;
    //totrhs05 = 0;
    //totrhs1 = 0;
    //totrhs2 = 0;
    //totrhs5 = 0;
    //totrhs10 = 0;
	//encrhs = 0;

	// Book output files and trees
	edm::Service<TFileService> fs;
	outTree = fs->make<TTree>("llpgtree","llpgtree");

	// Book histograms
	// hist1d[11] = fs->make<TH1D>("jetBcClTimeDiff", "jetBcClTimeDiff", jdtdiv, -1*jdtran, jdtran);
    //------ 1D Hists --------------------------------------------------------------------------

    //hist1d[0] = fs->make<TH1D>("wth_v_acc_all","Width v Acc E > 0 GeV;Width/2 [ns];acc",16,0,16);
    //hist1d[1] = fs->make<TH1D>("wth_v_acc_p5","Width v Acc E > 0.5 GeV;Width/2 [ns];acc",16,0,16);
    //hist1d[2] = fs->make<TH1D>("wth_v_acc_1","Width v Acc E > 1 GeV;Width/2 [ns];acc",16,0,16);
    //hist1d[3] = fs->make<TH1D>("wth_v_acc_2","Width v Acc E > 2 GeV;Width/2 [ns];acc",16,0,16);
    //hist1d[4] = fs->make<TH1D>("wth_v_acc_5","Width v Acc E > 5 GeV;Width/2 [ns];acc",16,0,16);
    //hist1d[5] = fs->make<TH1D>("wth_v_acc_10","Width v Acc E > 10 GeV;Width/2 [ns];acc",16,0,16);


	//example//hist1d[0] = new TH1D("run",addstr(fHTitle,"Run;Run").c_str(),88000,275000,363000);	

	//------ 2D Hists --------------------------------------------------------------------------

	//hist2d[0] = fs->make<TH2D>("amp_v_e_unrh","UnCaliRH Amp v Energy;amplitude;energy [GeV]",4000,0,400,200,0,20);
    //hist2d[1] = fs->make<TH2D>("amp_v_e_rh","RH Amp v Energy;amplitude;energy [GeV]",4000,0,400,200,0,20);
    //hist2d[2] = fs->make<TH2D>("untime_v_actdif","untime_v_actdif;nocorrtime [ns]; actdiff [ns]",500,-25,25,500,-25,25);//untime, actdif );
    //hist2d[3] = fs->make<TH2D>("untime_v_encdif","untime_v_encdif;nocorrtime [ns]; encdiff [ns]",500,-25,25,500,-25,25);//untime, encdif );
    //hist2d[4] = fs->make<TH2D>("untime_v_resd","untime_v_resd;nocorrtime [ns]; residual [ns]",500,-25,25,5000,-25,25);//untime, reso );
    //hist2d[5] = fs->make<TH2D>("encdif_v_amp","encdif_v_amp;encdiff [ns]; amplitude",500,-25,25,400,0,400);
    //hist2d[6] = fs->make<TH2D>("actdif_v_amp","encdif_v_amp;actdiff [ns]; amplitude",500,-25,25,400,0,400);
    //hist2d[7] = fs->make<TH2D>("encdif_v_energy","encdif_v_energy;encdiff [ns]; energy [GeV]",500,-25,25,200,0,100);

	std::cout << "Histograms Booked" << std::endl;

	// Create output Tree branches -----------------------------

	// Run, Lumi, Event info
	outTree->Branch("run", &run);
	outTree->Branch("lumi", &lumi);
	outTree->Branch("event", &event, "event/l");

    outTree->Branch("rhCaliID", &rhCaliID);
    outTree->Branch("rhCaliEnergy", &rhCaliEnergy);
    outTree->Branch("rhCaliRtTime", &rhCaliRtTime);
    outTree->Branch("rhCaliCCTime", &rhCaliCCTime);

    outTree->Branch("resRhID", &resRhID);
    outTree->Branch("resAmp", &resAmp);
    outTree->Branch("resE", &resE);
    outTree->Branch("resRtTime", &resRtTime);
    outTree->Branch("resCCTime", &resCCTime);
    outTree->Branch("resTOF", &resTOF);

    if( doTwoTier ){

        outTree->Branch("unrhJitter", &unrhJitter);
        outTree->Branch("unrhNonJitter", &unrhNonJitter);
        outTree->Branch("unrhEncNonJitter", &unrhEncNonJitter);
        outTree->Branch("unrhEnergy", &unrhEnergy);

    }//<<>>if( doTwoTier )

	if( doDiag ){

    	outTree->Branch("rhID", &rhID);
    	outTree->Branch("rhRtTime", &rhRtTime);
        outTree->Branch("rhCCTime", &rhCCTime);
        //outTree->Branch("rhTOF", &rhTOF);
        outTree->Branch("rhEnergy", &rhEnergy);
        //outTree->Branch("rhAmp", &rhAmp);

	    outTree->Branch("rhRtisOOT", &rhRtisOOT);
        outTree->Branch("rhCCisOOT", &rhCCisOOT);
	    //outTree->Branch("rhisGood", &rhisGood);
	    outTree->Branch("rhisWeird", &rhisWeird);
	    outTree->Branch("rhisDiWeird", &rhisDiWeird);
	    outTree->Branch("rhSwCross", &rhSwCross);
        //outTree->Branch("rhisGS6", &rhisGS6);
        //outTree->Branch("rhisGS1", &rhisGS1);
        //outTree->Branch("rhadcToGeV", &rhadcToGeV);
        //outTree->Branch("rhpedrms12", &rhpedrms12);

        outTree->Branch("phoEnergy", &phoEnergy);
        outTree->Branch("phoRhIds", &phoRhIds);
        outTree->Branch("phoPt", &phoPt);
        outTree->Branch("phoEta", &phoEta);
        outTree->Branch("phoPhi", &phoPhi);	
	    outTree->Branch("phoHadOverEM", &phoHadOverEM);
	    outTree->Branch("phoSigmaIEtaIEta", &phoSigmaIEtaIEta);
        outTree->Branch("phoCov2IEtaIEta", &phoCov2IEtaIEta);
        outTree->Branch("phoCov2IEtaIPhi", &phoCov2IEtaIPhi);
        outTree->Branch("phoCov2IPhiIPhi", &phoCov2IPhiIPhi);

	    outTree->Branch("phoEcalRHSumEtConeDR04", &phoEcalRHSumEtConeDR04);
	    outTree->Branch("phoHcalTwrSumEtConeDR04", &phoHcalTwrSumEtConeDR04);
	    outTree->Branch("phoTrkSumPtSolidConeDR04", &phoTrkSumPtSolidConeDR04);
	    outTree->Branch("phoTrkSumPtHollowConeDR04", &phoTrkSumPtHollowConeDR04);
		outTree->Branch("phoR9", &phoR9);
        outTree->Branch("phoSelType", &phoSelType);

        outTree->Branch("phoDiMass", &phoDiMass);
        outTree->Branch("phoDiAngle", &phoDiAngle);
        outTree->Branch("phoDiDr", &phoDiDr);
        outTree->Branch("phoDiPhi", &phoDiPhi);
        outTree->Branch("phoDiEta", &phoDiEta);

	}//<<>>if( doDiag )

}//>>>>void GammaResTool::beginJob()


// ------------ method called once each job just after ending the event loop	------------
void GammaResTool::endJob(){

	//hist1d[0]->Scale(1/totrhs0);
    //hist1d[1]->Scale(1/totrhs05);
    //hist1d[2]->Scale(1/totrhs1);
    //hist1d[3]->Scale(1/totrhs2);
    //hist1d[4]->Scale(1/totrhs5);
    //hist1d[5]->Scale(1/totrhs10);

	if( nEvents == 0 ) nEvents = 1;
	//if( DEBUG ) 
	std::cout << " Found Global % " << gloMatches/nEvents << " and Local % " << locMatches/nEvents << " of " << nEvents << std::endl;
	//std::cout << " Range Acceptance : " << encrhs << " / " << totrhs << " = " << encrhs/totrhs << std::endl;

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
