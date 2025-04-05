//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//`
//
// Original Author:  Jack W King III
//         Created:  Fri, 20 Dec 2024
//
//////////////////////////////////////////////////////////////////////

#include "ecal_config/KUCMSTimeCalibration.hh"
// ------------------------------------------- main function ------------------------------------------------------------
int main ( int argc, char *argv[] ){

    ///////////////////////////////////////////////////////////////////
    // this section should be moved to run file --- DONE --------------
    // no need to "recreate" cali maps every time, they will be saved -
    ///////////////////////////////////////////////////////////////////
    // standard usage for applying calibrations :
    // KUCMSTimeCalibration theCali;
    // theCali.setTag( tag ); 
    // 	-- or to use default skip this step 
    // 	-- only needs to be done once unless you want to change the tag to use
    // time = time - theCali.getCalibration( rechitDetID, run );
    //
    // to get smeared time :
    // theCali.setSmearTag( tag ); 
    //  -- or to use default skip this step
    //  -- only needs to be done once unless you want to change the tag to use
    // time = theCali.getSmearedTime( time, amplitude );
    ///////////////////////////////////////////////////////////////////

	// Cali Tags : Tags for calibrations to use
	std::string r2EOY( "EG_EOY_MINI" );
	std::string r2Fall17AOD( "RunIIFall17DRPremix" );
    std::string r2Fall17MINIAOD( "RunIIFall17MiniAODv2" );

	// IOV tags :  defines IOV periods to make calibration maps
    // for mc use "mc" for TTIov and XIov

    //-----//////////  making tt cali file :
    //std::string r2ulTag( "r2ul" );
    std::string r2ulTag( "r2ultt" );
    std::string mctag( "mc" );
    //std::string xiovtag( "prompt" );
    std::string xiovtag( "r2ulx" );

    std::string inputfilename( "kucmsTimeCaliR24FCCvRtTFile.txt");
    //std::string inputfilename( "kucmsTimeCaliTestFile.txt" ); // MET_AOD_R17_FULL 
    //std::string inputfilename( "kucmsTimeCaliR17File.txt" ); // EG_EOY_MINI PD

	std::string eosdir("root://cmseos.fnal.gov//store/user/jaking/");// input parameter!
    //std::string eosdir("root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/");

	//std::string indir("ecalTiming/gammares_llpana_pd/");// input paramter !  
    //std::string indir("/ecalTiming/gammares_llpana_mc/");// input paramter ! 
    //std::string indir("/kuncali/gammares_cali_mc/");
    //std::string indir("/kuncali/gammares_cali/");
    //std::string indir("/ecalTiming/gammares_cali_mc/");
    std::string indir("ecalTiming/gammares_ECAL_CC_HCAL_DI-v3/");
    
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//using true - true will start a new TClifile  - all histograms will be lost
    //////KUCMSTimeCalibration theCali( true, true ); // make a new TClifile and keep TCalifile open
	//KUCMSTimeCalibration theCali( true, true );
	// if Tcalifile kept open be sure to make and save tcalifile ( saving calihist will close tcalifile )

	////KUCMSTimeCalibration theCali; // open tfile in read only then close after setup
	//KUCMSTimeCalibration theCali( true );
	//theCali.SetEosDir(eosdir);
	//theCali.SetInDir(indir);

    //-----//////////  making tt cali  :
    //theCali.SetupIovMaps();
	// for PD
    //theCali.setTTIov( r2ulTag );
    //theCali.setXIov( xiovtag );
	// for MC
    //theCali.setTTIov( mctag );
    //theCali.setXIov( mctag );
    //theCali.makeTTCaliMapEGR( inputfilename, true ); // true == run only subset of events
    //theCali.makeTTCaliMapEGR( inputfilename );
    //theCali.makeCaliHists();
    //theCali.SaveCaliHists();
    //theCali.SaveTTRunFile();

    //-----//////////  making xtal cali :
    //theCali.SetupIovMaps();
    // for PD
    //theCali.setTTIov( r2ulTag );
    //theCali.setXIov( xiovtag );
    // for MC
    //theCali.setTTIov( mctag );
    //theCali.setXIov( mctag );
    //theCali.makeXCaliMapEGR( inputfilename, true ); // true == run only subset of events
    //theCali.makeXCaliMapEGR( inputfilename );
    //theCali.makeCaliHists();
    //theCali.SaveCaliHists();
    //theCali.SaveCaliRunFile();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // 2d resolution maps  
	//KUCMSTimeCalibration theCali; 
    //theCali.SetEosDir(eosdir);
    //theCali.SetInDir(indir);
    //  - 2d res plot  -- defaults to : false, true, false, "" ( scale, usecali, smear, name ext ) 
    //theCali.plot2dResolutionEGR( inputfilename );

	// 2d resolution maps smeared/uncali
    //theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_nocali" );// input, scale, usecali
    //theCali.setSmearTag("EG300202_DYF17");
    //theCali.plot2dResolutionEGR( inputfilename, true, true, true, "_smeared" );// input, scale, usecali, smear, name ext

	// make resolution paramters for iov periods
    //theCali.doResTimeFits();
	//theCali.SaveCaliRunFile();	

    // make resolution paramters for smearing
    //std::string sourceTag = "DYJetsToLLRunIIFall17AOD_0_X_ZEE_Data_Hist";
	//std::string run = "304476";
	//std::string destTag = "EG_EOY_MINI_" + run + "_X_ZEE_Data_Hist";
	//std::string smearTag = "EG" + run + "_DYF17";
    //theCali.makeSmearTag( sourceTag, destTag, smearTag );
	//theCali.SaveSmearFile();

	// standard resolution 
	//std::string histName = "DYJetsToLLRunIIFall17AOD_0_X_ZEE_Data_Hist_Smeared_v4";
    //std::string histName = "EG_EOY_MINI_304476_X_SRO_Data_Hist_NoCali_nocali_v4";
    //theCali.load2DResHist( histName );
	//theCali.doResTimeFit( histName );

	// extended range 2D
	//KUCMSTimeCalibration theCali;
    //theCali.SetEosDir(eosdir);
    //theCali.SetInDir(indir);
    //theCali.setLowEnergy( true );
    //theCali.setUseEffEnergy( true );
    //theCali.setUseEffEnergy( false );
    //theCali.SetXBinStr( "VARIABLE 5 10 15 20 25 30 40 50" ); // xle
    //theCali.SetXBinStr( "VARIABLE 5 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 1200 1800" ); // : xa
    //theCali.SetXBinStr( "VARIABLE 5 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 1800" ); // xj : justin
    //theCali.SetXBinStr( "VARIABLE 5 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 1200" ); // xb : * w/LE
	//theCali.SetXBinStr( "VARIABLE 1.0 2.0 3.0 4.0 5.0 7.5 10.0 12.5 15.0 17.5 20.0 22.5 25.0 30.0 40.0 60.0 120.0" ); // xea
    ////theCali.SetYBinStr( "CONSTANT 1800 -9 9" );
    //theCali.SetYBinStr( "CONSTANT 720 -9 9" );
	// justin profile y bins
    //theCali.SetYBinStr( "CONSTANT 240 -12 12" );
    //theCali.SetYBinStr( "CONSTANT 180 -9 9" );
    //theCali.SetYBinStr( "CONSTANT 180 -6 6" ); // *
    //theCali.SetYBinStr( "CONSTANT 180 -3 3" );
    //theCali.plot2dResolutionEGR( inputfilename, true, true, false, "_xea_pm12b240_v2" );// scale, cali, smear
    //theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_xea_pm9b720_v2" );// scale, cali, smear
    //theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_xa_pm9b720_r324fcc" );// scale, cali, smear

    // extended range resfit
    KUCMSTimeCalibration theCali;
	//std::string histName = "ResMap_302031_302393__ZEE_Data_Hist_std_v2";
	//std::string histName = "ResMap_303832_304616__DRO_Data_Hist_xb_pm9b180_v1";
	//std::string histName = "ResMap_303832_304616__DRO_Data_Hist_xea_pm9b180_v1";
	//std::string histName = "ResMap_303832_304616__SRO_Data_Hist_xea_pm9b720_v1";
	//std::string histName = "ResMap_305044_305081__ZEE_Data_Hist_xea_pm9b720_v1";
	//std::string histName = "ResMap_302031_302393__DRO_Data_Hist_xea_pm12b240_v2";
	std::string histName = "ResMap_0_999999__SRO_Data_Hist_NoCali_xa_pm9b720_r324fcc";
	theCali.load2DResHist( histName );
	theCali.setLowEnergy( true );
	//theCali.SetXBinStr( "VARIABLE 0.2 0.5 1 2 5 10 15 20 25 30 40 50" ); // xle
	theCali.SetXBinStr( "VARIABLE 5 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 1200 1800" ); // : xa
    //theCali.SetXBinStr( "VARIABLE 5 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 1200" ); // xb : * w/LE
    //theCali.SetXBinStr( "VARIABLE 1.0 2.0 3.0 4.0 5.0 7.5 10.0 12.5 15.0 17.5 20.0 22.5 25.0 30.0 40.0 60.0 120.0" ); // xea
    //theCali.SetXBinStr( "VARIABLE 5 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 1800" ); // xj
    //theCali.SetYBinStr( "CONSTANT 180 -6 6" );
    //theCali.SetYBinStr( "CONSTANT 180 -3 3" );
    theCali.SetYBinStr( "CONSTANT 720 -9 9" );
	//theCali.SetYBinStr( "CONSTANT 180 -9 9" );
	//theCali.SetYBinStr( "CONSTANT 240 -12 12" );
	theCali.doResTimeFit( histName );

	// plotting of mean time by run with calibraton : filename, start run, end run, usecali
	//theCali.plotMeanRunTimeEGR( inputfilename, 303838, 304796 );
    //theCali.plotMeanRunTimeEGR( inputfilename, 296399, 307554, true );
	//theCali.makeTTDiffMaps();// make trigger tower diffrence maps
	
    /////////////////////////////////////////////////////////////////////
    //-------------------------------------------------------------------
	
	std::cout << " -- Thats All Folks !!!!!!!!! " << std::endl;

    return 1;

}//<<>>int main ( int argc, char *argv[] )
