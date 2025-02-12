//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//`
//
// Original Author:  Jack W King III
//         Created:  Fri, 20 Dec 2024
//
//////////////////////////////////////////////////////////////////////

#include "KUCMSTimeCalibration.hh"
// ------------------------------------------- main function ------------------------------------------------------------
int main ( int argc, char *argv[] ){

	KUCMSTimeCalibration theCali;

    ///////////////////////////////////////////////////////////////////
    // this section should be moved to run file --- DONE --------------
    // no need to "recreate" cali maps every time, they will be saved -
    ///////////////////////////////////////////////////////////////////

	// Cali Tags : Tags for calibrations to use
	std::string r2EOY( "EG_EOY_MINI" );
	std::string r2Fall17AOD( "RunIIFall17DRPremix" );
    std::string r2Fall17MINIAOD( "RunIIFall17MiniAODv2" );

	// IOV tags :  defines IOV periods to make calibration maps
    // for mc use "mc" for TTIov and XIov

    //-----//////////  making tt cali file :
    std::string r2ulTag( "r2ul" );
    std::string mctag( "mc" );
    std::string xiovtag( "prompt" );

    //std::string inputfilename( "kucmsTimeCaliTestFile.txt" ); // MET_AOD_R17_FULL 
    //std::string inputfilename( "kucmsTimeCaliR17File.txt" ); // EG_EOY_MINI PD
    std::string inputfilename( "kucmsTimeCaliRunIIFall17File.txt" ); // RunIIFall17 MC

	//std::string eosdir("root://cmseos.fnal.gov//store/user/jaking/");// input parameter!
    std::string eosdir("root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/");
	//std::string indir("/ecalTiming/gammares_llpana_pd/");// input paramter !  
    //std::string indir("/ecalTiming/gammares_llpana_mc/");// input paramter ! 
    std::string indir("/kuncali/gammares_cali_mc/");
    //std::string indir("/ecalTiming/gammares_cali_mc/");

	theCali.SetEosDir(eosdir);
	theCali.SetInDir(indir);


    //-----//////////  making tt cali  :
	//theCali.LoadCaliHists( true ); // bool stayOpen
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
    //theCali.LoadCaliHists( true ); // bool stayOpen
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

    // make smear maps : 
    //  - 2d res plot
    //theCali.LoadCaliHists( true );
	//theCali.plot2dResolutionEGR( inputfilename );
    //theCali.SaveCaliHists();


	//theCali.LoadCaliHists( true );
    //theCali.plot2dResolutionEGR( inputfilename, false, false );// input, small, usecali
	//theCali.setTargetTag("EG_EOY_MINI");	
    //theCali.setTargetRun(300202);
    //theCali.setTargetRun(304476);
    //theCali.setSourceTag("DYJetsToLLRunIIFall17AOD");
    //theCali.setSourceRun(0);
    //theCali.plot2dResolutionEGR( inputfilename, false, true, true, "_v2" );// input, small, usecali, smear
	//theCali.SaveCaliHists();

	// make resolution paramters for smearing
    //theCali.doResTimeFits();
	theCali.LoadCaliHists( true );
    theCali.load2DResHist( "DYJetsToLLRunIIFall17AOD_0_X_ZEE_Data_Hist_Smeared_v2" );
	theCali.doResTimeFit( "DYJetsToLLRunIIFall17AOD_0_X_ZEE_Data_Hist_Smeared_v2" );
	theCali.SaveCaliHists();

	// plotting of mean time by run with calibraton : filename, start run, end run, usecali
	//theCali.plotMeanRunTimeEGR( inputfilename, 303838, 304796 );
    //theCali.plotMeanRunTimeEGR( inputfilename, 296399, 307554, true );
	//theCali.makeTTDiffMaps();// make trigger tower diffrence maps
	
    /////////////////////////////////////////////////////////////////////
    //-------------------------------------------------------------------
	
	std::cout << " -- Thats All Folks !!!!!!!!! " << std::endl;

    return 1;

}//<<>>int main ( int argc, char *argv[] )
