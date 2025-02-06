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
	std::string r2Fall17( "RunIIFall17DRPremix" );

	// IOV tags :  defines IOV periods to make calibration maps
    // for mc use "mc" for TTIov and XIov

    //-----//////////  making tt cali file :
    std::string r2ulTag( "r2ul" );
    std::string mctag( "mc" );
    std::string xiovtag( "prompt" );

    //std::string inputfilename( "kucmsTimeCaliTestFile.txt" ); // MET_AOD_R17_FULL 
    //std::string inputfilename( "kucmsTimeCaliR17File.txt" ); // EG_EOY_MINI PD
    std::string inputfilename( "kucmsTimeCaliRunIIFall17DRPremixFile.txt" ); // RunIIFall17DRPremix MC

	std::string eosdir("root://cmseos.fnal.gov//store/user/jaking/");// input parameter!
    //std::string eosdir("root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/");
	//std::string indir("/ecalTiming/gammares_llpana_pd/");// input paramter !  
    //std::string indir("/ecalTiming/gammares_llpana_mc/");// input paramter ! 
    //std::string indir("/kuncali/gammares_cali/");
    std::string indir("/ecalTiming/gammares_cali_mc/");

	theCali.SetEosDir(eosdir);
	theCali.SetInDir(indir);

    //-----//////////  making tt cali  :
    theCali.SetupIovMaps();
	// for PD
    //theCali.setTTIov( r2ulTag );
    //theCali.setXIov( xiovtag );
	// for MC
    //theCali.setTTIov( mctag );
    //theCali.setXIov( mctag );
    //theCali.makeTTCaliMapEGR( inputfilename, true ); // true == run only subset of events
    //theCali.makeTTCaliMapEGR( inputfilename );

    //-----//////////  making xtal cali :
    //theCali.SetupIovMaps();
    // for PD
    //theCali.setTTIov( r2ulTag );
    //theCali.setXIov( xiovtag );
    // for MC
    //theCali.setTTIov( mctag );
    //theCali.setXIov( mctag );
    //theCali.makeXCaliMapEGR( inputfilename, true );
    //theCali.makeXCaliMapEGR( inputfilename );

    // make smear maps : 
    //  - 2d res plot
	theCali.plot2dResolutionEGR( inputfilename );
    //theCali.plot2dResolutionEGR( inputfilename, false, false );// input, small, usecali

	// make resolution paramters for smearing
    //theCali.doResTimeFits();
	//theCali.doResTimeFits( false, true );// doLocal, doNoCali

	// plotting of mean time by run with calibraton : filename, start run, end run, usecali
	//theCali.plotMeanRunTimeEGR( inputfilename, 303838, 304796 );
    //theCali.plotMeanRunTimeEGR( inputfilename, 296399, 307554, true );
	//theCali.makeTTDiffMaps();// make trigger tower diffrence maps
	
    /////////////////////////////////////////////////////////////////////
    //-------------------------------------------------------------------
	
	std::cout << " -- Thats All Folks !!!!!!!!! " << std::endl;

    return 1;

}//<<>>int main ( int argc, char *argv[] )
