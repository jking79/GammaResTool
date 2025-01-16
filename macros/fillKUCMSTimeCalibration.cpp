//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
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

    // for mc use "mc" for TTIov and XIov

    //-----//////////  making tt cali file :
    std::string r2ulTag( "r2ul" );
    std::string mctag( "mc" );
    std::string xiovtag( "prompt" );
    std::string inputfilename( "kucmsTimeCaliTestFile.txt" );

	//-----////////// prep for making cali maps :
    //theCali.SetupIovMaps();
    //-----//////////  making tt cali  :
    //theCali.setTTIov( r2ulTag );
    //theCali.makeTTCaliMapEGR( inputfilename );
    //-----//////////  making xtal cali :
    //theCali.setXIov( xiovtag );
    //theCali.makeXCaliMapEGR( inputfilename );

    // make smear maps : 
    //  - 2d res plot
    //theCali.setXIov( xiovtag );
	//theCali.plot2dResolutionEGR( inputfilename );
	
    /////////////////////////////////////////////////////////////////////
    //-------------------------------------------------------------------
	
	std::cout << " -- Thats All Folks !!!!!!!!! " << std::endl;

    return 1;

}//<<>>int main ( int argc, char *argv[] )
