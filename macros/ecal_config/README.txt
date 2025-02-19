/////////////////////////////////////////////////////
//
// 
// KUCMS Calibration and Smearing
// KUCMSTimeCalibration.hh
// Author:  Jack W King III 
// Created:  Wed, 18 Sept 2024 19:19:35 GMT
//
//  ecal_config/
//	KUCMSEcalDetIDFunctions.hh
//	KUCMSHelperBaseClass.hh
//	KUCMSRootHelperBaseClass.hh
//	KUCMSTimeCalibration.hh
//	README.txt
//	UL2016_runlumi.txt
//	UL2017_runlumi.txt
//	UL2018_runlumi.txt
//	caliHistsTFile.root
//	caliRunConfig.txt
//	caliSmearConfig.txt
//	caliTTConfig.txt
//	fullinfo_detids_EB.txt
//	fullinfo_detids_EE.txt
//	fullinfo_v2_detids_EB.txt
//	fullinfo_v2_detids_EE.txt
//	howto.txt
//	reducedinfo_detids.txt
//	rhid_i12_list.txt
//	rhid_info_list.txt
//
//	fillKUCMSTimeCalibration.cpp ( used to set up calibration&smearing info )
//
////////////////////////////////////////////////////

Usage :

declare class in intial setup( do once )

KUCMSTimeCalibration theCali;

To use calibration :

	Set tag to use : ( default EG_EOY_MINI )
		theCali.setTag("EG_EOY_MINI");
	To get calibration :
		calibration = theCali.getCalibration( rechitID, run );
	To apply calibration subtract from rechit time :
		calibratedtime = rechittime - calibration;
			
		- or -

    Set tag to use : ( default EG_EOY_MINI )
        theCali.setTag("EG_EOY_MINI");
	To get calibrated time :
		calibratedtime = theCali.getCalibTime( rechittime, rechitID, run );


To use smearing :

    Set tag to use : ( default EG300202_DYF17 )
        theCali.setSmearTag("EG300202_DYF17");
	To get smeared time :
		smearedTime = theCali.getSmearedTime( rechitTime, rechitAmplitude );


`





