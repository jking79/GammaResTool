import FWCore.ParameterSet.Config as cms
from Configuration.ProcessModifiers.gpu_cff import gpu

# Calo geometry service model
#
# removed by tommaso
#
#ECAL conditions
#  include "CalibCalorimetry/EcalTrivialCondModules/data/EcalTrivialCondRetriever.cfi"
#
#TPG condition needed by ecalRecHit producer if TT recovery is ON
#from RecoLocalCalo.EcalRecProducers.ecalRecHitTPGConditions_cff import *
#ECAL reconstruction
#from RecoLocalCalo.EcalRecProducers.ecalGlobalUncalibRecHit_cfi import *
#from RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi import *
#from RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi import *
#from RecoLocalCalo.EcalRecProducers.ecalPreshowerRecHit_cfi import *
from RecoLocalCalo.EcalRecProducers.ecalDetIdToBeRecovered_cfi import *
#from RecoLocalCalo.EcalRecProducers.ecalCompactTrigPrim_cfi import *
#from RecoLocalCalo.EcalRecProducers.ecalTPSkim_cfi import *
#from RecoLocalCalo.EcalRecProducers.ecalDetailedTimeRecHit_cfi import *
from RecoLocalCalo.Configuration.ecalLocalRecoSequence_cff import *

from GammaResTool.GammaResTool.jwk_ku_ecalRecHit_cff import *
from GammaResTool.GammaResTool.jwk_ku_cc_ecalRecHit_cff import *
from GammaResTool.GammaResTool.jwk_ku_ecalMultiFitUncalRecHit_cff import *

ku_ecalUncalibRecHitSequence = cms.Sequence(ecalMultiFitUncalibRecHitBase*kuEcalMultiFitUncalibRecHit*ecalDetIdToBeRecovered)

ku_multi_ecalUncalibRecHitSequence = cms.Sequence(ecalMultiFitUncalibRecHitBase*
                                        kuEcalMultiFitUncalibRecHit*
                                        kuWtEcalMultiFitUncalibRecHit*
                                        kuCCEcalMultiFitUncalibRecHit*
                                        ecalDetIdToBeRecovered)

ku_lhc_ecalUncalibRecHitSequence  = cms.Sequence(#ecalMultiFitUncalibRecHitBase*
                                        kuEcalMultiFitUncalibRecHit*
                                        kuWtEcalMultiFitUncalibRecHit*
                                        kuCCEcalMultiFitUncalibRecHit*
                                        ecalDetIdToBeRecovered)

ku_reduced_multi_ecalUncalibRecHitSequence = cms.Sequence(#ecalMultiFitUncalibRecHitBase*
                                        kuEcalMultiFitUncalibRecHit*
                                        kuCCEcalMultiFitUncalibRecHit*
                                        ecalDetIdToBeRecovered)

ku_cc_gt_ecalUncalibRecHitSequence = cms.Sequence(#ecalMultiFitUncalibRecHitBase*
                                        #kuEcalMultiFitUncalibRecHit*
                                        kuCCEcalMultiFitUncalibRecHit*
                                        ecalDetIdToBeRecovered)

ku_cc_native_gt_ecalUncalibRecHitSequence = cms.Sequence(#ecalMultiFitUncalibRecHitBase*
                                        #kuEcalMultiFitUncalibRecHit*
                                        kuCCNativeEcalMultiFitUncalibRecHit*
                                        ecalDetIdToBeRecovered)

ku_reduced_nomulti_ecalUncalibRecHitSequence = cms.Sequence(#ecalMultiFitUncalibRecHitBase*
                                        kuEcalMultiFitUncalibRecHit*
                                        #kuCCEcalMultiFitUncalibRecHit*
                                        ecalDetIdToBeRecovered)

kucc_only_ecalUncalibRecHitSequence = cms.Sequence(kuCCEcalMultiFitUncalibRecHit*ecalDetIdToBeRecovered)

kuEcalLHCRecHit = ecalRecHitBase.clone(
	EErechitCollection = cms.string('kuRecHitsEE'),
	EEuncalibRecHitCollection = cms.InputTag("kuEcalMultiFitUncalibRecHit","kuEcalUncalibRecHitsEE"),
	EBuncalibRecHitCollection = cms.InputTag("kuEcalMultiFitUncalibRecHit","kuEcalUncalibRecHitsEB"),
	EBrechitCollection = cms.string('kuRecHitsEB'),
   skipTimeCalib = cms.bool(False),
   # below params for LHCInfo plots form badder
   killDeadChannels = cms.bool( True ),
   recoverEBVFE = cms.bool( False ),
   recoverEEVFE = cms.bool( False ),
   recoverEBFE = cms.bool( False ),
   recoverEEFE = cms.bool( False ),
   recoverEEIsolatedChannels = cms.bool( False ),
   recoverEBIsolatedChannels = cms.bool( False ),
	)

kuEcalRecHit = ecalRecHitBase.clone(
        EErechitCollection = cms.string('kuRecHitsEE'),
        EEuncalibRecHitCollection = cms.InputTag("kuEcalMultiFitUncalibRecHit","kuEcalUncalibRecHitsEE"),
        EBuncalibRecHitCollection = cms.InputTag("kuEcalMultiFitUncalibRecHit","kuEcalUncalibRecHitsEB"),
        EBrechitCollection = cms.string('kuRecHitsEB'),
        )

kuStcEcalRecHit = ecalRecHitBase.clone(
        EErechitCollection = cms.string('kuStcRecHitsEE'),
        EEuncalibRecHitCollection = cms.InputTag("kuEcalMultiFitUncalibRecHit","kuEcalUncalibRecHitsEE"),
        EBuncalibRecHitCollection = cms.InputTag("kuEcalMultiFitUncalibRecHit","kuEcalUncalibRecHitsEB"),
        EBrechitCollection = cms.string('kuStcRecHitsEB'),
	    skipTimeCalib = cms.bool(True),
        )

kuStcEcalLHCRecHit = ecalRecHitBase.clone(
	EErechitCollection = cms.string('kuStcRecHitsEE'),
    EEuncalibRecHitCollection = cms.InputTag("kuEcalMultiFitUncalibRecHit","kuEcalUncalibRecHitsEE"),
    EBuncalibRecHitCollection = cms.InputTag("kuEcalMultiFitUncalibRecHit","kuEcalUncalibRecHitsEB"),
    EBrechitCollection = cms.string('kuStcRecHitsEB'),
    skipTimeCalib = cms.bool(True),
    # below params for LHCInfo plots from badder
    killDeadChannels = cms.bool( True ),
    recoverEBVFE = cms.bool( False ),
    recoverEEVFE = cms.bool( False ),
    recoverEBFE = cms.bool( False ),
    recoverEEFE = cms.bool( False ),
    recoverEEIsolatedChannels = cms.bool( False ),
    recoverEBIsolatedChannels = cms.bool( False ),
        )

kuWtEcalRecHit = ecalRecHitBase.clone(
        EErechitCollection = cms.string('kuWtRecHitsEE'),
        EEuncalibRecHitCollection = cms.InputTag("kuWtEcalMultiFitUncalibRecHit","kuWtEcalUncalibRecHitsEE"),
        EBuncalibRecHitCollection = cms.InputTag("kuWtEcalMultiFitUncalibRecHit","kuWtEcalUncalibRecHitsEB"),
        EBrechitCollection = cms.string('kuWtRecHitsEB'),
        )

kuWtStcEcalRecHit = ecalRecHitBase.clone(
        EErechitCollection = cms.string('kuWtStcRecHitsEE'),
        EEuncalibRecHitCollection = cms.InputTag("kuWtEcalMultiFitUncalibRecHit","kuWtEcalUncalibRecHitsEE"),
        EBuncalibRecHitCollection = cms.InputTag("kuWtEcalMultiFitUncalibRecHit","kuWtEcalUncalibRecHitsEB"),
        EBrechitCollection = cms.string('kuWtStcRecHitsEB'),
        skipTimeCalib = cms.bool(True),
        )

kuCCStcEcalRecHit = ecalRecHitBase.clone(
        EErechitCollection = cms.string('kuCCStcRecHitsEE'),
        EEuncalibRecHitCollection = cms.InputTag("kuCCEcalMultiFitUncalibRecHit","kuCCEcalUncalibRecHitsEE"),
        EBuncalibRecHitCollection = cms.InputTag("kuCCEcalMultiFitUncalibRecHit","kuCCEcalUncalibRecHitsEB"),
        EBrechitCollection = cms.string('kuCCStcRecHitsEB'),
        skipTimeCalib = cms.bool(True),
        )

kuCCEcalRecHit = ecalRecHitBase.clone(
        EErechitCollection = cms.string('kuCCRecHitsEE'),
        EEuncalibRecHitCollection = cms.InputTag("kuCCEcalMultiFitUncalibRecHit","kuCCEcalUncalibRecHitsEE"),
        EBuncalibRecHitCollection = cms.InputTag("kuCCEcalMultiFitUncalibRecHit","kuCCEcalUncalibRecHitsEB"),
        EBrechitCollection = cms.string('kuCCRecHitsEB'),
        skipTimeCalib = cms.bool(True),
        )

kuCCNativeEcalRecHit = ecalRecHitCCBase.clone(
        EErechitCollection = cms.string('kuCCRecHitsEE'),
        EEuncalibRecHitCollection = cms.InputTag("kuCCNativeEcalMultiFitUncalibRecHit","EcalUncalibRecHitsEE"),
        EBuncalibRecHitCollection = cms.InputTag("kuCCNativeEcalMultiFitUncalibRecHit","EcalUncalibRecHitsEB"),
        EBrechitCollection = cms.string('kuCCRecHitsEB'),
        #skipTimeCalib = cms.bool(True),
        )


kuCCStcEcalLHCRecHit = ecalRecHitBase.clone(
        EErechitCollection = cms.string('kuCCStcRecHitsEE'),
        EEuncalibRecHitCollection = cms.InputTag("kuCCEcalMultiFitUncalibRecHit","kuCCEcalUncalibRecHitsEE"),
        EBuncalibRecHitCollection = cms.InputTag("kuCCEcalMultiFitUncalibRecHit","kuCCEcalUncalibRecHitsEB"),
        EBrechitCollection = cms.string('kuCCStcRecHitsEB'),
        skipTimeCalib = cms.bool(True),
        # below params for LHCInfo plots from badder
        killDeadChannels = cms.bool( False ),
        recoverEBVFE = cms.bool( False ),
        recoverEEVFE = cms.bool( False ),
        recoverEBFE = cms.bool( False ),
        recoverEEFE = cms.bool( False ),
        recoverEEIsolatedChannels = cms.bool( False ),
        recoverEBIsolatedChannels = cms.bool( False ),
        )

ku_ecalRecHitSequence        = cms.Sequence(ecalRecHitBase*
					    kuEcalRecHit
                                            #ecalCompactTrigPrim*
                                            #ecalTPSkim+
                                            #ecalPreshowerRecHit
				            )

ku_min_ecalRecHitSequence        = cms.Sequence(ecalRecHitBase*kuEcalRecHit)

kucc_only_ecalRecHitSequence        = cms.Sequence(kuCCStcEcalRecHit)

ku_multi_ecalRecHitSequence        = cms.Sequence(kuEcalRecHit*
					        #kuStcEcalRecHit*
                                                kuWtStcEcalRecHit*
                                                kuCCStcEcalRecHit
					       )

ku_lhc_ecalRecHitSequence        = cms.Sequence(kuEcalLHCRecHit*
                                                  #kuStcEcalRecHit*
                                                  kuWtStcEcalRecHit*
                                                  kuCCStcEcalLHCRecHit
                                               )

ku_reduced_multi_ecalRecHitSequence        = cms.Sequence(#kuEcalRecHit*
                                                  kuStcEcalRecHit*
                                                  #kuWtStcEcalRecHit*
                                                  kuCCStcEcalRecHit
                                               )

ku_cc_gt_ecalRecHitSequence        = cms.Sequence(#kuEcalRecHit*
                                                  kuCCEcalRecHit
                                                  #kuCCStcEcalRecHit
                                                  #kuWtStcEcalRecHit*
                                                  #kuCCStcEcalRecHit
                                               )

ku_cc_native_gt_ecalRecHitSequence        = cms.Sequence(#kuEcalRecHit*
                                                  #kuCCEcalRecHit
                                                  kuCCNativeEcalRecHit
                                                  #kuCCStcEcalRecHit
                                                  #kuWtStcEcalRecHit*
                                                  #kuCCStcEcalRecHit
                                               )


ku_reduced_flipped_ecalRecHitSequence     = cms.Sequence(#kuEcalRecHit*
                                                  kuStcEcalRecHit*
                                                  #kuWtStcEcalRecHit*
                                                  kuCCEcalRecHit
                                               )


ku_spike_multi_ecalRecHitSequence        = cms.Sequence(#kuEcalRecHit*
                                                  kuStcEcalLHCRecHit*
						  #kuStcEcalRecHit*
                                                  kuCCStcEcalLHCRecHit
                                               )

ku_spike_nomulti_ecalRecHitSequence        = cms.Sequence(#kuEcalRecHit*
                                                  kuStcEcalLHCRecHit
                                                  #kuStcEcalRecHit*
                                                  #kuCCStcEcalLHCRecHit
                                               )

# full sequences
ku_ecalLocalRecoSequence     	= cms.Sequence(ku_ecalUncalibRecHitSequence*ku_ecalRecHitSequence)

ku_min_ecalLocalRecoSequence    = cms.Sequence(ku_ecalUncalibRecHitSequence*kucc_only_ecalRecHitSequence)

kucc_only_ecalLocalRecoSequence  	= cms.Sequence(kucc_only_ecalUncalibRecHitSequence*kucc_only_ecalRecHitSequence)

ku_multi_ecalLocalRecoSequence   = cms.Sequence(ku_multi_ecalUncalibRecHitSequence*ku_multi_ecalRecHitSequence)

ku_lhc_ecalLocalRecoSequence   = cms.Sequence(ku_lhc_ecalUncalibRecHitSequence*ku_lhc_ecalRecHitSequence)

ku_reduced_multi_ecalLocalRecoSequence   = cms.Sequence(ku_reduced_multi_ecalUncalibRecHitSequence*ku_reduced_multi_ecalRecHitSequence)

ku_spike_multi_ecalLocalRecoSequence   = cms.Sequence(ku_reduced_multi_ecalUncalibRecHitSequence*ku_spike_multi_ecalRecHitSequence)

ku_reduced_flipped_ecalLocalRecoSequence   = cms.Sequence(ku_reduced_multi_ecalUncalibRecHitSequence*ku_reduced_flipped_ecalRecHitSequence)

ku_spike_nomulti_ecalLocalRecoSequence   = cms.Sequence(ku_reduced_nomulti_ecalUncalibRecHitSequence*ku_spike_nomulti_ecalRecHitSequence)

ku_cc_gt_ecalLocalRecoSequence   = cms.Sequence(ku_cc_gt_ecalUncalibRecHitSequence*ku_cc_gt_ecalRecHitSequence)

ku_cc_native_ecalLocalRecoSequence = cms.Sequence(ku_cc_native_gt_ecalUncalibRecHitSequence*ku_cc_native_gt_ecalRecHitSequence)

#from RecoLocalCalo.EcalRecProducers.ecalDetailedTimeRecHit_cfi import *
#_phase2_timing_ecalRecHitSequence = cms.Sequence( ku_ecalRecHitSequence.copy() + ecalDetailedTimeRecHit )
#from Configuration.Eras.Modifier_phase2_timing_cff import phase2_timing
#phase2_timing.toReplaceWith( ku_ecalRecHitSequence, _phase2_timing_ecalRecHitSequence )
#
#_fastSim_ecalRecHitSequence = ecalRecHitSequence.copyAndExclude([ecalCompactTrigPrim,ecalTPSkim])
#_fastSim_ecalUncalibRecHitSequence = ecalUncalibRecHitSequence.copyAndExclude([ecalDetIdToBeRecovered])
#from Configuration.Eras.Modifier_fastSim_cff import fastSim
#fastSim.toReplaceWith(ecalRecHitSequence, _fastSim_ecalRecHitSequence)
#fastSim.toReplaceWith(ecalUncalibRecHitSequence, _fastSim_ecalUncalibRecHitSequence)


