# GammaResTool

CMSSW install for Timing:

cmsrel CMSSW_13_0_7 ( or current CMSSW production version )
cd CMSSW_13_0_7/src 
cmsenv 
git cms-init

( make sure the analysis package is cloned into its own folder inside scr/  example( timing/timing/…  ) )
( yes ..  its required by CMSSW )
( initial scram b -j8 must be in src/ ?  after maybe scram-ed  from plugins/ or test/ …. )


mkdir GammaResTool
cd GammaResTool/
git clone https://github.com/jking79/GammaResTool.git
cd ../
scram b -j 8



