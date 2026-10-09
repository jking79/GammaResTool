# GammaResTool

CMS ECAL timing reconstruction, calibration, smearing, and analysis tools.

This repository contains tools developed for CMS ECAL timing studies, including
timing reconstruction, calibration, and smearing utilities.

## CMSSW Installation

### CMSSW Versions

The appropriate CMSSW release depends on the analysis and data-taking period.

Previously used configurations:

- LLPana / KUCMSNtuple: `CMSSW_13_3_0` / `CMSSW_13_3_3`
- EGammaRes timing studies: `CMSSW_14_0_11`
- Original timing development: `CMSSW_13_0_0` / `CMSSW_13_0_7`

Use the CMSSW release recommended for the analysis or production campaign.

### 1. Initialize the CMS Environment

Log in to an appropriate CMS computing environment (e.g., Fermilab LPC).

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
```

Set the SCRAM architecture appropriate for the CMSSW release and operating system.

For example, on an AlmaLinux 9 system using a CMSSW release built with GCC 13:

```bash
export SCRAM_ARCH=el9_amd64_gcc13
```

The appropriate architecture may differ for older CMSSW releases.

### 2. Create the CMSSW Release

Set the desired CMSSW version.

Example:

```bash
export CMSSW_RELEASE=CMSSW_14_0_11
```

Create and initialize the release area:

```bash
cmsrel ${CMSSW_RELEASE}
cd ${CMSSW_RELEASE}/src
cmsenv
```

Initialize the CMSSW Git repository:

```bash
git cms-init
```

The `git cms-init` step is required when working with or modifying CMSSW
packages. It is not required simply to run an unmodified CMSSW release.

### 3. Install GammaResTool

CMSSW analysis packages should be placed inside a package directory
under the CMSSW `src/` directory.

For GammaResTool, use the following directory structure:

```text
CMSSW_X_Y_Z/
└── src/
    └── GammaResTool/
        └── GammaResTool/
            ├── plugins/
            ├── test/
            ├── macros/
            └── ...
```

From the CMSSW `src/` directory:

```bash
mkdir -p GammaResTool
cd GammaResTool

git clone https://github.com/jking79/GammaResTool.git

cd ..
```

This produces the required nested package structure:

```text
src/GammaResTool/GammaResTool/
```

The nested directory structure follows the CMSSW package convention
`Subsystem/Package`.

### 4. Compile

From the CMSSW `src/` directory:

```bash
scram b -j 8
```

The initial build should be performed from the `src/` directory.

Subsequent builds may also be performed from the `src/` directory.
Building from individual package directories may work, but building
from `src/` is the recommended procedure.

### 5. Verify the Environment

Check the active CMSSW environment:

```bash
echo $CMSSW_VERSION
echo $CMSSW_BASE
echo $SCRAM_ARCH
```

The output should correspond to the selected CMSSW release and architecture.

To reactivate an existing CMSSW installation in a new shell:

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh

cd /path/to/CMSSW_X_Y_Z/src
cmsenv
```

---

## Historical ECAL Timing Development Configurations

The following configurations were used during previous ECAL timing
development efforts.

These instructions are retained for reference and reproducibility.

They are not required for a standard GammaResTool installation unless
the corresponding modified CMSSW reconstruction packages are needed.

### Original ECAL Timing Development

Original development based on the `nminafra/cmssw` repository.

From an initialized CMSSW `src/` directory:

```bash
git remote add NM https://github.com/nminafra/cmssw.git
git fetch NM

git checkout ecalTiming_10_2_5_v1

git cms-addpkg RecoLocalCalo/EcalRecProducers
```

**Note:** This configuration corresponds to an older CMSSW development
environment and should not be assumed compatible with current releases.

### Timing Analysis Package

The standalone Timing analysis package was developed in:

https://github.com/jking79/Timing

Installation:

```bash
cd $CMSSW_BASE/src

mkdir -p Timing
cd Timing

git clone https://github.com/jking79/Timing.git

cd Timing
git checkout jwk_dev

cd $CMSSW_BASE/src

scram b -j 8
```

The `jwk_dev` branch contains the corresponding development configuration.

Compatibility with the selected CMSSW release should be verified before use.

### ECAL Timing CC Development (CMSSW 13.0.0)

This configuration was used for the ECAL timing reconstruction
development based on the CC implementation.

CMSSW release:

```text
CMSSW_13_0_0
```

Related development release:

```text
CMSSW_13_1_X_2023-04-14-1100
```

Initialize the CMSSW environment:

```bash
cmsrel CMSSW_13_0_0

cd CMSSW_13_0_0/src

cmsenv

git cms-init
```

Add the development repository:

```bash
git remote add jwk https://github.com/jking79/cmssw.git

git fetch jwk
```

Check out the timing development branch:

```bash
git checkout --track jwk/ecalTiming_cc_deployment_13_0_0
```

Add the required CMSSW packages:

```bash
git cms-addpkg RecoLocalCalo/EcalRecProducers

git cms-addpkg DataFormats/EcalRecHit
```

Build the modified CMSSW environment:

```bash
scram b -j 8
```

**Important:** This configuration modifies CMSSW reconstruction packages
and is specific to the corresponding development release.

Do not apply these modifications to newer CMSSW releases without
verifying compatibility.

---

## KUCMSTimeCalibration

The `KUCMSTimeCalibration` class provides ECAL rechit timing calibration
and smearing utilities.

Location:

```text
GammaResTool/macros/ecal_config/
```

Primary class:

```text
KUCMSTimeCalibration.hh
```

Associated files:

```text
ecal_config/
├── KUCMSEcalDetIDFunctions.hh
├── KUCMSHelperBaseClass.hh
├── KUCMSRootHelperBaseClass.hh
├── KUCMSTimeCalibration.hh
├── README.txt
├── UL2016_runlumi.txt
├── UL2017_runlumi.txt
├── UL2018_runlumi.txt
├── caliHistsTFile.root
├── caliRunConfig.txt
├── caliSmearConfig.txt
├── caliTTConfig.txt
├── fullinfo_detids_EB.txt
├── fullinfo_detids_EE.txt
├── fullinfo_v2_detids_EB.txt
├── fullinfo_v2_detids_EE.txt
├── howto.txt
├── reducedinfo_detids.txt
├── rhid_i12_list.txt
└── rhid_info_list.txt
```

Additional utility:

```text
fillKUCMSTimeCalibration.cpp
```

This utility is used to initialize calibration and smearing information.

### Initializing the Calibration Class

Declare the calibration class once during initialization:

```cpp
KUCMSTimeCalibration theCali;
```

### Applying ECAL Timing Calibration

Select the calibration tag.

Default:

```cpp
theCali.setTag("EG_EOY_MINI");
```

Retrieve the calibration correction for a rechit:

```cpp
calibration = theCali.getCalibration(rechitID, run);
```

Apply the correction by subtracting it from the reconstructed rechit time:

```cpp
calibratedTime = rechitTime - calibration;
```

Alternatively, the calibrated time can be obtained directly:

```cpp
calibratedTime = theCali.getCalibTime(
    rechitTime,
    rechitID,
    run
);
```

Both methods provide the calibrated rechit time.

### Applying ECAL Timing Smearing

Select the smearing configuration.

Default:

```cpp
theCali.setSmearTag("EG300202_DYF17");
```

Apply the smearing:

```cpp
smearedTime = theCali.getSmearedTime(
    rechitTime,
    rechitAmplitude
);
```

The smearing configuration determines the timing resolution model
applied to the rechit.

---

## Notes

- Use the CMSSW release appropriate for the target dataset and analysis.
- Verify SCRAM architecture compatibility before creating a release.
- The standard GammaResTool installation does not require checking out
  modified CMSSW reconstruction packages.
- Historical timing reconstruction modifications are retained above
  for reference.
- When using modified ECAL reconstruction packages, confirm that the
  implementation is compatible with the selected CMSSW release.
- Calibration and smearing tags should be selected according to the
  intended dataset and timing reconstruction configuration.



