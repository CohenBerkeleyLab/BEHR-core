# BEHR Core
## Core code for the Berkeley High Resolution OMI NO2 Retrieval

This repository contains the core algorithm for the BEHR retrieval, written in Matlab.
It depends on several other repositories from the CohenBerkeleyLab GitHub. These can
be obtained automatically by using the `setup.sh` shell script provided in the root
directory. This will clone the other repositories as sibling folders to this one.
(This is a bash script which would be run on Mac or Linux machines in the Terminal
by changing to the BEHR-core directory and running the command `./setup.sh`. It should
also work in the Git Bash shell on Windows, but that has not been tested.)

Once those are cloned, you will also need to run the Matlab function `BEHR_initial_setup.m`
in the BEHR-core-utils repository. That will generate a `behr_paths.m` class in 
BEHR-core-utils/Utils/Constants, which contains all the necessary code and data paths for
BEHR to run. You will likely need to edit this file to point to the proper directories for
OMI SP, OMPIXCOR, MODIS, GLOBE, and WRF data. Once done, running `behr_paths.ValidatePaths()`
in Matlab should return True.

## Fair use policy
By using this code in your research you agree to the following terms in addition to the terms of reuse given in the license:
  1. Any publication that uses this code will include a citation to this repository. Cite as: 
     * Laughner, J.L., Zhu, Q., and Cohen, R.C. (2018, May 18). CohenBerkeleyLab/BEHR-core: Version 3.0B (Version v3.0B). Zenodo. http://doi.org/10.5281/zenodo.1249671
  1. Only the master branch is considered stable. All other branches are under development, subject to change,
     and are not recommended for scientific use.
  1. We do our best to ensure that the master branch is bug-free and scientifically sound. However, we cannot test all
     possible use cases. The user is ultimately responsible for ensuring that any results obtained using this code are
     scientifically accurate.
  1. If your research uses a branch other than master,
     please notify us as soon as possible that you intend to publish a manuscript using unpublished features of
     this code. If the publication competes with one of our own, we may ask that you delay publishing until we
     submit our manuscript.
