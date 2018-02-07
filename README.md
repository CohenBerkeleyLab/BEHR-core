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
