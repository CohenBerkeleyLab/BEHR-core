#!/bin/bash -i
# This script will figure out what SP .mat files we need to 
# keep BEHR up to date. It assumes that the record is continuous up to the last
# file is continuous, and that it only needs to run forward.  For simplicity, it
# will run up to the last day that MODIS albedo data is available for, since
# that release should be the limiting factor.
#
# Needs the environmental variables "SPDIR" (directory where the OMI_SP .mat
# files are kept), "MODDIR" (root directory for MODIS files - it will assume
# that the MCD43C3 product is in the MCD43C3 subdirectory within there), and
# MATRUNDIR (directory where it will generate files related to running MATLAB)
# Also relies on alias startmatlab being defined to open a terminal window
# Matlab session (no GUI).
#
# Should be set to run in crontab about a day after all the downloading is
# handled, that way we can be pretty sure all the data is in.

source ~/.bashrc

# Debugging level, set higher to print more information
DEBUG=1
# How many days back in time to look for OMI_SP .mat files. Must be < 0
# as it is tested against offset back in time.
stopoffset=-90

if [[ -z $SPDIR ]]
then
    echo "ERROR run_behr_main.sh: env. var. SPDIR unset"
    automessage.sh "ERROR run_behr_main.sh" "Env. variable SPDIR unset"
    exit 1
elif [[ -z $BEHRDIR ]]
then
    echo "ERROR run_behr_main.sh: env. var. BEHRDIR unset"
    automessage.sh "ERROR run_behr_main.sh" "Env. variable BEHRDIR unset"
    exit 1
elif [[ -z $AMFTOOLSDIR ]]
then
    echo "ERROR run_behr_main.sh: env. var. AMFTOOLSDIR unset"
    automessage.sh "ERROR run_behr_main.sh" "Env. variable AMFTOOLSDIR unset"
    exit 1
elif [[ -z $NO2PROFDIR ]]
then
    echo "ERROR run_behr_main.sh: env. var. NO2PROFDIR unset"
    automessage.sh "ERROR run_behr_main.sh" "Env. variable NO2PROFDIR unset"
    exit 1
elif [[ -z $MATRUNDIR ]]
then
    echo "ERROR run_read_omno2.sh: env. var. MATRUNDIR unset"
    automessage.sh "ERROR run_read_omno2" "Env. variable MATRUNDIR unset"
    exit 1
elif [[ $DEBUG -gt 0 ]]
then
    echo -e "SPDIR=${SPDIR}\nMODDIR=${MODDIR}\nMATRUNDIR=${MATRUNDIR}"
fi

# Check that the file server directories exist.
if [[ ! -d $SPDIR ]]
then
    echo "ERROR run_behr_main.sh: $SPDIR does not exist. Is the file server mounted?"
    automessage.sh "ERROR run_behr_main.sh" "$SPDIR does not exist. Is the file server mounted?"
    exit 1
fi
if [[ ! -d $BEHRDIR ]]
then
    echo "ERROR run_behr_main.sh: $BEHRDIR does not exist. Is the file server mounted?"
    automessage.sh "ERROR run_behr_main.sh" "$BEHRDIR does not exist. Is the file server mounted?"
    exit 1
fi
if [[ ! -d $NO2PROFDIR ]]
then
    echo "ERROR run_behr_main.sh: $NO2PROFDIR does not exist. Is the file server mounted?"
    automessage.sh "ERROR run_behr_main.sh" "$NO2PROFDIR does not exist. Is the file server mounted?"
    exit 1
fi

# Find the last existing OMI_BEHR_YYYYMMDD.mat file
offset=0
foundit=false
while true
do
    startdate=$(date -d "${offset} days" +'%Y%m%d')
    testfile="${BEHRDIR}/OMI_BEHR_*${startdate}.mat"

    if [[ $DEBUG -gt 1 ]]; then echo "Checking for $testfile"; fi

    for f in $testfile
    do
        if [[ -f $f ]]
        then
            if [[ $DEBUG -gt 0 ]]; then echo "Found $f"; fi
            offset=$((offset+1))
            startdate=$(date -d "${offset} days" +'%Y-%m-%d')
            foundit=true
            break
        fi
    done
    if $foundit
    then
        break
    elif [[ $offset -lt $stopoffset ]]
    then
        automessage.sh "run_behr_main.m failed" "No OMI_BEHR files found within $((-stopoffset)) days."
        exit 1
    else
        offset=$((offset - 1))
    fi
done

# Find the last existing OMI_SP_YYYYMMDD.mat file - this will be the end date
offset=0
foundit=false
while true
do
    enddate=$(date -d "${offset} days" +'%Y%m%d')
    testfile="${SPDIR}/OMI_SP_*${enddate}.mat"

    if [[ $DEBUG -gt 1 ]]; then echo "Checking for $testfile"; fi

    for f in $testfile
    do  
        if [[ -f $f ]]
        then
            if [[ $DEBUG -gt 0 ]]; then echo "Found $f"; fi
            enddate=$(date -d "${offset} days" +'%Y-%m-%d')
            foundit=true
            break
        fi  
    done
    if $foundit
    then
        break
    elif [[ $offset -lt $stopoffset ]]
    then
        automessage.sh "run_behr_main.m failed" "No OMI_SP files found within $((-stopoffset)) days."
        exit 1
    else
        offset=$((offset - 1))
    fi  
done


echo "Start: $startdate End: $enddate"

# I'm cheating and using the fact that "onCluster" will cause the read_omno2 MATLAB function
# to read the various directories from global variables, this way I don't have to worry
# about them being changed when I pull the BEHR git repo.

echo "warning('off', 'all'); global DEBUG_LEVEL; DEBUG_LEVEL=1;" >${MATRUNDIR}/runscript_behr.m 
echo "addpath(genpath('${HOME}/Documents/MATLAB/BEHR'))" >> ${MATRUNDIR}/runscript_behr.m 
echo "addpath(genpath('${HOME}/Documents/MATLAB/Classes'))" >> ${MATRUNDIR}/runscript_behr.m 
echo "addpath(genpath('${HOME}/Documents/MATLAB/Utils'))" >> ${MATRUNDIR}/runscript_behr.m 
echo "global onCluster; onCluster = true;" >> ${MATRUNDIR}/runscript_behr.m
echo "global numThreads; numThreads = 1;" >> ${MATRUNDIR}/runscript_behr.m
echo "global sp_mat_dir; sp_mat_dir = '$SPDIR'" >> ${MATRUNDIR}/runscript_behr.m
echo "global behr_mat_dir; behr_mat_dir = '$BEHRDIR'" >> ${MATRUNDIR}/runscript_behr.m
echo "global amf_tools_path; amf_tools_path = '$AMFTOOLSDIR'" >> ${MATRUNDIR}/runscript_behr.m
echo "global no2_profile_path; no2_profile_path = '$NO2PROFDIR'" >> ${MATRUNDIR}/runscript_behr.m

echo "BEHR_main('${startdate}', '${enddate}'); exit(0)" >> ${MATRUNDIR}/runscript_behr.m

startmatlab -r "run('${MATRUNDIR}/runscript_behr.m')" > "${MATRUNDIR}/mat-behr.log"

matexit=$?
if [[ $matexit -ne 0 ]]
then
    automessage.sh "MATLAB: BEHR_main.m failed" -f "$MATRUNDIR/mat-behr.log"
else
    automessage.sh "MATLAB: BEHR_main.m succeeded" -f "$MATRUNDIR/mat-behr.log"
fi
exit 0
