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

logfile="run-pub.log"

# Debugging level, set higher to print more information
DEBUG=1
# How many days back in time to look for OMI_SP .mat files. Must be < 0
# as it is tested against offset back in time.
stopoffset=-90

if [[ -z $BEHRDIR ]]
then
    echo "ERROR run_publishing.sh: env. var. BEHRDIR unset"
    automessage.sh "ERROR run_publishing.sh" "Env. variable BEHRDIR unset"
    exit 1
elif [[ -z $HDFSAVEDIR ]]
then
    echo "ERROR run_publishing.sh: env. var. HDFSAVEDIR unset"
    automessage.sh "ERROR run_publishing.sh" "Env. variable HDFSAVEDIR unset"
    exit 1
elif [[ -z $MATRUNDIR ]]
then
    echo "ERROR run_publishing.sh: env. var. MATRUNDIR unset"
    automessage.sh "ERROR run_publishing.sh" "Env. variable MATRUNDIR unset"
    exit 1
elif [[ $DEBUG -gt 0 ]]
then
    echo -e "SPDIR=${SPDIR}\nMODDIR=${MODDIR}\nMATRUNDIR=${MATRUNDIR}"
fi

# Check that the file server directories exist
if [[ ! -d $BEHRDIR ]]
then
    echo "ERROR run_publishing.sh: $BEHRDIR does not exist. Is the file server mounted?"
    automessage.sh "ERROR run_publishing.sh" "$BEHRDIR does not exist. Is the file server mounted?"
    exit 1
fi
if [[ ! -d $HDFSAVEDIR ]]
then
    echo "ERROR run_publishing.sh: $HDFSAVEDIR does not exist. Is the file server mounted?"
    automessage.sh "ERROR run_publishing.sh" "$HDFSAVEDIR does not exist. Is the file server mounted?"
    exit 1
fi
if [[ ! -d $MATRUNDIR ]]
then
    echo "ERROR run_publishing.sh: $MATRUNDIR does not exist."
    automessage.sh "ERROR run_publishing.sh" "$MATRUNDIR does not exist."
    exit 1
fi

# See if any matlab processes are running, if so, wait until they are done.
nproc=$(ps aux | grep -i matlab | wc -l)
safety=0
echo "Checking for running matlab processes" > $MATRUNDIR/$logfile
while [[ $nproc -gt 1 ]]
do
    echo "$(ps aux | grep -i matlab)"
    echo "$nproc other matlab processing active, waiting" >> $MATRUNDIR/$logfile
    sleep 1800 #wait half an hour
    nproc=$(ps aux | grep -i matlab | wc -l)
    safety=$((safety+1))
    if [[ $safety -gt 48 ]]; then
        automessage.sh "Problem with run_publishing.sh" "Waiting more than a day for other Matlab processes to finish. Please check on status of download computer."
    fi
done
echo "No matlab processes detected. Executing." >> $MATRUNDIR/$logfile

# Find the last existing HDF file
offset=0
foundit=false
while true
do
    startdate=$(date -d "${offset} days" +'%Y%m%d')
    testfile="${HDFSAVEDIR}/behr_hdf/OMI_BEHR_*${startdate}.hdf"

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
        automessage.sh "run_publishing.m failed" "No HDF files found within $((-stopoffset)) days"
        exit 1
    else
        offset=$((offset - 1))
    fi
done

# Find the last BEHR file - this will be the end date
# Find the last existing OMI_BEHR_YYYYMMDD.mat file
offset=0
foundit=false
while true
do
    enddate=$(date -d "${offset} days" +'%Y%m%d')
    testfile="${BEHRDIR}/OMI_BEHR_*${enddate}.mat"

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
        automessage.sh "run_publishing.m failed" "No OMI_BEHR files found within $((-stopoffset)) days."
        exit 1
    else
        offset=$((offset - 1))
    fi  
done

echo "Start: $startdate End: $enddate"

# I'm cheating and using the fact that "onCluster" will cause the read_omno2 MATLAB function
# to read the various directories from global variables, this way I don't have to worry
# about them being changed when I pull the BEHR git repo.

echo "warning('off', 'all'); global DEBUG_LEVEL; DEBUG_LEVEL=1;" >${MATRUNDIR}/runscript_pub.m 
echo "addpath(genpath('${HOME}/Documents/MATLAB/BEHR'))" >> ${MATRUNDIR}/runscript_pub.m 
echo "addpath(genpath('${HOME}/Documents/MATLAB/Classes'))" >> ${MATRUNDIR}/runscript_pub.m 
echo "addpath(genpath('${HOME}/Documents/MATLAB/Utils'))" >> ${MATRUNDIR}/runscript_pub.m 
echo "global onCluster; onCluster = true;" >> ${MATRUNDIR}/runscript_pub.m
echo "global numThreads; numThreads = 1;" >> ${MATRUNDIR}/runscript_pub.m
echo "global mat_file_dir; mat_file_dir = '$BEHRDIR'" >> ${MATRUNDIR}/runscript_pub.m
# Native HDF files
echo "global save_dir; save_dir = '${HDFSAVEDIR}/behr_hdf'" >> ${MATRUNDIR}/runscript_pub.m
echo "BEHR_publishing_v2('hdf','native',{},'${startdate}', '${enddate}'); " >> ${MATRUNDIR}/runscript_pub.m
# Native text files
echo "global save_dir; save_dir = '${HDFSAVEDIR}/behr_txt'" >> ${MATRUNDIR}/runscript_pub.m
echo "BEHR_publishing_v2('txt','native',{},'${startdate}', '${enddate}'); " >> ${MATRUNDIR}/runscript_pub.m
# Gridded HDF files
echo "global save_dir; save_dir = '${HDFSAVEDIR}/behr_regridded_hdf'" >> ${MATRUNDIR}/runscript_pub.m
echo "BEHR_publishing_v2('hdf','gridded',{},'${startdate}', '${enddate}');" >> ${MATRUNDIR}/runscript_pub.m
# Ensure that the website's record of the version number is up-to-date.
echo "update_behr_version_website; exit(0)" >> ${MATRUNDIR}/runscript_pub.m

startmatlab -r "run('${MATRUNDIR}/runscript_pub.m')" > "${MATRUNDIR}/mat-pub.log"

matexit=$?
if [[ $matexit -ne 0 ]]
then
    automessage.sh "MATLAB: run_publishing.m failed" -f "$MATRUNDIR/mat-pub.log"
else
    automessage.sh "MATLAB: run_publishing.m succeeded" -f "$MATRUNDIR/mat-pub.log"
fi
exit 0
