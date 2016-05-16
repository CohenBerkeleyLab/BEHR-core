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

DEBUG=1

if [[ -z $SPDIR ]]
then
    echo "ERROR run_read_omno2.sh: env. var. SPDIR unset"
    automessage.sh "ERROR run_read_omno2.sh" "Env. variable SPDIR unset"
    exit 1
elif [[ -z $MODDIR ]]
then
    echo "ERROR run_read_omno2.sh: env. var. MODDIR unset"
    automessage.sh "ERROR run_read_omno2.sh" "Env. variable MODDIR unset"
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

# Check that the SPDIR and MODDIR folders exist, if not, the file server
# likely isn't mounted
if [[ ! -d $SPDIR ]];
then
    echo "ERROR run_read_omno2.sh: $SPDIR is not a directory. Is the file server mounted?"
    automessage.sh "ERROR run_read_omno2.sh" "$SPDIR does not exist"
    exit 1
fi
if [[ ! -d $MODDIR ]];
then
    echo "ERROR run_read_omno2.sh: $MODDIR is not a directory. Is the file server mounted?"
    automessage.sh "ERROR run_read_omno2.sh" "$MODDIR does not exist"
    exit 1
fi


# Find the last existing OMI_SP_YYYYMMDD.mat file
offset=0
while true
do
    startdate=$(date -d "${offset} days" +'%Y%m%d')
    testfile="${SPDIR}/OMI_SP_${startdate}.mat"

    if [[ $DEBUG -gt 1 ]]; then echo "Checking for $testfile"; fi

    if [[ -f $testfile ]]
    then
        if [[ $DEBUG -gt 0 ]]; then echo "Found $testfile"; fi
        offset=$((offset+1))
        startdate=$(date -d "${offset} days" +'%Y-%m-%d')
        break
    else
        offset=$((offset - 1))
    fi
done

# Find the last MCD43C3 file
offset=0
while true
do
    y=$(date -d "${offset} days" +'%Y')
    doy=$(date -d "${offset} days" +'%j')
    fpat="MCD43C3.A${y}${doy}.*.hdf"
    testfile="${MODDIR}/MCD43C3/${y}/${fpat}"

    if [[ $DEBUG -gt 1 ]]; then echo "Checking for $testfile"; fi

    ls $testfile >& /dev/null
    if [[ $? -eq 0 ]]
    then
        if [[ $DEBUG -gt 0 ]]; then echo "Found $testfile"; fi
        enddate=$(date -d "${offset} days" +'%Y-%m-%d')
        break
    else
        offset=$((offset-1))
    fi
done

echo "Start: $startdate End: $enddate"

# I'm cheating and using the fact that "onCluster" will cause the read_omno2 MATLAB function
# to read the various directories from global variables, this way I don't have to worry
# about them being changed when I pull the BEHR git repo.

echo "warning('off', 'all'); global DEBUG_LEVEL; DEBUG_LEVEL=1;" >${MATRUNDIR}/runscript.m 
echo "addpath(genpath('${HOME}/Documents/MATLAB/BEHR'))" >> ${MATRUNDIR}/runscript.m 
echo "addpath(genpath('${HOME}/Documents/MATLAB/Classes'))" >> ${MATRUNDIR}/runscript.m 
echo "addpath(genpath('${HOME}/Documents/MATLAB/Utils'))" >> ${MATRUNDIR}/runscript.m 
echo "global onCluster; onCluster = true;" >> ${MATRUNDIR}/runscript.m
echo "global numThreads; numThreads = 1;" >> ${MATRUNDIR}/runscript.m
echo "global sp_mat_dir; sp_mat_dir = '/mnt/sat/SAT/BEHR/SP_Files_2014'" >> ${MATRUNDIR}/runscript.m
echo "global omi_he5_dir; omi_he5_dir = '/mnt/sat/SAT/OMI/OMNO2'" >> ${MATRUNDIR}/runscript.m
echo "global modis_myd06_dir; modis_myd06_dir = '/mnt/sat/SAT/MODIS/MYD06_L2'" >> ${MATRUNDIR}/runscript.m
echo "global modis_mcd43_dir; modis_mcd43_dir = '/mnt/sat/SAT/MODIS/MCD43C3'" >> ${MATRUNDIR}/runscript.m
echo "global globe_dir; globe_dir = '/mnt/sat/SAT/BEHR/GLOBE_Database'" >> ${MATRUNDIR}/runscript.m

echo "read_omno2_v_aug2012('${startdate}', '${enddate}'); exit(0)" >> ${MATRUNDIR}/runscript.m

startmatlab -r "run('${MATRUNDIR}/runscript.m')" > "${MATRUNDIR}/mat.log"

matexit=$?
if [[ $matexit -ne 0 ]]
then
    automessage.sh "MATLAB: read_omno2_v_aug2012.m failed" -f "$MATRUNDIR/mat.log"
else
    automessage.sh "MATLAB: read_omno2_v_aug2012.m succeeded" -f "$MATRUNDIR/mat.log"
fi
exit 0
