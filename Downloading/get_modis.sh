#!/bin/bash
#
# This function will contact the MODIS LAADS Web FTP server at ftp://ladsweb.nascom.nasa.gov/ and
# identify MYD06 and MCD43C1 files from the past three months that have not already been downloaded
# and retrieve them, putting them in the proper directory.  This should be run weekly (using cron).
#
# The MODIS root directory on the file server should be defined in the env. variable MODDIR. (This
# directory should contain subfolders for the various MODIS products.)
#
# Josh Laughner <joshlaugh5@gmail.com> 7 Aug 2015
# The automodis.py script must be there as well for this to work

DEBUG=1

source ~/.bashrc

fsource=${BASH_SOURCE[0]}
if [[ $(uname) == "Darwin" ]]; then
    thisfile=$(readlink $fsource)
    thisdir="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    scriptdir="${thisdir}/${thisfile}"
else
    thisfile=$(readlink -f $fsource)
    scriptdir=$(dirname $thisfile)
fi
# Define root MODIS directory on the file server from the env. var MODDIR.
# if empty, error.
if [[ $MODDIR == '' ]]
then
    echo "ERROR - GET_MODIS.SH: Env. variable MODDIR undefined"
    exit 1
fi

# Remote directories have the pattern ftp://ladsweb.nascom.nasa.gov/allData/
# <collection>/<product>/<year>/<day-of-year>/<data-files>. This is the root
# remote dir
REMOTEDIR="ftp://ladsweb.nascom.nasa.gov/allData"

## =============== ##
##     MCD43C1     ##
## =============== ##


# Retrieve the MCD43C1 data. There's only one file produced every 8 days, so we
# only need to look for the lone .hdf file in each day's folder.
for offset in `seq -90 -1`
do
    if [[ $(uname) == "Darwin" ]]
    then
        y=$(date -v ${offset}d +'%Y')
        doy=$(date -v ${offset}d +'%j')
    else
        y=$(date -d "${offset} days" +'%Y')
        doy=$(date -d "${offset} days" +'%j')
    fi
    fullLocalDir="${MODDIR}/MCD43C1/${y}/"
    fullRemDir="${REMOTEDIR}/6/MCD43C1/${y}/${doy}/"

    if [[ ! -d $fullLocalDir ]]
    then
        echo "Creating $fullLocalDir"
        mkdir -p $fullLocalDir
    fi

    cd $fullLocalDir
    if [[ $DEBUG -gt 0 ]]; then echo "Current directory is $(pwd)"; fi

    # This will retrieve a list of files in the remote directory to the .listing
    # file, which we can then parse to see if the remote files are present
    # locally. If wget returns a non-zero exit status, most likely the remote
    # directory did not exist b/c either that day has not been processed yet, or
    # it's not one of the days that MCD43C1 exists for. Yes, we could define
    # what days those should be, but this is safe against changes to the
    # processing schedule.

    if [[ $DEBUG -gt 0 ]]; then echo "Getting MCD43C1 file list for $y $doy"; fi

    wget -q --spider --no-remove-listing $fullRemDir
#    wget --spider --no-remove-listing $fullRemDir
    if [[ $? -ne 0 ]]; then continue; fi

    fname=$(grep '.hdf' .listing | awk '{print $9'})
    fname=$(echo $fname | tr -d '\r') # remove random carriage returns
    if [[ $DEBUG -gt 0 ]]; then echo "file = $fname"; fi
    # Only download if the file does not exist locally. If downloading, do not
    # include the header (ftp address) or any directories - just put the file
    # here.
    if [[ ! -f $fname ]]; 
    then 
        echo "Retrieving $fname"
        wget -q -nH -nd $fullRemDir$fname
    fi
done

## =============== ##
##      MYD06      ##
## =============== ##

# Cloud data requires a little more work. There are almost 300 granules per day;
# obviously we don't want to download all of them if we don't need them (and we
# don't). It's much more convinient to select which ones based on their lat/lon
# coordinates, so we will interface with the MODIS web services using SOAPpy.
# This will call the automodis.py Python script which will do this (the
# automodis.py script must be in the same folder as this one to work properly).
# This python script will return a list of URLs that fit the specified lon/lat
# criteria; this will then check if those files have already been downloaded and
# if not will call wget to retrieve the file in the proper folder.

if [[ $(uname) == "Darwin" ]]
then
    # Mac OSX implementation of date
    startdate=$(date -v -90d +'%Y-%m-%d 00:00:00')
    enddate=$(date +'%Y-%m-%d 23:59:59')
else
    # Ubuntu (hopefully standard Linux type) implementation
    startdate=$(date -d "-90 days" +'%Y-%m-%d 00:00:00')
    enddate=$(date +'%Y-%m-%d 23:59:59')
fi

# lon/lat bounds are specified as default values in the Python script
# only retrieve granules with daytime information

echo "Getting file list..."
#filelist=$(python ${scriptdir}/automodis.py --products MYD06_L2 --startTime "$startdate" --endTime "$enddate" --dayNightBoth 'DB')
python ${scriptdir}/automodis.py --products MYD06_L2 --startTime "$startdate" --endTime "$enddate" --dayNightBoth 'DB'
echo "Done."
filelist=$(cat $MATRUNDIR/modis_urls.txt)
for f in $filelist
do
    fname=$(basename $f)
    y=${fname:10:4}
    doy=${fname:14:3}

    fullLocalDir="${MODDIR}/MYD06_L2/${y}"

    if [[ ! -d $fullLocalDir ]]
    then
        mkdir -p $fullLocalDir
    fi

    cd $fullLocalDir

    if [[ ! -f $fname ]]
    then
        echo "Retrieving $fname"
        wget -q -nH -nd $f
    fi
done

exit 0
