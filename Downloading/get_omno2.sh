#!/bin/bash
# Get OMNO2 data, once the machine request has succeeded and the download notice
# has been placed in ~/Documents/dl_notices. That will be handled by a root
# crontab, since the download notice will be put in /home/nasa and owned by user
# nasa.
#
# It may take up to a half an hour for the download notice to be present, so
# this should run about a minute after the root script that requests the
# download notice, and will check once per minute to see if the root script has
# moved the download notice to the download staging folder.
#
# Basically all this will need to do is move the download notice into the
# download_staging folder on the file server and check if the specified file
# does not exist before downloading it. It will also ignore files ending in
# .xml, since these are unnecessary metadata files. The files will be downloaded
# to download_staging, then sorted into proper year/month folders.

# OMNO2DIR must be an environmental variable set pointing to the
# download_staging folder.  I made it and environmental variable so that this
# script can be moved from computer to computer w/o needing to be edited.

source ~/.bashrc

if [[ -z $OMNO2DIR ]]
then
    echo "ERROR get_omno2.sh: OMNO2DIR is not defined"
    exit 1
fi

cd "$OMNO2DIR"

# The root script will name the download notices DL_YYYYMMDD, so this is the
# file name we need to look for.
dlf="DN_$(date +'%Y%m%d')"

safety=0
while true
do
    if [[ -f $dlf ]]
    then
        echo "Using download notice $dlf"
        break
    elif [[ $safety -gt 360 ]]
    then
        automessage.sh "get_omno2.sh failed" "After six hours, the download notice $dlf has not been found"
        exit 1
    else
        echo "download notice $dlf not found, waiting 60 sec..."
        safety=$((safety+1))
        sleep 60
    fi
done

# Once it exists, loop through the files listed in it. Skip any ending in XML
# (unneeded metadata) for the rest, test that they do not exist first, then
# download.  Note that this will be running in the download_staging folder, so
# .. points to a directory containing subfolders by year, each of which is
# organized into monthly subfolders.

for f in $(cat $dlf)
do
    fname=$(basename $f)
    y=${fname:18:4}
    m=${fname:23:2}

    if [[ $fname == *.xml ]]
    then
        continue
    fi

    
    if [[ ! -d ../${y}/${m} ]]
    then
        mkdir -p ../${y}/${m}
    fi
    if [[ -f ../${y}/${m}/${fname} ]]
    then
        echo "File $fname exists"
        continue
    else
        wget -q -nH -nd $f
        mv $fname ../${y}/${m}/
    fi
    
done

exit 0
