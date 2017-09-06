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
# does not exist before downloading it. The files will be downloaded
# to download_staging, then sorted into proper year/month folders.

# OMNO2DIR must be an environmental variable set pointing to the
# download_staging folder.  I made it and environmental variable so that this
# script can be moved from computer to computer w/o needing to be edited.

source ~/.bashrc

# This requires 1 input: the dataset being downloaded (OMNO2 or OMPIXCOR)
if [[ -z $1 ]]; then
    (>&2 echo "$0 requires one input: OMNO2 or OMPIXCOR")
    exit 1
else
    dataset="$1"
fi

# This determines which staging directory will be used
if [[ $dataset == OMNO2 ]]; then
    staging_dir="$OMNO2DIR"
elif [[ $dataset == OMPIXCOR ]]; then
    staging_dir="$OMPIXCORDIR"
fi

if [[ -z $staging_dir ]]; then
    (>&2 echo "Dataset $dataset does not have a staging directory set via an environmental variable")
    exit 1
elif [[ ! -d $staging_dir ]]; then
    (>&2 echo "Staging directory $staging_dir does not exist")
    exit 1
fi

cd "$staging_dir"

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
        automessage.sh "$(basename $0) failed" "After six hours, the download notice $dlf has not been found"
        exit 1
    else
        echo "download notice $dlf not found, waiting 60 sec..."
        safety=$((safety+1))
        sleep 60
    fi
done

# Once it exists, loop through the files listed in it.
# Test that they do not exist first, then download.  Note that this will be 
# running in the download_staging folder, so .. points to a directory containing 
# subfolders by year, each of which is organized into monthly subfolders.

# This will ensure that the cookies file is always fresh so that if something got
# messed up before it'll be cleared
rm -f "$HOME/.urs_cookies"
touch "$HOME/.urs_cookies"

# Since we will get different file names for different products, we need to find
# the first date with format YYYYmMMDD in the name (though we only need the year
# and month). I'll include the underscore b/c the data time is preceeded by an
# underscore in the name, while the processing time is preceeded by a dash.
regex='_[0-9][0-9][0-9][0-9]m[0-9][0-9]'

for f in $(cat $dlf)
do
    fname=$(basename $f)
    if [[ $fname =~ $regex ]]; then
        fdate=${BASH_REMATCH[0]}
        y=${fdate:1:4}
        m=${fdate:6:2}
    else
        (>&2 echo "Could not find data date in the filename $fname")
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
        wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies -nH -nd $f
        mv $fname ../${y}/${m}/
    fi
    
done

exit 0
