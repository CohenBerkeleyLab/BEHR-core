#!/bin/bash
# Usage: ./order_omi.sh DATASET, where DATASET is either OMNO2 or OMPIXCOR.
#
# This script will need to be scheduled in the root crontab because it must be
# able to access the home directory for the nasa user in order to copy the download
# notice supplied by NASA. Basically, this script will order files from the last 120
# days, wait for the download notice to be delivered to /home/nasa/www, then copy
# it to either /home/nasa/OMNO2_dl_notices or /home/nasa/OMPIXCOR_dl_notices as a
# record (depending on which product it was told to request), then copy it to the
# directory specified by the environmental variable OMNO2DIR or OMPIXCORDIR (again,
# depending on which product it was told to request). These environmental variables
# need to be defined in the root's .bashrc file (in /root) for this script. (They
# also need to be set in the user josh's .bashrc for the get_omi.sh script that is
# called after this one.) Both of those variables should point to their respective
# "download_staging" folder on the Synology file server.

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




cd /home/nasa/www
# Remove any old notices
rm DN*.notify

# Construct the URL we will use to request the files we need. See
# http://disc.sci.gsfc.nasa.gov/additional/scienceTeam/s4pa_mri.shtml for
# information on the format. Note that %20 is used to represent a space in html.
syr=$(date -d "-120 days" +'%Y')
smn=$(date -d "-120 days" +'%m')
sday=$(date -d "-120 days" +'%d')

eyr=$(date +'%Y')
emn=$(date +'%m')
eday=$(date +'%d')

username="jlaughner"

#url="http://aurapar2.gesdisc.eosdis.nasa.gov/daac-bin/s4pa/s4pa_m2m_cgi.pl?user=${username}&dataset=OMNO2&version=003&startTime=${syr}-${smn}-${sday}%2000:00:00&endTime=${eyr}-${emn}-${eday}%2000:00:00"
#url="https://aurapar2.gesdisc.eosdis.nasa.gov/daac-bin/s4pa/s4pa_m2m_cgi.pl?user=${username}&dataset=OMNO2&version=003&startTime=${syr}-${smn}-${sday}%2000:00:00&endTime=${eyr}-${emn}-${eday}%2000:00:00"
url="https://aura.gesdisc.eosdis.nasa.gov/daac-bin/s4pa/s4pa_m2m_cgi.pl?user=${username}&dataset=${dataset}&version=003&startTime=${syr}-${smn}-${sday}%2000:00:00&endTime=${eyr}-${emn}-${eday}%2000:00:00"
echo $url

status=$(wget -O - "$url")
echo "Status = $status"

# The output of this wget command should contain the status "success"; if not,
# there's been a problem
if [[ ! $status =~ .*success.* ]]
then
    # send email about problem
    echo "Order failed. Contacting admin"
    msg="Order of $dataset failed. Message returned was:\n${status}"
    /usr/local/bin/automessage.sh "$dataset order failure" "$msg"
    exit 1
else
    echo "Order succeeded. Waiting for DN"
fi

timereq=$(date)

# Check for the delivery notification every 60 seconds
while true
do
    sleep 60

    # ls returns an exit code of 0 if a file matching DN*.notify is found.
    # $? retrieves the last exit code.
    ls DN*.notify >& /dev/null
    if [[ $? -eq 0 ]]
    then
        echo "DN received. Copying to necessary locations"
        new_fname="DN_$(date +'%Y%m%d')"
        cp DN*.notify ../${dataset}_dl_notices/${new_fname}
        cp DN*.notify ${staging_dir}/${new_fname}
        rm DN*.notify
        break
    fi

    # If nothing's come in for an hour, stop trying and email the admin
    timenow=$(date +'%H%M')
    if [[ $timenow -gt $(date -d "$timereq +1 hour" +'%H%M') ]]
    then
        automessage.sh "No delivery notice - $dataset" "$dataset DN requested at $timereq has not been delivered as of $(date)"
        exit 1
    fi

    echo "Still waiting for DN to be received"
done

exit 0
