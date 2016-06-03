#!/bin/bash

# Sends an email to adminstrator with the subject of the first argument
# and the message of the second
# Newlines (\n) do work, at least in the message

if [[ $# -lt 2 ]]
then
    echo "Usage: automessage.sh <subject> <message> OR automessage.sh <subject> -f <messagefile>"
    exit 1
fi

# user to contact
adminuser="jlaughner@berkeley.edu"

subject=''
msg=''
fromfile=0

while [[ $# > 0 ]]
do
key="$1"
    case $key in
        "-f")
        fromfile=1
        shift # shift the input arguments left by one
        ;;  
        *) # if not a flag, assume it's the subject, then message
        if [[ -z $subject ]]
        then
            subject=$key
        elif [[ -z $msg ]]
        then
            msg=$key
        else
            echo "Usage: automessage.sh <subject> <message> OR automessage.sh <subject> -f <messagefile>"
        fi
        shift
        ;;  
    esac
done


# use package sendEmail to send the mail via a gmail account set up just for this.
# You'll need to copy .gmailcredentials into your home directory and take ownership
# of it in order for this to work.

# Read in the credential, put them in an array such that ra[0] is the username
# and ra[1] the pw. This will keep the credentials as secure as I could (they
# should never show up as plaintext in the output to ps aux)
r=$(cat ~/.gmailcredentials)
ra=($r)

if [[ $fromfile -gt 0 ]]
then
    sendEmail -o tls=yes -f ${ra[0]} -t $adminuser -s smtp.gmail.com:587 -xu ${ra[0]} -xp ${ra[1]} -u $subject -o message-file=$msg
else
    sendEmail -o tls=yes -f ${ra[0]} -t $adminuser -s smtp.gmail.com:587 -xu ${ra[0]} -xp ${ra[1]} -u $subject -m $msg
fi
