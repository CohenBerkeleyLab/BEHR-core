#!/bin/bash

_mydir=$(dirname $0)
    
github_repos=('/volume1/share-sat/SAT/BEHR/MiscUtils.git' \
    '/volume1/share-sat/SAT/BEHR/BEHR_MatlabClasses_GitRepo.git' \
    '/volume2/share2/GROUP/GitRepos/MatlabPythonInterface.git')

berkeley_repos=('/volume1/share-sat/SAT/BEHR/MiscUtils.git' \
    '/volume1/share-sat/SAT/BEHR/BEHR_MatlabClasses_GitRepo.git' \
    '/volume2/share2/GROUP/GitRepos/MatlabPythonInterface.git')

manual_setup () {
    >&2 echo  "GIT does not appear to be installed on this computer.
If you are running this setup as a grad student on the BEHR project
at Berkeley, you should install GIT and properly clone BEHR from
the file server.

If you have obtained BEHR by downloading a release from GitHub, you'll
need to also download the following repositories:"

    for repo in ${github_repos[@]}; do
        >&2 echo "    $repo"
    done
}

github_setup () {
    # these need to change

    cd "$_mydir/.."

    for repo in ${github_repos[@]}; do
        git clone ${repo}
    done
}

berkeley_setup() {

    default_user=RCCohenLab
    echo "Enter the user name on the 128.32.208.13 server to use (empty will use $default_user)"
    read user

    if [[ -z $user ]]; then
        user=$default_user
    fi

    echo "You will need to enter your password for each rep (sorry, could not find a way around that)"
    read junk


    cd "$_mydir/.."

    for repo in ${berkeley_repos[@]}; do
        GIT_SSH_COMMAND="ssh -o PreferredAuthentications=password" git clone ${user}@128.32.208.13:${repo}
    done

}

# Is Git installed on the system?
isgit=$(which gitblah > /dev/null)
if [[ $isgit != 0 ]]; then
    manual_setup
    exit 0
fi

# Were we cloned from GitHub or the Berkeley file server?
origin_address="$(git remote get-url origin)"
if [[ $origin_address =~ 'github' ]]; then
    github_setup
    exit 0
else
    berkeley_setup
    exit 0
fi
