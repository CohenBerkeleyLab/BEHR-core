#!/bin/bash

_mydir=$(dirname $0)
_mydir=$(cd $_mydir; pwd -P)
clone_dir=$(dirname $_mydir) # by default, put the new repositories in the folder containing this one

# *************************************** #
# Setup variables needed in sub functions #
# *************************************** #

fork="CohenBerkeleyLab"

github_repos=("${fork}/BEHR-core.git" \
    "${fork}/BEHR-core-utils.git" \
    "${fork}/Matlab-Gen-Utils.git" \
    "${fork}/MatlabPythonInterface.git" )


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

print_usage() {
    cat <<EOU
Usage: $(basename $0) [--https]|[--ssh]|[--ssh=<alias>]
    Will clone the additional repositories that the BEHR code depends on. Code
    from GitHub can be cloned using either HTTPS or SSH protocols. Cloning with
    HTTPS can be done by anyone, but requires a username and password to push
    changes back. This is the default option. Cloning with SSH requires that you
    have a GitHub account that has access to the BEHR repositories configured 
    with an SSH key.

     Options:
        --https : force this to clone from GitHub using HTTPS protocol. 
        --ssh : force this to clone from GitHub using SSH protocol.
        --ssh=<alias> : clone using SSH protocol, substituting <alias> for 
            "github.com" in the URL. The default URL is "git@github.com:/${fork}/<repo>",
            but if you have an SSH alias setup that specifies which SSH key to use
            to access GitHub, you can specify that here. For example, if I had
            the following in my ~/.ssh/config file:

                Host lab-github
                HostName github.com
                User git
                IdentityFile ~/.ssh/id_rsa_lab_github
                IdentitiesOnly yes

            then calling

                ./$(basename $0) --ssh=lab-github

            would clone the repositories with "git@lab-github:/${fork}/..."

    The repositories will be cloned to:
EOU

    for repo in ${github_repos[@]}; do
        reponame=$(basename $repo)
        reponame=${reponame%.*}
        repopath="$clone_dir/$reponame"
        echo "        $reponame --> $repopath"
    done

    exit 0
}

github_setup () {
    cd "$clone_dir"

    for repo in ${github_repos[@]}; do
        git clone ${1}${repo}
    done
}


# *************** #
# Parse arguments #
# *************** #

force_ssh=false
force_https=false
ssh_alias="github.com"

for arg in "$@"; do
    case $arg in
        --ssh)
            force_ssh=true
            force_https=false
            ;;
        --https)
            force_ssh=false
            force_https=true
            ;;
        --ssh=*)
            ssh_alias="${arg#*=}"
            ;;
        *)
            # no need to explicitly specify help options, will print usage
            # whenever an unrecognized argument is given, including -h, --help
            print_usage
            ;;
    esac
    shift
done
            

# Is Git installed on the system?
which git > /dev/null
isgit=$?
if [[ $isgit != 0 ]]; then
    manual_setup
    exit 0
fi

# Were we cloned with SSH or HTTPS?
# Allow user to force using HTTPS or SSH with the options '--https' or '--ssh'
https_stem="https://github.com/"
ssh_stem="git@${ssh_alias}:"

if $force_https; then
    stem="$https_stem"
elif $force_ssh; then
    stem="$ssh_stem"
else
    remote=($(git remote -v | grep origin.*fetch))
    if [[ ${remote[1]} =~ '^https' ]]; then
        stem="$https_stem"
    else
        stem="$ssh_stem"
    fi
fi

github_setup "$stem"

