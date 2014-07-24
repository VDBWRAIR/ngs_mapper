#!/bin/bash

# Where python will be installed if needed
if [[ -z "$PYTHON_INSTALL_PREFIX" ]]
then
    PYTHON_INSTALL_PREFIX=/usr/local
fi

function ensure_pkg {
    if ! rpm -qa | grep -q "^${1}"
    then
        echo "Installing $1"
        yum install -y $1
    fi
}

function ensure_python {
    prefix=$1

    if [[ $(python -V 2>&1) =~ 2.7.[3456789] ]]
    then
        echo "Required version of python is already installed" 
        return 0
    else
        echo "Installing required version of python"
        ver='2.7.3'
        pushd $(mktemp -d) && wget --no-check-certificate https://www.python.org/ftp/python/${ver}/Python-${ver}.tgz -O- | tar xzf - && pushd Python-${ver} && ./configure --prefix $prefix && make && make install && popd && popd
        return $?
    fi
}

function install_pipeline {
    from=$1
    to=$2

    if [[ $from =~ .git$ ]] || [ -d ${from}/.git ]
    then
        echo "Installing pipeline into $to"
        cd $to
        pth=$(basename $from)
        pth=${pth/.git//}
        # Make sure no virtualenvironments are activated
        deactivate 1>/dev/null 2>&1
        if [ -d $pth ]
        then
            echo "Pipeline seems to already be installed at ${to}/${pth} so we will try to update"
            # Update and uninstall
            cd $pth
            # I guess stash will work for now in case there were changes(don't want to destroy something)
            git stash
            git pull
            ./uninstall.sh
        else
            git clone $from
            cd $pth
        fi

        # Ensure python is installed
        # I guess for now we will just hardcode it to go to /usr/local(which is the default anyways)
        # It only will be installed if we cannot detect Python 2.7.{3..9} anyways so if they already have python then no harm right?
        ensure_python ${PYTHON_INSTALL_PREFIX}
        # Ensure virtualenv is installed
        install_virtualenv ${to}/${pth} ${PYTHON_INSTALL_PREFIX}

        # Install
        ./install.sh
        _INSTALL_PATH=${to}/${pth}
        return $?
    else
        echo "$from is not a git repo"
        return 1
    fi
}

function install_virtualenv {
    pipeline_install_location=$1
    python_install_prefix=$2
    ${python_install_prefix}/bin/python -c "virtualenv" 2>/dev/null || pushd ${pipeline_install_location}/dependencies/virtualenv* && ${python_install_prefix}/bin/python setup.py install && popd
}

function install_system_packages {
    pkg_lst=$1
    echo "Ensuring all system package dependencies are installed"
    for pkg in ${pkg_lst[@]}
    do
        if ! ensure_pkg $pkg
        then
            return 1
        fi
    done
    return 0
}

function main {
    repo_pth=$1
    install_pth=$2
    pkg_lst=( libpng{,-devel} ncurses{,-devel} ImageMagick zlib{,-devel} freetype{,-devel} readline{,-devel} openssl{,-devel} gcc-c++ git )

    if [ ! $(id -u) -eq 0 ]
    then
        echo "Please run as root"
        return 1
    fi

    if [ -z $repo_pth ]
    then
        echo "Please specify the location of the miseqpipeline git repository"
        return 1
    fi

    install_system_packages $pkg_lst && install_pipeline $repo_pth $install_pth
    RET=$?

    if [ $RET -eq 0 ]
    then
        echo "The miseqpipeline is now installed into $install_pth. You can activate the pipeline by issuing the following command:"
        echo ". ${_INSTALL_PATH}/setup"
        return $RET
    else
        echo "There was an error during the installation."
        return $RET
    fi

}

main $@
