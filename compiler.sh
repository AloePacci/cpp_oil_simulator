#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo "Running $SCRIPT_DIR/compiler.sh"

cd $SCRIPT_DIR
if [ ! -f ./basic_sim/build ]; then
    cd basic_sim;
    mkdir build;
    cd ..;
fi
cd basic_sim/build

set -e
if [ "${1}" = "clean" ] ; then
    rm -r *
fi;

cmake ../.. -Wno-dev    
make

while true; do
    read -p "Do you wish to execute this program? [y/n]" yn
    case $yn in
        [Yy]* ) ./basic_sim/BasicSimulator; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

# ./basic_sim/build/basic_sim/BasicSimulator;
