#!/bin/bash

if [ -z "$1" ]
then
    echo 'Usage: view_profile.sh profile'
    exit 1
fi

pprof --text /usr/bin/Rscript "$1" | head -20
