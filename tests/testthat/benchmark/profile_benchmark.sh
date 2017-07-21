#!/bin/bash

commit_id=$(git rev-parse HEAD | perl -pe 's/(........).*/\1/')

CPUPROFILE="profile_${commit_id}.out" ./benchmark.sh
