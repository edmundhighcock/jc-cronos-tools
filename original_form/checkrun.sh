#!/bin/bash

qstat -u edhigh > qstattemp

compname=$(grep run"$1" qstattemp | awk '{print $8}' | cut -b 21-22)

ssh edhigh@privesaturne"$compname" "ls -rlt /scratch/rapsauve_edhigh/run"$1"*"

rm qstattemp
