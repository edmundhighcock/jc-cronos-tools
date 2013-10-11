#!/bin/bash

qstat -u $USER > qstattemp

compname=$(grep run"$1" qstattemp | awk '{print $8}' | cut -b 21-22)

ssh $USER@privesaturne"$compname" "ls -rlt /scratch/rapsauve_$USER/run"$1"*"

rm qstattemp
