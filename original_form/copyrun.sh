#!/bin/bash

qstat -u $USER > qstattemp

compname=$(grep run"$1" qstattemp | awk '{print $8}' | cut -b 21-22)

scp $USER@privesaturne"$compname":/scratch/rapsauve_$USER/run"$1"_*_"$2".mat $CRONOS_TEMP


