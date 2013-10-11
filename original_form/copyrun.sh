#!/bin/bash

qstat -u edhigh > qstattemp

compname=$(grep run"$1" qstattemp | awk '{print $8}' | cut -b 21-22)

scp edhigh@privesaturne"$compname":/scratch/rapsauve_edhigh/run"$1"_*_"$2".mat ~/IntegratedModelling/cronos_temp/


