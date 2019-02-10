#!/bin/bash

# https://github.com/andikleen/pmu-tools
#Simple script to upload a new version of HPCCG-RHT to Grid 5K and run them
#chmod +x executable, to give permission

echo This is a testing script
mpirun -np 1 HPCCG-WANG 50 50 25 5

echo Running again just to be sure
mpirun -np 1 HPCCG-WANG 50 50 25 3
