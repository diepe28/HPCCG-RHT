#!/bin/bash

# https://github.com/andikleen/pmu-tools
#Simple script to upload a new version of HPCCG-RHT to Grid 5K and run them
#chmod +x executable, to give permission

#sudo-g5k
#sudo apt-get install msr-tools
#sudo modprobe msr
#
#echo Printing values
#for cpu in {0..19}; do sudo /usr/sbin/rdmsr -p$cpu 0x1a0 -f 38:38; done
#
#echo Disabling turbo-mode
#for cpu in {0..19}; do sudo /usr/sbin/wrmsr -p$cpu 0x1a0 0x4000850089; do$
#
#echo Printing values
#for cpu in {0..19}; do sudo /usr/sbin/rdmsr -p$cpu 0x1a0 -f 38:38; done
#
#echo Done disabling turbo-mode

echo Running baseline
mpirun -np 1 HPCCG-WANG 200 200 800 10 > mybaseline-1rank.txt
mpirun -np 2 HPCCG-WANG 200 200 400 10 > mybaseline-2rank.txt
mpirun -np 4 HPCCG-WANG 200 200 200 10 > mybaseline-4rank.txt
mpirun -np 8 HPCCG-WANG 200 200 100 10 > mybaseline-8rank.txt
mpirun -np 16 HPCCG-WANG 200 200 50 10 > mybaseline-16rank.txt
mpirun -np 20 HPCCG-WANG 200 200 40 10 > mybaseline-20rank.txt
mpirun -np 40 HPCCG-WANG 200 200 20 10 > mybaseline-40rank.txt

