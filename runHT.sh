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

x=200
y=200
z=800
nRuns=10
nRanks=20
nCores=40

echo Running baseline
mpirun -np 1 HPCCG-WANG $x $y $z 800 $nRuns > baseline-1rank.txt
mpirun -np 2 HPCCG-WANG $x $y $z 400 $nRuns > baseline-2rank.txt
mpirun -np 4 HPCCG-WANG $x $y $z 200 $nRuns > baseline-4rank.txt
mpirun -np 8 HPCCG-WANG $x $y $z 100 $nRuns > baseline-8rank.txt
mpirun -np 16 HPCCG-WANG $x $y $z 50 $nRuns > baseline-16rank.txt
mpirun -np 20 HPCCG-WANG $x $y $z 40 $nRuns > baseline-20rank.txt
mpirun -np 40 HPCCG-WANG $x $y $z 20 $nRuns > baseline-40rank.txt

