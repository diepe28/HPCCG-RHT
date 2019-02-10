#!/bin/bash

# https://github.com/andikleen/pmu-tools
#Simple script to upload a new version of HPCCG-RHT to Grid 5K and run them
#chmod +x executable, to give permission

#sudo-g5k
#sudo apt-get install msr-tools
#sudo modprobe msr

#echo Printing values
#for cpu in {0..19}; do sudo /usr/sbin/rdmsr -p$cpu 0x1a0 -f 38:38; done

#echo Disabling turbo-mode
#for cpu in {0..19}; do sudo /usr/sbin/wrmsr -p$cpu 0x1a0 0x4000850089; do$

#echo Printing values
#for cpu in {0..19}; do sudo /usr/sbin/rdmsr -p$cpu 0x1a0 -f 38:38; done

#echo Done disabling turbo-mode

nRuns=10
nRanks=1
nCores=2
coresHT="0 20"
coresNOHT="0 2"

x=180
y=180
z=250

echo Running baseline 1 rank
mpirun -np $nRanks HPCCG-WANG $x $y $z $nRuns > 0-baseline-1rank.txt
mpirun -np $nRanks IMP-HPCCG-WANG $x $y $z $nRuns > 1-IMP-baseline-1rank.txt

echo Running replicated wang 1rank noHT
mpirun -np $nRanks HPCCG-WANG $x $y $z $nRuns $nCores $coresNOHT > 2-replicated-wang-1rank-noHT.txt
mpirun -np $nRanks IMP-HPCCG-WANG $x $y $z $nRuns $nCores $coresNOHT > 3-IMP-replicated-wang-1rank-noHT.txt

echo Running replicated wang 1 rank HT
mpirun -np $nRanks HPCCG-WANG $x $y $z $nRuns $nCores $coresHT > 4-replicated-wang-1rank-HT.txt
mpirun -np $nRanks IMP-HPCCG-WANG $x $y $z $nRuns $nCores $coresHT > 5-IMP-replicated-wang-1rank-HT.txt

echo Running replicated wang with var grouping 1 rank noHT
mpirun -np $nRanks HPCCG-WANG-VG $x $y $z $nRuns $nCores $coresNOHT > 6-replicated-wang-vg-1rank-noHT.txt
mpirun -np $nRanks IMP-HPCCG-WANG-VG $x $y $z $nRuns $nCores $coresNOHT > 7-IMP-replicated-wang-vg-1rank-noHT.txt

echo Running replicated wang with var grouping 1 rank HT
mpirun -np $nRanks HPCCG-WANG-VG $x $y $z $nRuns $nCores $coresHT > 8-replicated-wang-vg-1rank-HT.txt
mpirun -np $nRanks IMP-HPCCG-WANG-VG $x $y $z $nRuns $nCores $coresHT > 9-IMP-replicated-wang-vg-1rank-HT.txt

echo Running replicated wang just volatiles 1 rank noHT
mpirun -np $nRanks HPCCG-WANG-JV $x $y $z $nRuns $nCores $coresNOHT > 10-replicated-wang-jv-1rank-noHT.txt
mpirun -np $nRanks IMP-HPCCG-WANG-JV $x $y $z $nRuns $nCores $coresNOHT > 11-IMP-replicated-wang-jv-1rank-noHT.txt

echo Running replicated wang just volatiles 1 rank HT
mpirun -np $nRanks HPCCG-WANG-JV $x $y $z $nRuns $nCores $coresHT > 12-replicated-wang-jv-1rank-HT.txt
mpirun -np $nRanks IMP-HPCCG-WANG-JV $x $y $z $nRuns $nCores $coresHT > 13-IMP-replicated-wang-jv-1rank-HT.txt


coresNOHT="0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39"
coresHT="0 20 2 22 4 24 6 26 8 28 10 30 12 32 14 34 16 36 18 38 1 21 3 23 5 25 7 27 9 29 11 31 13 33 15 35 17 37 19 39"
nRanks=20
nCores=40
x=190
y=190
z=170

echo Running baseline 20,40 ranks
mpirun -np $nRanks HPCCG-WANG $x $y $z $nRuns > 14-mybaseline-20rank.txt

echo Running replicated wang 20 rank noHT
mpirun -np $nRanks HPCCG-WANG $x $y $z $nRuns $nCores $coresNOHT > 15-myreplicated-wang-20rank-noHT.txt
mpirun -np $nRanks IMP-HPCCG-WANG $x $y $z $nRuns $nCores $coresNOHT > 16-IMP-myreplicated-wang-20rank-noHT.txt

echo Running replicated wang 20 rank HT
mpirun -np $nRanks HPCCG-WANG $x $y $z $nRuns $nCores $coresHT > 17-myreplicated-wang-20rank-HT.txt
mpirun -np $nRanks IMP-HPCCG-WANG $x $y $z $nRuns $nCores $coresHT > 18-IMP-myreplicated-wang-20rank-HT.txt

echo Running replicated wang with var grouping 20 rank noHT
mpirun -np $nRanks HPCCG-WANG-VG $x $y $z $nRuns $nCores $coresNOHT > 19-myreplicated-wang-vg-20rank-noHT.txt
mpirun -np $nRanks IMP-HPCCG-WANG-VG $x $y $z $nRuns $nCores $coresNOHT > 20-IMP-myreplicated-wang-vg-20rank-noHT.txt

echo Running replicated wang with var grouping 20 rank HT
mpirun -np $nRanks HPCCG-WANG-VG $x $y $z $nRuns $nCores $coresHT > 21-myreplicated-wang-vg-20rank-HT.txt
mpirun -np $nRanks IMP-HPCCG-WANG-VG $x $y $z $nRuns $nCores $coresHT > 22-IMP-myreplicated-wang-vg-20rank-HT.txt

echo Running replicated wang just volatiles grouping 20 ranks noHT
mpirun -np $nRanks HPCCG-WANG-JV $x $y $z $nRuns $nCores $coresNOHT > 23-myreplicated-wang-jv-20rank-noHT.txt
mpirun -np $nRanks IMP-HPCCG-WANG-JV $x $y $z $nRuns $nCores $coresNOHT > 24-IMP-myreplicated-wang-jv-20rank-noHT.txt

echo Running replicated wang just volatiles grouping 20 ranks HT
mpirun -np $nRanks HPCCG-WANG-JV $x $y $z $nRuns $nCores $coresHT > 25-myreplicated-wang-jv-20rank-HT.txt
mpirun -np $nRanks IMP-HPCCG-WANG-JV $x $y $z $nRuns $nCores $coresHT > 26-IMP-myreplicated-wang-jv-20rank-HT.txt
