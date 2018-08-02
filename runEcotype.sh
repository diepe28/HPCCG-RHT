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

echo Running baseline 1 rank
mpirun -np 1 HPCCG-WANG 150 150 250 5 > 0-baseline-1rank.txt
mpirun -np 1 IMP-HPCCG-WANG 150 150 250 5 > 1-IMP-baseline-1rank.txt

echo Running replicated wang 1rank noHT
mpirun -np 1 HPCCG-WANG 150 150 250 5 2 0 2 > 2-replicated-wang-1rank-noHT.txt
mpirun -np 1 IMP-HPCCG-WANG 150 150 250 5 2 0 2 > 3-IMP-replicated-wang-1rank-noHT.txt

echo Running replicated wang 1 rank HT
mpirun -np 1 HPCCG-WANG 150 150 250 5 2 0 20 > 4-replicated-wang-1rank-HT.txt
mpirun -np 1 IMP-HPCCG-WANG 150 150 250 5 2 0 20 > 5-IMP-replicated-wang-1rank-HT.txt

echo Running replicated wang with var grouping 1 rank noHT
mpirun -np 1 HPCCG-WANG-VG 150 150 250 5 2 0 2 > 6-replicated-wang-vg-1rank-noHT.txt
mpirun -np 1 IMP-HPCCG-WANG-VG 150 150 250 5 2 0 2 > 7-IMP-replicated-wang-vg-1rank-noHT.txt

echo Running replicated wang with var grouping 1 rank HT
mpirun -np 1 HPCCG-WANG-VG 150 150 250 5 2 0 20 > 8-replicated-wang-vg-1rank-HT.txt
mpirun -np 1 IMP-HPCCG-WANG-VG 150 150 250 5 2 0 20 > 9-IMP-replicated-wang-vg-1rank-HT.txt

echo Running replicated wang just volatiles 1 rank noHT
mpirun -np 1 HPCCG-WANG-JV 150 150 250 5 2 0 2 > 10-replicated-wang-jv-1rank-noHT.txt
mpirun -np 1 IMP-HPCCG-WANG-JV 150 150 250 5 2 0 2 > 11-IMP-replicated-wang-jv-1rank-noHT.txt

echo Running replicated wang just volatiles 1 rank HT
mpirun -np 1 HPCCG-WANG-JV 150 150 250 5 2 0 20 > 12-replicated-wang-jv-1rank-HT.txt
mpirun -np 1 IMP-HPCCG-WANG-JV 150 150 250 5 2 0 20 > 13-IMP-replicated-wang-jv-1rank-HT.txt

echo Running baseline 20,40 ranks
mpirun -np 20 HPCCG-WANG 180 180 160 5 > 14-mybaseline-20rank.txt
mpirun -np 40 HPCCG-WANG 180 180 160 5 > 15-mybaseline-40rank.txt

echo Running replicated wang 20 rank noHT
mpirun -np 20 HPCCG-WANG 180 180 160 5 40 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 > 16-myreplicated-wang-20rank-noHT.txt
mpirun -np 20 IMP-HPCCG-WANG 180 180 160 5 40 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 > 17-IMP-myreplicated-wang-20rank-noHT.txt

echo Running replicated wang 20 rank HT
mpirun -np 20 HPCCG-WANG 180 180 160 5 40 0 20 2 22 4 24 6 26 8 28 10 30 12 32 14 34 16 36 18 38 1 21 3 23 5 25 7 27 9 29 11 31 13 33 15 35 17 37 19 39 > 18-myreplicated-wang-20rank-HT.txt
mpirun -np 20 IMP-HPCCG-WANG 180 180 160 5 40 0 20 2 22 4 24 6 26 8 28 10 30 12 32 14 34 16 36 18 38 1 21 3 23 5 25 7 27 9 29 11 31 13 33 15 35 17 37 19 39 > 19-IMP-myreplicated-wang-20rank-HT.txt

echo Running replicated wang with var grouping 20 rank noHT
mpirun -np 20 HPCCG-WANG-VG 180 180 160 5 40 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 > 20-myreplicated-wang-vg-20rank-noHT.txt
mpirun -np 20 IMP-HPCCG-WANG-VG 180 180 160 5 40 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 > 21-IMP-myreplicated-wang-vg-20rank-noHT.txt

echo Running replicated wang with var grouping 20 rank HT
mpirun -np 20 HPCCG-WANG-VG 180 180 160 5 40 0 20 2 22 4 24 6 26 8 28 10 30 12 32 14 34 16 36 18 38 1 21 3 23 5 25 7 27 9 29 11 31 13 33 15 35 17 37 19 39 > 22-myreplicated-wang-vg-20rank-HT.txt
mpirun -np 20 IMP-HPCCG-WANG-VG 180 180 160 5 40 0 20 2 22 4 24 6 26 8 28 10 30 12 32 14 34 16 36 18 38 1 21 3 23 5 25 7 27 9 29 11 31 13 33 15 35 17 37 19 39 > 23-IMP-myreplicated-wang-vg-20rank-HT.txt

echo Running replicated wang just volatiles grouping 20 ranks noHT
mpirun -np 20 HPCCG-WANG-JV 180 180 160 5 40 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 > 24-myreplicated-wang-jv-20rank-noHT.txt
mpirun -np 20 IMP-HPCCG-WANG-JV 180 180 160 5 40 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 > 25-IMP-myreplicated-wang-jv-20rank-noHT.txt

echo Running replicated wang just volatiles grouping 20 ranks HT
mpirun -np 20 HPCCG-WANG-JV 180 180 160 5 40 0 20 2 22 4 24 6 26 8 28 10 30 12 32 14 34 16 36 18 38 1 21 3 23 5 25 7 27 9 29 11 31 13 33 15 35 17 37 19 39 > 26-myreplicated-wang-jv-20rank-HT.txt
mpirun -np 20 IMP-HPCCG-WANG-JV 180 180 160 5 40 0 20 2 22 4 24 6 26 8 28 10 30 12 32 14 34 16 36 18 38 1 21 3 23 5 25 7 27 9 29 11 31 13 33 15 35 17 37 19 39 > 27-IMP-myreplicated-wang-jv-20rank-HT.txt
