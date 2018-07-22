#!/bin/bash

# https://github.com/andikleen/pmu-tools
#Simple script to upload a new version of HPCCG-RHT to Grid 5K and run them
#chmod +x executable, to give permission

echo Running baseline 1 rank
mpirun -np 1 HPCCG-WANG 150 150 250 5 > baseline-1rank.txt

echo Running replicated wang 1rank noHT
mpirun -np 1 HPCCG-WANG 150 150 250 5 2 0 2 > replicated-wang-1rank-noHT.txt
mpirun -np 1 IMP-HPCCG-WANG 150 150 250 5 2 0 2 > imp-replicated-wang-1rank-noHT.txt

echo Running replicated wang 1 rank HT
mpirun -np 1 HPCCG-WANG 150 150 250 5 2 0 16 > replicated-wang-1rank-HT.txt  
mpirun -np 1 IMP-HPCCG-WANG 150 150 250 5 2 0 16 > IMP-replicated-wang-1rank-HT.txt  

echo Running replicated wang with var grouping 1 rank noHT
mpirun -np 1 HPCCG-WANG-VG 150 150 250 5 2 0 2 > replicated-wang-vg-1rank-noHT.txt
mpirun -np 1 IMP-HPCCG-WANG-VG 150 150 250 5 2 0 2 > IMP-replicated-wang-vg-1rank-noHT.txt

echo Running replicated wang with var grouping 1 rank HT
mpirun -np 1 HPCCG-WANG-VG 150 150 250 5 2 0 16 > replicated-wang-vg-1rank-HT.txt 
mpirun -np 1 IMP-HPCCG-WANG-VG 150 150 250 5 2 0 16 > IMP-replicated-wang-vg-1rank-HT.txt 

echo Running replicated wang just volatiles 1rank noHT
mpirun -np 1 HPCCG-WANG-JV 150 150 250 5 2 0 2 > replicated-wang-jv-1rank-noHT.txt   
mpirun -np 1 IMP-HPCCG-WANG-JV 150 150 250 5 2 0 2 > IMP-replicated-wang-jv-1rank-noHT.txt   

echo Running replicated wang just volatiles 1rank HT
mpirun -np 1 HPCCG-WANG-JV 150 150 250 5 2 0 16 > replicated-wang-jv-1rank-HT.txt 
mpirun -np 1 IMP-HPCCG-WANG-JV 150 150 250 5 2 0 16 > replicated-wang-jv-1rank-HT.txt 

echo Running baseline 16 ranks
mpirun -np 16 HPCCG-WANG 150 150 110 5 > mybaseline-16rank.txt

echo Running replicated wang 16 rank noHT
mpirun -np 16 HPCCG-WANG 150 150 110 5 32 0 2 4 6 8 10 12 14 1 3 5 7 9 11 13 15 16 18 20 22 24 26 28 30 17 19 21 23 25 27 29 31 > myreplicated-wang-16rank-noHT.txt
mpirun -np 16 IMP-HPCCG-WANG 150 150 110 5 32 0 2 4 6 8 10 12 14 1 3 5 7 9 11 13 15 16 18 20 22 24 26 28 30 17 19 21 23 25 27 29 31 > IMP-myreplicated-wang-16rank-noHT.txt

echo Running replicated wang 16 rank HT
mpirun -np 16 HPCCG-WANG 150 150 110 5 32 0 16 2 18 4 20 6 22 8 24 10 26 12 28 14 30 1 17 3 19 5 21 7 23 9 25 11 27 13 29 15 31 > myreplicated-wang-16rank-HT.txt
mpirun -np 16 IMP-HPCCG-WANG 150 150 110 5 32 0 16 2 18 4 20 6 22 8 24 10 26 12 28 14 30 1 17 3 19 5 21 7 23 9 25 11 27 13 29 15 31 > IMP-myreplicated-wang-16rank-HT.txt

echo Running replicated wang with var grouping 16 rank noHT
mpirun -np 16 HPCCG-WANG-VG 150 150 110 5 32 0 2 4 6 8 10 12 14 1 3 5 7 9 11 13 15 16 18 20 22 24 26 28 30 17 19 21 23 25 27 29 31 > myreplicated-wang-vg-16rank-noHT.txt
mpirun -np 16 IMP-HPCCG-WANG-VG 150 150 110 5 32 0 2 4 6 8 10 12 14 1 3 5 7 9 11 13 15 16 18 20 22 24 26 28 30 17 19 21 23 25 27 29 31 > IMP-myreplicated-wang-vg-16rank-noHT.txt

echo Running replicated wang with var grouping 16 rank HT
mpirun -np 16 HPCCG-WANG-VG 150 150 110 5 32 0 16 2 18 4 20 6 22 8 24 10 26 12 28 14 30 1 17 3 19 5 21 7 23 9 25 11 27 13 29 15 31 > myreplicated-wang-vg-16rank-HT.txt
mpirun -np 16 IMP-HPCCG-WANG-VG 150 150 110 5 32 0 16 2 18 4 20 6 22 8 24 10 26 12 28 14 30 1 17 3 19 5 21 7 23 9 25 11 27 13 29 15 31 > IMP-myreplicated-wang-vg-16rank-HT.txt

echo Running replicated wang just volatiles grouping 16 rank noHT
mpirun -np 16 HPCCG-WANG-JV 150 150 110 5 32 0 2 4 6 8 10 12 14 1 3 5 7 9 11 13 15 16 18 20 22 24 26 28 30 17 19 21 23 25 27 29 31 > myreplicated-wang-jv-16rank-noHT.txt  
mpirun -np 16 IMP-HPCCG-WANG-JV 150 150 110 5 32 0 2 4 6 8 10 12 14 1 3 5 7 9 11 13 15 16 18 20 22 24 26 28 30 17 19 21 23 25 27 29 31 > IMP-myreplicated-wang-jv-16rank-noHT.txt  

echo Running replicated wang just volatiles grouping 16 rank HT
mpirun -np 16 HPCCG-WANG-JV 150 150 110 5 32 0 16 2 18 4 20 6 22 8 24 10 26 12 28 14 30 1 17 3 19 5 21 7 23 9 25 11 27 13 29 15 31 > myreplicated-wang-jv-16rank-HT.txt
mpirun -np 16 IMP-HPCCG-WANG-JV 150 150 110 5 32 0 16 2 18 4 20 6 22 8 24 10 26 12 28 14 30 1 17 3 19 5 21 7 23 9 25 11 27 13 29 15 31 > IMP-myreplicated-wang-jv-16rank-HT.txt
