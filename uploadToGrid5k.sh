#!/bin/bash

# https://github.com/andikleen/pmu-tools
#Simple script to upload a new version of HPCCG-RHT to Grid 5K and run them
#chmod +x executable, to give permission

#Example of use: ./uploadToGrid5K.sh
#It assumes is at the root of /miniXyce_RHT repo

folder="${PWD##*/}"
newFolder=$folder-Clean

cd ..
echo Folder Name: $folder $newFolder

cp -a $folder/ ./$newFolder
rm -rf $newFolder/.git
rm -rf $newFolder/.idea
rm -rf $newFolder/cmake-build-debug
mkdir $newFolder/cmake-build-debug
tar -czvf $newFolder.tar.gz $newFolder
rm -r -f $newFolder
echo "zip file created"

#echo "Copying files to Nancy..."
#scp $newFolder.tar.gz dperez@access.grid5000.fr:nancy/public
#echo "Copying files to Nantes..."
#scp $newFolder.tar.gz dperez@access.grid5000.fr:nantes/public
echo "Copying files to Lyon..."
scp $newFolder.tar.gz dperez@access.grid5000.fr:lyon/public
echo "Files copied to Grid5K Storage"
echo "Removing zip file"
rm HPCCG-RHT-Clean.tar.gz
echo "Success!!"

# remove previous folders
#rm -f -r nantes/public/HPCCG-RHT-Clean/ && rm -f -r nancy/public/HPCCG-RHT-Clean/

# Inside a node (if tar.gz was copied into public/ with appriate file structure)
# If name of folder is miniXyce_RHT (my machine)
#cd public/ && tar -xzvf HPCCG-RHT-Clean.tar.gz && rm HPCCG-RHT-Clean.tar.gz && cd HPCCG-RHT-Clean/cmake-build-debug/ && cmake .. && make

# If name of folder is miniXyce_RHT-master (directly from git)
#cd public/ && tar -xzvf HPCCG-RHT-Clean-master.tar.gz && rm HPCCG-RHT-Clean-master.tar.gz && cd HPCCG-RHT-Clean/cmake-build-debug/ && cmake .. && make

#nancy
#oarsub -p "cluster='graphite'" -I -l nodes=1,walltime=5

#nantes
#oarsub -p "cluster='econome'" -I -l nodes=1,walltime=5

#nantes
#oarsub -p "cluster='ecotype'" -I -l nodes=1,walltime=5

#lyon
#oarsub -p "cluster='nova'" -I -l nodes=1,walltime=5

#ssh dperez@access.grid5000.fr

# These are the commands for lyon/nova,

#Lyon - Nova 0,16 2,18 4,20... are hyperThreads (lstopo result)
#Core0       Core1     Core2        Core3     Core4       Core5    Core6      Core7
#0,16        2,18      4,20         6,22      8,24        10,26    12,28      14,30

#Core8       Core9     Core10       Core11    Core12      Core13   Core14     Core15
#1,17        3,19      5,21         7,23      9,25        11,27    13,29      15,31

# this as the commands, must copy and paste
if [ 1 -eq 0 ]; then

#After first 3 parameters (problem size) the parameters are:
#numExecutions numThreads ListOfCoresWhereThreadPinning... Examples:

#... 5 2 0 2 =
#... 5 executions, 2 threads, Leading in core 0, Trailing in core 2
#... 5 32 0 16 2 18 4 20 6 22 8 24 10 26 12 28 14 30 1 17 3 19 5 21 7 23 9 25 11 27 13 29 15 31 =
#... 5 executions, 32 threads, L1 core 0, T1 core 16, L2 core 2, T2 core 18, ..., L16 core 15, T16 core 31

# Not replicated app with different mpi ranks (to see HT results)
mpirun -np 1 HPCCG-RHT 80 80 800 5 > hpccg-norep-1rank.txt            &&
mpirun -np 16 HPCCG-RHT 80 80 50 5 > hpccg-norep-20rank.txt            &&
mpirun -np 32 HPCCG-RHT 80 80 25 5 > hpccg-norep-40rank.txt            

# 1 rank test baselines vs Our Best Approach (NewLimit + varGrouping 8)
mpirun -np 1 HPCCG-AC 80 80 800 5 2 0 2 > hpccg-rep-ac-1rank-noHT.txt   &&
mpirun -np 1 HPCCG-AC 80 80 800 5 2 0 16 > hpccg-rep-ac-1rank-HT.txt  &&

mpirun -np 1 HPCCG-SRMT 80 80 800 5 2 0 2 > hpccg-rep-srmt-1rank-noHT.txt   &&
mpirun -np 1 HPCCG-SRMT 80 80 800 5 2 0 16 > hpccg-rep-srmt-1rank-HT.txt  &&

mpirun -np 1 HPCCG-RHT 80 80 800 5 2 0 2 > hpccg-rep-rht-1rank-noHT.txt   &&
mpirun -np 1 HPCCG-RHT 80 80 800 5 2 0 16 > hpccg-rep-rht-1rank-HT.txt  

# 16 ranks test baselines vs Our Best Approach (NewLimit + varGrouping 8)
mpirun -np 16 HPCCG-AC 80 80 50 5 32 0 2 4 6 8 10 12 14 1 3 5 7 9 11 13 15 16 18 20 22 24 26 28 30 17 19 21 23 25 27 29 31 > hpccg-rep-ac-20rank-noHT.txt            &&
mpirun -np 16 HPCCG-AC 80 80 50 5 32 0 16 2 18 4 20 6 22 8 24 10 26 12 28 14 30 1 17 3 19 5 21 7 23 9 25 11 27 13 29 15 31 > hpccg-rep-ac-20rank-HT.txt            &&

mpirun -np 16 HPCCG-SRMT 80 80 50 5 32 0 2 4 6 8 10 12 14 1 3 5 7 9 11 13 15 16 18 20 22 24 26 28 30 17 19 21 23 25 27 29 31 > hpccg-rep-srmt-20rank-noHT.txt            &&
mpirun -np 16 HPCCG-SRMT 80 80 50 5 32 0 16 2 18 4 20 6 22 8 24 10 26 12 28 14 30 1 17 3 19 5 21 7 23 9 25 11 27 13 29 15 31 > hpccg-rep-srmt-20rank-HT.txt            &&

mpirun -np 16 HPCCG-RHT 80 80 50 5 32 0 2 4 6 8 10 12 14 1 3 5 7 9 11 13 15 16 18 20 22 24 26 28 30 17 19 21 23 25 27 29 31 > hpccg-rep-rht-20rank-noHT.txt            &&
mpirun -np 16 HPCCG-RHT 80 80 50 5 32 0 16 2 18 4 20 6 22 8 24 10 26 12 28 14 30 1 17 3 19 5 21 7 23 9 25 11 27 13 29 15 31 > hpccg-rep-rht-20rank-HT.txt

fi




