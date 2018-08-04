#!/bin/bash

# https://github.com/andikleen/pmu-tools
#Simple script to upload a new version of HPCCG-RHT to Grid 5K and run them
#chmod +x executable, to give permission

#Example of use: ./myScriptHPCCG.sh HPCCG-RHT g5k
#Example of use: ./myScriptHPCCG.sh HPCCG-RHT cenat

#rm -r -f nancy/public/HPCCG-RHT-Clean/
#rm nancy/public/HPCCG-RHT-Clean.tar.gz 
#rm -r -f lille/public/HPCCG-RH T-Clean/
#rm lille/public/HPCCG-RHT-Clean.tar.gz 

folder="$1"
newFolder="$1"-Clean
cluster="$2" #cenat | g5k

if [ "$2" == "cenat" ]
then
    cluster="dperez@cluster.cenat.ac.cr"
else
    cluster="dperez@access.grid5000.fr"
fi

echo Cluster http: $cluster

echo Folder Name: $folder $newFolder

cp -a $folder/ ./$newFolder
rm -rf $newFolder/.git
rm -rf $newFolder/.idea

#copying scripts to tempFolder
mkdir $newFolder/tempFolder
cp -a $newFolder/runNova.sh ./$newFolder/tempFolder
cp -a $newFolder/runEcotype.sh ./$newFolder/tempFolder
cp -a $newFolder/runEcotype-testHT.sh ./$newFolder/tempFolder
cp -a $newFolder/scriptTest.sh ./$newFolder/tempFolder

rm -rf $newFolder/cmake-build-debug

#renaming tempFolder to cmake-build-debug
mv $newFolder/tempFolder/ $newFolder/cmake-build-debug/

tar -czvf $newFolder.tar.gz $newFolder
rm -r -f $newFolder
echo "zip file created"

#echo "Copying files to Nancy..."
#scp $newFolder.tar.gz dperez@access.grid5000.fr:nancy/public

if [ "$2" == "cenat" ]
then
    echo "Copying files to Kabre..."
    scp $newFolder.tar.gz $cluster:~/public/workspace
    echo "Files copied to Kabre Storage"
else
    echo "Copying files to Nantes..."
    scp $newFolder.tar.gz $cluster:nantes/public
    echo "Copying files to Lyon..."
    scp $newFolder.tar.gz $cluster:lyon/public
    echo "Files copied to Grid5K Storage"
fi    
echo "Removing zip file"
rm HPCCG-RHT-Clean.tar.gz
echo "Success!!"

#rm -f -r nantes/public/HPCCG-RHT-Clean/ && rm -f -r nancy/public/HPCCG-RHT-Clean/

# Inside a node (if tar.gz was copied into public/ with appriate file structure)
#cd public/ && tar -xzvf HPCCG-RHT-Clean.tar.gz && rm HPCCG-RHT-Clean.tar.gz && cd HPCCG-RHT-Clean/cmake-build-debug/ && cmake .. && make

#Install newer gcc
#sudo-g5k && sudo -H /bin/bash 
#sudo echo "deb http://ftp.us.debian.org/debian testing main contrib non-free" > /etc/apt/sources.list.d/preferences.list &&
#sudo echo "Package: *
#Pin: release a=testing
#Pin-Priority: 100" > /etc/apt/preferences.d/preferences.list && exit
#sudo apt-get update && sudo apt-get install -t testing g++

#nancy
#oarsub -p "cluster='graphite'" -I -l nodes=1,walltime=5

#nantes
#oarsub -p "cluster='econome'" -I -l nodes=1,walltime=5

#nantes
#oarsub -p "cluster='ecotype'" -I -l nodes=1,walltime=5

#lyon
#oarsub -p "cluster='nova'" -I -l nodes=1,walltime=5

#kill a job (check job id at grid5k website)
#oardel 12345

#connect to a job with id 12345
#oarsub -C 12345

#execute a script file
#oarsub -p "cluster='nova'" -l nodes=1,walltime=4 "/home/dperez/public/HPCCG-RHT-Clean/cmake-build-debug/runNova.sh"

#oarsub -p "cluster='ecotype'" -l nodes=1,walltime=7 "/home/dperez/public/HPCCG-RHT-Clean/cmake-build-debug/runEcotype.sh"
#oarsub -p "cluster='ecotype'" -l nodes=1,walltime=7 "/home/dperez/public/HPCCG-RHT-Clean/cmake-build-debug/runEcotype-testHT.sh"
#ssh dperez@access.grid5000.fr

#extend a jobtime
#oarwalltime 12345 +1:30

#ssh dperez@cluster.cenat.ac.cr
#module load gcc/7.2.0 && module load cmake/3.12.0-rc2 && export CXX=/opt/compilers/gcc-7.2.0/bin/g++ && export CC=/opt/compilers/gcc-7.2.0/bin/gcc && 
#tar -xzvf HPCCG-RHT-Clean.tar.gz && rm HPCCG-RHT-Clean.tar.gz && cd HPCCG-RHT-Clean/cmake-build-debug/ && cmake .. && make
#qsub hpccg.pbs
#qdel jobId

#Usuario: dperez
#Clave: lnD8RhsF

