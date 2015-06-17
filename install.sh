#!/bin/bash
cd ${0%/*}
workingDir=$(pwd)
projectTitle="EAMPA Code"
binDir="bin"
binFile="eampa.x "
buildFiles="src/kinds.f90 src/msubs.f90 src/constants.f90 \
src/maths.f90 src/mMaths.f90 src/general.f90 src/units.f90 src/globals.f90  \
src/initialise.f90 \
src/loadData.f90 src/output.f90 \
src/eamGen.f90 \
src/readInput.f90 src/readEAM.f90 src/readConfig.f90 \
src/neighbourList.f90 src/preCalc.f90 \
src/calcEAM.f90 src/eval.f90 src/opti.f90 \
src/finalise.f90 \
src/eampa.f90"
commandLineA="mpif90 -g -Wall \
-fbounds-check -fcheck=all -mtune=native "
commandLineB="mpif90 -O3 -g -Wall \
-fbounds-check -fcheck=all -mtune=native "
commandLineC=" &> logs/build.log"
commandLineD=" &> logs/build_activity.log"
#----------------------------------------------------------------------------------
# Prepare directories
binPath=$binDir"/"$binFile
binPathLong=$workingDir"/"$binDir"/"$binFile
binPathOnlyLong=$workingDir"/"$binDir
dataDirDefault=$workingDir"/data"
dataTar=$workingDir"/data/activityData.tar.gz"
examplesDir=$workingDir"/examples"
#----------------------------------------------------------------------------------
echo "=============================================================="
echo "Installing "$projectTitle
echo "University of Birmingham"
echo "=============================================================="
echo "Working directory: "$workingDir
echo " "
echo "User Options:"
echo "------------------------------"
# User choices:
echo "Check code formatting and layout? (Default: N) [Y/N]"
read checkCode
echo "Remove unused variables? (Default: N) [Y/N]"
read checkUnused
echo "Level 3 Optimise? (Default: N) [Y/N]"
read optimise
echo "Data directory: (Default: $dataDirDefault)"
read dataDir
echo "Current working directory:"
pwd
# Set default data directory
if [[ $dataDir == "" || $dataDir == " " ]] 
then
  dataDir=$dataDirDefault
fi
mkdir -p $dataDir
echo "Make data directory "$dataDir
echo "Make bin directory"
mkdir -p bin
# Archive destination
dataTarExtract=$dataDir"/data.tar.gz"
# Update compile time in globals.f90
python make.py 0
# Check for unused variables
if [[ $checkUnused == "Y" || $checkUnused == "y" ]]  
then
  mkdir -p logs
  echo "Checking unused variables.";
  commandLine=$commandLineA" -o "$binPath$buildFiles$commandLineC
  eval $commandLine
  sleep 1 
  python make.py 3
else
  echo "Not checking unused variables.";
fi
# Check layout and formatting
if [[ $checkCode == "Y" || $checkCode == "y" ]]  
then
  echo "Checking formatting and layout.";
  python make.py 1
else
  echo "Not checking formatting and layout.";
fi
# Compile source
echo "Compile source"
mkdir -p bin
if [[ $optimise == "Y" || $optimise == "y" ]]  
then 
  commandLine=$commandLineB" -o "$binPath$buildFiles$commandLineD
  eval $commandLine
else 
  commandLine=$commandLineA" -o "$binPath$buildFiles$commandLineD
  eval $commandLine
fi
sleep 1 
rm *.mod
# Print output file details to user
echo "Bin file details:"
ls -l $binPathLong
# add export to profile so user can run activity.x
exportLine="export PATH=\"\$PATH:"$binPathOnlyLong"\""
profileFile=$HOME"/.bash_profile"
touch $profileFile
if grep -q "$binPathOnlyLong" "$profileFile"; 
then
  echo $profileFile" already has the bin directory path"
else
  echo "Adding "$binPathOnlyLong" to "$profileFile
  echo $exportLine >> $profileFile
fi
source $profileFile
# Extract DATA here
# e.g. 
#if [[ $dataTar != $dataTarExtract ]]
#then
#  cp $dataTar $dataTarExtract
#fi  
#cd $dataDir
#tar xzf data.tar.gz 
#if [[ $dataTar != $dataTarExtract ]]
#then
#  rm $dataTarExtract
#fi
echo "Installation complete"
