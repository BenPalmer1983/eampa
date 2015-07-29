#!/bin/bash
cd ${0%/*}
workingDir=$(pwd)
binDir="bin"
binFile="eampaDev.x "
buildFiles="src/kinds.f90 src/types.f90 \
src/msubs.f90 src/constants.f90 \
src/maths.f90 src/mMaths.f90 src/general.f90 src/units.f90 \
src/plot.f90 src/typeM.f90 \
src/globals.f90  \
src/initialise.f90 \
src/loadData.f90 src/output.f90 \
src/eamGen.f90 \
src/readInput.f90 src/readEAM.f90 \
src/readConfig.f90 src/bpConfig.f90 \
src/neighbourList.f90 src/neighbourListBP.f90 \
src/preCalc.f90 \
src/calcEAM.f90 src/bpCalcEAM.f90 \
src/eval.f90 src/evalBP.f90 \
src/opti.f90 \
src/finalise.f90 \
src/eampa.f90"
commandLineA="mpif90 -g -Wall -Wno-unused-function \
-fbounds-check -fcheck=all -mtune=native "
commandLineB="mpif90 -O3 -g -Wall -Wno-unused-function \
-fbounds-check -fcheck=all -mtune=native "
commandLineC=" &> logs/build.log"
#----------------------------------------------------------------------------------
binPath=$binDir"/"$binFile
binPathLong=$workingDir"/"$binDir"/"$binFile
binPathOnlyLong=$workingDir"/"$binDir
#----------------------------------------------------------------------------------
# Make section:
echo "Check code formatting and layout? (Default: N) [Y/N]"
#read checkCode
echo "Remove unused variables? (Default: N) [Y/N]"
#read checkUnused
echo "Level 3 Optimise? (Default: N) [Y/N]"
#read optimise
echo "Make link to /bin? (Default: N) [Y/N]"
#read linkToBin
echo "Current working directory:"
pwd
echo "Make bin directory"
mkdir -p bin
#-------------------------
# 
python make.py 0
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
if [[ $checkCode == "Y" || $checkCode == "y" ]]  
then
  echo "Checking formatting and layout.";
  python make.py 1
else
  echo "Not checking formatting and layout.";
fi
echo "Compile source"
mkdir -p bin
if [[ $optimise == "Y" || $optimise == "y" ]]  
then 
  commandLine=$commandLineB" -o "$binPath$buildFiles
  eval $commandLine
else 
  commandLine=$commandLineA" -o "$binPath$buildFiles
  eval $commandLine
fi
sleep 1 
rm *.mod
if [[ $linkToBin == "Y" || $linkToBin == "y" ]]  
then
  commandRemove="sudo rm /bin/"$binFile
  eval $commandRemove
  commandLink="sudo ln -s "$workingDir"/"$binPath" /bin"
  echo "Making link: "$commandLink
  eval $commandLink
fi
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
sleep 1 
reloadBashProfile="source "$profileFile
eval $reloadBashProfile