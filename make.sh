#!/bin/bash
cd ${0%/*}
workingDir=$(pwd)
# Directories
srcDir=$workingDir"/src"
modDir=$workingDir"/mod"
binDir=$workingDir"/bin"
mkdir -p $modDir
mkdir -p $binDir
# -----------Library Start--------
libDir="/code/lib/libBP"
modDirLib=$libDir"/mod"
libDirLib=$libDir"/lib"
libLink="-L"$libDirLib" -lBP"
#Make included library
#"/code/lib/makeLib.sh"
# -----------Library End--------
binFile="eampaDev.x "
fortLine="mpif90 -g \
-fbounds-check -mtune=native " # -Wno-unused-function -fcheck=all  -Wall -O3 
# Build files

buildFiles="types.f90 \
msubs.f90 \
typesM.f90 \
globals.f90  \
initialise.f90 \
loadData.f90 output.f90 \
eamGen.f90 \
readInput.f90 readEAM.f90 makePotential.f90 \
readConfig.f90 bpConfig.f90 \
neighbourList.f90 neighbourListBP.f90 neighbourListRelax.f90 \
preCalc.f90 \
calcEAM.f90 bpCalcEAM.f90 relaxCalcEAM.f90  \
eval.f90 evalBP.f90 relax.f90 \
opti.f90 \
finalise.f90 \
eampa.f90"

# Update compile time in globals.f90
#python make.py 0
# Compile
cd $srcDir
commandLine=$fortLine" -J "$modDir"  -I"$modDirLib" "
commandLine=$commandLine" "$buildFiles" "
commandLine=$commandLine" "$libLink" "
commandLine=$commandLine" -o "$binDir"/"$binFile
echo "----------------------------------------------------------------------------------"
echo "Compile with Fortran:"
echo " "
echo $commandLine
echo "----------------------------------------------------------------------------------"
eval $commandLine
# add export to profile so user can run activity.x
exportLine="export PATH=\"\$PATH:"$binDir"\""
profileFile=$HOME"/.bash_profile"
touch $profileFile
if grep -q "$binDir" "$profileFile"; 
then
  echo $profileFile" already has the bin directory path"
else
  echo "Adding "$binDir" to "$profileFile
  echo $exportLine >> $profileFile
fi