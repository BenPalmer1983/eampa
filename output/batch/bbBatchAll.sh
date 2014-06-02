#!/bin/bash
#MOAB -l walltime=00:10:00
#MOAB -j oe
#MOAB -q bbtest
cd "$HOME/pwscf/test\"
module load apps/intel/v2013.0.079
module load apps/QE/v5.0.2
mpirun pw.x < pwscfbatch1201.in > pwscfbatch1201.out
mpirun pw.x < pwscfbatch1202.in > pwscfbatch1202.out
mpirun pw.x < pwscfbatch1203.in > pwscfbatch1203.out
mpirun pw.x < pwscfbatch1204.in > pwscfbatch1204.out
mpirun pw.x < pwscfbatch1205.in > pwscfbatch1205.out
mpirun pw.x < pwscfbatch1206.in > pwscfbatch1206.out
mpirun pw.x < pwscfbatch1207.in > pwscfbatch1207.out
mpirun pw.x < pwscfbatch1208.in > pwscfbatch1208.out
mpirun pw.x < pwscfbatch1209.in > pwscfbatch1209.out
mpirun pw.x < pwscfbatch1210.in > pwscfbatch1210.out
mpirun pw.x < pwscfbatch1211.in > pwscfbatch1211.out
mpirun pw.x < pwscfbatch1212.in > pwscfbatch1212.out
mpirun pw.x < pwscfbatch1213.in > pwscfbatch1213.out
mpirun pw.x < pwscfbatch1214.in > pwscfbatch1214.out
mpirun pw.x < pwscfbatch1215.in > pwscfbatch1215.out
mpirun pw.x < pwscfbatch1216.in > pwscfbatch1216.out
mpirun pw.x < pwscfbatch1217.in > pwscfbatch1217.out
mpirun pw.x < pwscfbatch1218.in > pwscfbatch1218.out
mpirun pw.x < pwscfbatch1219.in > pwscfbatch1219.out
mpirun pw.x < pwscfbatch1220.in > pwscfbatch1220.out
mpirun pw.x < pwscfbatch1221.in > pwscfbatch1221.out
mpirun pw.x < pwscfbatch1222.in > pwscfbatch1222.out
mpirun pw.x < pwscfbatch1223.in > pwscfbatch1223.out
