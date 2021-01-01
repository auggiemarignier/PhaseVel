#!/bin/bash

ICUT=700
lmax=30
INDIR=hvh_all
INFILE=${INDIR}/hvh.000S179.asc.rwt3_std4
XYZINFILE=inputs/xyzinfile
OUTDIR=outputs
OUTPREFIX=${OUTDIR}/hvh.000S179.asc.rwt3_std4
JUNKDIR=junk

echo -- start time --

##
## This simple script calculates a global phase velocity map for fundamental Rayleigh waves with a period of about 200s (equivalent to the mode 0S61) using measurements from van Heijst and Woodhouse, 1997. The input data file is: hvh.000S061.asc.rwt3_std4

echo Build A matrix in spherical harmonics
hphasea ${INFILE} ${OUTPREFIX}a -lmax ${lmax} > ${JUNKDIR}/junka

echo Build AtA matrix 
hphaseata ${OUTPREFIX}a ${OUTPREFIX}ata 0 0 -ifd 1 > ${JUNKDIR}/junkata

echo Build Atd matrix
hphaseatd ${INFILE} ${OUTPREFIX}a ${OUTPREFIX}atd 0 0 -ifpriorm 0 -ifd 1 > ${JUNKDIR}/junkatd

echo Carry out eigenvalue-eigenvector decomposition of AtA matrix using damping in file out_apriori0
hphasedecom ${OUTPREFIX}ata ${INDIR}/out_apriori30  ${OUTPREFIX}evc > ${JUNKDIR}/junkevc

echo Obtain output phase velocity model
hphaseinv ${OUTPREFIX}evc ${OUTPREFIX}atd ${OUTPREFIX}inv ${OUTPREFIX}ata ${OUTPREFIX}dta -c ${ICUT} -ifpriorm 0  -iferror 0 > ${JUNKDIR}/junkinv

## Output model file in spherical harmonics is: ${OUTPREFIX}inv
echo Convert spherical harmonic model coefficients into dc/c for geographical coordinates

echo ${OUTPREFIX}inv > ${XYZINFILE}
echo ${OUTPREFIX}inv.xyz >> ${XYZINFILE}
echo 1 >> ${XYZINFILE}
echo 1 >> ${XYZINFILE}
echo 0 ${lmax} >> ${XYZINFILE}
echo -180 >> ${XYZINFILE}

raw2xyz_jr < ${XYZINFILE}

echo -- end time --
