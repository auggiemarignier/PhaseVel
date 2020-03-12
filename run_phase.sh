#!/bin/bash

INDIR=inputs
INFILE=${INDIR}/hvh.000S061.asc.rwt3_std4
XYZINFILE=${INDIR}/xyzinfile
OUTDIR=outputs
OUTPREFIX=${OUTDIR}/hvh.000S061.asc.rwt3_std4
JUNKDIR=junk

echo -- start time --

##
## This simple script calculates a global phase velocity map for fundamental Rayleigh waves with a period of about 200s (equivalent to the mode 0S61) using measurements from van Heijst and Woodhouse, 1997. The input data file is: hvh.000S061.asc.rwt3_std4

echo Build A matrix in spherical harmonics
hphasea ${INFILE} ${OUTPREFIX}a -lmax 20 > ${JUNKDIR}/junka

echo Build AtA matrix 
hphaseata ${OUTPREFIX}a ${OUTPREFIX}ata 0 0 > ${JUNKDIR}/junkata

echo Build Atd matrix
hphaseatd ${INFILE} ${OUTPREFIX}a ${OUTPREFIX}atd 0 0 -ifpriorm 0 > ${JUNKDIR}/junkatd

echo Build Atd matrix with crustal corrections  
## The file s_use000S061.gridj2.sph2 contains the crustal corrections.
hphaseatd ${INFILE} ${OUTPREFIX}a ${OUTPREFIX}atd_c  0 0 -ifpriorm 1 -ifile ${INDIR}/s_use000S061.gridj2.sph2 > ${JUNKDIR}/junkatdc

echo Carry out eigenvalue-eigenvector decomposition of AtA matrix using damping in file out_apriori0
hphasedecom ${OUTPREFIX}ata ${INDIR}/out_apriori0  ${OUTPREFIX}evc > ${JUNKDIR}/junkevc

echo Obtain output phase velocity model
hphaseinv ${OUTPREFIX}evc ${OUTPREFIX}atd_c ${OUTPREFIX}inv_c ${OUTPREFIX}ata ${OUTPREFIX}dta -c 300  -ifpriorm 0  -iferror 0 > ${JUNKDIR}/junkinv

## Output model file in spherical harmonics is: ${OUTPREFIX}inv
echo Convert spherical harmonic model coefficients into dc/c for geographical coordinates

echo ${OUTPREFIX}inv_c > ${XYZINFILE}
echo ${OUTPREFIX}inv_c.xyz >> ${XYZINFILE}
echo 1 >> ${XYZINFILE}
echo 1 >> ${XYZINFILE}
echo 0 20 >> ${XYZINFILE}
echo -180 >> ${XYZINFILE}

raw2xyz_jr < ${XYZINFILE}

echo -- end time --
