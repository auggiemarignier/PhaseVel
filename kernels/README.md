# Phase Velocity Sensitivity Kernels

Calculates the phase velocity sensitivity kernels using eigenfunctions from [MINEOS](https://github.com/geodynamics/mineos).

## Installation

Build system is configured by `cmake`.  This will download and build MINEOS (without installing the binaries), and then build the kernels scripts.

From the `kernels` directory, run

```bash
mkdir build
cd build
cmake ..
make
```

This will produce the `kernels` executable in the `build` directory.  This has only been tested on Mac OS X using `gfortran 10.2.0` as the fortran compiler.

## Usage

From the `build` directory, run

```bash
./kernels <path/to/earth/model> <output/directory>
```

The input earth model file should be of the same format as used by `minos_bran` in MINEOS.  The program will output the same outputs as `minos_bran` and `eigcon`, as well as ASCII files with the eigenfunctions and sensitivity kernels.