# Phase Velocity Sensitivity Kernels

Calculates the phase velocity sensitivity kernels using eigenfunctions from [MINEOS](https://github.com/geodynamics/mineos).

## Installation

Build system is configured by [cmake](https://cmake.org/).  This will download and build MINEOS (without installing the binaries), and then build the kernels scripts.

From the `kernels` directory, run

```bash
mkdir build
cd build
cmake ..
make
```

This will produce the `kernels` executable in the `build` directory.  This has only been tested on Mac OS using `gfortran 10.2.0` and `clang 12.0.5` as the fortran compiler.

## Usage

From the `build` directory, run

```bash
./kernels <path/to/earth/model> <output/directory>
```

The input earth model file should be of the same format as used by `minos_bran` in MINEOS.  The program will output the same outputs as `minos_bran` and `eigcon`, as well as ASCII files with the eigenfunctions and sensitivity kernels.

The sensitivity kernels are saved in `<output/directory>/kernelsasc`, with each file named `X.N.L.ASC` where `X` is the mode type (e.g. `S` for spheroidal mode), `N` is the overtone (e.g. 0 for the fundamental overtone), and `L` is the angular oder.  The files columns are `radius, Kkappa, Kmu, Kalpha, Kbeta`.

A simple python plotting script is given, that takes in the `<output/directory>` as a command line argument and plots a few of the `Kbeta` kernels.

```bash
python src/plot_kernels.py <output/directory>
```
