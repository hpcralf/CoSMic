# CoSMic

## Description 
A calibration-microsimulation approach to reduce uncertainty for
policy decisions on non-pharmacological interventions in the COVID-19
pandemic.
The package implements an age-structured spatial microsimulation model that
extends the Susceptible-Exposed-Infectious-Recovered (SEIR) framework.
Using an optimization approach based on subnational trends in the number of
intensive care patients, it is able to calibrate the model to the ongoing
spread of the epidemic and tries to estimate how the NPIs have affected it.
Based on these estimates the model can provide national and sub-national
forecasts for trends in the number of ICU patients and other indicators
under different scenarios regarding NPIs.

## Installation

### R-version
Since the CoSMic R-version is currently not developed up to the state of a 
proper R-package please install it using
```R
library(devtools)
devtools::install_github("hpcralf/CoSMic",build_manual=TRUE,build_vignettes=TRUE)
```

Alternatively you can check out the code with

```
git clone https://github.com/hpcralf/CoSMic.git
```
or download the source code zip-file from https://github.com/hpcralf/CoSMic/archive/main.zip.

Once having checked out the code or extracted the source-package, change to the code directory
(CoSMic in case of git-chekcout or CoSMic-main in case of zip-file download), start a R session
and execute

```R
library(devtools)
devtools::install(build_manual=TRUE,build_vignettes=TRUE)
```
### Fortran-version
Even though the Fortran model version is usable as a stand alone program it is meaningfull to
also install the R-version since the input files as well as the default set of input data can
be generated automatically using the function `convert.Rp.to.Fp` included in the R-version of
the model.

To build and install the the Fortran-version of the CoSMic model the GNU configure and build system is 
system is used. To build the Fortran-version a working MPI-Fortran-compiler is needed.

First checkout the current master branch from github

```Fortran
git clone https://github.com/hpcralf/CoSMic.git
```

Change into the CoSMic directory, create a build directory, save the CoSMic basepath to an 
environment variable and change into the build directory.

```Fortran
cd CoSMic
mkdir build
export COSMIC_PATH=$PWD
cd build
```

Set the `FC` and `FCFLAGS` environment variables to select the installed MPI-compiler and activate
OpenMP support. In case you are using a standard MPI-installation like OpenMPI along with the GNU
Fortran compiler gfortran on a Linux based system this will work like stated below.

```Fortran
export FC=mpif90
export FCFLAGS="-fopenmp -O3"
../Fortran/configure --prefix=${COSMIC_PATH}/install
```

Now build and install the Fortran-version and extend the `PATH` environment variable.

```Fortran
make
make install
export PATH=$COSMIC_PATH/install/bin:$PATH
```

## Usage
Please see

* https://github.com/hpcralf/CoSMic/blob/main/aux/CoSMic_0.11.1.0000.pdf
* https://github.com/hpcralf/CoSMic/blob/main/aux/CoSMic-Basic-Usage.html
