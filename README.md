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
Since CoSMic is currently not developed up to the state of a proper R-package please install
it using
```{r}
library(devtools)
devtools::install_github("hpcralf/CoSMic")
```

Alternatively you can check out the code with

```
git clone https://github.com/hpcralf/CoSMic.git
```
or download the source code zip-file from https://github.com/hpcralf/CoSMic/archive/main.zip.

Once having checked out the code or extracted the source-package change to the code directory
(CoSMic in case of git-chekcout or CoSMic-main case of zip-file download), start a R session
and execute

```R
library(devtools)
devtools::install()
```

## Usage
Please see

* https://github.com/hpcralf/CoSMic/blob/main/aux/CoSMic_0.11.1.0000.pdf
* https://github.com/hpcralf/CoSMic/blob/main/aux/CoSMic-Basic-Usage.html