# tsunami

[![Build Status](https://travis-ci.org/fortran-in-action/tsunami.svg?branch=master)](https://travis-ci.org/fortran-in-action/tsunami)
[![GitHub issues](https://img.shields.io/github/issues/fortran-in-action/tsunami.svg)](https://github.com/fortran-in-action/tsunami/issues)

A shallow water equations solver. Companion running example 
for the upcoming book Fortran in Action by Manning Publications.

## Getting started

### Getting the code

```
git clone https://github.com/fortran-in-action/tsunami
```

### Switching to different tags/snapshots

TODO

### Building the model

Dependencies:

* `gfortran`
* `cmake`

```
mkdir build
cd build
cmake ..
make
```

The executable will be built in the `build/bin` directory.

### Running the model

Inside the `build/bin` directory, type:

```
./tsunami
```

### Plotting the results

TODO
