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

### Switching between different tags

The code is organized in tags so you can easily go back and forth 
different stages of the app as it is taught in the book.
For example, below command will take you to the tag `2c`:

```
git checkout 2c
```

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

All the plotting code is located in the `plotting` directory.
You will need a Python build with numpy and matplotlib packages installed.
