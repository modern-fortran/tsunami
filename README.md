# tsunami

[![GitHub issues](https://img.shields.io/github/issues/modern-fortran/tsunami.svg)](https://github.com/modern-fortran/tsunami/issues)

A shallow water equations solver. 
Companion running example for [Modern Fortran: Building Efficient Parallel Applications](https://www.manning.com/books/modern-fortran?a_aid=modernfortran&a_bid=2dc4d442).

## Getting started

### Get the code

```
git clone https://github.com/modern-fortran/tsunami
```

### Build the code

```
cd tsunami
make
```

### Set up the Python environment for graphics

```
python3 -m venv venv
source venv/bin/activate
pip install -U pip
pip install -U -r requirements.txt
```

## Organization

* Chapter 2: Linear advection, single program
* Chapter 3: Linear advection, program with procedures
* Chapter 4: Non-linear gravity wave, incorporates modules
* Chapter 7: Parallel implementation with coarrays
* Chapter 8: Incorporating derived types
* Chapter 10: Incorporating overloaded operators
* Chapter 11: Teams and events
