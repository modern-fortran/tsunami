# tsunami

[![GitHub issues](https://img.shields.io/github/issues/modern-fortran/tsunami.svg)](https://github.com/modern-fortran/tsunami/issues)

A shallow water equations solver. 
Companion running example for [Modern Fortran: Building Efficient Parallel Applications](https://www.manning.com/books/modern-fortran?a_aid=modernfortran&a_bid=2dc4d442).

## Organization

* **Chapter 2**: We implement our first working solver which solves for linear advection of a wave in one dimension.
First encounter with declaration, variables, loops, branches, arithmetic, and writing to console.
* **Chapter 3**: We refactor our program from Chapter 2 to use procedures -- a function and a subroutine.
* **Chapter 4**: We refactor our program from Chapter 3 to define the procedures in external modules. 
We use this opportunity to augment the simulator to solve for non-linear gravity waves.
* **Chapter 7**: We parallelize the program from Chapter 4 using coarrays and observe the speed up.
* **Chapter 8**: We refactor our program from Chapter 7 to model our physical quantities (water height and velocity)
using a derived type.
* **Chapter 9**: In this chapter, we focus on the `Field` derived type and implement the common arithmetic operations
as type-bound methods. 
* **Chapter 10**: We continue the work on the code from Chapter 9, and overload the assignment operator to 
automatically synchronize the data across parallel images on every field assignment.
* **Chapter 12**: In the final chapter, we revisit the parallel code from Chapter 10 and explore how Fortran 2018
Teams, Events, and Collectives can be used for some more advanced parallel patterns.

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
