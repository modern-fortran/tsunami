# tsunami

[![Build Status](https://travis-ci.org/fortran-in-action/tsunami.svg?branch=master)](https://travis-ci.org/modern-fortran/tsunami)
[![GitHub issues](https://img.shields.io/github/issues/fortran-in-action/tsunami.svg)](https://github.com/modern-fortran/tsunami/issues)

A shallow water equations solver. Companion running example 
for [Modern Fortran: Building Efficient Parallel Applications](https://www.manning.com/books/modern-fortran?a_aid=modernfortran&a_bid=2dc4d442).

## Try it in the cloud

Click on the button below to create a free cloud server for up to 2 hours:

[![Dply](https://dply.co/b.svg)](https://dply.co/b/lHbdM5wp)

You will need a Github account with an SSH key associated with it.

The cloud-init script creates a `fortran` user, installs all the 
prerequisites such as the GNU Fortran compiler, OpenMPI, and OpenCoarrays.
It also clones the listings and tsunami repositories in the home directory.

Log in with ssh:

```
ssh fortran@ip-address
```

**Note: it may take up to several minutes before you can log in
due to server configuring and installing packages.**

This is a passwordless log-in that will use the SSH key that you
select in the server configuration step.
Once you log in as `fortran` you will have the following tools available:

* `gfortran`: The GNU Fortran compiler; use it to compile serial programs.
* `caf`: The coarray compile script; use it to compile parallel coarray programs.
* `cafrun`: The coarray execute script; use it to run parallel coarray programs.
* `mpif90`: The MPI compile script; use it to compile parallel MPI programs.
* `mpiexec`: The MPI execute script; use it to run parallel MPI programs.

## Getting the code

To run the model on your local machine, get the code using git:

```
git clone https://github.com/fortran-in-action/tsunami
```

### Switching between different tags

The code is organized in tags so you can easily go back and forth 
different stages of the app as it is taught in the book.
For example, below command will take you to the tag `3b`:

```
git checkout 3b
```

## Building the model

### System dependencies

On Debian-based systems:

```
sudo apt install gfortran cmake make
```

On Redhat-based systems:

```
sudo yum install gfortran cmake make
```

or

```
sudo dnf install gfortran cmake make
```

### Building OpenCoarrays

To build OpenCoarrays, follow the instruction in Appendix A of the book.
OpenCoarrays are required to build the parallel version of the model.

### Building tsunami

```
mkdir build
cd build
FC=caf cmake ..
make
```

The executable will be built in the `build/bin` directory.

## Running the model

Inside the `build/bin` directory, type:

```
./tsunami
```

## Plotting the results

All the plotting code is located in the `plotting` directory.
You will need a Python build (either 2.7 or 3.x is fine) 
with numpy and matplotlib packages installed.

If you're setting up Python from scratch, I recomment using Python 3
and creating a dedicated virtual environment. 
For example, in `tsunami/plotting` directory:

```
python3 -m venv venv            # create new python env
source venv/bin/activate        # activate python env
pip install -U pip              # upgrade installer
pip install -r requirements.txt # install dependencies
```
