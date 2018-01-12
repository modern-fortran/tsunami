[![Build Status](https://travis-ci.org/fortran-in-action/tsunami.svg?branch=master)](https://travis-ci.org/fortran-in-action/tsunami)
[![GitHub issues](https://img.shields.io/github/issues/fortran-in-action/tsunami.svg)](https://github.com/fortran-in-action/tsunami/issues)


# tsunami


A shallow water equations solver. Companion running example 
for the upcoming book Fortran in Action by Manning Publications.

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
For example, below command will take you to the tag `2c`:

```
git checkout 2c
```

### Building the model

Dependencies:

```sh
apt install gfortran cmake make
```

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
