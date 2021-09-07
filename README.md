# Examples for QSF library

## Introduction
This repository hosts example projects for the QSF header-only C++ library. 

QSF is a project aimed at providing *efficient* and *general* numerical tools for the field of Strong Field Laser Physics.

It's modular design based extensively on C++ templates and modern C++17 features, allows for creating clear and understandable projects keeping the **physics** first.

The examples include three-electron reduced-dimensionality models (`examples/nitrogen-3e`) and some simpler projects.

## Requirements

* compiler with C++17 support 
  > clang, gcc 9.1 or higher recommended
* cmake 3.8
* mpi library
  > e.g. openmpi, mpich
* fftw3 library
  > note: intel-provided fftw libraries might cause problems at the moment

## Usage

Download the repository recursively and navigate to the directory using:

```bash
git clone https://github.com/masteranza/QSF-examples.git --recurse-submodules 
cd QSF-examples
```

configure the project using one of the provided scripts (depending on the compiler), e.g:

```bash
./configure-gcc11
```
run any of the example projects (e.g. template-ini) using the `run.sh` script:

```bash
./run template-ini
```

You can inspect the scripts together with the `CMakeLists.txt` file in case of problems or if you'd like to learn more.

## Starting a new project

You can use any of the existing projects gathered in the `examples` folder as the starting point for your own projects.

It is advised to keep the directory name and the main `.cpp` filename identical.

There are three ways to provide configurations to a project:
- command line options (handled by `cxxopts` library, see nitrogen-3e)
- ini file (called `project.ini`) located next to the source file
- hardcoded in the `.cpp` file

Examples of each can be found in the `examples directory`.

