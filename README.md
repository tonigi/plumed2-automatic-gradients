Automatic Gradient Computation for Collective Variables in PLUMED 2
========================================

This repository contains example code from the paper *How to
Differentiate Collective Variables in Free Energy Codes:
Computer-Algebra Code Generation and Automatic Differentiation*,
illustrating two approaches to automated gradient computation for
collective variables in PLUMED.

It is a fork of the PLUMED 2 repository (www.plumed.org) taken at
release v2.4.0.

The new code is contained in submodules (directories) named
`src/curvature_codegen` (code generation approach from symbolic
expressions by SymPy) and `src/curvature_autodiff` (code
differentiation approach by the Stan Math library). You may use the
supplied C++ files as templates to implement your own CVs. The modules
can be enabled or disabled independently.

Example files with regression tests are provided in the directories
`regtest/curvature_codegen` and `regtest/curvature_autodiff`
respectively.  The rest of PLUMED 2 repository is unchanged.

To test, on most common machines the following instructions should get
you started. After extracting the distribution:

    ./configure 
    make -j4

You can enable and disable the modules independently, as follows:

    ./configure --enable-modules=+curvature_codegen:-curvature_autodiff



Approach 1 - Symbolic differentiation with code generation
--------------------

The notebook generating the "core" functions calculating the gradient
is in src/curvature_codegen/sympy_codegen directory. To regenerate the
code, execute the `CurvatureCodegen.ipynb` file (you will need Sympy,
available from www.sympy.org; the easiest way to install it is via
Conda).

To test:

    cd regtest/curvature_codegen/rt-m2		# Or any other of the examples
    ../../../src/lib/plumed driver --plumed plumed.dat --ixyz spiral.xyz

The above test calculates the radius of curvature at several consecutive
triplets of atoms along a spiral (see COLVAR).



Approach 2 - Automatic code differentiation
--------------------

Building the `curvature_autodiff` source files requires the Stan Math
library (distributed with the source). Depending on your system, you
may need to adjust Makefile paths.

To test:

    ./configure 
    make -j4
    cd regtest/curvature_autodiff/rt-1        # Or any other of the examples
    ../../../src/lib/plumed driver --plumed plumed.dat --ixyz spiral.xyz

(Note: if you are using GIT, the Stan Math library,
http://mc-stan.org/users/interfaces/math, is referenced as a
submodule; you may need to clone this repository with the --recursive
option.)






