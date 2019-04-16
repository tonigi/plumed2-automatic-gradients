Automatic Gradient Computation for Collective Variables in PLUMED 2
========================================

This repository contains example code from the paper *How to
Differentiate Collective Variables in Free Energy Codes:
Computer-Algebra Code Generation and Automatic Differentiation*,
[10.1016/j.cpc.2018.02.017](http://dx.doi.org/10.1016/j.cpc.2018.02.017),
illustrating two approaches to automated gradient computation for
collective variables in PLUMED.

It is a fork of the PLUMED 2 repository (www.plumed.org) taken at
release v2.4.0.

Depending on how you obtained this archive, you will find either the
full distribution, or the new modules with their regtests only. The
new code is contained in submodules (directories) named
`src/curvature_codegen` (code generation approach from symbolic
expressions by SymPy) and `src/curvature_autodiff` (code
differentiation approach with the Stan Math library). Example files
with regression tests are provided in the directories
`regtest/curvature_codegen` and `regtest/curvature_autodiff`
respectively.  The rest of PLUMED 2 repository is unchanged.


If you are using GIT, the Stan Math library is referenced as a
submodule: it will be automatically downloaded if you clone the
repository with the --recursive option, or (after clone) 

	git submodule update --init --recursive

To test, on most common machines the following instructions should get
you started. After extracting the distribution:

    ./configure 
    make -j4


You may use the supplied C++ files as templates to implement your own
CVs. The modules can be enabled or disabled independently, as follows:

    ./configure --enable-modules=+curvature_codegen:-curvature_autodiff




Approach 1 - Symbolic differentiation with code generation
--------------------

The notebook generating the "core" functions calculating the gradient
is in src/curvature_codegen/sympy_codegen directory. To regenerate the
code, execute the `CurvatureCodegen.ipynb` file (you will need Sympy,
available from www.sympy.org; the easiest way to install it is via
Conda).

The example CV defines the keywords `CURVATURE_CODEGEN` and
`CURVATURE_MULTICOLVAR_CODEGEN`.

To test:

    cd regtest/curvature_codegen/rt-m2		# Or any other of the examples
    ../../../src/lib/plumed driver --plumed plumed.dat --ixyz spiral.xyz

The above test calculates the radius of curvature at several consecutive
triplets of atoms along a spiral (see COLVAR).



Approach 2 - Automatic code differentiation
--------------------

Building the `curvature_autodiff` source files requires the Stan Math
library (available at https://github.com/stan-dev/math; tested with
release 2.16). Depending on your system, you may need to adjust
Makefile paths.

The example CV defines the keyword `CURVATURE_AUTODIFF`.

To test:

    cd regtest/curvature_autodiff/rt-1        # Or any other of the examples
    ../../../src/lib/plumed driver --plumed plumed.dat --ixyz spiral.xyz


Note that `regtest/curvature_autodiff/rt-2` invokes
`CURVATURE_CODEGEN` for comparison.

