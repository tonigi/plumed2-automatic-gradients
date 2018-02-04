Automatic reverse-mode differentiation for collective variables
--------------

Toni Giorgino


This directory shows how to use the
[Stan Math library](http://mc-stan.org/users/interfaces/math) to
implement auto-differentiating collective variable functions.  The
`STANMATHDEMO` colvar can be copied, modified and used as a template.

Automatic reverse-mode code differentiation is one of the two
approaches described in the paper

> T. Giorgino. Automatic differentiation of collective variables in
> PLUMED2. CPC (submitted)

...and complements it. The other approach, automatic symbolic
differentiation, is shown in the `curvature_codegen` module.

Note that the Stan Math library is included here as a git
submodule. Make sure you have submodules enabled. The easiest way is
to use the `git clone --recursive` when you first clone PLUMED2's
repository.

