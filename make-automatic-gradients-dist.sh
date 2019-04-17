#!/bin/bash

outzip=plumed2-automatic-gradients-dist.zip
# rm $outzip
# git submodule deinit --all
# make clean
git archive -o $outzip --prefix=plumed2-automatic-gradients-PATCH/ HEAD src/curvature* regtest/curvature* README.md
