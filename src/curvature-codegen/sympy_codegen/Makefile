CFLAGS=-O3
LOADLIBES=-lm

default: test_curvature

# --------

# To regenerate code:
# curvature.c curvature.h: CurvatureCodegen.ipynb
#     jupyter nbconvert --execute CurvatureCodegen.ipynb


test_curvature: test_curvature.o curvature_codegen.o

