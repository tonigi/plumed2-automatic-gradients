CFLAGS=-O3
LOADLIBES=-lm

default: test_curvature_plain test_curvature_cse

clean:
	-rm *.o

# --------

# To regenerate code:
# curvature.c curvature.h: CurvatureCodegen.ipynb
#     jupyter nbconvert --execute CurvatureCodegen.ipynb


test_curvature_plain: test_curvature.o curvature_codegen.o
	$(CC) $(LOADLIBES) -o $@ $^

test_curvature_cse: test_curvature.o curvature_codegen_cse.o
	$(CC) $(LOADLIBES) -o $@ $^



