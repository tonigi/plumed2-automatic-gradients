CFLAGS=-O3
LOADLIBES=-lm

default: test_curvature_plain test_curvature_cse

# To regenerate code:
codegen:
	jupyter nbconvert --execute --ExecutePreprocessor.timeout=600 CurvatureCodegen.ipynb


clean:
	-rm *.o


test_curvature_plain: test_curvature.o curvature_codegen.o
	$(CC) $(LOADLIBES) -o $@ $^

test_curvature_cse: test_curvature.o curvature_codegen_cse.o
	$(CC) $(LOADLIBES) -o $@ $^


# This is one way to query hardware counters for double precision
# operations performed by a command (Linux only). The event code may
# differ by CPU architecture. See also ...
#
# http://www.bnikolic.co.uk/blog/hpc-howto-measure-flops.html
# http://icl.cs.utk.edu/projects/papi/wiki/PAPITopics:SandyFlops#Floating_Point_Flavors

results: default
	-rm *.results
	for i in `seq 10`; do perf stat -o plain.results --append -e r538010 -e r531010 -e r530211 ./test_curvature_plain; done 
	for i in `seq 10`; do perf stat -o cse.results --append   -e r538010 -e r531010 -e r530211 ./test_curvature_cse; done 


