OPT = -O3

#OPT = -O3 -m64 -static
#OPT = -Wall -g

CFLAGS = -Wall -pedantic -g
CC = clang
LIBS = -lm -lOpenCL

all: target

target: structure


#valgrind: OPT = -g -O1
#valgrind: CFLAGS += -g
#valgrind: clean
#valgrind: target

#prof: CFLAGS += -pg
prof: CC = gcc
prof: target


#debug: OPT = -g -O1
#debug: CFLAGS += -g
#debug: clean
#debug: target

UPDATES = ForwardAndBackward.o UpdateQ.o UpdateZ.o UpdateP.o UpdateLocPrior.o UpdateAlpha.o UpdateGeno.o UpdateEpsilon.o UpdateFst.o UpdateLambda.o

structure: structure.o params.o datain.o output.o ran.o mymath.o debug.o $(UPDATES) Init.o Util.o Kernels.o randGen.o CalcLike.o
	$(CC) -o structure structure.o CalcLike.o params.o datain.o output.o ran.o debug.o mymath.o $(UPDATES) Init.o Util.o Kernels.o randGen.o $(OPT) $(LIBS) $(CFLAGS)

#STRAT: STRAT.o params.o datain.o ran.o mymath.o
#	$(CC) -o STRAT STRAT.o params.o datain.o ran.o mymath.o $(OPT) $(LIBS)

#STRAT.o: STRAT.c
#	$(CC) -c STRAT.c  $(OPT) $(CFLAGS)

structure.o: structure.c
	$(CC) -c structure.c $(OPT) $(CFLAGS)

debug.o: debug.c
	$(CC) -c debug.c $(OPT) $(CFLAGS)

CalcLike.o: CalcLike.c
	$(CC) -c CalcLike.c $(OPT) $(CFLAGS)

output.o: output.c
	$(CC) -c output.c $(OPT) $(CFLAGS)

datain.o: datain.c
	$(CC) -c datain.c $(OPT) $(CFLAGS)

params.o: params.c
	$(CC) -c params.c $(OPT) $(CFLAGS)

ran.o: ran.c
	$(CC) -c ran.c $(OPT) $(CFLAGS)

mymath.o: mymath.c
	$(CC) -c mymath.c $(OPT) $(CFLAGS)

randGen.o: randGen.c
	$(CC) -c randGen.c $(OPT) $(CFLAGS)

UpdateQ.o: Updates/UpdateQ.c
	$(CC) -c Updates/UpdateQ.c $(OPT) $(CFLAGS)

UpdateZ.o: Updates/UpdateZ.c
	$(CC) -c Updates/UpdateZ.c $(OPT) $(CFLAGS)

UpdateP.o: Updates/UpdateP.c
	$(CC) -c Updates/UpdateP.c $(OPT) $(CFLAGS)

UpdateLocPrior.o: Updates/UpdateLocPrior.c
	$(CC) -c Updates/UpdateLocPrior.c $(OPT) $(CFLAGS)

UpdateAlpha.o: Updates/UpdateAlpha.c
	$(CC) -c Updates/UpdateAlpha.c $(OPT) $(CFLAGS)

UpdateGeno.o: Updates/UpdateGeno.c
	$(CC) -c Updates/UpdateGeno.c $(OPT) $(CFLAGS)

UpdateEpsilon.o: Updates/UpdateEpsilon.c
	$(CC) -c Updates/UpdateEpsilon.c $(OPT) $(CFLAGS)

UpdateFst.o: Updates/UpdateFst.c
	$(CC) -c Updates/UpdateFst.c $(OPT) $(CFLAGS)

UpdateLambda.o: Updates/UpdateLambda.c
	$(CC) -c Updates/UpdateLambda.c $(OPT) $(CFLAGS)

ForwardAndBackward.o: Updates/ForwardAndBackward.c
	$(CC) -c Updates/ForwardAndBackward.c $(OPT) $(CFLAGS)

Init.o: Init.c
	$(CC) -c Init.c $(OPT) $(CFLAGS)

Util.o: Util.c
	$(CC) -c Util.c $(OPT) $(CFLAGS)

Kernels.o: Kernels.c
	$(CC) -c Kernels.c $(OPT) $(CFLAGS)

KernelTest: Kernels.c
	$(CC) -o KernelTest Kernels.c  -Wall -pedantic -lOpenCL

clean:
	@rm -f *.o structure

test: structure
	CL_LOG_ERRORS=stdout /usr/bin/time --verbose ./structure -m testdata/mainparams
errtest: structure
	CL_LOG_ERRORS=stdout /usr/bin/time --verbose ./structure -m testdata/mainparams -D 1407409831
errtest2: structure
	CL_LOG_ERRORS=stdout /usr/bin/time --verbose ./structure -m testdata/mainparams -D 1407500105
errtest3: structure
	CL_LOG_ERRORS=stdout /usr/bin/time --verbose ./structure -m testdata/mainparams -D 1407500453
errtest4: structure
	CL_LOG_ERRORS=stdout /usr/bin/time --verbose ./structure -m testdata/mainparams -D 1407501403
goodtest: structure
	CL_LOG_ERRORS=stdout /usr/bin/time --verbose ./structure -m testdata/mainparams -D 1407420643
cerrtest: structure
	CL_LOG_ERRORS=stdout /usr/bin/time --verbose ./cpustructure -m testdata/mainparams -D 1407409831
wcachetest: structure
	CUDA_CACHE_DISABLE=0 /usr/bin/time ./structure -m testdata/mainparams

bigtest:
	/usr/bin/time --verbose ./structure -m bigdata/mainparams
cbigtest:
	/usr/bin/time --verbose ./cpustructure -m bigdata/mainparams

callgrind: prof structure
	CUDA_CACHE_DISABLE=0 valgrind --tool=callgrind ./structure -m testdata/mainparams -i testdata/testdata1.str -o testdata1.out

cleanTest:
	rm -f testdata1.out_f

gpuprof: prof structure
	ltpv "./structure -m testdata/mainparams -i testdata/testdata1.str -o testdata1.out"

beautify:
	astyle -n --style=kr -xC79 -c -j -Y -xw -R "*.c" "*.h" "*.cl"
