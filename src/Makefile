OPT = -O3
#OPT = -O3 -m64 -static
#OPT = -Wall -g 

CFLAGS = -Wall -pedantic
CC = gcc
LIBS = -lm -lOpenCL

all: target

target: structure 


#valgrind: OPT = -g -O1
#valgrind: CFLAGS += -g
#valgrind: clean
#valgrind: target

prof: CFLAGS += -g
prof: CFLAGS += -pg
prof: target

#debug: OPT = -g -O1
#debug: CFLAGS += -g
#debug: clean
#debug: target

UPDATES = ForwardAndBackward.o UpdateQ.o UpdateZ.o UpdateP.o UpdateLocPrior.o UpdateAlpha.o UpdateGeno.o UpdateEpsilon.o UpdateFst.o UpdateLambda.o

structure: structure.o params.o datain.o output.o ran.o mymath.o  $(UPDATES) Init.o Util.o Kernels.o
	$(CC) -o structure structure.o params.o datain.o output.o ran.o  mymath.o $(UPDATES) Init.o Util.o Kernels.o $(OPT) $(LIBS) $(CFLAGS)

#STRAT: STRAT.o params.o datain.o ran.o mymath.o
#	$(CC) -o STRAT STRAT.o params.o datain.o ran.o mymath.o $(OPT) $(LIBS)

#STRAT.o: STRAT.c
#	$(CC) -c STRAT.c  $(OPT) $(CFLAGS)

structure.o: structure.c
	$(CC) -c structure.c $(OPT) $(CFLAGS)

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

KernelsTest: Kernels.c
	$(CC) -o Kernels Kernels.c  -Wall -pedantic -lOpenCL

clean:
	@rm -f *.o structure

test: structure
	./structure -m testdata/mainparams -i testdata/testdata1.str -o testdata1.out

callgrind: prof structure
	valgrind --tool=callgrind ./structure -m testdata/mainparams -i testdata/testdata1.str -o testdata1.out

cleanTest:
	rm -f testdata1.out_f

