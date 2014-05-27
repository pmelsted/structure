OPT = -O3
#OPT = -O3 -m64 -static
#OPT = -Wall -g 

CFLAGS = -Wall -pedantic
CC = gcc
LIBS = -lm

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

structure: structure.o params.o datain.o output.o ran.o mymath.o ForwardAndBackward.o UpdateQ.o UpdateZ.o UpdateP.o UpdateLocPrior.o 
	$(CC) -o structure structure.o params.o datain.o output.o ran.o  mymath.o ForwardAndBackward.o UpdateQ.o UpdateZ.o UpdateP.o UpdateLocPrior.o $(OPT) $(LIBS) $(CFLAGS)

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

UpdateQ.o: UpdateQ.c
	$(CC) -c UpdateQ.c $(OPT) $(CFLAGS)

UpdateZ.o: UpdateZ.c
	$(CC) -c UpdateZ.c $(OPT) $(CFLAGS)

UpdateP.o: UpdateP.c
	$(CC) -c UpdateP.c $(OPT) $(CFLAGS)

UpdateLocPrior.o: UpdateLocPrior.c
	$(CC) -c UpdateLocPrior.c $(OPT) $(CFLAGS)

ForwardAndBackward.o: ForwardAndBackward.c
	$(CC) -c ForwardAndBackward.c $(OPT) $(CFLAGS)
clean:
	@rm -f *.o structure

test: structure
	./structure -m testdata/mainparams -i testdata/testdata1.str -o testdata1.out

callgrind: prof structure
	valgrind --tool=callgrind ./structure -m testdata/mainparams -i testdata/testdata1.str -o testdata1.out

cleanTest:
	rm -f testdata1.out_f

