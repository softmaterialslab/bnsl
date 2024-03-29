# This is a makefile.
# Use option -p in CC for profiling with gprof

PROG = postprocessor
BIN = ../../bin/
bindir = bin
OBJ = main.o functions.o
CCF = CC

# nanoHUB flags. 
nanoHUBCC = g++ -O3 -g -Wall -std=c++11
nanoHUBLFLAG = -lm -L${BOOST_LIBDIR} -lboost_program_options -lboost_filesystem -lboost_system
nanoHUBCFLAG = -I${BOOST_INCDIR}
nanoHUBOFLAG = -o
# BigRed2 flags. 
BigRed2CC = CC -O3 -g -std=c++11
#BigRed2LFLAG = -lgsl -lgslcblas -lm -lboost_program_options -lboost_filesystem -lboost_system
BigRed2LFLAG = -L/N/soft/cle7/boost/gnu/1.73.0/lib -L/N/soft/cle7/gsl/2.6/lib -lgsl -lgslcblas -lm -lboost_filesystem -lboost_system -lboost_program_options
BigRed2CFLAG = -c
BigRed2OFLAG = -o
# General purpose flags.
CC = g++ -O3 -g -Wall -std=c++11
LFLAG = -lm -L${BOOST_LIBDIR} -lboost_program_options -lboost_filesystem -lboost_system
CFLAG = -c
OFLAG = -o

all: $(PROG)

install: all
	@echo "Installing $(PROG) into $(bindir) directory"; cp -f $(PROG) $(BIN)

nanoHUB-install:
	. /etc/environ.sh; use -e -r boost-1.62.0-mpich2-1.3-gnu-4.7.2; make CCF=nanoHUB all
	@echo "Installing $(PROG) into $(bindir) directory on nanohub"; cp -f $(PROG) $(BIN)

cluster-install:
	module load PrgEnv-gnu && module load boost/gnu; make CCF=BigRed2 all
	@echo "Installing $(PROG) into $(bindir) directory on the computing cluster"; cp -f $(PROG) $(BIN)

$(PROG) : $(OBJ)
ifeq ($(CCF),BigRed2)	
	$(BigRed2CC) $(BigRed2OFLAG) $(PROG) $(OBJ) $(BigRed2LFLAG)
%.o : %.cpp
	$(BigRed2CC) -c $(BigRed2CFLAG) $< -o $@
else ifeq ($(CCF),nanoHUB)
	$(nanoHUBCC) $(nanoHUBOFLAG) $(PROG) $(OBJ) $(nanoHUBLFLAG)
%.o : %.cpp
	$(nanoHUBCC) -c $(nanoHUBCFLAG) $< -o $@
else
	$(CC) $(OFLAG) $(PROG) $(OBJ) $(LFLAG)
%.o : %.cpp
	$(CC) -c $(CFLAG) $< -o $@	
endif

.PHONY: all install nanoHUB-install cluster-install
