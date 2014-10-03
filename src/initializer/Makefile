###### Makefile for the initializer code
##
## FFTW libraries need to be explicitely linked.
##	Code uses FFTW 2 or FFTW 3version, and 
##      loading modules 
##      sets up FFTW_HOME environment variable.
##
## FLAGS:
##      -DLONG_INTEGER
##         If you want to use long integers -- this 
##         may be needed for particle IDs
##
##      -DDOUBLE_REAL
##         If you want positions and velocities in 
##         double precision
##
##      -DNOFFTW
##         For tests that do not require FFT, 
##         above flag can be passed. The code can 
##         than be run on a machine without FFTW.
##
##      -DFFTW2
##         If FFTW version 2 is being used, pass this 
##         flag to the compiler
##
##      -DFFTW3
##         If FFTW version 3 is being used, pass this 
##         flag to the compiler
##
##      -DTESTING 
##         Non-production flag -- makes code perform 
##         several tests, most importantly
##         it checks for reality of the density field.
##
##      -DH5PART
##         For output in HDF5 format
##
##      -DPARALLEL_IO
##         For using parallel IO in HDF5

default: all

OBJDIR = ${HACC_OBJDIR}

include ../include.mk
include ../dfft/include.mk
include ../halo_finder/include.mk

PROGLIST = init init2 libinit.a
PROGS = $(addprefix $(OBJDIR)/,$(PROGLIST))

CPP = ${HACC_MPI_CXX}
CC = ${HACC_MPI_CC}

#CPPFLAGS = -O3 -funroll-loops -DFFTW3 -DTESTING
#CFLAGS = -O3 -funroll-loops -std=c99
CPPFLAGS += -DFFTW3 -DTESTING
CFLAGS += 

CPPFLAGS += ${FFTW_CXXFLAGS}
CPPFLAGS += ${DFFT_CXXFLAGS}

#for initial P(k) check
CPPFLAGS += ${HF_CXXFLAGS}

ifndef FFTW_HOME
FFTW_HOME = $(shell dirname $(shell dirname $(shell which fftw-wisdom)))
endif
#LDFLAGS += -L${FFTW_HOME}/lib -Xlinker -rpath -Xlinker ${FFTW_HOME}/lib
LDFLAGS += -L${FFTW_HOME}/lib
FFTW_FLAGS = -I$(FFTW_HOME)/include

LDFLAGS += -L$(HF_HOME)/$(OBJDIR)
LDFLAGS += -ldfft -lbigchunk

## For FFTW2 uncomment this line:
##FFTWLIBS= -lfftw_mpi  -lrfftw_mpi -lrfftw -lfftw
## For FFTW3 uncomment this line:
ifeq ("${FFTW_WRAPPER}", "essl")
CFLAGS += -DESSL_FFTW=1
CPPFLAGS += -DESSL_FFTW=1
XTRAFFTWLIBS += -L${FFTW_HOME}/bglib64 -L${ESSL_LIBDIR} -L${XLF_LIBDIR}
XTRAFFTWLIBS += -lfftw3_esslbg -lesslsmpbg -lxlfmath -lxlf90_r
else
ifeq ("${FFTW_THREADING}", "omp")
CFLAGS += -DFFTW3_THREADS=1
CPPFLAGS += -DFFTW3_THREADS=1
FFTWLIBS = -lfftw3_mpi -lfftw3_omp -lfftw3
else
FFTWLIBS = -lfftw3_mpi -lfftw3
endif
endif

## For HDF5 support uncomment these:
##LDFLAGS = -L$(FFTW_HOME)/lib -L$(H5Part)/src -L($HDF5HOME)/lib
##FFTW_FLAGS = -I$(FFTW_HOME)/include -I($H5Part)/src -I($HDF5HOME)/include
##FFTWLIBS = -lfftw3_mpi -lfftw3 -lhdf5 -lH5Part -lz

CFLAGS += ${HACC_MPI_CFLAGS}
CPPFLAGS += ${HACC_MPI_CXXFLAGS}
LDFLAGS += ${HACC_MPI_LDFLAGS}

LDFLAGS += ${XTRAFFTWLIBS}

########## SHOULD NOT CHANGE ANYTHING BELOW ##############################

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: %.c | $(OBJDIR)
	${HACC_MPI_CC} ${CFLAGS} -c -o $@ $<

$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	${HACC_MPI_CXX} ${CPPFLAGS} -c -o $@ $<

CPP_SOURCES = Cosmology.cpp MT_Random.cpp DataBase.cpp \
              Parallelization.cpp Initializer.cpp

C_SOURCES = PerfMon.c

INCLUDE = Cosmology.h DataBase.h MT_Random.h Parallelization.h \
          PerfMon.h TypesAndDefs.h Initializer.h \
          Basedata.h 

CPP_OBJLIST = $(CPP_SOURCES:.cpp=.o)
CPP_OBJ := $(addprefix $(OBJDIR)/,$(CPP_OBJLIST))

C_OBJLIST = $(C_SOURCES:.c=.o)
C_OBJ := $(addprefix $(OBJDIR)/,$(C_OBJLIST))


all: $(C_SOURCES) $(CPP_SOURCES) $(PROGS)
	rm -rf $(OBJDIR)/*.o

tarball:
	tar -zcf source.tar.gz $(C_SOURCES) $(CPP_SOURCES) \
                               $(INCLUDE) README Makefile indat cmb.tf


$(OBJDIR)/main.o: main.cpp $(HF_HOME)/Definition.h | $(OBJDIR)
	$(CPP) -c $(HF_TYPE_FLAGS) -I$(HF_HOME) $(CPPFLAGS) $(FFTW_FLAGS) $< -o $@

$(OBJDIR)/init: $(C_OBJ) $(CPP_OBJ) $(OBJDIR)/main.o
	make -C ${DFFT}
	make -C $(HF_HOME) $(OBJDIR)/libbigchunk.a
	$(CPP) $(CPPFLAGS) $(C_OBJ) $(CPP_OBJ) $(OBJDIR)/main.o -o $@ ${DFFT_LDFLAGS} $(FFTWLIBS) $(LDFLAGS)

$(OBJDIR)/init2.o: init2.cpp $(HF_HOME)/Definition.h $(HF_HOME)/InitialExchange.h | $(OBJDIR)
	$(CPP) -c $(HF_TYPE_FLAGS) -I$(HF_HOME) $(CPPFLAGS) $(FFTW_FLAGS) $< -o $@

$(OBJDIR)/init2: $(C_OBJ) $(CPP_OBJ) $(HF_HOME)/Definition.h $(HF_HOME)/InitialExchange.h $(OBJDIR)/init2.o
	make -C ${DFFT}
	make -C $(HF_HOME) $(OBJDIR)/libparticle.a $(OBJDIR)/libpartition.a $(OBJDIR)/libbigchunk.a $(OBJDIR)/libbigchunk.a
	$(CPP) $(CPPFLAGS) $(C_OBJ) $(CPP_OBJ) $(OBJDIR)/init2.o -o $@ ${DFFT_LDFLAGS} $(FFTWLIBS) -L$(HF_HOME)/$(OBJDIR) -lparticle -lpartition $(LDFLAGS)

#.cpp.o: $(INCLUDE)
#	$(CPP) -c $(CPPFLAGS) $(FFTW_FLAGS) $< -o $@
#
#.c.o: $(INCLUDE)
#	$(CC) -c $(CFLAGS) $(FFTW_FLAGS) $< -o $@

$(OBJDIR)/libinit.a: $(OBJDIR)/libinit.a($(C_OBJ) $(CPP_OBJ))
	ranlib $@

