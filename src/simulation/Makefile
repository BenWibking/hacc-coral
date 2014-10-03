OBJDIR = ${HACC_OBJDIR}

TARGETS += $(OBJDIR)/libsimulation.a

all: ${TARGETS}

include ../include.mk
include ../initializer/include.mk
include ../halo_finder/include.mk
include ../dfft/include.mk

CFLAGS += ${HACC_MPI_CFLAGS}
CXXFLAGS += ${HACC_MPI_CXXFLAGS}
LDFLAGS += ${HACC_MPI_LDFLAGS}

CXXFLAGS += ${MPI_COMPILE_FLAGS} ${FFTW_CXXFLAGS}

ifeq (${MPI_NAME},mvapich2)
	CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
endif
CXXFLAGS += ${DFFT_CXXFLAGS}

CFLAGS += ${INIT_CFLAGS}
CXXFLAGS += ${INIT_CXXFLAGS}

ifndef FFTW_HOME
FFTW_HOME = $(shell dirname $(shell dirname $(shell which fftw-wisdom)))
endif
LDFLAGS  += -L${FFTW_HOME}/lib -Xlinker -rpath -Xlinker ${FFTW_HOME}/lib

ifeq ("${FFTW_MAJOR_VERSION}", "3")
CFLAGS += -DFFTW3=1
CXXFLAGS += -DFFTW3=1
endif

CFLAGS += ${HF_CFLAGS}
CXXFLAGS += ${HF_CXXFLAGS}

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: %.c | $(OBJDIR)
	${HACC_MPI_CC} ${CFLAGS} -c -o $@ $<

$(OBJDIR)/%.o: %.cxx | $(OBJDIR)
	${HACC_MPI_CXX} ${CXXFLAGS} -c -o $@ $<

SOURCES += Domain.cxx
SOURCES += MC3Options.cxx
SOURCES += SimpleTimings.cxx
SOURCES += lightcone.cxx
SOURCES += skewers.cxx
SOURCES += spline.cxx
SOURCES += TimeStepper.cxx
SOURCES += RestartIO.cxx
OBJLIST = $(SOURCES:.cxx=.o)
OBJECTS := $(addprefix $(OBJDIR)/,$(OBJLIST))
$(OBJDIR)/libsimulation.a: $(OBJDIR)/libsimulation.a($(OBJECTS))
	ranlib $@

$(OBJDIR)/cic_check: cic_check.c
	${HACC_CC} -o $@ cic_check.c

$(OBJDIR)/cic_check2: cic_check2.c
	${HACC_CC} -o $@ cic_check2.c
