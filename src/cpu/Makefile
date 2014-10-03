OBJDIR = ${HACC_OBJDIR}

PROGLIST = hacc_tpm #hacc_pm hacc_tpm hacc_pk
PROGS = $(addprefix $(OBJDIR)/,$(PROGLIST))

all: ${PROGS}

include ../include.mk
include ../simulation/include.mk
include ../dfft/include.mk
include ../halo_finder/include.mk
include ../initializer/include.mk

#CXXFLAGS += ${MPI_COMPILE_FLAGS} 
CXXFLAGS += ${FFTW_CXXFLAGS}

ifeq (${MPI_NAME},mvapich2)
CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
endif
CXXFLAGS += ${DFFT_CXXFLAGS}

CFLAGS += ${INIT_CFLAGS}
CXXFLAGS += ${INIT_CXXFLAGS}
LDFLAGS += ${INIT_LDFLAGS} -linit

CFLAGS += ${HF_CFLAGS}
CXXFLAGS += ${HF_CXXFLAGS}
LDFLAGS += ${HF_LDFLAGS} -lparticle -lgrid -lpartition -lBHForceTree -lhalofinder -lhalotime -lbigchunk

CFLAGS += ${SIM_CFLAGS}
CXXFLAGS += ${SIM_CXXFLAGS}
LDFLAGS += ${SIM_LDFLAGS} -lsimulation

ifndef FFTW_HOME
FFTW_HOME = $(shell dirname $(shell dirname $(shell which fftw-wisdom)))
endif
#LDFLAGS  += -L${FFTW_HOME}/lib -Xlinker -rpath -Xlinker ${FFTW_HOME}/lib
LDFLAGS += -L${FFTW_HOME}/lib

ifeq ("${FFTW_MAJOR_VERSION}", "2")
LDFLAGS += -lfftw_mpi -lfftw
else
CFLAGS += -DFFTW3=1
CXXFLAGS += -DFFTW3=1
ifeq ("${FFTW_WRAPPER}", "essl")
CFLAGS += -DESSL_FFTW=1
CXXFLAGS += -DESSL_FFTW=1
XTRALDFLAGS += -L${FFTW_HOME}/bglib64 -L${ESSL_LIBDIR} -L${XLF_LIBDIR}
XTRALDFLAGS += -lfftw3_esslbg -lesslsmpbg -lxlfmath -lxlf90_r
else
ifeq ("${FFTW_THREADING}", "omp")
CFLAGS += -DFFTW3_THREADS=1
CXXFLAGS += -DFFTW3_THREADS=1
LDFLAGS += -lfftw3_mpi -lfftw3_omp -lfftw3
else
LDFLAGS += -lfftw3_mpi -lfftw3
endif
endif
endif
LDFLAGS  += ${DFFT_LDFLAGS} -ldfft

CFLAGS += ${HACC_MPI_CFLAGS}
CXXFLAGS += ${HACC_MPI_CXXFLAGS}
LDFLAGS += ${HACC_MPI_LDFLAGS}
LDFLAGS += ${XTRALDFLAGS}

CFLAGS += -I.
CXXFLAGS += -I.

#uncomment these for cosmo output from a single rank (0)
#CFLAGS += -DSINGLE_RANK_OUTPUT
#CXXFLAGS += -DSINGLE_RANK_OUTPUT


$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: %.c | $(OBJDIR)
	${HACC_MPI_CC} ${CFLAGS} -c -o $@ $<

$(OBJDIR)/%.o: %.cxx | $(OBJDIR)
	${HACC_MPI_CXX} ${CXXFLAGS} -c -o $@ $<


$(HF_HOME)/$(OBJDIR)/%:
	make -C $(HF_HOME) $(OBJDIR)/$(@F)

$(INIT_HOME)/$(OBJDIR)/%:
	make -C $(INIT_HOME) $(OBJDIR)/$(@F)

$(SIM_HOME)/$(OBJDIR)/%:
	make -C $(SIM_HOME) $(OBJDIR)/$(@F)

$(DFFT)/$(OBJDIR)/%:
	make -C $(DFFT) $(OBJDIR)/$(@F)


driver_hires-local.cxx: ${SIM_HOME}/driver_hires.cxx
	ln -s $< $@

driver_pm-local.cxx: ${SIM_HOME}/driver_pm.cxx
	ln -s $< $@

driver_pk-local.cxx: ${SIM_HOME}/driver_pk.cxx
	ln -s $< $@

mc3-local.cxx: ${SIM_HOME}/mc3.cxx
	ln -s $< $@

MC3Extras-local.cxx: ${SIM_HOME}/MC3Extras.cxx
	ln -s $< $@


#$(OBJDIR)/hacc_pm: $(OBJDIR)/driver_pm-local.o $(OBJDIR)/mc3-local.o $(OBJDIR)/MC3Extras-local.o $(OBJDIR)/Particles.o $(SIM_HOME)/$(OBJDIR)/libsimulation.a $(INIT_HOME)/$(OBJDIR)/libinit.a $(HF_HOME)/$(OBJDIR)/libparticle.a $(HF_HOME)/$(OBJDIR)/libgrid.a $(HF_HOME)/$(OBJDIR)/libhalofinder.a $(HF_HOME)/$(OBJDIR)/libbhforcetree.a $(HF_HOME)/$(OBJDIR)/libpartition.a $(HF_HOME)/$(OBJDIR)/libhalotime.a $(HF_HOME)/$(OBJDIR)/libbigchunk.a $(DFFT)/$(OBJDIR)/libdfft.a
#	${HACC_MPI_CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

$(OBJDIR)/hacc_pm: $(OBJDIR)/driver_pm-local.o $(OBJDIR)/mc3-local.o $(OBJDIR)/MC3Extras-local.o $(OBJDIR)/Particles.o $(SIM_HOME)/$(OBJDIR)/libsimulation.a $(INIT_HOME)/$(OBJDIR)/libinit.a $(HF_HOME)/$(OBJDIR)/libparticle.a $(HF_HOME)/$(OBJDIR)/libgrid.a $(HF_HOME)/$(OBJDIR)/libhalofinder.a $(HF_HOME)/$(OBJDIR)/libBHForceTree.a $(HF_HOME)/$(OBJDIR)/libpartition.a $(HF_HOME)/$(OBJDIR)/libhalotime.a $(HF_HOME)/$(OBJDIR)/libbigchunk.a $(DFFT)/$(OBJDIR)/libdfft.a
	${HACC_MPI_CXX} ${CXXFLAGS} -o $@ $(OBJDIR)/driver_pm-local.o $(OBJDIR)/mc3-local.o $(OBJDIR)/MC3Extras-local.o $(OBJDIR)/Particles.o ${LDFLAGS}


$(OBJDIR)/hacc_pk: $(OBJDIR)/driver_pk-local.o $(OBJDIR)/mc3-local.o $(OBJDIR)/MC3Extras-local.o $(OBJDIR)/Particles.o $(SIM_HOME)/$(OBJDIR)/libsimulation.a $(INIT_HOME)/$(OBJDIR)/libinit.a $(HF_HOME)/$(OBJDIR)/libparticle.a $(HF_HOME)/$(OBJDIR)/libgrid.a $(HF_HOME)/$(OBJDIR)/libhalofinder.a $(HF_HOME)/$(OBJDIR)/libBHForceTree.a $(HF_HOME)/$(OBJDIR)/libpartition.a $(HF_HOME)/$(OBJDIR)/libhalotime.a $(HF_HOME)/$(OBJDIR)/libbigchunk.a $(DFFT)/$(OBJDIR)/libdfft.a
	${HACC_MPI_CXX} ${CXXFLAGS} -o $@ $(OBJDIR)/driver_pk-local.o $(OBJDIR)/mc3-local.o $(OBJDIR)/MC3Extras-local.o $(OBJDIR)/Particles.o ${LDFLAGS}


#$(OBJDIR)/hacc_tpm: $(OBJDIR)/driver_hires-local.o $(OBJDIR)/mc3-local.o $(OBJDIR)/MC3Extras-local.o $(OBJDIR)/Particles.o $(SIM_HOME)/$(OBJDIR)/libsimulation.a $(INIT_HOME)/$(OBJDIR)/libinit.a $(HF_HOME)/$(OBJDIR)/libparticle.a $(HF_HOME)/$(OBJDIR)/libgrid.a $(HF_HOME)/$(OBJDIR)/libhalofinder.a $(HF_HOME)/$(OBJDIR)/libbhforcetree.a $(HF_HOME)/$(OBJDIR)/libpartition.a $(HF_HOME)/$(OBJDIR)/libhalotime.a $(HF_HOME)/$(OBJDIR)/libbigchunk.a $(DFFT)/$(OBJDIR)/libdfft.a
#	${HACC_MPI_CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

$(OBJDIR)/hacc_tpm: $(OBJDIR)/driver_hires-local.o $(OBJDIR)/mc3-local.o $(OBJDIR)/MC3Extras-local.o $(OBJDIR)/Particles.o $(SIM_HOME)/$(OBJDIR)/libsimulation.a $(INIT_HOME)/$(OBJDIR)/libinit.a $(HF_HOME)/$(OBJDIR)/libparticle.a $(HF_HOME)/$(OBJDIR)/libgrid.a $(HF_HOME)/$(OBJDIR)/libhalofinder.a $(HF_HOME)/$(OBJDIR)/libBHForceTree.a $(HF_HOME)/$(OBJDIR)/libpartition.a $(HF_HOME)/$(OBJDIR)/libhalotime.a $(HF_HOME)/$(OBJDIR)/libbigchunk.a $(DFFT)/$(OBJDIR)/libdfft.a
	${HACC_MPI_CXX} ${CXXFLAGS} -o $@ $(OBJDIR)/driver_hires-local.o $(OBJDIR)/mc3-local.o $(OBJDIR)/MC3Extras-local.o $(OBJDIR)/Particles.o ${LDFLAGS}


clean::
	rm -rf $(OBJDIR); rm -f *~ driver_hires-local.cxx driver_pm-local.cxx mc3-local.cxx MC3Extras-local.cxx driver_pk-local.cxx viz_rho-local.cxx
