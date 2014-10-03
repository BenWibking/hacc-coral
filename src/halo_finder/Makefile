OBJDIR = ${HACC_OBJDIR}
#SOBJDIR = ${HACC_OBJDIR}_serial


PROGS    = \
	$(OBJDIR)/HaloFinder \
	$(OBJDIR)/HaloTest \
	$(OBJDIR)/HaloTestP \
	$(OBJDIR)/SubHaloTest \
	$(OBJDIR)/ForceTestP \
	$(OBJDIR)/BinHalos \
	$(OBJDIR)/OneToN \
	$(OBJDIR)/NToN \
	$(OBJDIR)/Gadget2ToCosmo \
	$(OBJDIR)/CosmoToGadget2 \
	$(OBJDIR)/GridTestP \
	$(OBJDIR)/InitTestP \
	$(OBJDIR)/ForceSpeed \
	$(OBJDIR)/ForceRandom \
	$(OBJDIR)/ForceTreeTest


LIBS	= \
	$(OBJDIR)/libpartition.a \
	$(OBJDIR)/libparticle.a \
	$(OBJDIR)/libgrid.a \
	$(OBJDIR)/libhalofinder.a \
	$(OBJDIR)/libBHForceTree.a \
	$(OBJDIR)/libhalotime.a \
	$(OBJDIR)/libbigchunk.a


# Probably need to fix dependencies on some of the header files.
HDR= \
	RCBForceTree.h \
	RCOForceTree.h \
	BHForceTree.h \
	BHTree.h \
	ChainingMesh.h \
	CosmoHalo.h \
	CosmoHaloFinder.h \
	CosmoHaloFinderP.h \
	Definition.h \
	FOFHaloProperties.h \
	ForceLaw.h \
	ForceTree.h \
	HaloCenterFinder.h \
	InitialExchange.h \
	Message.h \
	Partition.h \
	ParticleDistribute.h \
	ParticleExchange.h \
	SODHalo.h \
	SubHaloFinderDev.h \
	SubHaloFinder.h \
	TestExchange.h \
	Timer.h \
	Timings.h


all: ${LIBS}
libs: ${LIBS}
progs: ${PROGS}

include include.mk
include ../include.mk
include ../dfft/include.mk
CXXFLAGS += ${DFFT_CXXFLAGS}

CFLAGS += ${HACC_MPI_CFLAGS}
CXXFLAGS += ${HACC_MPI_CXXFLAGS}
LDFLAGS += ${HACC_MPI_LDFLAGS}

#CXXFLAGS += ${MPI_COMPILE_FLAGS} 
CXXFLAGS += ${HF_TYPE_FLAGS}
CXXFLAGS += ${HF_WARNING}
#CXXFLAGS += -g
#CXXFLAGS += -DDEBUG
CXXFLAGS += -I.
#LDFLAGS  = ${MPI_LD_FLAGS} -lmpi_cxx

CFLAGS += -I../dfft
CXXFLAGS += -I../dfft


$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: %.c | $(OBJDIR)
	${HACC_MPI_CC} ${CFLAGS} -c -o $@ $<

$(OBJDIR)/%.o: %.cxx | $(OBJDIR)
	${HACC_MPI_CXX} ${CXXFLAGS} -c -o $@ $<


#$(SOBJDIR):
#	mkdir -p $(SOBJDIR)
#
#$(SOBJDIR)/%.o: %.c | $(SOBJDIR)
#	${HACC_CC} -c -o $@ $<
#
#$(SOBJDIR)/%.o: %.cxx | $(SOBJDIR)
#	${HACC_CXX} -c -o $@ $<


#$(OBJDIR)/ForceLaw.o: ForceLaw.h


dims-local.c: ../dfft/dims.c
	cp -f ../dfft/dims.c dims-local.c


$(OBJDIR)/ForceSpeed: $(OBJDIR)/ForceSpeed.o $(OBJDIR)/ForceLaw.o
	${HACC_MPI_CXX} -o $@ $^ ${HACC_MPI_LDFLAGS}


$(OBJDIR)/ForceRandom: $(OBJDIR)/ForceRandom.o $(OBJDIR)/ForceLaw.o
	${HACC_MPI_CXX} -o $@ $^ ${HACC_MPI_LDFLAGS}


$(OBJDIR)/HaloTest.o: CosmoHaloFinder.h
$(OBJDIR)/HaloTest: $(OBJDIR)/HaloTest.o $(OBJDIR)/CosmoHaloFinder.o
	${HACC_MPI_CXX} -o $@ ${CXXFLAGS} $^ ${LDFLAGS}


$(OBJDIR)/HaloTestP: $(OBJDIR)/HaloTestP.o $(OBJDIR)/SubHaloFinder.o $(OBJDIR)/HaloFinderInput.o $(OBJDIR)/libhalofinder.a $(OBJDIR)/libparticle.a $(OBJDIR)/libhalotime.a $(OBJDIR)/libpartition.a $(OBJDIR)/libbigchunk.a $(OBJDIR)/dims-local.o
	${HACC_MPI_CXX} -o $@ $^ ${CXXFLAGS} ${LDFLAGS}


$(OBJDIR)/HaloFinder: $(OBJDIR)/HaloFinder.o $(OBJDIR)/SubHaloFinder.o $(OBJDIR)/HaloFinderInput.o $(OBJDIR)/libhalofinder.a $(OBJDIR)/libparticle.a $(OBJDIR)/libhalotime.a $(OBJDIR)/libpartition.a $(OBJDIR)/libbigchunk.a $(OBJDIR)/dims-local.o
	${HACC_MPI_CXX} -o $@ $^ ${CXXFLAGS} ${LDFLAGS}


$(OBJDIR)/SubHaloTest: $(OBJDIR)/SubHaloTest.o $(OBJDIR)/SubHaloFinderDev.o $(OBJDIR)/HaloFinderInput.o $(OBJDIR)/libhalofinder.a $(OBJDIR)/libparticle.a $(OBJDIR)/libpartition.a $(OBJDIR)/libhalotime.a $(OBJDIR)/libbigchunk.a $(OBJDIR)/dims-local.o
	${HACC_MPI_CXX} -o $@ $^ ${CXXFLAGS} ${LDFLAGS}


$(OBJDIR)/ForceTestP: $(OBJDIR)/ForceTestP.o $(OBJDIR)/ForceTree.o $(OBJDIR)/HaloFinderInput.o $(OBJDIR)/libparticle.a $(OBJDIR)/libBHForceTree.a $(OBJDIR)/libpartition.a $(OBJDIR)/libhalotime.a $(OBJDIR)/libbigchunk.a $(OBJDIR)/dims-local.o
	${HACC_MPI_CXX} -o $@ $^ ${CXXFLAGS} ${LDFLAGS}

$(OBJDIR)/ForceTreeTest: $(OBJDIR)/ForceTreeTest.o $(OBJDIR)/libparticle.a $(OBJDIR)/libBHForceTree.a $(OBJDIR)/libpartition.a $(OBJDIR)/libhalotime.a $(OBJDIR)/libbigchunk.a $(OBJDIR)/dims-local.o
	${HACC_MPI_CXX} -o $@ $^ ${CXXFLAGS} ${LDFLAGS}

$(OBJDIR)/GridTestP: $(OBJDIR)/GridTestP.o $(OBJDIR)/libpartition.a $(OBJDIR)/libgrid.a $(OBJDIR)/libbigchunk.a $(OBJDIR)/dims-local.o
	${HACC_MPI_CXX} -o $@ $^ ${CXXFLAGS} ${LDFLAGS}


$(OBJDIR)/InitTestP: $(OBJDIR)/InitTestP.o $(OBJDIR)/TestExchange.o $(OBJDIR)/libparticle.a $(OBJDIR)/libhalotime.a $(OBJDIR)/libpartition.a $(OBJDIR)/libbigchunk.a $(OBJDIR)/dims-local.o
	${HACC_MPI_CXX} -o $@ $^ ${CXXFLAGS} ${LDFLAGS}


$(OBJDIR)/OneToN: $(OBJDIR)/OneToN.o
	${HACC_MPI_CXX} -o $@ $^ ${CXXFLAGS} 


$(OBJDIR)/NToN: $(OBJDIR)/NToN.o
	${HACC_MPI_CXX} -o $@ $^ ${CXXFLAGS} 


$(OBJDIR)/Gadget2ToCosmo: $(OBJDIR)/Gadget2ToCosmo.o
	${HACC_MPI_CXX} -o $@ $^ ${CXXFLAGS} 


$(OBJDIR)/CosmoToGadget2: $(OBJDIR)/CosmoToGadget2.o
	${HACC_MPI_CXX} -o $@ $^ ${CXXFLAGS} 


$(OBJDIR)/BinHalos: $(OBJDIR)/BinHalos.o
	${HACC_MPI_CXX} -o $@ $^ ${CXXFLAGS} ${LDFLAGS}


$(OBJDIR)/libpartition.a: $(OBJDIR)/libpartition.a($(OBJDIR)/Partition.o)
	ranlib $@


PARTICLE_SOURCES += ParticleDistribute.cxx
PARTICLE_SOURCES += ParticleExchange.cxx
PARTICLE_SOURCES += InitialExchange.cxx
PARTICLE_SOURCES += Message.cxx
PARTICLE_OBJLIST = $(PARTICLE_SOURCES:.cxx=.o)
PARTICLE_OBJECTS := $(addprefix $(OBJDIR)/,$(PARTICLE_OBJLIST))
$(OBJDIR)/libparticle.a: $(OBJDIR)/libparticle.a($(PARTICLE_OBJECTS))
	ranlib $@


$(OBJDIR)/libgrid.a: $(OBJDIR)/libgrid.a($(OBJDIR)/GridExchange.o)
	ranlib $@


HALO_SOURCES += CosmoHaloFinder.cxx
HALO_SOURCES += CosmoHaloFinderP.cxx
HALO_SOURCES += FOFHaloProperties.cxx
HALO_SOURCES += HaloCenterFinder.cxx
HALO_SOURCES += SODHalo.cxx
HALO_SOURCES += BHTree.cxx
HALO_SOURCES += ChainingMesh.cxx
HALO_OBJLIST = $(HALO_SOURCES:.cxx=.o)
HALO_OBJECTS := $(addprefix $(OBJDIR)/,$(HALO_OBJLIST))
$(OBJDIR)/libhalofinder.a: $(OBJDIR)/libhalofinder.a($(HALO_OBJECTS))
	ranlib $@


HT_SOURCES += Timings.cxx
HT_SOURCES += Timer.cxx
HT_OBJLIST = $(HT_SOURCES:.cxx=.o)
HT_OBJECTS := $(addprefix $(OBJDIR)/,$(HT_OBJLIST))
$(OBJDIR)/libhalotime.a: $(OBJDIR)/libhalotime.a($(HT_OBJECTS))
	ranlib $@


FORCE_SOURCES += ForceLaw.cxx
FORCE_SOURCES += BHForceTree.cxx
FORCE_SOURCES += RCBForceTree.cxx
FORCE_SOURCES += RCOForceTree.cxx
FORCE_SOURCES_C += BGQStep16.c
FORCE_SOURCES_C += BGQCM.c
FORCE_OBJLIST = $(FORCE_SOURCES:.cxx=.o)
FORCE_OBJLIST_C = $(FORCE_SOURCES_C:.c=.o)
FORCE_OBJECTS := $(addprefix $(OBJDIR)/,$(FORCE_OBJLIST) $(FORCE_OBJLIST_C))
$(OBJDIR)/libBHForceTree.a: $(OBJDIR)/libBHForceTree.a($(FORCE_OBJECTS) $(FORCE_OBJLIST_C))
	ranlib $@

$(OBJDIR)/libbigchunk.a: $(OBJDIR)/libbigchunk.a($(OBJDIR)/bigchunk.o)
	ranlib $@
