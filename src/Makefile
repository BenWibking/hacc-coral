#SUBDIRS += cell
SUBDIRS += cpu
SUBDIRS += dfft
SUBDIRS += halo_finder
SUBDIRS += initializer
SUBDIRS += simulation
#SUBDIRS += util

clean::
	rm -f *~
	-for d in $(SUBDIRS); do (cd $$d; make clean ); done
