FFTW_CFLAGS = -I${FFTW_HOME}/include
FFTW_CXXFLAGS = -I${FFTW_HOME}/include

clean::
	rm -rf $(OBJDIR); rm -f *~ *-local.*
