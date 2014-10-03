#!/bin/sh

ROOT=${HOME}/install/essl_fftw
LIBDIR32=${ROOT}/bglib
export LIBDIR32

LIBDIR64=${ROOT}/bglib64
export LIBDIR64

INCLUDEDIR=${ROOT}/include
export INCLUDEDIR

FFTW3ROOT=/soft/libraries/essl/5.1.1-0.beta/essl/5.1/FFTW3
export FFTW3ROOT

CC=bgxlc_r
export CC

INCLUDE='-I/soft/libraries/essl/5.1.1-0.beta/essl/5.1/include -I/soft/libraries/essl/5.1.1-0.beta/essl/5.1/FFTW3/include'
export INCLUDE

make -e -f ${FFTW3ROOT}/src/Makefile install
