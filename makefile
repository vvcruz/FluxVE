# FLUXVE COMPILATION MAKEFILE
#
# Vinicius Vaz da Cruz
#

PROG=fluxve
CC=gcc
FC=ifort -nofor_main
#FC=gfortran
#ifort flags
FFLAGS=-traceback 
#gfortran flags
#FFLAGS=-O3
CFLAGS=-O3 

all: main

main: src/main.o src/rdinput.o src/readallspl.o src/read.o src/mtrxdiag.o src/lapackdiag.o src/splinesurf.o src/readpot.o src/jacobi.o src/readspl.o src/Jflux.o src/morse_vib.o src/fourier.o
	$(FC) $(FFLAGS) -o $(PROG) src/main.o src/rdinput.o src/readallspl.o src/read.o src/readpot.o src/mtrxdiag.o src/lapackdiag.o src/splinesurf.o src/jacobi.o src/readspl.o src/Jflux.o src/morse_vib.o src/fourier.o -lfftw3 -lgsl -lgslcblas

main.o: src/main.c src/readf.h
	$(CC) $(CFLAGS) -c src/main.c

rdinput.o: src/rdinput.c
	$(CC) $(CFLAGS) -c src/rdinput.c

Jflux.o: src/Jflux.c src/Jflux.h
	$(CC) $(CFLAGS) -c src/Jflux.c -lgsl -lgslcblas -lm

read.o: src/read.f src/readf.h
	$(FC) $(FFLAGS) -c src/read.f

readpot.o: src/readpot.f
	$(FC) $(FFLAGS) -c src/readpot.f

readspl.o: src/readspl.f
	$(FC) $(FFLAGS) -c src/readspl.f

mtrxdiag.o: src/mtrxdiag.f 
	$(FC) $(FFLAGS) -c src/mtrxdiag.f

lapackdiag.o: src/lapackdiag.f
	$(FC) $(FFLAGS) -c src/lapackdiag.f

splinesurf.o: src/splinesurf.f src/splinesurf.h
	$(FC) $(FFLAGS) -c src/splinesurf.f

jacobi.o: src/jacobi.c src/jacobi.h
	$(CC) $(CFLAGS) -c src/jacobi.c

morse_vib.o: src/morse_vib.c src/morse_vib.h
	$(CC) $(CFLAGS) -c src/morse_vib.c

fourier.o: src/fourier.c src/fourier.h
	$(CC) $(CFLAGS) -c src/fourier.c -lfftw3 -lm

clean:
	rm *~ src/*~ src/*.o
