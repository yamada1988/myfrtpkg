FC=gfortran
LIBS=-L/usr/lib64 -lfftw3 -I/usr/include

all: testfftw.out

testfftw.out:test_fftw.f90
	$(FC) $(LIBS) $<  -o $@

clean:
	rm -f testfftw.out
