# Makefile for the mex-version of the limited memory bundle method
#
# put your compile flags here
FF = gfortran -O2
# insert your mex command here
MEX = /usr/local/matlab2007/bin/mex
# insert your fortran library here
# (full file name if you use a static library,
# -l<library name> if you use a dynamically linked library)
FORTRAN_LIB = -lgfortran
#FORTRAN_LIB = /usr/lib/gcc/i486-linux-gnu/4.2.4/libgfortran.a

all: lmbmex.o lmbm.o lmsub.o matcal.o lmbm_driver

lmbmex.o: lmbmex.f
	$(FF) -c lmbmex.f

lmbm.o: lmbm.f
	$(FF) -c lmbm.f

lmsub.o: lmsub.f
	$(FF) -c lmsub.f

matcal.o: matcal.f
	$(FF) -c matcal.f

lmbm_driver: lmbm.o lmbmex.o lmsub.o matcal.o
	$(MEX) lmbm_driver.c objfunc.c lmbm.o lmbmex.o lmsub.o matcal.o $(FORTRAN_LIB)

clean:
	rm -f lmbm_driver.mexglx *.o
