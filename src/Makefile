SHELL=/bin/bash

SRCS = rc2d_pm.f90 rc2d.f90 rc2d_util.f90 rc2d_fileio.f90 \
		rc2d_prs.f90 rc2d_vector.f90 rc2d_wake.f90 rc2d_params.f90

.SUFFIXES: .o .cpp .f90
OBJS = $(SRCS:.f90=.o)

FC  =	gfortran
#FC  =	ifort
CMD =	rc2d

FFLAGS =	-O3 -cpp -D_SPH -D_WINDMILL -ljson -L../json-fortran/src/lib -I../json-fortran/src/include # gfortran
#FFLAGS =	-O3 -cpp -D_SPH -D_WINDMILL -ipo -no-prec-div -fp-model fast=2 -xHost -free -qopenmp \ # ifort
#		-ljson -L../json-fortran/src/lib -I../json-fortran/src/include
LDFLAGS = 
LIBS   =

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS) -o $(@) $(OBJS) $(LDFLAGS) $(LIBS)

.f90.o:
	$(FC) $(FFLAGS) -c $<

clean:
	$(RM) $(OBJS) $(CMD) *.mod

