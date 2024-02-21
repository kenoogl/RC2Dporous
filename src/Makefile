SHELL=/bin/bash

SRCS = rc2d_array.f90 rc2d_param.f90 rc2d_vector.f90 rc2d_fileio.f90  \
		rc2d_prs.f90 rc2d_util.f90 rc2d_pm.f90 rc2d.f90  

.SUFFIXES: .o .cpp .f90
OBJS = $(SRCS:.f90=.o)

FC  =	gfortran
CMD =	rc2d

FFLAGS =  -O3 #-ipo -O3 -no-prec-div -fp-model fast=2 -xHost -free -qopenmp
LDFLAGS = 
LIBS   =

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS) -o $(@) $(OBJS) $(LDFLAGS) $(LIBS)

.f90.o:
	$(FC) $(FFLAGS) -c $<

clean:
	$(RM) $(OBJS) $(CMD) *.mod

