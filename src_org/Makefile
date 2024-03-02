SHELL=/bin/zsh

SRCS = 2dporous-wd-change.f

.SUFFIXES: .o .cpp .f90 .f
OBJS = $(SRCS:.f=.o)

FC  =	ifort
CMD =	2dp

FFLAGS =  -ipo -O3 -no-prec-div -fp-model fast=2 -xHost -fixed
LDFLAGS =
LIBS   =

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS) -o $(@) $(OBJS) $(LDFLAGS) $(LIBS)

.f.o:
	$(FC) $(FFLAGS) -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

