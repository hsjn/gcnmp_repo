
OBJECTS   = drive.o equdif.o desys.o com.o
#These libs are important to get linkup:
#LIBRARIES      = /opt/intel/Compiler/11.0/081/mkl/lib/em64t/libmkl_lapack.a  \
#                 /opt/intel/Compiler/11.0/081/mkl/lib/em64t/libmkl_em64t.a  \
#                 /opt/intel/Compiler/11.0/081/mkl/lib/em64t/libmkl_intel_lp64.a\
#                 /opt/intel/Compiler/11.0/081/lib/intel64/libguide.a
LIBRARIES	=
#mac has framework accelerate:
#which is supposed to give lapack access
#gfortran -o eigtest -framework accelerate eigtest.f95

COMP      = gfortran
#COMP	= ifort  on m4400 
OPTION = $(DEBUG)
#OPTION = $(OPTIMIZE)
FC = $(COMP)
CXX = gcc
INC = ./
COMPILE   = $(FC) -I$(INC) -L$(LIBS)  -M$(MSG) $(OMP) $(AUTP) -c \
                 $(OPTION) $(LISTING) $(SPECIAL)
#LINK = -Wl,-noinhibit-exec -Wl,-Map,mapfilename\
#       -L/opt/intel/Compiler/11.0/081/mkl/lib/em64t -lmkl_em64t -lmkl_intel_lp64
#
LOAD      = $(FC)   -o $(NAME) $(LINK) 
PROTECT   = chmod 755
DELETE    = rm -f
.SUFFIXES:	
NAME	=	grav_sling

all: $(NAME)
	@echo
	@echo ">>  `date '+%a %d-%h-%y %r'`  `pwd`  `uname -mns`  $(LOGNAME)"
	@echo
	@echo make of $(NAME) in `pwd` on `hostname` completed
$(NAME):           $(OBJECTS)  $(LIBRARIES)
	$(LOAD)    $(OBJECTS)   $(LIBRARIES) $(FTN90LIB)
	$(PROTECT) $(NAME)



drive.o:	drive.f90 equdif.o com.o
		$(COMPILE) drive.f90

equdif.o:	 equdif.f90
		$(COMPILE) equdif.f90
desys.o:	desys.f90 com.o
		$(COMPILE) desys.f90
com.o:	com.f90
		$(COMPILE) com.f90
clean:
	$(DELETE) $(NAME) *.a *.o *.mod
