NAME = bocal

F90FILES = *.f90

OFILES = $(F90FILES:.f90=.o)

#FC = ifort
FC = gfortran
#FLAGS = -O3 -fopenmp
#FLAGS = -O0  -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interfaces -warn interfaces -I/usr/include
# FLAGS = -O0 -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -I/usr/include
FLAGS = -O2 -g  
# linkagem de bibliotecas ldd(comando linux para listar as bibliotecas)
#LIBS = -lnetcdf -lnetcdff -L/home/schiavo/libs/lib -I/home/schiavo/libs/include -lcgns 
#--------------------------------------------------------------------------------------

%.o : %.f90
	$(FC) $(FLAGS) -c $*.f90 $(LIBS)

$(NAME) : $(OFILES)
	$(FC) $(FLAGS) $(OFILES) -o $(NAME) $(LIBS)

clear :
	clear; rm -f *.o *.mod $(NAME) *.txt *genmod*

run :
	./$(NAME)
