NAME = bocal

#MODFILES = modules/vars.f90
F90FILES = vars.f90 output_routines.f90 functions.f90 mesh.f90 jacobian_eta.f90 jacobian_ksi.f90 metric_terms.f90 euler_explicit.f90 implicit_beam_warming.f90 create_block_tridiagonal.f90 thomas_block_tridiagonal.f90 thomas_block_tridiagonal_v0.f90 thomas_block_tridiagonal_v1.f90 blktriad.f90 inv.f90 initial_cond_curv.f90 fluxes_curvilinear.f90 vector_curvilinear.f90 boundary_conditions_curv.f90 residue.f90 artificial_dissipation_D4.f90 artificial_dissipation_D2.f90 inverse.f90 projeto1.f90

#OFILES = $(MODFILES:modules/*.f90=.o) $(F90FILES:.f90=.o)
OFILES = $(F90FILES:.f90=.o)

#FC = ifort
FC = gfortran
#FLAGS = -O3 -fopenmp
#FLAGS = -O0  -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interfaces -warn interfaces -I/usr/include
#FLAGS = -O0 -g -fpe:0 -traceback -check all -check bounds -check uninit -ftrapuv -debug all -I/usr/include
#FLAGS = -O0 -g -mkl
# gfortran fin bad operations (be careful!!) and NaN
#FLAGS = -g -fbacktrace -ffpe-trap=zero,overflow,underflow
#
FLAGS = -g 
# linkagem de bibliotecas ldd(comando linux para listar as bibliotecas)
LIBS = -lblas -llapack
#LIBS = -lnetcdf -lnetcdff -L/home/schiavo/libs/lib -I/home/schiavo/libs/include -lcgns 
#--------------------------------------------------------------------------------------

%.o : %.f90
	$(FC) $(FLAGS) -c $*.f90 $(LIBS)

$(NAME) : $(OFILES)
	$(FC) $(FLAGS) $(OFILES) -o $(NAME) $(LIBS)
	./$(NAME)

clear :
	clear; rm -f *.o *.mod $(NAME) *.txt *genmod*

clear_all:
	clear; rm -f *.o *.mod *.dat $(NAME) *.txt *genmod*
