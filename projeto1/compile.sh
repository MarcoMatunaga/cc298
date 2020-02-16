#!/bin/bash

clear
sh cleanbuild.sh

#===================================================================
# directories to look for Fortran source files
#===================================================================

ndir=1

srcdir[0]="./"

#===================================================================
# define flags for different compilers
#===================================================================

#ifort_flags_dbg=" -g -O0 -fpe:0 -traceback -implicitnone -warn declarations -warn unused -warn ignore_loc -warn truncated_source -check all -check bounds -check uninit -ftrapuv -gen-interface -warn interfaces -debug all -static_mpi -check_mpi -fp-speculation=off -init=zero -init=arrays -init=snan"
#ifort_flags_opt=" -O3 -xHost -static -fast -assume buffered_io"

gfortran_flags_dbg="-O0 -g -fbacktrace -fcheck=all -fcheck=bounds -Wall" 

gfortran_flags_opt="-O0 -g -fbacktrace -fcheck=all -fcheck=bounds -Wall"

#===================================================================
# define necessary libs and linker flags (LDFLAGS)
#===================================================================

# LAPACK library for gradients (intel mkl).
#LIBS="-mkl"

# LAPACK library for gradients (General Linux implementation).
LIBS="-lblas -llapack"

#===================================================================
# Ask user which compiler and compilation mode to use
#===================================================================

# Compiler.
# FC = 0 - ifort
# FC = 1 - gfortran
# FC = 2 - MPI intel
# FC = 3 - MPI gfortran
#
# Otimization.
# MODE = 0 - Debug
# MODE = 1 - Optimized

echo "Setting MPI intel and Optimized compilation parameters..."

FC=1
MODE=0

#===================================================================
# process template for mkmf
#===================================================================

# every compiled object and temp files are written on a
# temp folder and later removed

if   [ $FC == 0 ]
then
  compiler="ifort"

  if [ $MODE == 0 ]
  then
    fflags=${ifort_flags_dbg}
  elif [ $MODE == 1 ]
  then
    fflags=${ifort_flags_opt}
  fi

elif [ $FC == 1 ]
then
  compiler="gfortran"

  if [ $MODE == 0 ]
  then
    fflags=${gfortran_flags_dbg}
  elif [ $MODE == 1 ]
  then
    fflags=${gfortran_flags_opt}
  fi

elif [ $FC == 2 ]
then
  compiler="mpif90"

  if [ $MODE == 0 ]
  then
    fflags=${ifort_flags_dbg}
  elif [ $MODE == 1 ]
  then
    fflags=${ifort_flags_opt}
  fi

elif [ $FC == 3 ]
then
  compiler="mpif90"

  if [ $MODE == 0 ]
  then
    fflags=${ifort_flags_dbg}
  elif [ $MODE == 1 ]
  then
    fflags=${ifort_flags_opt}
  fi

fi

echo \
"FC = ${compiler}
LD = ${compiler}
FFLAGS = ${fflags}
LDFLAGS = ${LIBS}
" > mkmf_template

#===================================================================
# begin compilation
#===================================================================

# initialize the directories list
dirlist=${srcdir[0]}

for ((  i = 1 ;  i < $ndir;  i++  ))
do
  dirlist="${dirlist} ${srcdir[i]}"
done

echo ''
echo 'Looking for Fortran sources in: '${dirlist}

echo ''
echo '-----------------------------------------------------'
echo 'compiling...'
echo '-----------------------------------------------------'

# Create Makefile with mkmf and the template file
# NOTE: if there is not a mkmf in the system use the line below
#perl mkmf -x -t mkmf_template ${dirlist}
# NOTE: if you have a mkmf in your system use the line below
mkmf -x -t mkmf_template ${dirlist}
echo 'done.'

