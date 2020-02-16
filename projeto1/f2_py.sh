#!/bin/bash

f2py -L/usr/lib/ -lblas -llapack -c mod_vars.f90 mod_output_routines.f90 mod_functions.f90 mod_diagonalization.f90 mod_vars_sw.f90 mod_flux_vector_splitting.f90 mod_vars_proj3.f90 src*.f90 -m bocal
