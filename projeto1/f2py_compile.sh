#!/bin/bash

clear
sh cleanbuild.sh

#===================================================================
# directories to look for Fortran source files
#===================================================================

ndir=1
srcdir[0]="./"

#===================================================================
# define flags
#===================================================================

compiler="f2py" 
fflags=" "
libs="-lblas -llapack"
module="-m bocal"

#===================================================================
# create the template for the mkmf
#===================================================================

echo \
"
FC=${compiler}
LD=${compiler}
FFLAGS="${libs}${fflags}"
MODULE=${module}
" > mkmf_template

mkmf -x -t mkmf_template ${dirlist}
