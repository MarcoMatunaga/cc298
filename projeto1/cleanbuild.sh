#!/bin/bash
# based on the cleanbuild.sh from BRU3D
rmfiles ()
{ 
 for ARG in $*
 do 
   if [ -f $ARG ]
   then 
   rm -f $ARG
   fi
 done
}

# files to be deleted
rmfiles *.out mkmf_template *.a *.o *.mod *.dat *__genmod.f90
