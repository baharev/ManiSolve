#!/bin/bash
set -e
rm -f *.o *.a
gcc -c -fPIC -ggdb3 -O2    -march=native -std=c99 -Wall -Wextra fd05.c
gcc -c -fPIC -ggdb3 -Ofast -march=native -std=c99 -Wall -Wextra randpoint.c
gfortran -c -fPIC -ggdb3 -O2 -march=native va27d.f 
gfortran -c -fPIC -ggdb3 -O2 -march=native ddeps.f 
gfortran -c -fPIC -ggdb3 -Ofast -march=native ddot.f 
ar -q libva27.a *.o
rm -f *.o