#!/bin/bash
g++ -c planeWaves.C -o planeWaves.o -O3 -std=c++11
g++ -c UEGHamiltonian.C -o UEGHamiltonian.o -O3 -std=c++11
g++ -c UEG_MAIN_binarytest.C  -o FCIQMC_MAIN.o -O3 -std=c++11
g++ planeWaves.o UEGHamiltonian.o FCIQMC_MAIN.o -o FCIQMC_MAINexec -O3 -std=c++11
exit 0
