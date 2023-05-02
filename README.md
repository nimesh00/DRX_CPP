# DRX_CA
C++ codebase of cellular Automata based Dynamic Recrystallization (Under Development)

The source files are in the src folder. "main.cpp" is one of the working implementataions of the paper "Simulation of discontinuous dynamic recrystallization in pure Cu using a probabilistic cellular automaton" by Hallberg, Håkan; Wallin, Mathias; Ristinmaa, Matti.

This codebase has been with me for some years now so it might be difficult for me to resolve any issues that may occur while running this program but i will be happy to assist in any way I can.

## TO RUN THE SIMULATION

Most of the parameters are defined in "helpers.h" file in include/ folder. Change them accordingly and run $cmake .. && make$ in the build folder (I would suggest you start from an empty build folder). For any new executables you might want to make, add them to the CMakeLists.txt file similarly to how it is done for "main" executable.

## DEPENDENCIES
In case you want to enable plotting (live or post simulation plots of stress vs strain) then you will need gnuplot program installed for that. Apart from that you will need a C++ compiler like gcc to compile the program. C++ versions beyond C++14 are not supported so make sure you select the correct version while compiling.
