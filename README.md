The Ising model in C++
------------------------
There are one .cpp file, Isingmodel.cpp, and one .h file, random.h.

random.h
--------
Includes the function random, which generates random numbers and returns a uniform deviate between 0.0 and 1.0.

Isingmodel.cpp
--------------
Main program to explore temperature-dependent behaviour in ferromagnets with an Ising model. We use the Markov Chain Monte Carlo approach to sample spin configurations and computes the mean energy, mean magnetization, the heat capacity and the susceptibility. The parallelization using OpenMP to estimate the speed-up factor and to numerically estimate the critical temperature at which our 2D system undergoes a phase transition.

Build: g++-12 -O3 Isingmodel.cpp -forenmp -std=c+11 -I /usr/loca/opt/armadillo/include -L /usr/loca/opt/armadillo/lib -larmadillo  -o Isingmodel.exe  
Build: export OMP_NUM_THREADS = 4 (only fot he third part)
Run: ./Isingmodel.exe
In order to execute it correctly we must specify in which part we are 1, 2 or 3:
1 to only execute the quantities as a function of temperatures, for a given temperature range.
2 to only execute the quantities as a function of Monte Carlo cycles, for a given mcs range.
3 The same as in 1 but with parallelization

i.e. For the first part; Run: ./Isingmodel.exe 1

PenningTrap.py
--------------
Python script that reads the datas obtained in the C++ and generates plots of the the analytical results and the values obtained by using Markov Chain Monte Carlo for 2 x 2 lattice with periodic boundary conditions. The analytical results were obtained by functions. Also, it generates the plots of relative errors for the mean energy, magnetization, susceptibility and heat capacity. Normalized histograms of mean energy generated samples for two different temperatures for orderend and unordered initial states and the quantities versus temperature for some lattices are also plot there. Finally, the critical temperature versus each lattice has been represented.


Run command: python3 Isingmodel.py

