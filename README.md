Time dependent Schrödinger equation for a double slit in a box
------------------------
There are one .cpp file, Project5.cpp, and one .hpp file, Project5.hpp.

Project5.hpp
--------------
Includes the functions:
- k: converts the indices of the matrix into a vector index.
- mat_construct: constructs the matrices A and B following the Crank-Nicolson method
- mat_fill: fills the matrices diagonals and subdiagonals
- u_next: obtains the wave function in the next time step
- set_up_u_0: sets up the initial wave function
- V: sets up the slit barrier 



Project5.cpp
--------------
Main program to to simulate the two-dimensional time-dependent Schrödinger equation, and use it to study a double-slit-in-a-box setup.

Build: g++ -O3 Project5.cpp -std=c++14 -larmadillo -o Project5.exe
Run: ./Project5.exe

PenningTrap.py
--------------
Python script that reads the datas obtained in the C++ and generates plots of the the analytical results and the values obtained by using Markov Chain Monte Carlo for 2 x 2 lattice with periodic boundary conditions. The analytical results were obtained by functions. Also, it generates the plots of relative errors for the mean energy, magnetization, susceptibility and heat capacity. Normalized histograms of mean energy generated samples for two different temperatures for orderend and unordered initial states and the quantities versus temperature for some lattices are also plot there. Finally, the critical temperature versus each lattice has been represented.


Run command: python3 Isingmodel.py

