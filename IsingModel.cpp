/*

 IsingModel.cpp

    Created by Elena Muñoz, David Martínez, Antonio Gómez and Alejandro Carballido.

    Code based on Morten's structure in: Morten Hjorth-Jensen, "Computational Physics Lecture Notes Fall 2015"

Description of the program:
 
 Main program to explore temperature-dependent behaviour in ferromagnets with an Ising model. We use the Markov Chain Monte Carlo approach to sample spin configurations and computes the mean energy, mean magnetization, the heat capacity and the susceptibility. The parallelization using OpenMP to estimate the speed-up factor and to numerically estimate the critical temperature at which our 2D system undergoes a phase transition.
 
--> Instructions to compile:
 
 g++-12 -O3 Isingmodel.cpp -forenmp -std=c+11 -I /usr/loca/opt/armadillo/include -L /usr/loca/opt/armadillo/lib -larmadillo -o Isingmodel.exe
 
--> Instructions to run:
 
 ./Isingmodel.exe output_file_name part
 
*/

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <assert.h>
#include <cmath>
#include <armadillo>
#include <math.h>
#include "omp.h"  // Header file with OpenMP
#include "random.h" // Header file with the functions to generate random numbers


using namespace std;
using namespace arma;


ofstream ofile;

// Inline function to obtain the position in the lattice of the spin's neighbours, with periodic boundary conditions
// -Arguments:
//      int i: Position of the spin where you are computing the energy
//      int limit: Length of the lattice (n_spins)
//      int add: Integer that takes the value -1 for the neigbour to the left and +1 for the one inmediately to the right

inline int periodic(int i, int limit, int add) {
    return (i+limit+add) % (limit);
} // end function periodic


// Function to read data from screen
// -Arguments:
//      int& n_spins: Length of the lattice
//      int& mcs: Number of MonteCarlo cycles
//      double& initial_temp: Initial temperature
//      double& final_temp: Final temperature
//      double& temp_step: Length of the steps between different temperatures


void read_input(int& n_spins, int& mcs, double& initial_temp, double& final_temp, double& temp_step)
{
    n_spins=40;
    mcs=100000;
    initial_temp=2.25;
    final_temp=2.36;
    temp_step=0.01;
} // end function read_input


// Function to initialise energy, spin matrix and magnetization
// -Arguments:
//      int& n_spins: Length of the lattice
//      double& temp: Temperature of the state. When the temperature is below 1.5 J/kb the initial spins start ordered and when it's over that temperature it keeps the same configuration as the last temperature computed.
//      mat& spin_matrix: Armadillo matrix that stores the orientation of the spins (+1/-1)
//      double& E: Energy of the state; before calling initialize it has to be set to 0
//      double& M: Magnetization of the state; before calling initialize it has to be set to 0
//      int& random: Variable for setting the "level of randomness" of the initialization; when it equals 1 it orders the initial spins and when set to 2 it initialize a disordered lattice. If it equals another number it takes into account the temperature for setting the initial lattice
//      long& seed: number for initializing a random sequence (it has to be negative)

void initialize(int n_spins, double temp, mat& spin_matrix, double& E, double& M, int random, long& seed)
{
    // Setup spin matrix
    
    if (random==1) temp=1.49; // Set temp<1.5 for every temperature when asking for ordered lattices (random==1)
    if (random==2){ // For random initialization flipping random spins
        for(int y =0; y < n_spins; y++) {
            for (int x= 0; x < n_spins; x++){
                // Find random position
                int ix = (int) (ran(&seed)*(double)n_spins);
                int iy = (int) (ran(&seed)*(double)n_spins);
                if (ran(&seed) <= 0.5) { // Random selection for flipping with the same probability for flipping or not
                    spin_matrix(iy,ix) *= -1;
                }
            }
        }
    }
    else{
        
        for(int y =0; y < n_spins; y++) {
            for (int x= 0; x < n_spins; x++){
                if (temp < 1.5) spin_matrix(y,x) = 1; // lattice ordered only for low temperatures
            }
        }
    }
    
    // setup initial energy and magnetization
    
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            M += (double) spin_matrix(y,x);
            E -= (double) spin_matrix(y,x)*(spin_matrix(periodic(y,n_spins,-1),x) + spin_matrix(y,periodic(x,n_spins,-1)));
        }
    }
} // end function initialise
    



// Function for Metropolis algorithm. If the probability of a new random state is bigger than the current state it changes the state and computes the new energy and magnetization
// -Arguments:
//      int& n_spins: Length of the lattice
//      long& seed: number for initializing a random sequence (it has to be negative)
//      mat& spin_matrix: Armadillo matrix that stores the orientation of the spins (+1/-1)
//      double& E: Energy of the state; before calling initialize it has to be set to 0
//      double& M: Magnetization of the state; before calling initialize it has to be set to 0
//      vec exp_term: Armadillo vector that contains the values of the possible probabilities depending on the energy

void Metropolis(int n_spins, long& seed, mat& spin_matrix, double& E, double& M, vec exp_term)
{
    // loop over all spins
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            // Find random position
            int ix = (int) (ran(&seed)*(double)n_spins);
            int iy = (int) (ran(&seed)*(double)n_spins);
            int deltaE = 2*spin_matrix(iy,ix)*
            (spin_matrix(iy,periodic(ix,n_spins,-1))+
             spin_matrix(periodic(iy,n_spins,-1),ix) +
             spin_matrix(iy,periodic(ix,n_spins,1)) +
             spin_matrix(periodic(iy,n_spins,1),ix));
            // Here we perform the Metropolis test
            if (ran(&seed) <= exp_term(deltaE+8) ) {
                spin_matrix(iy,ix) *= -1; // flip one spin and accept new spin config
                // update energy and magnetization
                M += (double) 2*spin_matrix(iy,ix);
                E += (double) deltaE;
            }
        }
    }
} // end of Metropolis sampling over spins





// Function that prints to file the results of the calculations
// -Arguments:
//      int& n_spins: Length of the lattice
//      int& mcs: Number of MonteCarlo cycles
//      double& temp: Temperature of the state.
//      vec& average: Armadillo vector that contains the sum of the values of the energy for each state, energy squared, magnetization and magnetization squared.
void output(int n_spins, int mcs, double temp, vec& average, int part)
{
    // Compute the normalized results
    double norm = 1/((double) (mcs)); // divided by total number of cycles
    double Eaverage = average(0)*norm;
    double E2average = average(1)*norm;
    double Maverage = average(2)*norm;
    double M2average = average(3)*norm;
    double Mabsaverage = average(4)*norm;
    // all expectation values are per spin, divide by 1/n_spins/n_spins
    double Evariance = (E2average- Eaverage*Eaverage)/n_spins/n_spins;
    double Mvariance = (M2average - Maverage*Maverage)/n_spins/n_spins;
    double M2variance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;
    
    if (part==2)
    {
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << setw(15) << setprecision(3) << temp;
        ofile << setw(15) << setprecision(1) << mcs;
        ofile << setw(15) << setprecision(8) << Eaverage/n_spins/n_spins;
        ofile << setw(15) << setprecision(8) << Mabsaverage/n_spins/n_spins<< endl;;
    }
    
    else
    {
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << setw(15) << setprecision(3) << temp;
        ofile << setw(15) << setprecision(8) << Eaverage/n_spins/n_spins; //Expectation value of the energy
        ofile << setw(15) << setprecision(8) << Evariance/temp/temp; //Expectation value of the heat capacity
        ofile << setw(15) << setprecision(8) << Mabsaverage/n_spins/n_spins; //Expectation value of the magnetization
        ofile << setw(15) << setprecision(8) << M2variance/temp << endl; //Expectation value of the energy
        
    }

} // end output function


// main program
int main(int argc, const char* argv[])
{
    string outfilename;
    int part;
    long seed;
    int n_spins, mcs;
    double initial_temp, final_temp, temp_step;
    
    // Read in output file and part, abort if there are too few command-line arguments
    
    if( argc <= 1 ){
        cout << "Error: " << argv[0] <<
        " read also output file and select a part of the program on same line" << endl;
        exit(1);
    }
    else{
        outfilename=argv[1];
        part=atoi(argv[2]);
    }
    if (part==1){
    
        vec exp_term =vec(17);
        vec average=vec(5);
        double E, M;
        
        ofile.open(outfilename);
        // Read in initial values such as size of lattice, temp and cycles
        read_input(n_spins, mcs, initial_temp, final_temp, temp_step);
        
        //Create spin matrix
        mat spin_matrix = mat(n_spins, n_spins).fill(1);
        
        // every node has its own seed for the random numbers, this is important else
        // if one starts with the same seed, one ends with the same random numbers
        seed = -1; // random starting point
        
        
        
        // Start Monte Carlo sampling by looping over the different temperatures
        for ( double temp = initial_temp; temp <= final_temp; temp+=temp_step){
            // initialise energy and magnetization
            E = M = 0.;
            // setup array for possible energy changes
            for(int de =-8; de <= 8; de++) exp_term(de+8) = 0;
            for(int de =-8; de <= 8; de+=4) exp_term(de+8) = exp(-de/temp);
            // initialise array for expectation values
            for(int i = 0; i < 5; i++) average(i) = 0.;
            initialize(n_spins, temp, spin_matrix, E, M, 0, seed);
            // start Monte Carlo computation
            for (int cycles = 1; cycles <= mcs; cycles++){
                Metropolis(n_spins, seed, spin_matrix, E, M, exp_term);
                // update expectation values
                average(0) += E; average(1) += E*E;
                average(2) += M; average(3) += M*M; average(4) += fabs(M);
               
            }
            // print results
            output(n_spins, mcs, temp, average, part);
            
        }
        ofile.close(); // close output file
    }
    
    else if (part==2){
        
        double E, M;
        vec exp_term = vec(17);
        vec average = vec(5);
        
        ofile.open(outfilename);
        
        int n_spins, mcs;
        double initial_temp, final_temp, temp_step;

        read_input(n_spins, mcs, initial_temp, final_temp, temp_step);
        
        mat spin_matrix = mat(n_spins, n_spins).fill(1);
        
        seed = -1; // random starting point
        
        
        for ( double temp = initial_temp; temp <= final_temp; temp+=temp_step){
            // Start Monte Carlo sampling by looping over T first
            for (int burn_in_mcs = 100; burn_in_mcs<= 100000; burn_in_mcs+=100){
                mcs = burn_in_mcs;
                // initialise energy and magnetization
                E = M = 0.;
                // setup array for possible energy changes
                for(int de =-8; de <= 8; de++) exp_term(de+8) = 0;
                for(int de =-8; de <= 8; de+=4) exp_term(de+8) = exp(-de/temp);
                // initialise array for expectation values
                for(int i = 0; i < 5; i++) average(i) = 0.;
                initialize(n_spins, temp, spin_matrix, E, M, 1, seed); //ordered spins
                // start Monte Carlo computation
                for (int cycles = 1; cycles <= mcs; cycles++){
                    Metropolis(n_spins, seed, spin_matrix, E, M, exp_term);
                    // update expectation values
                    average(0) += E; average(1) += E*E;
                    average(2) += M; average(3) += M*M; average(4) += fabs(M);
                }
        
                // print results
                output(n_spins, mcs, temp, average, part);
            }
            for (int burn_in_mcs = 100; burn_in_mcs<= 100000; burn_in_mcs+=100){
                mcs = burn_in_mcs;
                // initialise energy and magnetization
                E = M = 0.;
                // setup array for possible energy changes
                for(int de =-8; de <= 8; de++) exp_term(de+8) = 0;
                for(int de =-8; de <= 8; de+=4) exp_term(de+8) = exp(-de/temp);
                // initialise array for expectation values
                for(int i = 0; i < 5; i++) average(i) = 0.;
                initialize(n_spins, temp, spin_matrix, E, M, 2, seed); //disordered spins
                // start Monte Carlo computation
                for (int cycles = 1; cycles <= mcs; cycles++){
                    Metropolis(n_spins, seed, spin_matrix, E, M, exp_term);
                    // update expectation values
                    average(0) += E; average(1) += E*E;
                    average(2) += M; average(3) += M*M; average(4) += fabs(M);
                }
                // print results
                output(n_spins, mcs, temp, average, part);
            }
        }
        ofile.close(); // close output file
    }
    
    
    
    else if (part==3){
    
        
        // Read in initial values such as size of lattice, temp and cycles
        read_input(n_spins, mcs, initial_temp, final_temp, temp_step);
        
        const int n_temp_max = (final_temp-initial_temp)/temp_step;
    
        ofile.open(outfilename);
        mat averages = mat(n_temp_max+1, 6); //create a matrix for saving the results of every thread
        
        #pragma omp parallel
        {
            mat spin_matrix = mat(n_spins,n_spins).fill(1);
            double temp, E, M;
            vec exp_term = vec(17);
        
            seed = -1; // random starting point
            
            #pragma omp for
            for (int n_temp = 0; n_temp <= n_temp_max; n_temp++)
            {
                
                temp = n_temp*temp_step+initial_temp;
                
                // initialise energy and magnetization
                E = M = 0.;
                
                // setup array for possible energy changes
                for(int de =-8; de <= 8; de++) exp_term(de+8) = 0;
                for(int de =-8; de <= 8; de+=4) exp_term(de+8) = exp(-de/temp);
                
                initialize(n_spins, temp, spin_matrix, E, M, 0, seed);
                
                // start Monte Carlo computation
                for (int cycles = 1; cycles <= mcs; cycles++){
                    Metropolis(n_spins, seed, spin_matrix, E, M, exp_term);
                    // update expectation values
                    averages(n_temp,0) += E; averages(n_temp,1) += E*E;
                    averages(n_temp,2) += M; averages(n_temp,3) += M*M;
                    averages(n_temp,4) += fabs(M);
                }
                //Also store temperature for later output
                averages(n_temp,5) += temp;
                
            }
            ofile.close(); // close output file
        }
    
        ofile.open(outfilename);
        
        // print results into the output file
        vec average = vec(5);
        
        for (int i=0; i<=n_temp_max; i++){
            for (int j=0; j<5; j++){
                average(j)=averages(i,j); //store in a vec because output only accepts vecs
            }
            output(n_spins, mcs, averages(i,5), average, part);
        }
        
        ofile.close();
        
    }
    return 0;
}
