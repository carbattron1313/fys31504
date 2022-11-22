#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 12:54:51 2022

@author: Elena Muñoz Rivas, Alejandro Carballido Mantecón, Antonio Gómez Garrido
and David Martínez Hernández.
"""

import matplotlib.pyplot as plt
import numpy as np
from typing import List
import math
from scipy.stats import linregress

# Ising model in two dimensions.
#Plot the analytical results and the values obtained by using Markov Chain
#Monte Carlo for 2 x 2 lattice with periodic boundary conditions.

#Load the values obtained by MCMC for T=1J/K_B:
Isingmodel1 = np.loadtxt("Isingmodel_prob4.txt")

#Import the values of the mean energy, the specific heat capacity, the
#susceptibility and the magnetization in Lists:
temperature = Isingmodel1[:,0]
Eaverage = Isingmodel1[:,1]
Evariance = Isingmodel1[:,2]
Maverage = Isingmodel1[:,3]
Mvariance = Isingmodel1[:,4]

#Number of spins, for 2 lattice (N = 2^2 = 4)
n_spins: int = 4

#Analytical solution:
#Function of the mean energy:
def energy(temperature: List[float])-> List[float]:
    """
    Calculte the mean energy of the system for two dimensions with boundary
    conditions and a lattice of 2.

    Args:
        temperature: the temperature of the system, for the exercise 4,
        T = 1J/k_b.
    Returns:
        the mean energy
    """
    J: int = 1
    energy_an: List[float] = []
    for temp in temperature:
        beta: float = 1/temp
        energy_an.append(- J*8*math.sinh(8*beta*J)/(math.cosh(8*beta*J)+3))
    return energy_an


#Function of the mean square energy:
def energy_2(temperature: List[float])-> List[float]:
    """
    Calculte the mean square energy of the system for two dimensions assuming a
    2x2 lattice with boundary conditions.

    Args:
        temperature: the temperature of the system, for the exercise 4,
        T = 1J/k_b.
    Returns:
        the mean square energy
    """
    J: int = 1
    energy_an_2: List[float] = []
    for temp in temperature:
        beta: float = 1/temp
        energy_an_2.append(J*64*math.cosh(8*beta*J)/(math.cosh(8*beta*J)+3))
    return energy_an_2


#Function of the mean magnetization:
def magnetization(temperature: List[float])-> List[float]:
    """
    Calculte the mean magnetization of the system for two dimensions assuming a
    2x2 lattice with boundary conditions.

    Args:
        temperature: the temperature of the system, for the exercise 4,
        T = 1J/k_b.
    Returns:
        the mean magnetization
    """
    J: int = 1
    magnetization_an: List[float] = []
    for temp in temperature:
        beta: float = 1/temp
        magnetization_an.append((4+2*math.exp(8*beta*J))/(math.cosh(8*beta*J)+3))
    return magnetization_an


#Function of the mean square magnetization:
def magnetization_2(temperature: List[float])-> List[float]:
    """
    Calculte the mean square magnetization of the system for two dimensions
    assuming a 2x2 lattice with boundary conditions.

    Args:
        temperature: the temperature of the system, for the exercise 4,
        T = 1J/k_b.
    Returns:
        the mean square magnetization
    """
    J: int = 1
    magnetization_an_2: List[float] = []
    for temp in temperature:
        beta: float = 1/temp
        magnetization_an_2.append(abs(8*(1+math.exp(8*beta*J))/(math.cosh(8*beta*J)+3)))
    return magnetization_an_2



#List of the analytical values obtained by the functions:
Eanalytical2: List[float] = energy_2(temperature)
Eanalytical: List[float] = energy(temperature)
Manalytical2: List[float] = magnetization_2(temperature)
Manalytical: List[float] = magnetization(temperature)

#Create Lists of the analytical values, needed for the for loop:
e_analytical: List[float] = []
m_analytical: List[float] = []
Cv: List[float] = []
Chi: List[float] = []


#Empty list to include the relative errors for the mean energy, magnetization,
#susceptibility and heat capacity:
err_rel_e: List[float] = []
err_rel_m: List[float] = []
err_rel_chi: List[float] = []
err_rel_cv: List[float] = []


#Loop used to obtained the values of the heat capacity, the susceptibility,
#the mean energy per spin and the magnetization per spin for different
#teperatures:
for i in range(len(temperature)):
    Cv.append( 1/(temperature[i]**2*n_spins) * (Eanalytical2[i] - Eanalytical[i]**2))
    Chi.append(1/(temperature[i]*n_spins) * (Manalytical2[i] - Manalytical[i]**2))
    e_analytical.append(Eanalytical[i]/n_spins)
    m_analytical.append(Manalytical[i]/n_spins)


#Loop used to obtained the relative error of the heat capacity, the
#susceptibility, the mean energy per spin and the magnetization per spin for
#different teperatures:
for i in range(len(temperature)):
    err_rel_cv.append(abs(Cv[i]-Evariance[i])/Cv[i])
    err_rel_chi.append(abs(Chi[i]-Mvariance[i])/Cv[i])
    err_rel_e.append(abs(e_analytical[i]-Eaverage[i])/Cv[i])
    err_rel_m.append(abs(m_analytical[i]-Maverage[i])/Cv[i])


#Do the log of  the relatives errors for each temperature:
err_rel_log: List[float]
err_rel_log_cv = np.log(err_rel_cv)
err_rel_log_chi = np.log(err_rel_chi)
err_rel_log_e = np.log(err_rel_e)
err_rel_log_m = np.log(err_rel_m)



#Plot the experimental (MCMC method) and  the analytical mean energy per spin
#versus temperature:
plt.plot(temperature,Eaverage, label='$\epsilon_{exp}$', color = "darkgreen")
plt.plot(temperature,e_analytical, label='$\epsilon_{analytical}$' , color = "blue")
plt.xlabel("T (J$/k_B$)")
plt.ylabel("$\epsilon$ (J)")
plt.grid(True)
plt.legend()
plt.savefig('Isingmodel_energy1.pdf')

plt.show()



#Plot the experimental (MCMC method) and  the analytical heat capacity
#versus temperature:
plt.plot(temperature,Evariance, label='$C_{v(analytical)}$', color = "purple")
plt.plot(temperature,Cv, label='$C_{v(analytical)}$', color = "darkorange")
plt.xlabel("T(J$/k_B$)")
plt.ylabel("$C_v$ ($J k_B$)")
plt.grid(True)
plt.legend()
plt.savefig('Isingmodel_cv1.pdf')

plt.show()


#Plot the experimental (MCMC method) and  the analytical mean magnetization
#versus temperature:
plt.plot(temperature,Maverage, label='$m_{exp}$', color = "darkgreen")
plt.plot(temperature,m_analytical, label='$m_{analytical}$', color = "blue")
plt.xlabel("T (J/$k_B$)")
plt.ylabel("$m$")
plt.grid(True)
plt.legend()
plt.savefig('Isingmodel_mag1.pdf')
plt.show()



#Plot the experimental (MCMC method) and  the analytical mean susceptibility per spin
#versus temperature:
plt.plot(temperature,Mvariance, label='$\chi_{exp}$',  color = "purple")
plt.plot(temperature, Chi, label='$\chi_{analytical}$', color = "darkorange")
plt.xlabel("T (J/$k_B$)")
plt.ylabel("$\chi$")
plt.grid(True)
plt.legend()
plt.savefig('Isingmodel_chi1.pdf')
plt.show()

#=============================================================================


#Plot the relatives errors as functions of temperature:
plt.plot(temperature, err_rel_e, label='$\epsilon_{errel}$', color = "darkgreen")
plt.plot(temperature,err_rel_m, label='$m_{errel}$' , color = "blue")
plt.plot(temperature,err_rel_cv, label='$C_{v errel}$',  color = "purple")
plt.plot(temperature, err_rel_chi, label='$\chi_{errel}$', color = "darkorange")
plt.xlabel("T (J$/k_B$)")
plt.ylabel("$r_{rel}$")
plt.grid(True)
plt.legend()
plt.savefig('Isingmodel_err_rel.pdf')

plt.show()



#Plot the log of the relatives errors as functions of temperature:
plt.plot(temperature, err_rel_log_e, label='$\epsilon_{errel}$', color = "darkgreen")
plt.plot(temperature,err_rel_log_m, label='$m_{errel}$' , color = "blue")
plt.plot(temperature,err_rel_log_cv, label='$C_{v errel}$',  color = "purple")
plt.plot(temperature, err_rel_log_chi, label='$\chi_{errel}$', color = "darkorange")
plt.xlabel("T (J$/k_B$)")
plt.ylabel("$log(r_{rel})$")
plt.grid(True)
plt.legend()
plt.savefig('Isingmodel_err_rel_log.pdf')

plt.show()



#=============================================================================

#Load the values obtained by MCMC for two temperatures and different mcs,
#starting from ordered (all spins pointing the same way) and unordered (random)
#initial states:
Isingmodel_prob5 = np.loadtxt("Isingmodel_prob5.txt")

#Import the values of temperature in List:
temperature2 = Isingmodel_prob5[:,0]


#Create empty List for the values of Monte Carlo cycles and the mean energies
# when the temperature is T=1 J/K_B and T=2.4 J/K_B starting from ordered (ord)
#and unordered (unord) initial states:
mcs_list: List[float] = []

energy_1_ord: List[float] = []
energy_2_ord: List[float] = []
energy_1_unord: List[float] = []
energy_2_unord: List[float] = []

magnetization_1_ord: List[float] = []
magnetization_2_ord: List[float] = []
magnetization_1_unord: List[float] = []
magnetization_2_unord: List[float] = []


#Loop for to create the List for different tempertures and initial states.
#In the .txt the values where mixed so we make ifs to separate them from
#ordered and unordered and from one temperature to another:
count: int = 0
mcs: int = 100000/100
for i in range(len(temperature2)):
        count += 1
        if count<=mcs:
            mcs_list.append(float(Isingmodel_prob5[i,1]))
            energy_1_ord.append(Isingmodel_prob5[i,2])
            magnetization_1_ord.append(Isingmodel_prob5[i,3])

        elif mcs< count<= 2*mcs:
            energy_1_unord.append(Isingmodel_prob5[i,2])
            magnetization_1_unord.append(Isingmodel_prob5[i,3])


        elif 2*mcs< count<=3*mcs:
            energy_2_ord.append(Isingmodel_prob5[i,2])
            magnetization_2_ord.append(Isingmodel_prob5[i,3])
        elif 3*mcs< count:
            energy_2_unord.append(Isingmodel_prob5[i,2])
            magnetization_2_unord.append(Isingmodel_prob5[i,3])


#Plot the mean energy when T = 1J/K_B:
plt.plot(mcs_list, energy_1_ord, label='$\epsilon_{ord1}$', color = "darkgreen")
plt.plot(mcs_list, energy_1_unord, label='$\epsilon_{unord1}$' , color = "blue")
plt.ylabel("$\epsilon (J)$")
plt.xlabel("mcs")
plt.grid(True)
plt.legend()
plt.savefig('Isingmodel_energy1_mcs.pdf')

plt.show()



#Plot the mean energy when T = 1J/K_B:
plt.plot(mcs_list[150:], energy_1_ord[150:], label='$\epsilon_{ord1}$', color = "darkgreen")
plt.plot(mcs_list[150:], energy_1_unord[150:], label='$\epsilon_{unord1}$' , color = "blue")
plt.ylabel("$\epsilon (J)$")
plt.xlabel("mcs")
plt.grid(True)
plt.legend()
plt.savefig('Isingmodel_energy1_mcs_cut.pdf')

plt.show()



#Plot the mean magnetization when T = 1J/K_B:
plt.plot(mcs_list, magnetization_1_ord, label='$m_{ord1}$',  color = "purple")
plt.plot(mcs_list, magnetization_1_unord, label='$m_{unord1}$', color = "darkorange")
plt.ylabel("$m$")
plt.xlabel("mcs")
plt.grid(True)
plt.legend()
plt.savefig('Isingmodel_mag1_mcs.pdf')

plt.show()


#Plot the mean magnetization when T = 1J/K_B:
plt.plot(mcs_list[150:], magnetization_1_ord[150:], label='$m_{ord1}$',  color = "purple")
plt.plot(mcs_list[150:], magnetization_1_unord[150:], label='$m_{unord1}$', color = "darkorange")
plt.ylabel("$m$")
plt.xlabel("mcs")
plt.grid(True)
plt.legend()
plt.savefig('Isingmodel_mag1_mcs_cut.pdf')

plt.show()


#Plot the mean energy when T = 2.4J/K_B:
plt.plot(mcs_list, energy_2_ord, label='$\epsilon_{ord2}$', color = "darkgreen")
plt.plot(mcs_list, energy_2_unord, label='$\epsilon_{unord2}$' , color = "blue")
plt.ylabel("$\epsilon$ (J)")
plt.xlabel("mcs")
plt.grid(True)
plt.legend()
plt.savefig('Isingmodel_energy2_mcs.pdf')

plt.show()

#Plot the mean magnetization when T = 2.4J/K_B:
plt.plot(mcs_list, magnetization_2_ord, label='$m_{ord2}$',  color = "purple")
plt.plot(mcs_list, magnetization_2_unord, label='$m_{unord2}$', color = "darkorange")
plt.ylabel("$m$")
plt.xlabel("mcs")
plt.grid(True)
plt.legend()
plt.savefig('Isingmodel_mag2_mcs.pdf')

plt.show()


#=============================================================================

#Normalized histograms of mean energy generated samples when T = 1 J/K_B and
#T = 2.4 J/K_B for orderend and unordered initial states:
energy_1_ord_norm: List[float] = []
for i in energy_1_ord:
    energy_1_ord_norm.append(i/np.sum(energy_1_ord))





(counts, bins, patches) = plt.hist(energy_1_ord_norm,  bins=100, density= False, edgecolor = 'black')
plt.xlabel("$\epsilon_{ord1}$ (J)")
plt.ylabel("$p_{\epsilon}(\epsilon; T)$")
plt.savefig('Isingmodel_ehistogram_1ord.pdf')

plt.show()


(counts, bins, patches) = plt.hist(energy_2_ord,  bins=100, density=False, edgecolor = 'black')
plt.xlabel("$\epsilon_{ord2}$ (J)")
plt.ylabel("$p_{\epsilon}(\epsilon; T)$")
plt.savefig('Isingmodel_ehistogram_2ord.pdf')

plt.show()


(counts, bins, patches) = plt.hist(energy_1_unord,  bins=100, density=False, edgecolor = 'black')
plt.xlabel("$\epsilon_{unord1}$ (J)")
plt.ylabel("$p_{\epsilon}(\epsilon; T)$")
plt.savefig('Isingmodel_ehistogram_1unord.pdf')

plt.show()


(counts, bins, patches) = plt.hist(energy_2_unord,  bins=100, density=False, edgecolor = 'black')
plt.xlabel("$\epsilon_{unord2}$ (J)")
plt.ylabel("$p_{\epsilon}(\epsilon; T)$")
plt.savefig('Isingmodel_ehistogram_2unord.pdf')

plt.show()


#=============================================================================

#Load the values obtained by MCMC for different lattices of size (L = 40, 60,
#80 and 100) for a range of temperatures T = [2.1,2.4] J/k_B with small steps
#(0.03). And the values when we have made a zoom in.
Isingmodel_prob7_100 = np.loadtxt("Isingmodel_prob7_100.txt")
Isingmodel_prob7_100_zoom = np.loadtxt("Isingmodel_prob7_100_zoom.txt")

#Import the values of the quantities in List (L=100):
temperature3 = list(Isingmodel_prob7_100[:6,0])
energy_3_100= list(Isingmodel_prob7_100[:6,1])
cv_3_100 = list(Isingmodel_prob7_100[:6,2])
magnetization_3_100 = list(Isingmodel_prob7_100[:6,3])
chi_3_100 = list(Isingmodel_prob7_100[:6,4])

temperature3_f = list(Isingmodel_prob7_100[9:,0])
energy_3_100_f= list(Isingmodel_prob7_100[9:,1])
cv_3_100_f = list(Isingmodel_prob7_100[9:,2])
magnetization_3_100_f = list(Isingmodel_prob7_100[9:,3])
chi_3_100_f = list(Isingmodel_prob7_100[9:,4])

temperature3_zoom = list(Isingmodel_prob7_100_zoom[:,0])
energy_3_100_zoom = list(Isingmodel_prob7_100_zoom[:,1])
cv_3_100_zoom = list(Isingmodel_prob7_100_zoom[:,2])
magnetization_3_100_zoom = list(Isingmodel_prob7_100_zoom[:,3])
chi_3_100_zoom = list(Isingmodel_prob7_100_zoom[:,4])

for i in range(len( energy_3_100_zoom)):
    temperature3.append(temperature3_zoom[i])
    energy_3_100.append(energy_3_100_zoom[i])
    cv_3_100.append(cv_3_100_zoom[i])
    magnetization_3_100.append(magnetization_3_100_zoom[i])
    chi_3_100.append(chi_3_100_zoom[i])

for i in range(len( energy_3_100_f)):
    temperature3.append(temperature3_f[i])
    energy_3_100.append(energy_3_100_f[i])
    cv_3_100.append(cv_3_100_f[i])
    magnetization_3_100.append(magnetization_3_100_f[i])
    chi_3_100.append(chi_3_100_f[i])
#=======================
Isingmodel_prob7_80 = np.loadtxt("Isingmodel_prob7_80.txt")
Isingmodel_prob7_80_zoom = np.loadtxt("Isingmodel_prob7_80_zoom.txt")

#Import the values of the quantities in List (L=80):
energy_3_80= list(Isingmodel_prob7_80[:6,1])
cv_3_80 = list(Isingmodel_prob7_80[:6,2])
magnetization_3_80 = list(Isingmodel_prob7_80[:6,3])
chi_3_80 = list(Isingmodel_prob7_80[:6,4])

energy_3_80_f= list(Isingmodel_prob7_80[9:,1])
cv_3_80_f = list(Isingmodel_prob7_80[9:,2])
magnetization_3_80_f = list(Isingmodel_prob7_80[9:,3])
chi_3_80_f = list(Isingmodel_prob7_80[9:,4])

energy_3_80_zoom = list(Isingmodel_prob7_80_zoom[:,1])
cv_3_80_zoom = list(Isingmodel_prob7_80_zoom[:,2])
magnetization_3_80_zoom = list(Isingmodel_prob7_80_zoom[:,3])
chi_3_80_zoom = list(Isingmodel_prob7_80_zoom[:,4])

for i in range(len( energy_3_80_zoom)):
    energy_3_80.append(energy_3_80_zoom[i])
    cv_3_80.append(cv_3_80_zoom[i])
    magnetization_3_80.append(magnetization_3_80_zoom[i])
    chi_3_80.append(chi_3_80_zoom[i])

for i in range(len( energy_3_80_f)):
    energy_3_80.append(energy_3_80_f[i])
    cv_3_80.append(cv_3_80_f[i])
    magnetization_3_80.append(magnetization_3_80_f[i])
    chi_3_80.append(chi_3_80_f[i])

#=======================
Isingmodel_prob7_60 = np.loadtxt("Isingmodel_prob7_60.txt")
Isingmodel_prob7_60_zoom = np.loadtxt("Isingmodel_prob7_60_zoom.txt")

#Import the values of the quantities in List (L=60):
energy_3_60= list(Isingmodel_prob7_60[:6,1])
cv_3_60 = list(Isingmodel_prob7_60[:6,2])
magnetization_3_60 = list(Isingmodel_prob7_60[:6,3])
chi_3_60 = list(Isingmodel_prob7_60[:6,4])

energy_3_60_f = list(Isingmodel_prob7_60[9:,1])
cv_3_60_f = list(Isingmodel_prob7_60[9:,2])
magnetization_3_60_f = list(Isingmodel_prob7_60[9:,3])
chi_3_60_f = list(Isingmodel_prob7_60[9:,4])

energy_3_60_zoom = list(Isingmodel_prob7_60_zoom[:,1])
cv_3_60_zoom = list(Isingmodel_prob7_60_zoom[:,2])
magnetization_3_60_zoom = list(Isingmodel_prob7_60_zoom[:,3])
chi_3_60_zoom = list(Isingmodel_prob7_60_zoom[:,4])

for i in range(len( energy_3_60_zoom)):
    energy_3_60.append(energy_3_60_zoom[i])
    cv_3_60.append(cv_3_60_zoom[i])
    magnetization_3_60.append(magnetization_3_60_zoom[i])
    chi_3_60.append(chi_3_60_zoom[i])

for i in range(len( energy_3_60_f)):
    energy_3_60.append(energy_3_60_f[i])
    cv_3_60.append(cv_3_60_f[i])
    magnetization_3_60.append(magnetization_3_60_f[i])
    chi_3_60.append(chi_3_60_f[i])

#=======================
Isingmodel_prob7_40 = np.loadtxt("Isingmodel_prob7_40.txt")
Isingmodel_prob7_40_zoom = np.loadtxt("Isingmodel_prob7_40_zoom.txt")

#Import the values of the quantities in List (L=40):
energy_3_40= list(Isingmodel_prob7_40[:6,1])
cv_3_40 = list(Isingmodel_prob7_40[:6,2])
magnetization_3_40 = list(Isingmodel_prob7_40[:6,3])
chi_3_40 = list(Isingmodel_prob7_40[:6,4])
energy_3_40_f= list(Isingmodel_prob7_40[9:,1])
cv_3_40_f = list(Isingmodel_prob7_40[9:,2])
magnetization_3_40_f = list(Isingmodel_prob7_40[9:,3])
chi_3_40_f = list(Isingmodel_prob7_40[9:,4])

energy_3_40_zoom = list(Isingmodel_prob7_40_zoom[:,1])
cv_3_40_zoom = list(Isingmodel_prob7_40_zoom[:,2])
magnetization_3_40_zoom = list(Isingmodel_prob7_40_zoom[:,3])
chi_3_40_zoom = list(Isingmodel_prob7_40_zoom[:,4])

for i in range(len( energy_3_40_zoom)):
    energy_3_40.append(energy_3_40_zoom[i])
    cv_3_40.append(cv_3_40_zoom[i])
    magnetization_3_40.append(magnetization_3_40_zoom[i])
    chi_3_40.append(chi_3_40_zoom[i])

for i in range(len( energy_3_40_f)):
    energy_3_40.append(energy_3_40_f[i])
    cv_3_40.append(cv_3_40_f[i])
    magnetization_3_40.append(magnetization_3_40_f[i])
    chi_3_40.append(chi_3_40_f[i])


#Plot the experimental (MCMC method) energies using parallelitation for a
#lattice of L=40:
plt.plot(temperature3,energy_3_40, label='$40 X 40$',  color = "purple")
plt.plot(temperature3, energy_3_60, label='$60 X 60$', color = "darkorange")
plt.plot(temperature3, energy_3_80, label='$80 X 80$', color = "darkgreen")
plt.plot(temperature3, energy_3_100, label='$100 X 100$', color = "blue")
plt.xlabel("T (J/$k_B$)")
plt.ylabel("$\epsilon$ (J)")
plt.grid(True)
plt.legend()
plt.savefig('Isingmodel_e_c.pdf')
plt.show()

#Plot the experimental (MCMC method) heat capacity using parallelitation for a
#lattice of L=100:
plt.plot(temperature3,cv_3_40, label='$40 X 40$',  color = "purple")
plt.plot(temperature3, cv_3_60, label='$60 X 60$', color = "darkorange")
plt.plot(temperature3, cv_3_80, label='$80 X 80$', color = "darkgreen")
plt.plot(temperature3, cv_3_100, label='$100 X 100$', color = "blue")
plt.xlabel("T (J/$k_B$)")
plt.ylabel("$C_V (J k_B)$")
plt.grid(True)
plt.legend()
plt.savefig('Isingmodel_cv_c.pdf')
plt.show()

#Plot the experimental (MCMC method) magnetization using parallelitation for a
#lattice of L=100:
plt.plot(temperature3, magnetization_3_40, label='$40 X 40$',  color = "purple")
plt.plot(temperature3, magnetization_3_60, label='$60 X 60$', color = "darkorange")
plt.plot(temperature3, magnetization_3_80, label='$80 X 80$', color = "darkgreen")
plt.plot(temperature3, magnetization_3_100, label='$100 X 100$', color = "blue")
plt.xlabel("T (J/$k_B$)")
plt.ylabel("$m$")
plt.grid(True)
plt.legend()
plt.savefig('Isingmodel_m_c.pdf')
plt.show()

#Plot the experimental (MCMC method) susceptibility using parallelitation for a
#lattice of L=100:
plt.plot(temperature3,chi_3_40, label='$40 X 40$',  color = "purple")
plt.plot(temperature3, chi_3_60, label='$60 X 60$', color = "darkorange")
plt.plot(temperature3, chi_3_80, label='$80 X 80$', color = "darkgreen")
plt.plot(temperature3, chi_3_100, label='$100 X 100$', color = "blue")
plt.xlabel("T (J/$k_B$)")
plt.ylabel("$\chi$")
plt.grid(True)
plt.legend()
plt.savefig('Isingmodel_chi_c.pdf')
plt.show()


#=============================================================================
#Creating the List with the values of T_c and each lattice L:
L_rev: List[float] = [1/100, 1/80, 1/60, 1/40]
temp_c: List[float] = [2.28, 2.29, 2.3,2.32]


# Get slope, intercept from linregress() to plot y' = intercept + slope*x
(slope, intercept, rvalue, pvalue, stderr) = linregress(L_rev,temp_c)
temp_c_pred: List[float] = []
# Plot linear regression line.
for i in L_rev:
    temp_c_pred.append(intercept + slope*(i))

plt.plot(L_rev,temp_c_pred, color="green", label="Fitted line")
plt.plot(L_rev,temp_c, ".", color="green", label="Fitted line")
plt.xlabel("$L**{-1}$")
plt.ylabel("$T_c$ (J/$k_B$)")
plt.legend()
plt.savefig('Isingmodel_chi_c.pdf')
plt.show()




