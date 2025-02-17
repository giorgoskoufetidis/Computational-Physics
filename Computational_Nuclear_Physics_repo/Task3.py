"""
@author: Giorgos Koufetidis
Computational Nuclear Physics 
Task 3 
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt 

# Define the species in the decay chain
species = ['244Ct', '240Cm', '236Pu', '232U', '228Th', '224Ra', '220Rn', '216Po', '212Pb']
# Define the half-life of each species in seconds
half_life= [1140, 2.33e6, 104400, 2.271e9, 5.992e8, 3197e6, 56, 0.15, 39600]

# Initial conditions: only 244Ct is present initially with 1000 units, all others are 0
N0 = [1000, 0, 0, 0, 0, 0, 0, 0, 0]

# Calculate the decay constants (λ) for each species
l = [np.log(2) / hl for hl in half_life]

# Define the system of differential equations for the decay chain
def deriv(N, t):
    dn_dt = [0] * len(half_life)  # Initialize the list with zeros
    dn_dt[0] = -l[0] * N[0]       # The decay rate of the first species
    for i in range(1, len(half_life)):
        dn_dt[i] = l[i-1] * N[i-1] - l[i] * N[i]  # The decay rate for each subsequent species
    return dn_dt

# Define the time points where the solution is computed (split into three segments)
t1 = np.linspace(0, 1e5, int(1e5))       # First segment: from 0 to 1e5 seconds
t2 = np.linspace(1e5, 1e8, int(1e3))     # Second segment: from 1e5 to 1e8 seconds
t3 = np.linspace(1e8, 1e13, int(1e5))    # Third segment: from 1e8 to 1e13 seconds

# Solve the differential equations for each time segment
N1 = odeint(deriv, N0, t1).T             # Solve for the first time segment
N2 = odeint(deriv, N1[:, -1], t2).T      # Use the end conditions of N1 as initial for N2
N3 = odeint(deriv, N2[:, -1], t3).T      # Use the end conditions of N2 as initial for N3

# Concatenate the time arrays and the solutions
t = np.concatenate((t1, t2[1:], t3[1:]))
N = np.concatenate((N1, N2[:, 1:], N3[:, 1:]), axis=1)

# Plot the concentrations for each species over time
for i in range(len(species)):
    plt.semilogx(t, N[i], label=species[i])

# Add title, labels, legend, and grid to the plot
plt.title("General decay chain")
plt.xlabel('t(sec)')
plt.legend(fontsize=7)
plt.grid()
plt.show()














"""
for i in range(len(half_life)):
    plt.semilogx(t1,N1[i],label = species[i])

    
    
plt.title("General decay chain")
plt.xlabel('t(sec)')
plt.legend()
plt.show()
for i in range(len(half_life)):
    plt.semilogx(t2,N2[i],label=species[i])
    
plt.title("General decay chain")
plt.xlabel('t(sec)')
plt.legend()
plt.show()
for i in range(len(half_life)):
    plt.semilogx(t3,N3[i],label=species[i])
    
plt.title("General decay chain")
plt.xlabel('t(sec)')
plt.legend()
plt.show()




def deriv(N, t):
     Return dX/dt for each of the species. 
    return (-l[0] * N[0],                              
             l[0] * N[0] - l[1] * N[1],   
             l[1] * N[1] - l[2] * N[2],                 
             l[2] * N[2] - l[3] * N[3],                 
             l[3] * N[3] - l[4] * N[4],
             l[4] * N[4] - l[5] * N[5],
             l[5] * N[5] - l[6] * N[6],
             l[6] * N[6] - l[7] * N[7],
             l[7] * N[7] - l[8] * N[8],     
            )


plt.semilogx(t,N[0],label='244Ct')
plt.semilogx(t,N[1],label = '240Cm')
plt.semilogx(t,N[2],label = '236Pu')
plt.semilogx(t,N[3],label = '232U')
plt.semilogx(t,N[4],label = '228Th')
plt.semilogx(t,N[5],label = '224Ra')
plt.semilogx(t,N[6],label = '220Rn')
plt.semilogx(t,N[7],label = '216Po')
plt.semilogx(t,N[8],label = '212Pb')
plt.title("General decay chain")
plt.xlabel('t(sec)')
plt.legend()


pylab.plot(t, N[0], label=r'$[\,\mathrm{^{212}Pb}]$', c='k', ls='--')
pylab.plot(t, N[1], label=r'$[\,\mathrm{^{212}Bi}]$', c='k', ls=':')
pylab.plot(t, N[4], label=r'$[\,\mathrm{^{208}Pb}]$', c='k')
pylab.plot(t, N[2]*1.e2, label=r'$[\,\mathrm{^{208}Tl}] \times 10^{2}$',
           c='gray', ls='--')
pylab.plot(t, N[3]*1.e5, label=r'$[\,\mathrm{^{212}Po}] \times 10^{5}$',
           c='gray', ls=':')

Css = (1-np.exp(-l[0]*t))     # all intermediates in steady-state
pylab.plot(t, Css, c='k', ls='-.', label=r'$[\,\mathrm{^{208}Pb}]$ (SSA)')

pylab.legend()
pylab.xlabel(r'$t\;/\mathrm{s}$')
pylab.ylabel(r'conc. /arb. units')
pylab.show()
"""