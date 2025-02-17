"""
@author: Giorgos Koufetidis
Computational Nuclear Physics 
Task 4 
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as smp

# Define constants
atom_radius = 6.69  # Radius of the atom in femtometers
surface_atom = 0.45  # Surface thickness parameter in femtometers
r_0 = 0.6777  # Parameter for normalization, not used in this code
E = [126, 183]  # Energies in MeV
Z = 79  # Atomic number
h_c = 197  # Planck's constant times speed of light in MeV*fm
e = 1.602e-19  # Elementary charge in Coulombs
e0 = 55.26349406 * 10**-3  # Permittivity constant in e^2*MeV^-1*fm^-1

# Define symbolic variables
t = smp.symbols('t', real=True)
q, r, R, a = smp.symbols('q, r, R, a', real=True, positive=True)

# Define the normalized radial density function for the nucleus
pr_Normal = 4 * smp.pi * r**2 / (1 + smp.exp((r - R) / a))

# Define a sample function x for integration purposes
x = smp.sin(t)

# Function to calculate normalization constant r0
def r_Normalization(R_, a_):
    integral_value = smp.Integral(pr_Normal.subs({a: a_, R: R_}), (r, 0, smp.oo)).evalf()
    return Z / integral_value

# Calculate integral of x * exp(I*q*r*cos(t)) over t from 0 to pi
I1 = smp.integrate(x * smp.exp(smp.I * q * r * smp.cos(t)), (t, 0, smp.pi)).simplify()
print(I1)

# Store the normalization constant r0
r0 = r_Normalization(atom_radius, surface_atom)

# Plotting the normalized charge density
radious = np.linspace(0, 10, 100)
def density_(radious):
    return r0 / (1 + np.exp((radious - atom_radius) / surface_atom))

plt.plot(radious, density_(radious))
plt.grid()
plt.xlabel("r(fm)")
plt.ylabel("ρ(r)")
plt.title("Electric Charge Density")
plt.show()

# Function to define the radial probability distribution
def pr(r, a, R, r0):
    return r0 / (1 + smp.exp((r - R) / a))

# Function to calculate the momentum transfer q
def q(E, t):
    t_rad = smp.pi * t / 180  # Convert angle t from degrees to radians
    return (2 * E / h_c) * smp.sin(t_rad / 2)

# Function to calculate the integral I1
def I1(q, R, a, r0):
    integral = smp.Integral(r * smp.sin(q * r) * pr(r, surface_atom, atom_radius, r0), (r, 0, 10))
    return 4 * smp.pi / q * integral

# Function to calculate the differential cross-section ds/dw
def ds_dw(E, t, R, a, r0):
    q_value = q(E, t)
    i1_value = I1(q_value, atom_radius, surface_atom, r0)**2
    return (E / (2 * smp.pi))**2 * (1 / h_c)**4 * 1 / q_value**4 * (1 / e0)**2 * i1_value

# Define the scattering angles theta
theta = np.linspace(20, 130, 90)

# Calculate and plot the differential cross-section for each energy
for i in range(len(E)):
    data = [ds_dw(E[i], t, atom_radius, surface_atom, r0) for t in theta]
    plt.semilogy(theta, data, "-.", label=f"E={E[i]}")
    plt.grid()
    plt.xlabel(r"$\theta$ (degrees)")
    plt.ylabel(r"$\frac{d\sigma}{d\Omega} ({fm}^2/steradian)$")
    plt.title(r"Differential Cross-Section $\frac{d\sigma}{d\Omega}$ vs scattering angle (θ)")
    plt.grid(True)
    plt.legend()
plt.show()





















"""

sol = 2 * smp.pi**2 * r_0 * a**2 * smp.cos(q * R) / smp.sinh(smp.pi * q * a)
theta_values = np.linspace(0,130,200)
ds_dw_list = []
for i in theta_values:
    ds_dw_list.append(ds_dw(q.subs(t_rad,i),6.63,0.45,r_0))
  

def ds_dwFinal(E,q):
    return (E/2*np.pi)**2*(1/h_c)**4*(1/q**4)*(1/e0)**2*ds_dw(q.subs(t_rad,i),6.63,0.45,r_0)**2

for i in theta_values:
    ds_dw_list.append(ds_dwFinal(E,q.subs(t_rad,i)))
plt.plot(theta_values,ds_dw_list)


ds_dw1 = (E / (2 * smp.pi))**2 * (1 / h_c)**4 * (1 / q**4) * (1/ e0)**2 * sol**2
ds_dw2 = ds_dw1.subs({R:atom_radious, a:serfice_atom })
ds_dw3 =smp.lambdify(t,ds_dw2,'numpy')

plt.semilogy(theta_values,ds_dw3(theta_values))
plt.xlabel('θ (degrees)')
plt.ylabel('Differential Cross-section')
plt.title('Differential Cross-section vs. θ')
plt.grid(True)
plt.show()
"""