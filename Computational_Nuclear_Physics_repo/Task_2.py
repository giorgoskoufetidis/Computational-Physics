"""
Task 2 - Liquid-Drop Model
@author: Giorgos Koufetidis

This script simulates the liquid-drop model of the atomic nucleus, solving for mass numbers (A), neutron numbers (N), and proton numbers (Z) as functions of density (n). It uses the Fermi momentum and Coulomb terms to derive key nuclear properties.
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from prettytable import PrettyTable

# Constants
m_p = 938.272  # Proton mass in MeV/c^2
m_n = 939.565  # Neutron mass in MeV/c^2
av = 15.71511  # Volume term coefficient in MeV
ac = 0.71363   # Coulomb term coefficient in MeV
a_s = 17.53638 # Surface term coefficient in MeV
a_a = 23.37837 # Asymmetry term coefficient in MeV
cl = 3.40665 * 1e-3  # Coulomb logarithmic term coefficient in MeV
dm = m_n - m_p  # Mass difference between neutron and proton

# List of elements with atomic numbers and symbols
elements = [
    [1, 'H'], [2, 'He'], [3, 'Li'], [4, 'Be'], [5, 'B'],
    [6, 'C'], [7, 'N'], [8, 'O'], [9, 'F'], [10, 'Ne'],
    [11, 'Na'], [12, 'Mg'], [13, 'Al'], [14, 'Si'], [15, 'P'],
    [16, 'S'], [17, 'Cl'], [18, 'Ar'], [19, 'K'], [20, 'Ca'],
    [21, 'Sc'], [22, 'Ti'], [23, 'V'], [24, 'Cr'], [25, 'Mn'],
    [26, 'Fe'], [27, 'Co'], [28, 'Ni'], [29, 'Cu'], [30, 'Zn'],
    [31, 'Ga'], [32, 'Ge'], [33, 'As'], [34, 'Se'], [35, 'Br'],
    [36, 'Kr'], [37, 'Rb'], [38, 'Sr'], [39, 'Y'], [40, 'Zr'],
    [41, 'Nb'], [42, 'Mo'], [43, 'Tc'], [44, 'Ru'], [45, 'Rh'],
    [46, 'Pd'], [47, 'Ag'], [48, 'Cd'], [49, 'In'], [50, 'Sn'],
    [51, 'Sb']
]

# Functions
def pf(n):
    """Calculate the Fermi momentum for a given nucleon number density (n)."""
    pi_sq = np.pi**2
    return np.power(3 * n * pi_sq, 1 / 3)

def a_bar(n, cl):
    """Calculate the modified Coulomb term."""
    return ac - cl * pf(n)

def solution(initial_guess, n, cl):
    """
    Solve the system of equations for the liquid-drop model.
    
    Args:
        initial_guess (list): Initial guesses for the root-finding algorithm.
        n (float): Nucleon number density in fm^-3.
        cl (float): Coulomb logarithmic term coefficient.
    
    Returns:
        tuple: x1, y1, A, Z, N (solutions to the equations).
    """
    def equations_iso(vars):
        x, y = vars
        eq1 = -a_s / x**2 + 2 * a_bar(n, cl) * x * y**2
        eq2 = -dm + 2 * a_bar(n, cl) * x**2 * y - 4 * a_a * (1 - 2 * y) + np.power(y, 1 / 3) * pf(n)
        return [eq1, eq2]

    x1, y1 = fsolve(equations_iso, initial_guess)
    A = round(x1**3)
    Z = round(A * y1)
    N = A - Z
    return x1, y1, A, Z, N

def linear_approximation(n):
    """Provide a linear approximation for the initial guess."""
    x = 3.90610 + 0.03023 * pf(n)
    y = 0.45405 - 0.00419 * pf(n)
    return x, y

def find_element(Z):
    """Find the element symbol by atomic number."""
    for element in elements:
        if element[0] == Z:
            return element[1]
    return "Unknown"

# Main Function
def main():
    """Main function to compute and visualize the liquid-drop model."""
    n_values = np.logspace(-2.2, 3.34, 30)  # Nucleon number densities
    initial_guess = [3, 1]  # Initial guesses for solving equations
    
    # Lists to store results
    A_list, Z_list, N_list = [], [], []
    x1_list, y1_list, pf_list = [], [], []
    found_elements = []

    for n in n_values:
        x, y = linear_approximation(n)
        x1, y1, A, Z, N = solution(initial_guess, n, cl)
        A_list.append(A)
        Z_list.append(Z)
        N_list.append(N)
        x1_list.append(round(x1, 3))
        y1_list.append(round(y1, 3))
        pf_list.append(round(pf(n), 2))
        found_elements.append(find_element(Z))

    # Display results in a table
    table = PrettyTable()
    table.add_column("Pf (MeV)", pf_list)
    table.add_column("x", x1_list)
    table.add_column("y", y1_list)
    table.add_column("A", A_list)
    table.add_column("N", N_list)
    table.add_column("Z", Z_list)
    table.add_column("Element", found_elements)
    print(table)

    # Plot A, N, and Z vs Pf
    plt.figure()
    plt.plot(pf_list, A_list, label='A')
    plt.plot(pf_list, N_list, label='N')
    plt.plot(pf_list, Z_list, label='Z')
    plt.title("A, Z, N vs Pf (MeV)")
    plt.xlabel("Pf (MeV)")
    plt.legend()
    plt.grid()
    plt.show()

    # Plot A, N, and Z vs n
    plt.figure()
    plt.semilogx(n_values, A_list, label='A')
    plt.semilogx(n_values, N_list, label='N')
    plt.semilogx(n_values, Z_list, label='Z')
    plt.title(r'$A, Z, N \ \mathrm{vs} \ n \ (\mathrm{fm}^{-3})$')
    plt.xlabel(r'$n \ (\mathrm{fm}^{-3})$')
    plt.legend()
    plt.grid()
    plt.show()

    # Plot x and y vs Pf
    plt.figure()
    plt.plot(pf_list, x1_list, label='x')
    plt.plot(pf_list, y1_list, label='y')
    plt.title("x and y vs Pf (MeV)")
    plt.xlabel("Pf (MeV)")
    plt.ylabel("x, y")
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    main()
