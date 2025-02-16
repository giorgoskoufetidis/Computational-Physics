"""
Task 1 - Energy and Length Computation
@author: Giorgos Koufetidis
This script computes energy and length values for given potential parameters.
It uses numerical methods to solve physics-based equations and generates tabular and graphical results.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import bisect


# Constants
E = 2.224  # Binding energy in MeV
m_c = 939  # Mass constant in MeV
hbar_c = 197.327  # Reduced Planck constant in MeV·fm
gamma = 0.326  # Constant used in wavefunction calculation


def cot(x):
    """Calculate the cotangent of x."""
    return 1 / np.tan(x)


def calculate_values(pair):
    """
    Calculate L, L_square, V0, and ξ1 for a given V·L^2 pair.

    Args:
        pair (float): V·L^2 value.

    Returns:
        tuple: L (fm), L_square, V0 (MeV), ξ1.
    """
    constant = (2 * m_c / hbar_c ** 2) * pair
    valid_upper_bound = np.sqrt(constant)  # Maximum valid x value

    # Ensure the domain for bisect is valid
    if valid_upper_bound <= 1:
        raise ValueError(f"No valid domain for pair={pair}.")

    lower_bound = 1  # Fixed lower bound
    upper_bound = min(3, valid_upper_bound)  # Dynamic upper bound based on `pair`

    # Define the function to find roots
    def f(x):
        try:
            sqrt_term = np.sqrt(constant - x**2)
            return sqrt_term + x * cot(x)
        except ValueError:
            return np.nan

    # Ensure function values are valid at bounds
    if np.isnan(f(lower_bound)) or np.isnan(f(upper_bound)):
        raise ValueError(f"Invalid function values for pair={pair}. Check bounds.")

    # Find root using bisection method
    try:
        root_bisection = bisect(f, lower_bound, upper_bound)
    except ValueError as e:
        raise ValueError(f"Bisection failed for pair={pair}: {e}")

    j1 = root_bisection
    L = np.sqrt(pair - (hbar_c ** 2) / (2 * m_c) * j1**2) / np.sqrt(E)
    L_square = pair / E - (hbar_c ** 2) / (2 * m_c) * j1**2 / E
    V0 = pair / L**2
    return round(L, 3), round(L_square, 3), round(V0, 3), round(j1, 3)


def compute_wavefunctions(L, V0, k):
    """
    Compute wavefunction coefficients A and B.

    Args:
        L (float): Length in fm.
        V0 (float): Potential in MeV.
        k (float): Wave number.

    Returns:
        tuple: Coefficients A and B.
    """
    ratio_A_B = np.exp(-gamma * L) / np.sin(k * L)
    integral_a = quad(lambda x: np.sin(k * x)**2, 0, L)[0]
    integral_b = quad(lambda x: np.exp(-2 * gamma * x), L, np.inf)[0]
    z = integral_a * ratio_A_B**2 + integral_b
    B = round(np.sqrt(1 / z), 3)
    A = round(ratio_A_B * B, 3)
    return A, B


def plot_n_r(A, B, L, k, V0):
    """
    Plot n(r) for given A, B, L, k, and V0.

    Args:
        A (float): Coefficient A.
        B (float): Coefficient B.
        L (float): Length in fm.
        k (float): Wave number.
        V0 (float): Potential in MeV.
    """
    r1 = np.linspace(0, L, 1000)
    r2 = np.linspace(L, 4 * L, 1000)
    y1 = A * np.sin(k * r1)
    y2 = B * np.exp(-gamma * r2)

    plt.plot(r1, y1, label="n(r) for 0<r<L")
    plt.plot(r2, y2, label="n(r) for r>L")
    plt.legend()
    plt.ylim(-round(A, 1) - 1, round(A, 1) + 1)
    plt.title(f"n(r) for L={round(L, 2)} and V0={round(V0, 2)}")
    plt.axhline(y=0, color='black')
    plt.xlabel("r (fm)")
    plt.ylabel("n(r)")
    plt.grid()
    plt.show()


def main():
    """Main function to calculate and plot results."""
    pairs = [60, 100, 150, 200, 250, 300, 350, 400, 450]
    print(f"{'V·L^2 (MeV·fm)':<15} {'ξ1':<10} {'L (fm)':<10} {'V0 (MeV)':<10} {'k (fm^-1)':<10} {'A':<10} {'B':<10}")

    for pair in pairs:
        try:
            L, _, V0, j1 = calculate_values(pair)
            k = np.sqrt((2 * m_c * (V0 - E)) / hbar_c**2)
            A, B = compute_wavefunctions(L, V0, k)
            print(f"{pair:<15} {j1:<10} {L:<10} {V0:<10} {round(k, 3):<10} {A:<10} {B:<10}")
            plot_n_r(A, B, L, k, V0)
            
        except ValueError as e:
            print(f"Error for pair={pair}: {e}")


if __name__ == "__main__":
    main()
