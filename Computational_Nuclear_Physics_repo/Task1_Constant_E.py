"""
Task 1 - Energy-Length Calculation with Constant E
@author: Giorgos Koufetidis

This script calculates energy (E), length (L), and potential (V0) using randomized inputs for V0 and L, while treating E as a constant. It performs error analysis, linear regression, and visualization.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import bisect
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import random

# Constants
m_c = 939  # Mass constant in MeV
hbar_c = 197  # Reduced Planck constant in MeV·fm
E_real = 2.224  # Constant binding energy in MeV
gamma = 0.326  # Constant for wavefunction computation


def cot(x):
    """Calculate the cotangent of x."""
    return 1 / np.tan(x)


def compute_root(V0, L):
    """
    Compute the root (ξ1) for a given V0 and L.

    Args:
        V0 (float): Potential in MeV.
        L (float): Length in fm.

    Returns:
        float: ξ1 root value.
    """
    constant = 2 * m_c * V0 * L**2 / hbar_c**2
    valid_upper_bound = np.sqrt(constant)  # Ensure square root is valid

    if valid_upper_bound <= 1:
        raise ValueError(f"No valid domain for V0={V0}, L={L}.")

    # Define the function for root finding
    def f(x):
        try:
            sqrt_term = np.sqrt(constant - x**2)
            return sqrt_term + x * cot(x)
        except ValueError:
            return np.nan

    # Solve using bisection
    root = bisect(f, 1, min(3, valid_upper_bound))
    return root


def compute_error(V0, L, j1, E_constant):
    """
    Compute the error between the calculated and real binding energy.

    Args:
        V0 (float): Potential in MeV.
        L (float): Length in fm.
        j1 (float): Root value ξ1.
        E_constant (float): Constant binding energy.

    Returns:
        float: Error percentage.
    """
    E_calculated = V0 - (hbar_c**2 * j1**2) / (2 * m_c * L**2)
    return abs((E_calculated - E_constant) / E_constant) * 100


def compute_wavefunctions(L, k):
    """
    Compute wavefunction coefficients A and B.

    Args:
        L (float): Length in fm.
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


def main():
    """Main function to calculate and analyze energy-length relationships."""
    random.seed(42)  # Seed for reproducibility
    results = []  # List to store results for table

    # Generate random V0 and L values and compute results
    for _ in range(10000):  # Perform multiple iterations
        V0 = random.uniform(20, 90)  # V0 range: 20-90 MeV
        L = random.uniform(1.2, 4)  # L range: 1.2-4 fm

        try:
            j1 = compute_root(V0, L)
            error = compute_error(V0, L, j1, E_real)

            if error < 1:  # Store only results with error < 1%
                k = np.sqrt(2 * m_c * (V0 - E_real) / hbar_c**2)
                A, B = compute_wavefunctions(L, k)
                results.append((V0, L, j1, E_real, error, k, A, B))

        except ValueError:
            continue  # Skip invalid cases

    # Prepare data for visualization
    V0_values, L_values = zip(*[(row[0], row[1]) for row in results])

    # Perform linear regression
    model = LinearRegression()
    model.fit(np.array(V0_values).reshape(-1, 1), L_values)
    L_pred = model.predict(np.array(V0_values).reshape(-1, 1))
    r_squared = r2_score(L_values, L_pred)

    # Print results summary
    print(f"{'V0 (MeV)':<10} {'L (fm)':<10} {'ξ1':<10} {'E (MeV)':<10} {'Error (%)':<10} {'k (fm^-1)':<10} {'A':<10} {'B':<10}")
    for row in results[:10]:  # Display first 10 results
        print(f"{row[0]:<10.3f} {row[1]:<10.3f} {row[2]:<10.3f} {row[3]:<10.3f} {row[4]:<10.3f} {row[5]:<10.3f} {row[6]:<10.3f} {row[7]:<10.3f}")

    # Plot results
    plt.scatter(V0_values, L_values, color='blue', label='Data')
    plt.plot(V0_values, L_pred, color='red', label=f'Linear Fit (R²={r_squared:.2f})')
    plt.xlabel('V0 (MeV)')
    plt.ylabel('L (fm)')
    plt.title('Relationship Between V0 and L')
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()
