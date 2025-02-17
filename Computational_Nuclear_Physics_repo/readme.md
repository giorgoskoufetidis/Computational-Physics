
# Energy-Length Calculation with Constant E

## Description

This script calculates energy $E$, length $L$, and potential $V_0$ using randomized inputs for $V_0$ and $L$, while treating $E$ as a constant. It performs error analysis, linear regression, and visualization to explore the relationship between potential and length.

## Features

- Generates random values for potential $V_0$ and length $L$.
- Computes the root $\xi_1$ using the bisection method.
- Calculates energy error percentage.
- Determines wavefunction coefficients ($A$ and $B$).
- Performs linear regression on $V_0$ and $L$.
- Visualizes the relationship between $V_0$ and $L$.

## Mathematical Formulation

### Root Calculation

The root $\xi_1$ is determined by solving:

$$
\sqrt{\frac{2mV_0L^2}{\hbar^2} - \xi^2} + \xi \cot(\xi) = 0
$$

using the bisection method. The constant term:

$$
C = \frac{2m_c V_0 L^2}{\hbar_c^2}
$$

is used to define the valid range for $\xi$.

### Energy Calculation

The calculated energy $E$ is given by:

$$
E = V_0 - \frac{\hbar_c^2 \xi_1^2}{2m_c L^2}
$$

The percentage error is computed as:

$$
\text{Error} = \left| \frac{E - E_{\text{real}}}{E_{\text{real}}} \right| \times 100
$$

### Wavefunction Coefficients

The coefficients $A$ and $B$ are determined based on the wavefunction continuity conditions. The ratio $A/B$ is given by:

$$
\frac{A}{B} = \frac{e^{-\gamma L}}{\sin(kL)}
$$

where $\gamma$ and $k$ are defined as:

$$
\gamma = \sqrt{\frac{2m_c E}{\hbar_c^2}}, \quad k = \sqrt{\frac{2m_c (V_0 - E)}{\hbar_c^2}}
$$

The normalization condition leads to the final expressions for $A$ and $B$.

## Execution

The main function generates random values for $V_0$ and $L$, calculates roots, errors, and wavefunction coefficients, and performs linear regression.

### Steps:

1. Randomly generate $V_0$ (20-90 MeV) and $L$ (1.2-4 fm).
2. Compute $\xi_1$ and check for a valid domain.
3. Calculate error percentage and filter results (error < 1%).
4. Compute wavefunction coefficients $A$ and $B$.
5. Perform linear regression on $V_0$ vs. $L$.
6. Visualize the results using matplotlib.


# Liquid-Drop Model

## Author
Giorgos Koufetidis

## Description
This script simulates the Liquid-Drop Model of the atomic nucleus, solving for mass numbers (\(A\)), neutron numbers (\(N\)), and proton numbers (\(Z\)) as functions of nucleon number density (\(n\)). It uses the Fermi momentum and Coulomb terms to derive key nuclear properties.

## Mathematical Formulation

### Fermi Momentum
The Fermi momentum \(p_f\) is calculated as:
$$
p_f = (3\pi^2 n)^{\frac{1}{3}}
$$
where \(n\) is the nucleon number density.

### Modified Coulomb Term
The modified Coulomb term is defined as:
$$
\bar{a}(n) = a_c - c_l \cdot p_f
$$
where \(a_c\) is the Coulomb coefficient and \(c_l\) is the logarithmic term coefficient.

### Nuclear Mass Equations
The system of equations used to find the equilibrium conditions is:

1. **Surface Energy Balance:**
   $$
   -\frac{a_s}{x^2} + 2\bar{a}(n)xy^2 = 0
   $$

2. **Proton-Neutron Balance:**
   $$
   -\Delta m + 2\bar{a}(n)x^2y - 4a_a(1 - 2y) + y^{\frac{1}{3}} p_f = 0
   $$

Here:
- \(x\) and \(y\) are dimensionless parameters,
- \(a_s\) is the surface term coefficient,
- \(a_a\) is the asymmetry energy coefficient,
- \(\Delta m\) is the neutron-proton mass difference.

### Computation of \(A\), \(Z\), and \(N\)
The solutions of the above equations yield:
$$
A = x^3, \quad Z = Ay, \quad N = A - Z
$$
where:
- \(A\) is the mass number,
- \(Z\) is the proton number,
- \(N\) is the neutron number.

## Execution
The main function of the script performs the following steps:
1. Generates logarithmically spaced nucleon number densities \(n\).
2. Uses a linear approximation for the initial guess.
3. Solves the system of equations for each \(n\) value.
4. Retrieves corresponding element symbols.
5. Displays the results in a table.
6. Visualizes trends with several plots.

