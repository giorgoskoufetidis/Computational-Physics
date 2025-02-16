# Energy-Length Calculation with Constant E

## Description

This script calculates energy (E), length (L), and potential (V0) using randomized inputs for V0 and L, while treating E as a constant. It performs error analysis, linear regression, and visualization to explore the relationship between potential and length.

## Features

- Generates random values for potential (V0) and length (L).
- Computes the root \( \  xi_1 \) using the bisection method.
- Calculates energy error percentage.
- Determines wavefunction coefficients (A and B).
- Performs linear regression on V0 and L.
- Visualizes the relationship between V0 and L.

## Mathematical Formulation

### Root Calculation

The root \( \xi_1 \) is determined by solving:

\[
\sqrt{\frac{2mV_0L^2}{\hbar^2} - \xi^2} + \xi \cot(\xi) = 0
\]

using the bisection method. The constant term:

\[
C = \frac{2m_c V_0 L^2}{\hbar_c^2}
\]

is used to define the valid range for \( \xi \).

### Energy Calculation

The calculated energy \( E \) is given by:

\[
E = V_0 - \frac{\hbar_c^2 \xi_1^2}{2m_c L^2}
\]

The percentage error is computed as:

\[
\text{Error} = \left| \frac{E - E_{\text{real}}}{E_{\text{real}}} \right| \times 100
\]

### Wavefunction Coefficients

The coefficients A and B are determined based on the wavefunction continuity conditions. The ratio \( A/B \) is given by:

\[
\frac{A}{B} = \frac{e^{-\gamma L}}{\sin(kL)}
\]

where \( \gamma \) and \( k \) are defined as:

\[
\gamma = \sqrt{\frac{2m_c E}{\hbar_c^2}}, \quad k = \sqrt{\frac{2m_c (V_0 - E)}{\hbar_c^2}}
\]

The normalization condition leads to the final expressions for A and B.

## Execution

The main function generates random values for V0 and L, calculates roots, errors, and wavefunction coefficients, and performs linear regression.

### Steps:

1. Randomly generate V0 (20-90 MeV) and L (1.2-4 fm).
2. Compute \( \xi_1 \) and check for valid domain.
3. Calculate error percentage and filter results (error < 1%).
4. Compute wavefunction coefficients A and B.
5. Perform linear regression on V0 vs. L.
6. Visualize the results using matplotlib.

To execute the script, run:

```bash
python script_name.py
