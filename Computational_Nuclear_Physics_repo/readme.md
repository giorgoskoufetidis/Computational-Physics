
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


# Liquid Drop Model of the Nucleus

## Description

This script implements the **Liquid Drop Model** to estimate the binding energy of atomic nuclei. The model is based on the semi-empirical mass formula (SEMF), which approximates the nuclear binding energy using five key terms.

## Features

- Computes nuclear binding energy using the **semi-empirical mass formula (SEMF)**.
- Considers **volume, surface, Coulomb, asymmetry, and pairing terms**.
- Supports input for **mass number $A$ and atomic number $Z$**.
- Calculates **binding energy per nucleon** for stability analysis.
- Plots **binding energy trends** across different nuclei.

## Mathematical Formulation

The **binding energy** $B(A, Z)$ is given by:

$$
B(A, Z) = a_v A - a_s A^{2/3} - a_c \frac{Z(Z-1)}{A^{1/3}} - a_a \frac{(A-2Z)^2}{A} + \delta(A, Z)
$$

where:

- $A$ = mass number (total number of nucleons).
- $Z$ = atomic number (number of protons).
- $a_v$, $a_s$, $a_c$, and $a_a$ are **empirical coefficients**.
- $\delta(A, Z)$ is the **pairing term**, defined as:

$$
\delta(A, Z) = \begin{cases}
    + a_p A^{-3/4} & \text{if } Z \text{ and } N \text{ are both even} \\ % & is your "\tab"-like command (it's a tab alignment character)
    0 & \text{if } A \text{ is odd} \\ % & is your "\tab"-like command (it's a tab alignment character)
    - a_p A^{-3/4} & \text{if } Z \text{ and } N \text{ are both odd.}
\end{cases}
$$



where $N = A - Z$ is the neutron number.

## Execution

The script calculates binding energy for various nuclei based on user input.

### Steps:

1. Input values for **mass number $A$** and **atomic number $Z$**.
2. Compute individual contributions to **binding energy**.
3. Apply the **pairing term** based on nuclear structure.
4. Compute the **binding energy per nucleon**:

$$
B/A = \frac{B(A, Z)}{A}
$$

5. Visualize trends of **binding energy vs. mass number**.

