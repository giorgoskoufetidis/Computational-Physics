
import numpy as np  

class EOSPath:
    """
    EOSPath is a class that represents an equation of state (EOS) for polytropic stars. 
    It allows for the computation of pressure, energy density, and other related quantities 
    based on a piecewise polytropic EOS.
    Attributes:
        segment_densities (np.ndarray): Array of densities that define the boundaries of the 
            polytropic segments.
        gammas (np.ndarray): Array of adiabatic indices (gamma) for each polytropic segment.
        K_values (np.ndarray): Array of proportionality constants (K) for each polytropic segment.
        epsilon_i (np.ndarray): Array of energy densities at the segment boundaries.
        P_i (np.ndarray): Array of pressures at the segment boundaries.
        name (str): Name of the EOS.
    Methods:
        __init__(segment_densities, gammas, K_base=1.0, name="EOS"):
            Initializes the EOSPath object with the given segment densities, adiabatic indices, 
            base proportionality constant, and name.
        _compute_K_epsilons(K_base):
            Computes the proportionality constants (K), energy densities (epsilon), and pressures 
            (P) for each polytropic segment based on the given base constant.
        pressure(rho):
            Computes the pressure for a given density (rho) based on the EOS.
        energy_density(rho):
            Computes the energy density for a given density (rho) based on the EOS.
        energy_density_from_pressure(P):
            Computes the energy density for a given pressure (P) by inverting the EOS.
        get_curve(rho_values):
            Computes the pressure and energy density curves for a given array of densities.
    """
    def __init__(self, segment_densities, gammas, K_base=1.0, name="EOS"):
        self.segment_densities = np.array(segment_densities)
        self.gammas = np.array(gammas)
        self.K_values, self.epsilon_i, self.P_i = self._compute_K_epsilons(K_base)
        self.name = name

    def _compute_K_epsilons(self, K_base):
        N = len(self.gammas)
        K_values = np.zeros(N)
        epsilon_i = np.zeros(N)
        P_i = np.zeros(N)

        rho = self.segment_densities
        gamma = self.gammas

        K_values[0] = K_base

        if gamma[0] == 1:
            P_i[0] = K_base * rho[0]
            epsilon_i[0] = rho[0]
        else:
            P_i[0] = K_base * rho[0] ** gamma[0]
            epsilon_i[0] = rho[0] + P_i[0] / (gamma[0] - 1)

        for i in range(1, N):
            rho_prev = rho[i-1]
            rho_curr = rho[i]
            gamma_prev = gamma[i-1]
            gamma_curr = gamma[i]
            K_prev = K_values[i-1]

            if gamma_prev == 1:
                P_prev = K_prev * rho_curr
                eps_prev = epsilon_i[i-1] * rho_curr / rho_prev + K_prev * rho_curr * np.log(rho_curr / rho_prev)
            else:
                P_prev = K_prev * rho_curr ** gamma_prev
                eps_prev = rho_curr + P_prev / (gamma_prev - 1)
                eps_prev += (epsilon_i[i-1] - rho_prev - K_prev * rho_prev**gamma_prev / (gamma_prev-1)) * (rho_curr / rho_prev)

            K_values[i] = P_prev / rho_curr**gamma_curr
            P_i[i] = P_prev
            epsilon_i[i] = eps_prev

        return K_values, epsilon_i, P_i

    def pressure(self, rho):
        idx = np.searchsorted(self.segment_densities, rho) - 1
        idx = np.clip(idx, 0, len(self.gammas) - 1)
        return self.K_values[idx] * rho ** self.gammas[idx]

    def energy_density(self, rho):
        rho = np.atleast_1d(rho)  # Always treat rho as array
        idx = np.searchsorted(self.segment_densities, rho) - 1
        idx = np.clip(idx, 0, len(self.gammas) - 1)

        gamma = self.gammas[idx]
        K = self.K_values[idx]
        rho_i = self.segment_densities[idx]
        eps_i = self.epsilon_i[idx]
        P_i = self.P_i[idx]

        result = np.zeros_like(rho)

        # Case where gamma == 1
        mask_gamma1 = (gamma == 1)
        result[mask_gamma1] = eps_i[mask_gamma1] * rho[mask_gamma1] / rho_i[mask_gamma1] + K[mask_gamma1] * rho[mask_gamma1] * np.log(rho[mask_gamma1] / rho_i[mask_gamma1])

        # Case where gamma != 1
        mask_not_gamma1 = ~mask_gamma1
        P = self.pressure(rho[mask_not_gamma1])
        result[mask_not_gamma1] = (eps_i[mask_not_gamma1] / rho_i[mask_not_gamma1] - P_i[mask_not_gamma1] / (eps_i[mask_not_gamma1] * (gamma[mask_not_gamma1] - 1))) * (P / K[mask_not_gamma1])**(1 / gamma[mask_not_gamma1]) + P / (gamma[mask_not_gamma1] - 1)

        if result.shape[0] == 1:
            return result[0]  # Return scalar if input was scalar
        return result


    def energy_density_from_pressure(self, P):
        for i in range(len(self.segment_densities)):
            gamma = self.gammas[i]
            K = self.K_values[i]
            if i == len(self.segment_densities) - 1 or P < self.P_i[i+1]:
                rho = (P / K)**(1 / gamma)
                return self.energy_density(rho)

        gamma = self.gammas[-1]
        K = self.K_values[-1]
        rho = (P / K)**(1 / gamma)
        return self.energy_density(rho)

    def get_curve(self, rho_values):
        P_values = self.pressure(np.array(rho_values))
        E_values = self.energy_density(np.array(rho_values))
        return np.array(E_values), np.array(P_values)
                