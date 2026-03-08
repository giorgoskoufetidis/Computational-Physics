import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from scipy.integrate import solve_ivp
import csv
import os
from itertools import product
from concurrent.futures import ProcessPoolExecutor, as_completed
from CRUST import CRUST
import sympy as sp
from sympy.utilities.lambdify import lambdify
# Constants
r_sat = 2.7e14  # g/cm^3
c = 2.998e8  # m/s
MeV_to_J = 1.60218e-13  # 1 MeV in Joules

class CRUST_Differential:

    def eoscsym(self):
        pp = sp.Symbol('pp')
        eosfunc = 0.00873 + 103.17338 * (1 - sp.exp(-pp / 0.38527)) + 7.34979 * (1 - sp.exp(-pp / 0.01211))
        eosfunc_diff = eosfunc.diff(pp)
        return lambdify(pp, eosfunc_diff, 'numpy')

    def eoscsym2(self):
        pp = sp.Symbol('pp')
        eosfunc = 0.00015 + 0.00203 * (1 - sp.exp(-pp / 344827.5)) + 0.10851 * (1 - sp.exp(-pp / 7692.3076))
        eosfunc_diff = eosfunc.diff(pp)
        return lambdify(pp, eosfunc_diff, 'numpy')

    def eoscsym3(self):
        pp = sp.Symbol('pp')
        eosfunc = 0.0000051 * (1 - sp.exp(-pp * 0.2373 * 10**10)) + 0.00014 * (1 - sp.exp(-pp * 0.4021 * 10**8))
        eosfunc_diff = eosfunc.diff(pp)
        return lambdify(pp, eosfunc_diff, 'numpy')

    def eoscsym4(self):
        pp = sp.Symbol('pp')
        log_pp = sp.log(pp, 10)
        eosfunc = 10 * (31.93753 + 10.82611 * log_pp + 1.29312 * log_pp**2 +
                        0.08014 * log_pp**3 + 0.00242 * log_pp**4 + 0.000028 * log_pp**5)
        eosfunc_diff = eosfunc.diff(pp)
        return lambdify(pp, eosfunc_diff, 'numpy')

    def crust_equation1(self, P):
        return self.eoscsym()(P)

    def crust_equation2(self, P):
        return self.eoscsym2()(P)

    def crust_equation3(self, P):
        return self.eoscsym3()(P)

    def crust_equation4(self, P):
        return self.eoscsym4()(P)

    def equation(self, P):
        if 9.34375e-5 <= P <= 0.184:
            return self.crust_equation1(P)
        elif 4.1725e-8 <= P < 9.34375e-5:
            return self.crust_equation2(P)
        elif 1.44875e-11 <= P < 4.1725e-8:
            return self.crust_equation3(P)
        elif P < 1.44875e-11:
            return self.crust_equation4(P)
        else:
            raise ValueError("Pressure out of bounds")

# Conversion function
def conv_to_MeV(value):
    result = value * 1e3  # to kg/m^3
    result = result * c ** 2  # to J/m^3
    result = result / MeV_to_J  # to MeV/m^3
    result = result * 1e-45  # to MeV/fm^3
    return result
# HLPS_2 function
def HLPS_2(P):
    return 172.858 * (1 - np.exp(-P / 22.8644)) + 2777.75 * (1 - np.exp(-P / 1909.97)) + 161.553
def HLPS_3(P):
    return 131.811 * (1 - np.exp(-P/4.41577)) + 924.143 * (1 - np.exp(- P/523.736)) + 81.5682

# EOSPath class
class EOSPath:
    def __init__(self, P0, E0, segment_densities, gammas, name="EOS"):
        self.segment_densities = np.array(segment_densities)
        self.gammas = np.array(gammas)
        self.P0 = P0
        self.E0 = E0
        self.name = name
        self.Pi = []
        self.Ei = []
        self.K_list = []
        self.linear_start = None
        self.linear_eps = None

    def calculate_K(self, P, rho, gamma):
        return P / rho ** gamma

    def calculate_P(self, rho, K, gamma):
        return K * rho ** gamma

    def calculate_E(self, P, Ki, gamma, rho_i, P_i, E_i):
        if gamma == 1:
            return (E_i / rho_i) * P / Ki + np.log(1 / rho_i) * P - P * np.log(Ki / P)
        else:
            return (E_i / rho_i - P_i / (rho_i * (gamma - 1))) * (P / Ki) ** (1 / gamma) + P / (gamma - 1)

    def dE_dP(self, P, Ki, gamma, rho_i, P_i, E_i):
        if gamma == 1:
            return E_i / (rho_i * Ki) + np.log(1 / rho_i) + 1 - np.log(Ki / P)
        else:
            return (1 / (Ki * gamma)) * (E_i / rho_i - P_i / (rho_i * (gamma - 1))) * (P / Ki) ** (1 / gamma - 1) + 1 / (gamma - 1)

    def EOS_Linear(self, P, P_transition, E_transition):
        
        return E_transition + (P - P_transition)
    
    def calculate_parameters(self):
        # Initialize storage lists
        Pi_vals = [self.P0]
        Ei_vals = [self.E0]
        Ki_vals = []

        n = len(self.gammas)

        # Calculate Pi, Ei, Ki arrays for each segment
        for i in range(n):
            gamma = self.gammas[i]
            rho_i = self.segment_densities[i]
            rho_next = self.segment_densities[i + 1]

            K_i = self.calculate_K(Pi_vals[i], rho_i, gamma)
            Ki_vals.append(K_i)

            P_next = self.calculate_P(rho_next, K_i, gamma)
            Ei_next = self.calculate_E(P_next, K_i, gamma, rho_i, Pi_vals[i], Ei_vals[i])

            Pi_vals.append(P_next)
            Ei_vals.append(Ei_next)

        # Fine pressure sampling for slope checking
        P_range = np.linspace(Pi_vals[0], Pi_vals[-1], 10_000)
        P_tr = None
        E_tr = None

        # Function to compute slope dE/dP for piecewise polytropic EOS at arbitrary P
        def slope_at_P(P):
            # Find which segment P belongs to
            for i in range(len(Pi_vals) - 1):
                if Pi_vals[i] <= P <= Pi_vals[i + 1]:
                    return self.dE_dP(P, Ki_vals[i], self.gammas[i], self.segment_densities[i], Pi_vals[i], Ei_vals[i])
            # If beyond known segments, assume linear EOS slope = 1
            if self.linear_start and P > self.linear_start:
                return 1
            return 1  # fallback

        # Find the first pressure where slope < 1 (causality violation)
        for P in P_range:
            slope = slope_at_P(P)
            if slope < 1:
                P_tr = P
                # Calculate corresponding energy at transition
                for i in range(len(Pi_vals) - 1):
                    if Pi_vals[i] <= P_tr <= Pi_vals[i + 1]:
                        E_tr = self.calculate_E(P_tr, Ki_vals[i], self.gammas[i], self.segment_densities[i], Pi_vals[i], Ei_vals[i])
                        break
                # Truncate lists at transition
                for idx, P_bound in enumerate(Pi_vals):
                    if P_bound >= P_tr:
                        Pi_vals = Pi_vals[:idx + 1]
                        Ei_vals = Ei_vals[:idx + 1]
                        Pi_vals[-1] = P_tr
                        Ei_vals[-1] = E_tr
                        break
                break

        # Save results to the object
        self.Pi = Pi_vals
        self.Ei = Ei_vals
        self.K_list = Ki_vals
        self.linear_start = P_tr
        self.linear_eps = E_tr
 
 
    def get_energy_from_pressure(self, P):
        for i in range(len(self.Pi) - 1):
            if self.Pi[i] <= P <= self.Pi[i + 1]:
                return self.calculate_E(P, self.K_list[i], self.gammas[i], self.segment_densities[i], self.Pi[i], self.Ei[i])
        if P > self.Pi[-1]:
            if self.linear_start is None or self.linear_eps is None:
                # fallback if linear start/eps not set
                return self.Ei[-1]  # or some other reasonable fallback
            return self.linear_eps + (P - self.linear_start)

crust = CRUST()  
crust_diff = CRUST_Differential()
# TOV equation
def tov_rhs(r, z, eos_object):
    M, P, y = z
    if P <= 0:
        return [0, 0, 0]
    elif  0.184 < P <= 1.722:
        epsilon = HLPS_2(P)
        dP_de = (HLPS_2(P + 1e-8) - HLPS_2(P - 1e-8)) / (2e-8)
    elif P > 1.722:
        eps_plus = eos_object.get_energy_from_pressure(P + 1e-8)
        eps_minus = eos_object.get_energy_from_pressure(P - 1e-8)

        if eps_plus is None or eps_minus is None:
            # fallback to avoid crash
            return [0, 0, 0]

        epsilon = eos_object.get_energy_from_pressure(P)
        dP_de = (eps_plus - eps_minus) / (2e-8)
    else:
        epsilon = crust.equation(P)
        crust_plus = crust.equation(P + 1e-11)
        crust_minus = crust.equation(P - 1e-11)
        dP_de = (crust_plus - crust_minus) / (2e-11)
        if epsilon is None or crust_plus is None or crust_minus is None:
            print(f"[Error] CRUST energy undefined at P={P:.5e}")
            return [0, 0, 0]

    F =  (1-1.474 * 11.2 * ( 10 **(-6)) * (r** 2) * (epsilon - P)) * ((1-2.948 * M / r) ** (-1))
    J = 1.474 * 11.2 * (10**(-6)) * (r**2) * (5 * epsilon + 9 * P + (epsilon + P) / (1 / dP_de)) *  ((1-2.948 * M / r) **(-1)) - 6 * \
    ((1-2.948 * M / r) ** (-1)) - 4 * ((1.474 ** 2) * (M ** 2)/(r ** 2)) *  ((1 + 11.2 *  (10 **(-6)) * (r** 3) * (P / M)) ** 2) * ((1-2.948 * M / r) ** (-2))
    dM_dr = 11.2e-6 * r ** 2 * epsilon
    dP_dr = -1.474 * (epsilon * M / r ** 2) * (1 + P / epsilon) * (1 + 11.2e-6 * r ** 3 * P / M) * (1 - 2.948 * M / r) ** (-1)
    dyrdt = (-y * y - y * F - J) / r
    return [dM_dr, dP_dr, dyrdt]

def compute_k2(beta, yR):
    k2 = (8 *(beta ** 5)/5) * ((1-2 * beta)** 2) * (2-yR +2 * beta * (yR - 1)) * (2 *
     beta * (6-3 * yR +3 * beta *  (5 * yR -8))+ 4 * beta ** 3 * (13-11 * yR + beta * (3 *
     yR - 2) + 2 * beta ** 2 * (1 + yR)) + 3 * ((1 - 2 * beta)** 2) * (2 - yR + 2 * beta * (yR - 1)) * np.log(1 - 2 * beta)) ** (-1)
    return k2

def stop_when_pressure_small(r, y):
    return y[1] - 1e-10
stop_when_pressure_small.terminal = True
stop_when_pressure_small.direction = -1
# Solve one TOV star
def solve_star(P_central, eos_object):
    z0 = [1e-12, P_central, 2.0]  # M, P, y
    r_span = (1e-6, 2e6)

    sol = solve_ivp(lambda r, z: tov_rhs(r, z, eos_object),
                    r_span,
                    z0,
                    method='RK45',
                    atol=1e-10,
                    rtol=1e-8,
                    events=stop_when_pressure_small,
                    max_step=0.01)

    if sol.status == 1:
        R_star = sol.t_events[0][0]
        M_star = sol.y_events[0][0][0]
        yR = sol.y_events[0][0][2]
        beta = 1.474 * M_star / R_star
        k2 = compute_k2(beta, yR)
        return M_star, R_star, P_central, k2, yR
    else:
        return None

def process_model(args):
    model_name, P0, E0, segment_densities, gammas, initial_pressures = args

    eos_object = EOSPath(P0, E0, segment_densities, gammas, name=model_name)
    eos_object.calculate_parameters()

    output_path = os.path.join("TOV_results_tidal_HLPS2", f"{model_name}_TOV.csv")
    has_data = False  # Track if anything succeeded

    # Create and open file before loop
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Mass", "Radius", "Pressure", "k2", "yR", "Type"])  # header

        for P_center in initial_pressures:
            sol = solve_star(P_center, eos_object)
            if sol is not None:
                M_star, R_star, P_center, k2, yR = sol
                writer.writerow([M_star, R_star, P_center, k2, yR, 0])
                has_data = True
                print(f"Model:{model_name} | Mass:{M_star:.4f} | Radius:{R_star:.4f} | Pc={P_center:.2f} | k2={k2:.4f} | yR={yR:.4f}")
            else:
                print(f"Warning: Model {model_name} Pc={P_center:.2f} did not converge.")

    if not has_data:
        print(f"⚠️ No successful TOV solution for {model_name}. File was created but is empty.")

    return model_name, 1 if has_data else 0


def get_segment_rhos(gamma_1, num_segments):
    ρ_sat_MeV = conv_to_MeV(r_sat)
    if gamma_1 <= 2:
        ρ_high = 16 * ρ_sat_MeV
    elif 2 < gamma_1 <= 3:
        ρ_high = 12 * ρ_sat_MeV
    else:
        ρ_high = 9 * ρ_sat_MeV

    log_ρ0 = np.log(ρ_sat_MeV)
    log_ρn = np.log(7.5 * ρ_sat_MeV)
    segment_rhos = np.exp(np.linspace(log_ρ0, log_ρn, num_segments + 1))
    segment_rhos[-1] = ρ_high
    return segment_rhos

if __name__ == "__main__":
    P_sat = 1.722
    # P_sat = 1.722 # in MeV/fm^3
    # P_sat = 2.816 # in MeV/fm^3
    E_sat = HLPS_2(P_sat)
    segments = 5
    gamma_options = [1, 2, 3, 4]

    all_gamma_paths = list(product(gamma_options, repeat=segments))
    eos_jobs = []
    for idx, gamma_path in enumerate(all_gamma_paths):
        gamma_1 = gamma_path[0]
        segment_rhos = get_segment_rhos(gamma_1, segments)
        model_name = f"EOS_{idx + 1}"
        eos = EOSPath(P_sat, E_sat, segment_rhos, gamma_path, name=model_name)
        eos.calculate_parameters()
        P_final = eos.Pi[-1]
        ic1 = np.arange(1.8, 5, 0.1)
        ic2 = np.arange(5, 3501 if P_final > 800 else P_final, 1)
        initial_pressures = np.concatenate((ic1, ic2))
        eos_jobs.append((model_name, P_sat, E_sat, segment_rhos, gamma_path, initial_pressures))

    os.makedirs("TOV_results_tidal_HLPS2", exist_ok=True)
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_model, job) for job in eos_jobs]
        for future in as_completed(futures):
            model_name, success = future.result()
            print(f"✅ Finished {model_name}" if success else f"❌ {model_name} had no results")
