# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 08:21:15 2024

@author: User
"""
import numpy as np

from CRUST import CRUST
from EOS import EOS
from scipy.integrate import solve_ivp
from functools import partial
from tqdm import tqdm
from multiprocessing import Pool

'''
# Example EOS function (from your thesis
crust = CRUST(1e-5)
print(crust.equation())


# Define TOV right-hand side equations
def tov_rhs(r, z):
    M = z[0]
    P = z[1]
    
    if P >= 0.184:
        eos = EOS(P)
        epsilon = eos.MDI_1()
    else:
        crust = CRUST(P)
        epsilon = crust.equation()
    if P <= 0:
        return [0, 0]  # Stop integration when pressure drops below zero
    dM_dr = 11.2e-6 * r**2 * epsilon
    dP_dr = - 1.474 * (epsilon * M / (r**2)) * (1 + P/epsilon) * (1 + 11.2 * (10**(-6)) * (r** 3) * P/M) * ((1 - 2.948 * M/r)** (-1))
    dzdr = [dM_dr,dP_dr]
    return dzdr

step = 0.0001
G = 6.674 * 10 ** (- 33)
n = 100
ic1 = np.arange(1.5,5,0.5)
ic2 = np.arange(5,10,1)
ic = np.concatenate((ic1,ic2),axis = None)

j = 0
for i in ic:
    z0 = [0.000000000001,i]
    P0 = 1.0
    rmin = 0.00000001
    rmax = 0.01
    Mf = np.array([])
    Pf = np.array([])
    R = np.array([])
    z0_old = np.array([])
    j += 1
    print(j)
    
    while P0 > 1e-12:
        res = solve_ivp(tov_rhs,(rmin,rmax),z0,method = 'LSODA',atol = 1e-12,rtol = 1e-8)
        z0_old = z0[1]
        z0[0] = res.y[0][~np.isnan(res.y[0])][-1]
        z0[1] = res.y[1][~np.isnan(res.y[1])][-1]
        if z0[0] < 0:
            break
        if z0_old == z0[1]:
            break
        
        rmin = res.t[-1]
        rmax = rmin + 0.001
        P0 = z0[1]
        Mf = np.append(Mf,res.y[0][~np.isnan(res.y[0])][-1])
        Pf = np.append(Pf,res.y[1][~np.isnan(res.y[1])][-1])
        R = np.append(R,res.t[~np.isnan(res.y[1])])
        if Pf[-1] < 0:
            idx = np.argwhere(Pf < 0)[0,0]
            Pf = np.delete(Pf,np.s_[idx::], 0)
            Mf = np.delete(Mf,np.s_[idx::], 0)
            R = np.delete(R,np.s_[idx::], 0)
        


'''
'''
# Define TOV right-hand side equations
def tov_rhs(r, z, index):
    M = z[0]
    P = z[1]
    if P >= 0.184:
        eos = EOS(P)
        epsilon = eos.models()[index]
    else:
        crust = CRUST(P)
        epsilon = crust.equation()
    if P <= 0:
        return [0, 0]  
    dM_dr = 11.2e-6 * r**2 * epsilon
    dP_dr = -1.474 * (epsilon * M / r**2) * (1 + P / epsilon) * (1 + 11.2e-6 * r**3 * P / M) * (1 - 2.948 * M / r)**(-1)
    return [dM_dr, dP_dr]

# EOS models
models = [
    "MDI_1", "MDI_2", "MDI_3", "MDI_4", "NLD", "HHJ_1", "HHJ_2", "SKa",
    "SkI4", "HLPS_2", "HLPS_3", "SCVBB", "WFF_1", "WFF_2", "PS", "W", "BGP", "BL_1",
    "BL_2", "DH", "APR_1"
]

# Initial conditions
ic1 = np.arange(1.5, 5, 0.1) #1.5,5,0.1
ic2 = np.arange(5, 1201, 1) # 5, 1201,1
init = np.concatenate((ic1, ic2), axis=None)

initial_conditions = {model: init for model in models}

# Process function for a single model
def process_initial_conditions(model_name, ic, index):
    mass_list = []
    radius_list = []
    pressure_list = []

    for i in ic:
        z0 = [1e-12, i]  # Initial [Mass, Pressure]
        P0 = i
        rmin = 1e-8
        rmax = 0.01

        Mf = np.array([])  # Mass
        Pf = np.array([])  # Pressure
        R = np.array([])   # Radius

        while P0 > 1e-12:
            rhs_with_index = partial(tov_rhs, index=index)  
            res = solve_ivp(rhs_with_index, (rmin, rmax), z0, method='LSODA', atol=1e-12, rtol=1e-8)
            if res.success:
                # Append results from this step
                R_step = res.t[~np.isnan(res.t)]
                Mf_step = res.y[0][~np.isnan(res.y[0])]
                Pf_step = res.y[1][~np.isnan(res.y[1])]

                if any(Pf_step < 0):
                    idx = np.argmax(Pf_step < 0)
                    R_step = R_step[:idx]
                    Mf_step = Mf_step[:idx]
                    Pf_step = Pf_step[:idx]

                R = np.append(R, R_step)
                Mf = np.append(Mf, Mf_step)
                Pf = np.append(Pf, Pf_step)

                z0 = [Mf[-1], Pf[-1]]
                rmin = R[-1]
                rmax = rmin + 0.001
                P0 = z0[1]
            else:
                break

        if len(R) > 0 and len(Mf) > 0:
            mass_list.append(Mf[-1])  # Final mass
            radius_list.append(R[-1])  # Final radius
            pressure_list.append(Pf[-1])

    return model_name, mass_list, radius_list, pressure_list

# Multiprocessing wrapper
def parallel_process_model(model_name, index):
    ic = initial_conditions[model_name]
    return process_initial_conditions(model_name, ic, index)

if __name__ == "__main__":
    with open("TOV_results_for_multiple_models.txt", "w") as f:
        with Pool(processes=2) as pool:
            # Progress bar for models
            for model_name, mass_list, radius_list, pressure_list in tqdm(
                pool.imap_unordered(
                    partial(parallel_process_model, index=0),  
                    models
                ),
                total=len(models),
                desc="Processing Models",
            ):
                f.write(f"# model: {model_name} \n")
                f.write("Mass Radius Pressure \n")
                for m, r, p in zip(mass_list, radius_list, pressure_list):
                    f.write(f"{m:.8e} {r:.8e} {p:.8e}\n")
'''


def tov_rhs(r, z, index=None):
    M = z[0]
    P = z[1]
    if P >= 0.184:
        eos = EOS(P)  
        epsilon = eos.models()[index]
    else:
        crust = CRUST(P)  
        epsilon = crust.equation()
    if P <= 0:
        return [0, 0]
    dM_dr = 11.2e-6 * r**2 * epsilon
    dP_dr = -1.474 * (epsilon * M / r**2) * (1 + P / epsilon) * (1 + 11.2e-6 * r**3 * P / M) * (1 - 2.948 * M / r)**(-1)
    return [dM_dr, dP_dr]


models = [
    "MDI_1", "MDI_2", "MDI_3", "MDI_4", "NLD", "HHJ_1", "HHJ_2", "SKa",
    "SkI4", "HLPS_2", "HLPS_3", "SCVBB", "WFF_1", "WFF_2", "PS", "W", "BGP", "BL_1",
    "BL_2", "DH", "APR_1"
]


ic1 = np.arange(1.5, 5, 0.1)  # 1.5, 5, 0.1
ic2 = np.arange(5, 1201, 1)   # 5, 1201, 1
init = np.concatenate((ic1, ic2), axis=None)

initial_conditions = {model: init for model in models}


def process_initial_conditions(model_name, ic, index):
    mass_list = []
    radius_list = []
    pressure_list = []

    for i in ic:
        z0 = [1e-12, i]  
        P0 = i
        rmin = 1e-8
        rmax = rmin + 0.01  

        Mf = np.array([])  
        Pf = np.array([])  
        R = np.array([])   

        while P0 > 1e-12:
            res = solve_ivp(lambda r, z: tov_rhs(r, z, index), (rmin, rmax), z0, method='LSODA', atol=1e-12, rtol=1e-8)
            if res.success:
                # Append results from this step
                R_step = res.t[~np.isnan(res.t)]
                Mf_step = res.y[0][~np.isnan(res.y[0])]
                Pf_step = res.y[1][~np.isnan(res.y[1])]

                if any(Pf_step < 0):
                    idx = np.argmax(Pf_step < 0)
                    R_step = R_step[:idx]
                    Mf_step = Mf_step[:idx]
                    Pf_step = Pf_step[:idx]

                R = np.append(R, R_step)
                Mf = np.append(Mf, Mf_step)
                Pf = np.append(Pf, Pf_step)

                if len(R) == 0 or len(Mf) == 0:
                    break

                z0 = [Mf[-1], Pf[-1]]
                rmin = R[-1]
                rmax = rmin + 0.001  
                P0 = z0[1]
            else:
                break

        if len(R) > 0 and len(Mf) > 0:
            mass_list.append(Mf[-1])  
            radius_list.append(R[-1])  
            pressure_list.append(Pf[-1])

    return model_name, mass_list, radius_list, pressure_list


def parallel_process_model(args):
    model_name, index = args  
    
    ic = initial_conditions[model_name]
    return process_initial_conditions(model_name, ic, index)

if __name__ == "__main__":
    with open("TOV_results_for_multiple_models.txt", "w") as f:
        with Pool(processes=8) as pool:
           
            model_indices = [[model_name, idx] for idx, model_name in enumerate(models)]

            # Progress bar for models
            for model_name, mass_list, radius_list, pressure_list in tqdm(
                pool.imap_unordered(parallel_process_model, model_indices),
                total=len(models),
                desc="Processing Models",
            ):
                f.write(f"# model: {model_name} \n")
                f.write("Mass Radius Pressure \n")
                if mass_list and radius_list and pressure_list:
                    for m, r, p in zip(mass_list, radius_list, pressure_list):
                        f.write(f"{m:.8e} {r:.8e} {p:.8e}\n")

