import numpy as np
from prettytable import PrettyTable
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

hbar_c = 197.326

class QEOS:
    def __init__(self, ms, delta, B):
        self.ms = ms
        self.delta = delta
        self.B = B

    def alpha(self):
        a = -self.ms**2 / 6 + 2 * self.delta**2 / 3
        return a

    def energy_density(self, P):
        a = self.alpha()
        m = np.sqrt(-3 * a + np.sqrt(4 / 3 * np.pi**2 * (self.B + P) * hbar_c**3 + 9 * a**2))
        epsilon = 3 * P + 4 * self.B - 9 * a * m**2 / (np.pi**2 * hbar_c**3)
        return epsilon


def stability_function(delta, ms):
    mn = 939
    Betta = -ms**2 * mn**2 / (12 * np.pi**2) + delta**2 * mn**2 / (3 * np.pi**2) + mn**4 / (108 * np.pi**2)
    return Betta / (197.326)**3

def stability_point(ms, mn, delta):
    """Calculate the upper limit of B based on Equation (10), with proper unit conversion."""
    conversion_factor = (197.326)**3
    limit = -ms**2 * mn**2 / (12 * np.pi**2) + delta**2 * mn**2 / (3 * np.pi**2) + mn**4 / (108 * np.pi**2)
    return limit / conversion_factor


delta = [50, 100, 150]
ms = 95  # MeV
B = list(range(60, 181, 5))
eos = {}
mn = 939  # MeV

i = 1
for b in B:
    for d in delta:
        stability_limit = stability_point(ms, mn, d)
        if b < stability_limit:
            qeos = QEOS(ms, d, b)
            eos[f"CLF_{i}"] = qeos
            i += 1


table = PrettyTable()
table.field_names = ["CLF", "B", "D", "ms"]
for clf, qeos in eos.items():
    table.add_row([clf, qeos.B, qeos.delta, qeos.ms])
print(table)
ms = np.linspace(0, 300, 1000)
plt.plot(stability_function(50, ms), ms, label="Δ = 50")
plt.plot(stability_function(100, ms), ms, label="Δ = 100")
plt.plot(stability_function(150, ms), ms, label="Δ = 150")
plt.legend()
for i in range(1, len(eos) + 1):
    plt.scatter(eos[f"CLF_{i}"].B, eos[f"CLF_{i}"].ms, label=f"CLF_{i}")

plt.xlabel(r"B(MeV/fm$^3$)")
plt.ylabel("Strange Quark Mass (MeV)")
plt.title("Stability Function and EOS Points")
plt.grid(True, linestyle="--", linewidth=0.5)
plt.show()

x = np.linspace(0, 1200, 1200)
y = eos["CLF_1"].energy_density(x)
plt.plot(x, y, label="CLF_1 EOS")
plt.xlabel(r"r(fm)")
plt.ylabel(r"$\epsilon(r)$ (MeV)")
plt.legend()
plt.show()


def tov_rhs(r, z):
    M = z[0]  
    P = z[1]  

    epsilon = eos["CLF_1"].energy_density(P)
    if P <= 0: 
        return [0, 0]

    # TOV equations
    dM_dr = 11.2e-6 * r**2 * epsilon
    dP_dr = -1.474 * (epsilon * M / r**2) * (1 + P / epsilon) * \
            (1 + 11.2e-6 * r**3 * P / M) * (1 - 2.948 * M / r)**(-1)

    return [dM_dr, dP_dr]


def solve_tov_clf1():
 
    with open("TOV_results_CLF_1.txt", "w") as f:
        f.write("# model: CLF_1\n")
        f.write("Mass Radius Pressure\n")

        initial_pressures = np.concatenate((np.arange(1.5, 5, 0.1), np.arange(5, 1201, 1)), axis=None)
        i = 0
        for P0 in initial_pressures:
            z0 = [1e-12, P0]  # Initial conditions: [Mass, Pressure]
            rmin = 1e-8
            rmax = 0.01

            while P0 > 1e-12:
                res = solve_ivp(tov_rhs, (rmin, rmax), z0, method='LSODA', atol=1e-12, rtol=1e-8)
                if res.success:
                    R_step = res.t[~np.isnan(res.t)]
                    Mf_step = res.y[0][~np.isnan(res.y[0])]
                    Pf_step = res.y[1][~np.isnan(res.y[1])]

                    # Stop if pressure becomes negative
                    if any(Pf_step < 0):
                        idx = np.argmax(Pf_step < 0)
                        R_step = R_step[:idx]
                        Mf_step = Mf_step[:idx]
                        Pf_step = Pf_step[:idx]

                    # Write results directly to file
                    for j in range(len(R_step)):
                        f.write(f"{Mf_step[j]:.8e} {R_step[j]:.8e} {Pf_step[j]:.8e}\n")

                    z0 = [Mf_step[-1], Pf_step[-1]]
                    rmin = R_step[-1]
                    rmax = rmin + 0.001
                    P0 = z0[1]
                else:
                    break

            i += 1
            print(f"Completed iteration {i}.")

if __name__ == "__main__":
    solve_tov_clf1()
