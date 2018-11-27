from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

# parameters
vsg1 = 1.202
vsg2 = 1
vsn1 = 0.856
vsn2 = 1
vsfr1 = 2.8
vsfr2 = 2.8
vex = 0
vsf = 0.6
va = 20
vi = 3.3
kdg = 1
kdn = 1
kdfr = 1
kdf = 0.09
Kag1 = 0.28
Kag2 = 0.55
Kan = 0.55
Kafr = 0.5
Kaf = 5
Kig = 2
Kin1 = 0.28
Kin2 = 2
Kifr = 0.5
Ka = 0.7
Ki = 0.7
Kd = 2
FpTristablility = 0.06
Fp0 = 0.03  # as seen in fig. 2
FPend = 0.11  # as seen in fig. 2


r = 3
s = 4
q = 4
u = 3
v = 4
w = 4
z = 4


# simulation time
start = 0  # min
end = 50 # max
plotPoints = 10
# timeGrid = np.linspace(start, end, plotPoints) #right now  not used


# returns the ODE as the operator d/dt(), where Fp is a constant (for phase space diagrams)
def dif_Fp_const(t, x):
    G = x[0]
    N = x[1]
    FR = x[2]
    ERK = x[3]
    dG = (vsg1 * ((ERK ** r) / ((Kag1 ** r) + (ERK ** r))) + vsg2 * ((G ** s) / ((Kag2 ** s) + (G ** s)))) * (
            (Kig ** q) / ((Kig ** q) + (N ** q))) - kdg * G
    dN = (vsn1 * ((Kin1 ** u) / ((Kin1 ** u) + (ERK ** u))) + vsn2 * ((N ** v) / ((Kan ** v) + (N ** v)))) * (
            (Kin2 ** w) / ((Kin2 ** w) + (G ** w))) - kdn * N
    dFR = vsfr1 * ((Kifr) / (Kifr + N)) + vsfr2 * ((G) / (Kafr + G)) - kdfr * FR
    dERK = va * FR * ((FpTristablility) / (Kd + FpTristablility)) * ((1 - ERK) / (Ka + 1 - ERK)) - vi * ((ERK) / (Ki + ERK))  # vin->vi ; typo?
    return [dG, dN, dFR, dERK]


# returns the ODE as the operator d/dt(), where Fp is a parameter (for bifurcation diagrams)
def dif_Fp_param(t, x):
    G = x[0]
    N = x[1]
    FR = x[2]
    ERK = x[3]
    Fp = x[4]
    dG = (vsg1 * ((ERK ** r) / ((Kag1 ** r) + (ERK ** r))) + vsg2 * ((G ** s) / ((Kag2 ** s) + (G ** s)))) * (
            (Kig ** q) / ((Kig ** q) + (N ** q))) - kdg * G
    dN = (vsn1 * ((Kin1 ** u) / ((Kin1 ** u) + (ERK ** u))) + vsn2 * ((N ** v) / ((Kan ** v) + (N ** v)))) * (
            (Kin2 ** w) / ((Kin2 ** w) + (G ** w))) - kdn * N
    dFR = vsfr1 * ((Kifr) / (Kifr + N)) + vsfr2 * ((G) / (Kafr + G)) - kdfr * FR
    dERK = va * FR * ((Fp) / (Kd + Fp)) * ((1 - ERK) / (Ka + 1 - ERK)) - vi * ((ERK) / (Ki + ERK))  # vin->vi ; typo?
    dFp = ((FPend-Fp0)/(50-start))*t  # is the maximal time
    return [dG, dN, dFR, dERK, dFp]

# solve with V0 as list (different starting points). Max N and G are hardcoded as 2.2
def varying_v0_solver(step_size):
    results_list = []
    step = step_size
    while step <= 2.2:
        results_list.append(solve_ivp(dif_Fp_const, (start, end), [step, 0.2, 0, 0]))  # order for v0: [G0,N0,FR0,ERK0]
        results_list.append(solve_ivp(dif_Fp_const, (start, end), [0.2, step, 0, 0]))  # order for v0: [G0,N0,FR0,ERK0]
        step += step_size
    return results_list


resultsPhaseSpace = varying_v0_solver(0.2)
resultBifurcation = solve_ivp(dif_Fp_param, (start, end), [0, 0, 0, 0, 0.03])  # order for v0: [G0,N0,FR0,ERK0,Fp])

G_data_bif = resultBifurcation.y[0]
N_data_bif = resultBifurcation.y[1]
Fp_data_bif = resultBifurcation.y[4]
t_data_bif = resultBifurcation.t

# square grid of 4 subplots
plt.figure(figsize=(9,6)) # generate slightly larger figure

# Gata6(Fgf4) plot - bifurcation diagram
plot1 = plt.subplot(2,2,1)
plot1.plot(Fp_data_bif, G_data_bif, label='upper left')
plot1.set_xlabel('Fgf4', fontsize=15)
plot1.set_ylabel('Gata6', fontsize=15)
# Nanog(Fgf4) plot - bifurcation diagram
plot2 = plt.subplot(2,2,2)
plot2.plot(Fp_data_bif, N_data_bif, label='upper right')
plot2.set_xlabel('Fgf4', fontsize=15)
plot2.set_ylabel('Nanog', fontsize=15)
# Nanog(Gata6) plot - phase space diagram
plot3 = plt.subplot(2,2,3)
for result in resultsPhaseSpace:
    G_data_ps = result.y[0]
    N_data_ps = result.y[1]
    plot3.plot(G_data_ps, N_data_ps, label='lower left')
    plot3.set_xlabel('Gata6', fontsize=15)
    plot3.set_ylabel('Nanog', fontsize=15)
# Fgf4(t) plot from the bifurcation plots
plot4 = plt.subplot(2,2,4)
plt.plot(t_data_bif, Fp_data_bif, label='lower right')
plot4.set_xlabel('Time', fontsize=15)
plot4.set_ylabel('Fgf4', fontsize=15)
plt.show()
