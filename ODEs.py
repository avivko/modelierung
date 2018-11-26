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
Kd = 3
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
end = 100  # max
plotPoints = 100
# timeGrid = np.linspace(start, end, plotPoints) #right now  not used


# returns the ODE as the operator d/dt()
def dif(t, x):
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
    dFp = Fp0+((FPend-Fp0)/(end-start))*t
    return [dG, dN, dFR, dERK, dFp, t]


# solve with V0 as list (different starting points). Max N and G are hardcoded as 2.2
def varying_v0_solver(step_size):
    results_list = []
    step = step_size
    while step <= 2.2:
        results_list.append(solve_ivp(dif, (start, end), [step, 0, 0, 0, 0, 0]))  # order for v0: [G0,N0,FR0,ERK0,Fp,t]
        results_list.append(solve_ivp(dif, (start, end), [0, step, 0, 0, 0, 0]))  # order for v0: [G0,N0,FR0,ERK0,Fp,t]
        step += step_size
    return results_list


results = varying_v0_solver(0.2)

# square grid of 4 subplots
plt.figure(figsize=(9,6)) # generate slightly larger figure

for result in results:
    G_data = result.y[0]
    N_data = result.y[1]
    FR_data = result.y[2]
    ERK_data = result.y[3]
    Fp_data = result.y[4]
    t_data = result.y[5]
    # access first 3 one with ax.
    plot1 = plt.subplot(2,2,1)
    plot1.plot(Fp_data, G_data, label='upper left')
    plot1.set_xlabel('Fgf4', fontsize=15)
    plot1.set_ylabel('Gata6', fontsize=15)
    plot2 = plt.subplot(2,2,2)
    plot2.plot(Fp_data, N_data, label='upper right')
    plot2.set_xlabel('Fgf4', fontsize=15)
    plot2.set_ylabel('Nanog', fontsize=15)
    plot3 = plt.subplot(2,2,3)
    plot3.plot(G_data, N_data, label='lower left')
    plot3.set_xlabel('Gata6', fontsize=15)
    plot3.set_ylabel('Nanog', fontsize=15)
    # access the last with plt. (linear)
    plot4 = plt.subplot(2,2,4)
    plt.plot(t_data, Fp_data, label='lower right')
    plot4.set_xlabel('Time', fontsize=15)
    plot4.set_ylabel('Fgf4', fontsize=15)
plt.show()