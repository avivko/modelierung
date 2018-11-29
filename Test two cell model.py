from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from time import gmtime, strftime


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
FpTristability = 0.06
Fp0 = 0.03  # as seen in fig. 2
FPend = 0.11  # as seen in fig. 2
FR0 = 2.8
ERK0 = 0.25


r = 3
s = 4
q = 4
u = 3
v = 4
w = 4
z = 4


# simulation time
start = 0  # min
end = 50  # max
plotPoints = 500
timeGrid = np.linspace(start, end, plotPoints)

#####two cell model:######







# Inital condintions
G10 = 0
G20 = 0
N20 = 0
N10 = 0
FR10 = 2.8
FR20 = 2.8
ERK20 = 0.25
ERK10 = 0.25
F0 = 0.066
gam = 0.03  # gamma deviation around the avergage extracellular concentration of Fgf4
def dif_twocell(t, x):
    G1 = x[0]
    G2 = x[1]
    N1 = x[2]
    N2 = x[3]
    FR1 = x[4]
    FR2 = x[5]
    ERK1 = x[6]
    ERK2 = x[7]
    Fs1 = x[8]
    Fs2 = x[9]
    Fp1 = x[10]
    Fp2 = x[11]
    F = x[12]

    dG1 = (vsg1 * ((ERK1 ** r) / ((Kag1 ** r) + (ERK1 ** r))) + vsg2 * ((G1 ** s) / ((Kag2 ** s) + (G1 ** s)))) * (
            (Kig ** q) / ((Kig ** q) + (N1 ** q))) - kdg * G1
    dG2 = (vsg1 * ((ERK2 ** r) / ((Kag1 ** r) + (ERK2 ** r))) + vsg2 * ((G2 ** s) / ((Kag2 ** s) + (G2 ** s)))) * (
            (Kig ** q) / ((Kig ** q) + (N2 ** q))) - kdg * G2
    dN1 = (vsn1 * ((Kin1 ** u) / ((Kin1 ** u) + (ERK1 ** u))) + vsn2 * ((N1 ** v) / ((Kan ** v) + (N1 ** v)))) * (
            (Kin2 ** w) / ((Kin2 ** w) + (G1 ** w))) - kdn * N1
    dN2 = (vsn1 * ((Kin1 ** u) / ((Kin1 ** u) + (ERK2 ** u))) + vsn2 * ((N2 ** v) / ((Kan ** v) + (N2 ** v)))) * (
            (Kin2 ** w) / ((Kin2 ** w) + (G2 ** w))) - kdn * N2
    dFR1 = vsfr1 * ((Kifr) / (Kifr + N1)) + vsfr2 * ((G1) / (Kafr + G1)) - kdfr * FR1
    dFR2 = vsfr1 * ((Kifr) / (Kifr + N2)) + vsfr2 * ((G2) / (Kafr + G2)) - kdfr * FR2
    dERK1 = va * FR1 * ((Fp1) / (Kd + Fp1)) * ((1 - ERK1) / (Ka + 1 - ERK1)) - vi * (
            (ERK1) / (Ki + ERK1))  # vin->vi ; typo?
    dERK2 = va * FR2 * ((Fp2) / (Kd + Fp2)) * ((1 - ERK2) / (Ka + 1 - ERK2)) - vi * (
            (ERK2) / (Ki + ERK2))
    dFs1 = vsf * ((N1 ** z) / ((Kaf ** z)) + (N1 ** z)) - kdf * Fs1 + vex
    dFs2 = vsf * ((N2 ** z) / ((Kaf ** z)) + (N2 ** z)) - kdf * Fs2 + vex
    dFp1 = (1 - gam) * F
    dFp2 = (1 + gam) * F
    dF = (dFp1 + dFp2)/(2)




    return [dG1, dG2, dN1, dN2, dFR1, dFR2, dERK1, dERK2, dFs1, dFs2, dFp1, dFp2, dF ]



#solving

resulttwocell = solve_ivp(dif_twocell,(start, end), [G10, G20, N10, N20, FR10, FR20, ERK10, ERK20, 0, 0, 0, 0, F0], t_eval=timeGrid)




#plotting
#plot A
ax = plt.axes()

ax.set_xlabel('Time', fontsize = 15)
ax.set_ylabel('Gata6, Nanog', fontsize = 15)
ax.set_ylim(bottom = 0, top = 2.2)
ax.set_xlim(left = 0, right = 50)
xdata = resulttwocell.t
ydata = resulttwocell.y[0]
plt.plot(xdata, ydata)
xdata = resulttwocell.t
ydata = resulttwocell.y[2]
plt.plot(xdata, ydata)
plt.show()

#plot B
ax = plt.axes()
ax.set_xlabel('Time', fontsize = 15)
ax.set_ylabel('Gata6, Nanog', fontsize = 15)
ax.set_ylim(bottom = 0, top = 2.2)
ax.set_xlim(left = 0, right = 50)
xdata = resulttwocell.t
ydata = resulttwocell.y[1]
plt.plot(xdata, ydata)
xdata = resulttwocell.t
ydata = resulttwocell.y[3]
plt.plot(xdata, ydata)
plt.show()




ax = plt.axes()
ax.set_xlabel('FGF4', fontsize = 15)
ax.set_ylabel('Time', fontsize = 15)
ax.set_ylim(bottom = 0, top = 0.075)
ax.set_xlim(left = 0, right = 50)
xdata = resulttwocell.t
ydata = resulttwocell.y[10]
plt.plot(xdata, ydata)
xdata = resulttwocell.t
ydata = resulttwocell.y[11]
plt.plot(xdata, ydata)
plt.show()


