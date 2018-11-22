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
Fp = 4 # Fgf4 concentration ?

r = 3
s = 4
q = 4
u = 3
v = 4
w = 4
z = 4

# simulation time
start = 0  # min
end = 100  # min
plotPoints = 1000
timeGrid = np.linspace(start, end, plotPoints)


# returns the ODE as the operator d/dt()
def dif(x, t):
    G, N, FR, ERK = x
    dG = (vsg1 * ((ERK ** r) / ((Kag1 ** r) + (ERK ** r))) + vsg2 * ((G ** s) / ((Kag2 ** s) + (G ** s)))) * ((Kig ** q) / ((Kig ** q)+(N ** q))) - kdg *G
    dN = (vsn1*((Kin1 ** u)/((Kin1 ** u)+(ERK ** u))) + vsn2 * ((N ** v)/((Kan ** v) + (N ** v)))) * ((Kin2 ** w)/((Kin2 ** w)+(G ** w))) - kdn * N
    dFR = vsfr1 * ((Kifr)/(Kifr + N)) + vsfr2 * ((G)/(Kafr + G)) - kdfr * FR
#    dERK = va * FR * ((Fp)/())
#   return [dG, dN, dFR, dERK]


#v0 = [G0,N0,FR0,ERK0]
#result = solve_ivp(f, y0, timeGrid)
