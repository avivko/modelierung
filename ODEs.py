from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from time import gmtime, strftime

timeAndDate = strftime("%Y-%m-%d %H:%M:%S", gmtime())

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
    dERK = va * FR * ((FpTristablility) / (Kd + FpTristablility)) * ((1 - ERK) / (Ka + 1 - ERK)) - vi * (
            (ERK) / (Ki + ERK))  # vin->vi ; typo?
    return [dG, dN, dFR, dERK]

# returns the ODE as the operator d/dt(), where Fp is a parameter is varied (for phase space movie)
def dif_Fp_var(t, x):
    G = x[0]
    N = x[1]
    FR = x[2]
    ERK = x[3]
    dG = (vsg1 * ((ERK ** r) / ((Kag1 ** r) + (ERK ** r))) + vsg2 * ((G ** s) / ((Kag2 ** s) + (G ** s)))) * (
            (Kig ** q) / ((Kig ** q) + (N ** q))) - kdg * G
    dN = (vsn1 * ((Kin1 ** u) / ((Kin1 ** u) + (ERK ** u))) + vsn2 * ((N ** v) / ((Kan ** v) + (N ** v)))) * (
            (Kin2 ** w) / ((Kin2 ** w) + (G ** w))) - kdn * N
    dFR = vsfr1 * ((Kifr) / (Kifr + N)) + vsfr2 * ((G) / (Kafr + G)) - kdfr * FR
    dERK = va * FR * (Fp / (Kd + Fp)) * ((1 - ERK) / (Ka + 1 - ERK)) - vi * (
            (ERK) / (Ki + ERK))  # vin->vi ; typo?
    return [dG, dN, dFR, dERK]


# returns the ODE as the operator d/dt(), where Fp is a parameter (for bifurcation diagrams)
def dif_Fp_param(t, x):  # not used right now !
    G = x[0]
    N = x[1]
    FR = x[2]
    ERK = x[3]
    dG = (vsg1 * ((ERK ** r) / ((Kag1 ** r) + (ERK ** r))) + vsg2 * ((G ** s) / ((Kag2 ** s) + (G ** s)))) * (
            (Kig ** q) / ((Kig ** q) + (N ** q))) - kdg * G
    dN = (vsn1 * ((Kin1 ** u) / ((Kin1 ** u) + (ERK ** u))) + vsn2 * ((N ** v) / ((Kan ** v) + (N ** v)))) * (
            (Kin2 ** w) / ((Kin2 ** w) + (G ** w))) - kdn * N
    dFR = vsfr1 * ((Kifr) / (Kifr + N)) + vsfr2 * ((G) / (Kafr + G)) - kdfr * FR
    dERK = va * FR * Fp(t) / (Kd + Fp(t)) * ((1 - ERK) / (Ka + 1 - ERK)) - vi * (
            (ERK) / (Ki + ERK))  # vin->vi ; typo?  # is the maximal time
    return [dG, dN, dFR, dERK]


# varying Fp for the bifurcation function.
def Fp(t):  # not used right now !
    return ((FPend - Fp0) / (end - start)) * t + Fp0


# solve with V0 as list (different starting points). Max N and G are hardcoded as 2.2
def varying_v0_solver(step_size,use_varying_Fp=0):
    results_list = []
    step = step_size
    if use_varying_Fp != 1:
        print('Fp constant')
        while step <= 1.6:
            results_list.append(solve_ivp(dif_Fp_const, (start, end), [step, 0.2, FR0, ERK0],
                                      t_eval=timeGrid))  # order for v0: [G0,N0,FR0,ERK0]
            results_list.append(solve_ivp(dif_Fp_const, (start, end), [0.2, step, FR0, ERK0],
                                      t_eval=timeGrid))  # order for v0: [G0,N0,FR0,ERK0]
            step += step_size
    elif use_varying_Fp == 1:
        print('Fp varied for video')
        while step <= 1.6:
            results_list.append(solve_ivp(dif_Fp_var, (start, end), [step, 0.2, FR0, ERK0],
                                      t_eval=timeGrid))  # order for v0: [G0,N0,FR0,ERK0]
            results_list.append(solve_ivp(dif_Fp_var, (start, end), [0.2, step, FR0, ERK0],
                                     t_eval=timeGrid))  # order for v0: [G0,N0,FR0,ERK0]
            step += step_size
    return results_list

def solver_for_Fp():
    value_list = []
    for t_value in timeGrid:
        value_list.append(Fp(t_value))
    return value_list


def four_plots():  # right now not used/isn't useful
    resultsPhaseSpace = varying_v0_solver(0.2)
    resultBifurcation = solve_ivp(dif_Fp_param, (start, end), [0, 0, FR0, ERK0],
                                  t_eval=timeGrid)  # order for v0: [G0,N0,FR0,ERK0,Fp])

    G_data_bif = resultBifurcation.y[0]
    N_data_bif = resultBifurcation.y[1]
    Fp_data_bif = solver_for_Fp()
    t_data_bif = resultBifurcation.t
    # square grid of 4 subplots
    plt.figure(figsize=(9, 6))  # generate slightly larger figure

    # Gata6(Fgf4) plot - bifurcation diagram
    plot1 = plt.subplot(2, 2, 1)
    plot1.plot(Fp_data_bif, G_data_bif, label='upper left')
    plot1.set_xlabel('Fgf4', fontsize=15)
    plot1.set_ylabel('Gata6', fontsize=15)

    # Nanog(Fgf4) plot - bifurcation diagram
    plot2 = plt.subplot(2, 2, 2)
    plot2.plot(Fp_data_bif, N_data_bif, label='upper right')
    plot2.set_xlabel('Fgf4', fontsize=15)
    plot2.set_ylabel('Nanog', fontsize=15)

    # Nanog(Gata6) plot - phase space diagram
    plot3 = plt.subplot(2, 2, 3)
    plot3.set_xlabel('Gata6', fontsize=15)
    plot3.set_ylabel('Nanog', fontsize=15)
    plot3.set_ylim(bottom=0, top=2.2)
    plot3.set_xlim(left=0, right=2.2)
    for result in resultsPhaseSpace:
        G_data_ps = result.y[0]
        N_data_ps = result.y[1]
        plot3.plot(G_data_ps, N_data_ps, label='lower left')

    # Fgf4(t) plot from the bifurcation plots
    plot4 = plt.subplot(2, 2, 4)
    plt.plot(t_data_bif, Fp_data_bif, label='lower right')
    plot4.set_xlabel('Time', fontsize=15)
    plot4.set_ylabel('Fgf4', fontsize=15)
    plt.show()


#make a phase space plot (plots and also saved picture)
def phase_space_only():
    resultsPhaseSpace = varying_v0_solver(0.2, 0)
    ax = plt.axes()
    ax.set_xlabel('Gata6', fontsize=15)
    ax.set_ylabel('Nanog', fontsize=15)
    ax.set_ylim(bottom=0, top=2.2)
    ax.set_xlim(left=0, right=2.2)
    for result in resultsPhaseSpace:
        G_data_ps = result.y[0]
        N_data_ps = result.y[1]
        plt.plot(G_data_ps, N_data_ps, label='lower left')
    plt.savefig('phaseplot{}'.format(timeAndDate), dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='png',
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None, metadata=None)
    plt.show()


#make a dynamic phase space plot with varying Fp (movie gets saved)
def phase_plot_animation():
    fig, ax = plt.subplots()
    line, = ax.plot([], [], lw=2)
    Fp_value_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
    xdata, ydata = [],[]

    def init():
        ax.set_xlabel('Gata6', fontsize=15)
        ax.set_ylabel('Nanog', fontsize=15)
        ax.set_ylim(bottom=0, top=2.2)
        ax.set_xlim(left=0, right=2.2)
        ax.grid()
        line.set_data([], [])
        Fp_value_text.set_text('')
        return line,

    def update(frame):
        global Fp
        Fp = frame
        xdata, ydata = [], []
        for result in varying_v0_solver(0.2, 1):
            xdata.append(result.y[0])
            ydata.append(result.y[1])
            line.set_data(xdata, ydata)
        Fp_value_text.set_text('Fp={0:.2f}'.format(Fp))
        print('current Fp value is', Fp)
        return line,

    animation = FuncAnimation(fig, update, frames=np.linspace(0, 0.11, 128), init_func=init, blit=True)
    animation.save('animation_{}.mp4'.format(timeAndDate), fps=10, extra_args=['-vcodec', 'libx264'])
    print('movie saved')


def main():
    phase_plot_animation()


if __name__ == '__main__':
    main()
