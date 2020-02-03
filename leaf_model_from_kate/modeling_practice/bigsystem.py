
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

import bigfuncs
import mgradient_funcs

# parameters to change transcription promotion per formulas below

# dY1 = betas[0] * (ws[0] * X1 > Ks[0]) - alphas[0] * Y1
# dX2 = betas[1] * (ws[1] * X1 + ws[2] * Y2 > Ks[1]) - alphas[1] * X2
# dY2 = betas[2] * (ws[3] * X2 > Ks[2]) - alphas[2] * Y2
#
# dZ1 = betas[3] * (ws[4] * X1 > Ks[3]) - betas[4] * (ws[5] * Y1 > Ks[4]) - alphas[3] * Z1
# dZ2 = betas[5] * (ws[6] * X2 > Ks[5]) - betas[6] * (ws[7] * Y2 > Ks[6]) - alphas[4] * Z2
# dZ3 = betas[7] * (ws[8] * X2 + ws[9] * Y2 > Ks[7] ) - alphas[5] * Z3

#         0    1   2   3    4   5    6   7   8   9
alphas = [1,   1,  1,  1,   1,  1]
betas =  [1,   1,  1,  1,   1,  1,   1,  1]
ws =     [.4, .4, .4, .4,  .4, .4,  .4, .4, .4, .4]
Ks =     [.2, .7, .2, .1, .3, .1, .3, .7]



# how long to run simulation, seconds
runtime = 8
# intervals between seconds
interval = 1000
num_time_points = runtime*interval

# array of time points
t = np.linspace(0, runtime, num_time_points)
dt = t[1]

# where to start the conc of Y1, X2, Y2, Z1-3
start_concs = [0, 0, 0, 0, 0, 0]

# when X stimulus, S_x starts (around 1)
# to get the same values every time, set the seed
np.random.seed(13)
S_x_start_time = np.random.random()

# initiate X1 concentration
# conc X saturates
Xst = 1
# rate dependency
alphx = 1

def main():

    # solve diffeq for FFL actors

    # Y1 = mgradient_funcs.defint(start_concs[0], bigfuncs.simpleActivation, num_time_points, dt, args=(betas[0], ws[0], X1, Ks[0], alphas[0]))
    # X2 = mgradient_funcs.defint(start_concs[1], bigfuncs.coherentFFL, num_time_points, dt, args=(betas[1], X1, Y1, ws[1], ws[2], Ks[1], alphas[1]))
    # Y2 = mgradient_funcs.defint(start_concs[2], bigfuncs.simpleActivation, num_time_points, dt,  args=(betas[2], ws[3], X2, Ks[2], alphas[2]))
    #
    # Z1 = mgradient_funcs.defint(start_concs[3], bigfuncs.incoherentFFL, num_time_points, dt, args=(betas[3], betas[4], X1, Y1, ws[4], ws[5], Ks[3], Ks[4], alphas[3]))
    # Z2 = mgradient_funcs.defint(start_concs[4], bigfuncs.incoherentFFL, num_time_points, dt, args=(betas[5], betas[6], X2, Y2, ws[6], ws[7], Ks[5], Ks[6], alphas[4]))
    # Z3 = mgradient_funcs.defint(start_concs[5], bigfuncs.coherentFFL, num_time_points, dt, args=(betas[7], X2, Y2, ws[8], ws[9], Ks[7], alphas[5]))
    # #

    # adapt to use odeint
    sol = odeint(bigfuncs.dYdt, start_concs, t, args=(betas, ws, Ks, alphas, Xst, S_x_start_time, alphx))

    #
    graphing_X1s = [bigfuncs.X1(Xst, S_x_start_time, alphx, time) for time in t]
    # make plots of dynamics of the system over time
    fig = plt.figure()
    plt.subplots_adjust(hspace = 0.5)
    plt.style.use(['dark_background', 'seaborn-talk'])
    fig.set_size_inches(8, 6)
    plt.rcParams['axes.facecolor'] = 'k'
    plt.rcParams['savefig.facecolor'] = 'k'
    plt.title('Protein Concentration of B. subtilis Interlocking FFLs')
    # fix weird buggy axes
    plt.axis('off')

    # plots x's
    ax = fig.add_subplot(3, 1, 1)
    ax.plot(t, graphing_X1s, color='darkseagreen')
    ax.plot(t, sol[:, 1], color='mediumorchid')
    plt.ylim(-.1, 1.1)
    plt.ylabel('X(t)')

    # plot y's
    ax = fig.add_subplot(3, 1, 2)
    ax.plot(t, sol[:, 0], color='lightsteelblue', linestyle='--')
    ax.plot(t, sol[:, 2], color='plum', linestyle='--')
    ax.set_ylim(-.1, 1.1)
    plt.ylabel('Y(t)')

    # plot z's
    ax = fig.add_subplot(3, 1, 3)
    ax.plot(t, sol[:, 3], color='darkseagreen')
    ax.plot(t, sol[:, 4], color='mediumorchid')
    ax.plot(t, sol[:, 5], color='darkmagenta')
    ax.set_ylim(-.1, 1.1)
    plt.ylabel('Z(t)')

    plt.xlabel('time (s)')

    plt.show()
    return None


if __name__ == '__main__':
    main()