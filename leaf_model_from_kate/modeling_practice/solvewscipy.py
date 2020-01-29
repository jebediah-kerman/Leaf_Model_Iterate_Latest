
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp

import scipyfuncs


# parameters per book
alpha = 1
alphapr = 1
w1 = w2 = w1pr = w2pr = 0.5
w3pr = 0.4
Ky = 0.01
Kz = 0.68
beta = 1
betapr = 1

# how long to run simulation, seconds
runtime = 8
# intervals between seconds
interval = 100
num_time_points = runtime*interval

# array of time points
t = np.linspace(0, runtime, num_time_points)

# where to start the conc of Y and Z
Y_start = 0
Z_start = 0

# initiate X1, X2 empty ranges

X1 = []
X2 = []

def main():
    # adjust X arrays based on simulation to run

    # sim A
    # X1.append([.5, 1])
    # X1.append([5, 8])
    # X2.append([3.5, 4])

    #sim b
    X1.append([.5,1])
    X1.append([4,6])
    X2.append([1.5,2])
    X2.append([7,8])

    # actually calculate X1 X2 for graphing
    X1_graphing = [scipyfuncs.X_pulsing(time, X1) for time in t]
    X2_graphing = [scipyfuncs.X_pulsing(time, X2) for time in t]

    sol = odeint(scipyfuncs.dYdt, [Y_start, Z_start], t, args=(X1, X2, w1, w2, Ky, alpha, beta, w1pr, w2pr, w3pr, Kz, alphapr, betapr))

    # make plots of dynamics of the system over time
    fig = plt.figure()
    plt.style.use(['dark_background', 'seaborn-talk'])
    fig.set_size_inches(8, 6)
    plt.rcParams['axes.facecolor'] = 'k'
    plt.rcParams['savefig.facecolor'] = 'k'
    plt.subplots_adjust(hspace = 0.5)
    ax = fig.add_subplot(4, 1, 1)
    ax.plot(t, X1_graphing, color='tab:blue')
    plt.title('Neuron Voltages C. elegans Multi-FFL')
    plt.ylabel('X1(t)')
    ax.set_ylim(-.1, 1.1)
    ax = fig.add_subplot(4, 1, 2)
    ax.plot(t, X2_graphing, color='tab:green')
    plt.ylabel('X2(t)')
    ax.set_ylim(-.1, 1.1)
    ax = fig.add_subplot(4, 1, 3)
    ax.plot(t, sol[:, 0], color='tab:orange')
    plt.ylabel('Y(t)')
    ax.axhline(0.45, color='w', linestyle='--')
    ax.set_ylim(-.1, 1.1)
    ax = fig.add_subplot(4, 1, 4)
    ax.plot(t, sol[:, 1], color='tab:pink')
    ax.set_ylim(-.1, 1.1)
    plt.ylabel('Z(t)')
    plt.xlabel('time (s)')
    plt.show()
    return None


if __name__ == '__main__':
    main()