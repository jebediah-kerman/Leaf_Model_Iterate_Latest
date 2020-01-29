
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines

import funcs
import mgradient_funcs

# x values to use per book
sim_a = False

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
t = 8
# intervals between seconds
interval = 1000

# from time parameters
dt = 1 / float(interval)
data_points = t*interval

# where to start the conc of Y and Z
Y_start = [0]
Z_start = [0]

# initiate X1, X2 voltages per book

X1 = np.zeros(data_points)
X2 = np.zeros(data_points)

def main():
    # adjust X arrays based on simulation to run
    if sim_a:
        X1[int(.5*interval):1*interval] = 1
        X1[5*interval:] = 1
        X2[int(3.5*interval):4*interval] = 1
    else:
        X1[int(.5*interval):1*interval] = 1
        X1[4*interval:6*interval] = 1
        X2[int(1.5*interval):2*interval] = 1
        X2[7*interval:] = 1

    # initialize y and z arrays
    # ys = np.zeros(data_points)
    # zs = np.zeros(data_points)

    # calculate Y and Z at each time based on X values
    # ys[0] = Y_start
    # zs[0] = Z_start

    ys = mgradient_funcs.defint(Y_start, funcs.dydt, data_points, dt, [X1, X2, w1, w2, Ky, alpha, beta])
    zs = mgradient_funcs.defint(Z_start, funcs.dzdt, data_points, dt, [X1, X2, ys, w1pr, w2pr, w3pr, Kz, alphapr, betapr])

    # for i in range(1,data_points):
        # ys[i] = ys[i-1] + funcs.dydt(X1[i-1], X2[i-1], w1, w2, Ky, alpha, beta, ys[i-1])*dt
        # zs[i] = zs[i-1] + funcs.dzdt(X1[i-1], X2[i-1], ys[i-1], w1pr, w2pr, w3pr, Kz, alphapr, betapr, zs[i-1])*dt

    # even values for t
    ts = np.arange(0,t, dt)

    # make plots of dynamics of the system over time
    fig = plt.figure()
    plt.subplots_adjust(hspace = 0.5)
    ax = fig.add_subplot(4, 1, 1)
    ax.plot(ts, X1, color='tab:blue')
    plt.title('Neuron Voltages C. elegans Multi-FFL')
    plt.ylabel('X1(t)')
    ax = fig.add_subplot(4, 1, 2)
    ax.plot(ts, X2, color='tab:green')
    plt.ylabel('X2(t)')
    ax = fig.add_subplot(4, 1, 3)
    ax.plot(ts, ys, color='tab:orange')
    plt.ylabel('Y(t)')
    ax.axhline(0.45, color='k', linestyle='--')
    ax = fig.add_subplot(4, 1, 4)
    ax.plot(ts, zs, color='tab:pink')
    plt.ylabel('Z(t)')
    plt.xlabel('time (s)')
    plt.show()
    return None


if __name__ == '__main__':
    main()