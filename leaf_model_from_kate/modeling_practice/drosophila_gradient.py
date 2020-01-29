import numpy as np
from scipy.integrate import odeint, ode
import matplotlib.pyplot as plt

import dgradientfuncs as f

# Simulation time parameters
# number of iterations to run
num_t = 100
# total time to run arbitrary units
T = 100000
# time values
dt = 1

# discrete xi from (0,l)
num_xi = 100
N = num_xi + 1
# length of region defined by sum of xi
L = 1
x = np.linspace(0, L, N)
dx = x[1]

# parameters
alphaC = 1
alphaI = .001

k = 1
dI = .5
dC = 1
dM = .001

# starting concentrations across synctium
I_start = 0.5
M_start = 0.5
C_start = 0

I = []
M = []
C = []

# define boundaries
ne_width = 10
# i is only expressed in the NE
I.append([I_start]*ne_width +[0] * (N - 2*ne_width) + [I_start] * ne_width)

# complexes to begin
C.append([C_start] * N)

# m is uniformly expressed
M.append([M_start] * N)
P = np.zeros(N)

# P is only expressed in the DR
P[int(ne_width+1):int(N-ne_width)] = .5

def main():
    # solve using the odespy solver
    solver = ode(f.dYdt).set_integrator('vode')
    solver.set_f_params(N, dI, dC, dM, k, alphaI, alphaC, P, dx)
    solver.set_initial_value([I[0], C[0], M[0]], T)

    t = []

    while solver.successful() and solver.t < T:
        solver.integrate(solver.t+dt)
        I.append(solver.y[0])
        C.append(solver.y[1])
        M.append(solver.y[2])
        t.append(solver.t)


    I_array = np.array(I)
    C_array = np.array(C)
    M_array = np.array(M)

    # test robustness of two models of morphogen gradients to canges (delta) in starting morphogen concentration
    # sol = odeint(f.dYdt, [startI, startC, startM], t, args=(N, dI, dC, dM, k, alphaI, alphaC, P, dx))
    #

    # ## compare values at steady state
    # print('The steady state value for xi=L at t=T is '+str(m_robust[num_t, 50]) +' per the numerical solution and '+
    #       str(mgradient_funcs.mst_power(m_concs[0], .5, beta, alpha)) + ' for the analytical solution for '
    #                                                                         'the robust method')
    #
    # print('For a threshold of  '+str(threshold) + ' the difference between Mo and Mo\' for the robust method is '+
    #       str(m_robust_at_threshold - m_delta_at_threshold))
    #
    # # note x and t are saved in m the reverse of what would make sense ... m(t,x)
    #
    fig = plt.figure()
    plt.subplots_adjust(hspace=0.5)
    plt.subplots_adjust(wspace=0.5)

    ax = fig.add_subplot(3, 1, 1)
    # at steady state (last sim) plot I, M across all x
    plt.title('Concentration of Protein M and it\'s Inhibitor I')
    ax.plot(x, I_array[num_t,:], color='palevioletred')
    ax.plot(x, M_array[num_t, :], color='cornflowerblue')
    ax.set_ylim(-0.05, 1.05)
    plt.ylabel('[Protein]')
    plt.xlabel('Position')

    ax = fig.add_subplot(3, 1, 2)
    # check that the midpoint concentrations are at steady state for I, C, M, P (should always be)
    plt.title('Midline is at Steady State')
    ax.plot(t, I_array[:, int(N/2)], color='palevioletred')
    ax.plot(t, C_array[:, int(N / 2)], color='mediumslateblue')
    ax.plot(t, M_array[:, int(N / 2)], color='cornflowerblue')
    ax.plot(t, P[int(N / 2)], color='mediumslateblue', linestyle='--')
    ax.set_ylim(-0.05, 1.05)
    plt.ylabel('[Protein]')
    plt.xlabel('Time')

    ax = fig.add_subplot(3, 1, 3)
    # check that the end concentrations are at steady state for I, C, M, P (should always be)
    plt.title('Left Edge is at Steady State')
    ax.plot(t, I_array[:, 0], color='palevioletred')
    ax.plot(t, C_array[:, 0], color='mediumslateblue')
    ax.plot(t, M_array[:, 0], color='cornflowerblue')
    ax.plot(t, P[0], color='mediumslateblue', linestyle='--')
    ax.set_ylim(-0.05, 1.05)
    plt.ylabel('[Protein]')
    plt.xlabel('Time')

    plt.show()

    return None

if __name__ == '__main__':
    main()