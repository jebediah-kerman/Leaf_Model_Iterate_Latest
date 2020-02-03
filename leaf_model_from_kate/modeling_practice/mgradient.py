import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

import mgradient_funcs

# Simulation time parameters
# number of iterations to run
num_t = 100
# total time to run arbitrary units
T = 100
# time values
t = np.linspace(0, T, num_t+1)
dt = t[1]

# discrete xi from (0,l)
num_xi = 100
N = num_xi + 1
# length of region defined by sum of xi
L = 1
x = np.linspace(0, L, N)

# other parameters for dmdt
beta = .01
alpha = 1
dx = x[1]

# initial conditions m(x,0)
m_concs = np.zeros(N)
m_concs[0] = mgradient_funcs.M_o(False)
m_concs[1:] = mgradient_funcs.M_i()

# test robustness with a change to starting concentration
m_concs_delta = m_concs.copy()
m_concs_delta[0] = mgradient_funcs.M_o(True)

def main():
    # iterate over each discrete xi and solve ode for time t

    # test robustness of two models of morphogen gradients to canges (delta) in starting morphogen concentration
    m_robust = odeint(mgradient_funcs.dMdt, m_concs, t, args=(N, beta, alpha, dx))
    m_robust_delta = odeint(mgradient_funcs.dMdt, m_concs_delta, t, args=(N, beta, alpha, dx))

    m_not_robust = odeint(mgradient_funcs.dMdt_not_robust, m_concs, t, args=(N, beta, alpha, dx))
    m_not_robust_delta = odeint(mgradient_funcs.dMdt_not_robust, m_concs_delta, t, args=(N, beta, alpha, dx))


    threshold = 0.01
    m_robust_at_threshold = mgradient_funcs.get_xi_at_threshold(threshold, m_robust, num_t, num_xi)
    m_delta_at_threshold = mgradient_funcs.get_xi_at_threshold(threshold, m_robust_delta, num_t, num_xi)
    m_not_at_threshold = mgradient_funcs.get_xi_at_threshold(threshold, m_not_robust, num_t, num_xi)
    mn_delta_at_threshold = mgradient_funcs.get_xi_at_threshold(threshold, m_not_robust_delta, num_t, num_xi)

    ## compare values at steady state
    print('The steady state value for xi=L at t=T is '+str(m_robust[num_t, 50]) +' per the numerical solution and '+
          str(mgradient_funcs.mst_power(m_concs[0], .5, beta, alpha)) + ' for the analytical solution for '
                                                                            'the robust method')
    print('The steady state value for xi=L at t=T is '+str(m_not_robust[num_t, 50]) +' per the numerical solution and '+
          str(mgradient_funcs.mst_simple(m_concs[0], .5, beta, alpha)) + ' for the analytical solution for '
                                                                            'the robust method')

    print('For a threshold of  '+str(threshold) + ' the difference between Mo and Mo\' for the robust method is '+
          str(m_robust_at_threshold - m_delta_at_threshold))
    print('For a threshold of  ' + str(threshold) + ' the difference between Mo and Mo\' for the not robust method is ' +
          str(m_not_at_threshold - mn_delta_at_threshold))
    # using kate's basic euler method
    # m_robust = mgradient_funcs.defint(m_concs, mgradient_funcs.dmdt, num_t+1, dt, [N, beta, alpha, dx])
    # m_robust_delta = mgradient_funcs.defint(m_concs_delta, mgradient_funcs.dmdt, num_t+1, dt, [N, beta, alpha, dx])
    #
    # m_not_robust = mgradient_funcs.defint(m_concs, mgradient_funcs.dmdt_not_robust, num_t+1, dt, [N, beta, alpha, dx])
    # m_not_robust_delta = mgradient_funcs.defint(m_concs_delta, mgradient_funcs.dmdt_not_robust, num_t + 1, dt, [N, beta, alpha, dx])

    # make plots egads
    # note x and t are saved in m the reverse of what would make sense ... m(t,x)

    # plots of m_robust over time and space
    fig = plt.figure()
    plt.style.use(['dark_background', 'seaborn-talk'])
    fig.set_size_inches(8, 6)
    plt.rcParams['axes.facecolor'] = 'k'
    plt.rcParams['savefig.facecolor'] = 'k'
    plt.subplots_adjust(hspace=0.5)
    plt.subplots_adjust(wspace=0.5)

    # robust plots
    fig.suptitle('[M]o in Robust Morphogen Gradient')


    ax = fig.add_subplot(3, 2, 1)
    plt.title('[M] at discrete t=')
    plt.ylabel('M(x,0)')
    # ax.plot(x, m_robust[0, :], color='tab:orange')
    plt.scatter(x, m_robust[0, :], c=cm.RdPu(m_robust[0, :]), edgecolor='none')
    ax.set_ylim(-0.05, 1.05)

    ax = fig.add_subplot(3, 2, 2)
    ax.plot(t, m_robust[:, 0], color='rebeccapurple')
    plt.title('[M] at discrete xi=')
    plt.ylabel('M(0,t)')
    ax.set_ylim(-0.05, 1.05)

    ax = fig.add_subplot(3, 2, 3)
    plt.ylabel('M(x,T/2)')
    # ax.plot(x, m_robust[int(num_t / 2), :], color='tab:orange')
    plt.scatter(x, m_not_robust[int(num_t / 2), :], c=cm.RdPu(m_not_robust[int(num_t / 2), :]), edgecolor='none')
    ax.set_ylim(-0.05, 1.05)

    ax = fig.add_subplot(3, 2, 4)
    ax.plot(t, m_robust[:,int(num_xi/2)], color='pink')
    plt.ylabel('M(L/2,t)')
    ax.set_ylim(-0.05, 1.05)

    ax = fig.add_subplot(3, 2, 5)
    plt.ylabel('M(x,T)')
    # ax.plot(x, m_robust[num_t, :], color='tab:orange')
    plt.scatter(x, m_robust[num_t, :], c=cm.RdPu(m_robust[num_t, :]), edgecolor='none')
    plt.xlabel('Location in Tissue (x)')
    ax.set_ylim(-0.05, 1.05)

    ax = fig.add_subplot(3, 2, 6)
    ax.plot(t, m_robust[:,num_xi], color='lavenderblush')
    plt.ylabel('M(L,t)')
    plt.xlabel('Time (t)')
    ax.set_ylim(-0.05, 1.05)

    plt.show()

    # make the same plots for delta
    # fig = plt.figure()
    # plt.style.use(['dark_background', 'seaborn-talk'])
    # fig.set_size_inches(8, 6)
    # plt.rcParams['axes.facecolor'] = 'k'
    # plt.rcParams['savefig.facecolor'] = 'k'
    # plt.subplots_adjust(hspace=0.5)
    # plt.subplots_adjust(wspace=0.5)
    #
    # fig.suptitle('[M] Delta in Robust Morphogen Gradient')
    # ax = fig.add_subplot(3, 2, 1)
    # plt.title('[M] at discrete t=')
    # plt.ylabel('M(x,0)')
    # ax.plot(x, m_robust_delta[0, :], color='tab:orange')
    # ax.set_ylim(-0.05, 1.05)
    #
    # ax = fig.add_subplot(3, 2, 2)
    # ax.plot(t, m_robust_delta[:, 0], color='tab:blue')
    # plt.title('[M] at discrete xi=')
    # plt.ylabel('M(0,t)')
    # ax.set_ylim(-0.05, 1.05)
    #
    # ax = fig.add_subplot(3, 2, 3)
    # plt.ylabel('M(x,T/2)')
    # ax.plot(x, m_robust_delta[int(num_t / 2), :], color='tab:orange')
    # ax.set_ylim(-0.05, 1.05)
    #
    # ax = fig.add_subplot(3, 2, 4)
    # ax.plot(t, m_robust_delta[:, int(num_xi / 2)], color='tab:blue')
    # plt.ylabel('M(L/2,t)')
    # ax.set_ylim(-0.05, 1.05)
    #
    # ax = fig.add_subplot(3, 2, 5)
    # plt.ylabel('M(x,T)')
    # ax.plot(x, m_robust_delta[num_t, :], color='tab:orange')
    # plt.xlabel('Location in Tissue (x)')
    # ax.set_ylim(-0.05, 1.05)
    #
    # ax = fig.add_subplot(3, 2, 6)
    # ax.plot(t, m_robust_delta[:, num_xi], color='tab:blue')
    # plt.ylabel('M(L,t)')
    # plt.xlabel('Time (t)')
    # ax.set_ylim(-0.05, 1.05)
    #
    #
    #
    # plt.show()

    # not robust plots
    fig = plt.figure()
    plt.style.use(['dark_background', 'seaborn-talk'])
    fig.set_size_inches(8, 6)
    plt.rcParams['axes.facecolor'] = 'k'
    plt.rcParams['savefig.facecolor'] = 'k'
    plt.subplots_adjust(hspace=0.5)
    plt.subplots_adjust(wspace=0.5)

    # make the same plots for the not robust solution
    ax = fig.add_subplot(3, 2, 1)
    fig.suptitle('[M]o in Simple Morphogen Gradient')
    plt.title('[M] at discrete t=')
    plt.ylabel('M(x,0)')
    plt.scatter(x, m_not_robust[0, :], c=cm.RdPu(m_not_robust[0, :]), edgecolor='none')
    # ax.plot(x, m_not_robust[0, :], color='tab:orange')
    ax.set_ylim(-0.05, 1.05)

    ax = fig.add_subplot(3, 2, 2)
    ax.plot(t, m_not_robust[:, 0], color='rebeccapurple')
    plt.title('[M] at discrete xi=')
    plt.ylabel('M(0,t)')
    ax.set_ylim(-0.05, 1.05)

    ax = fig.add_subplot(3, 2, 3)
    plt.ylabel('M(x,T/2)')
    # ax.plot(x, m_not_robust[int(num_t / 2), :], color='tab:orange')
    plt.scatter(x, m_not_robust[int(num_t / 2), :], c=cm.RdPu( m_not_robust[int(num_t / 2), :]), edgecolor='none')
    ax.set_ylim(-0.05, 1.05)

    ax = fig.add_subplot(3, 2, 4)
    ax.plot(t, m_not_robust[:, int(num_xi / 2)], color='pink')
    plt.ylabel('M(L/2,t)')
    ax.set_ylim(-0.05, 1.05)

    ax = fig.add_subplot(3, 2, 5)
    plt.ylabel('M(x,T)')
    # ax.plot(x, m_not_robust[num_t, :], color='tab:orange')
    plt.scatter(x, m_not_robust[num_t, :], c=cm.RdPu(m_not_robust[num_t, :]), edgecolor='none')
    plt.xlabel('Location in Tissue (x)')
    ax.set_ylim(-0.05, 1.05)

    ax = fig.add_subplot(3, 2, 6)
    ax.plot(t, m_not_robust[:, num_xi], color='lavenderblush')
    plt.ylabel('M(L,t)')
    plt.xlabel('Time (t)')
    ax.set_ylim(-0.05, 1.05)

    plt.show()

    # plots comparing m_robust and m_not_robust gradients at final time point
    fig = plt.figure()
    plt.style.use(['dark_background', 'seaborn-talk'])
    fig.set_size_inches(8, 6)
    plt.rcParams['axes.facecolor'] = 'k'
    plt.rcParams['savefig.facecolor'] = 'k'
    plt.subplots_adjust(hspace=0.5)
    plt.subplots_adjust(wspace=0.5)

    ax = fig.add_subplot(2, 1, 2)
    plt.title('Robust Gradient')
    plt.ylabel('M(x,T)')

    ax.set_yscale('log')
    ax.set_yticks([.01, .1, 1, 10])
    ax.set_ylim(.001, 10)

    ax.plot(x, m_robust[num_t, :], color='burlywood')
    ax.plot(x, m_robust_delta[num_t, :], color='papayawhip', linestyle='--')
    plt.xlabel('Location in Tissue (x)')

    ax = fig.add_subplot(2, 1, 1)
    plt.title('Not Robust Gradient')
    plt.ylabel('M(x,T)')
    ax.set_yscale('log')
    ax.set_yticks([.01, .1, 1, 10])
    ax.set_ylim(.001, 10)

    ax.plot(x, m_not_robust[num_t, :], color='indianred')
    ax.plot(x, m_not_robust_delta[num_t, :], color='lightcoral', linestyle='--')
    plt.xlabel('Location in Tissue (x)')


    plt.show()

    # make 3D figure
    # robust plot
    fig = plt.figure()
    plt.style.use(['dark_background', 'seaborn-talk'])
    fig.set_size_inches(8, 6)
    plt.rcParams['axes.facecolor'] = 'k'
    plt.rcParams['savefig.facecolor'] = 'k'
    ax = fig.add_subplot(111, projection='3d')

    SX, ST = np.meshgrid(x, t)
    ax.plot_surface(SX, ST, m_robust, cmap='RdPu')
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('M(x,t)')
    ax.set_title('[M] in a Robust Morphogen Gradient')

    plt.show()

    # not robust plot
    fig = plt.figure()
    plt.style.use(['dark_background', 'seaborn-talk'])
    fig.set_size_inches(8, 6)
    plt.rcParams['axes.facecolor'] = 'k'
    plt.rcParams['savefig.facecolor'] = 'k'
    ax = fig.add_subplot(111, projection='3d')

    SX, ST = np.meshgrid(x, t)
    ax.plot_surface(SX, ST, m_not_robust, cmap='RdPu')
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('M(x,t)')
    ax.set_title('[M] in Simple Morphogen Gradient')

    plt.show()

    return None

if __name__ == '__main__':
    main()