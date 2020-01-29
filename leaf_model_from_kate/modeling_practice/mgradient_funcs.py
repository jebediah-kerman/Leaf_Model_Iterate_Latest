
import numpy as np

def M_o(delta):
    m_o = 10

    if delta:
        m_o = 5

    return m_o


def M_i():
    return 0

def dsdt():
    return 0

# to solve using odeint
def dMdt(m, t, N, beta, alpha, dx):
    dmdt = np.zeros(N)

    dmdt[0] = dsdt()

    for i in range(1, N-1):
        dmdt[i] = (beta / dx ** 2) * (m[i+1] - 2*m[i] + m[i-1]) - alpha*m[i] ** 2

    dmdt[-1] = (2 * beta / dx ** 2) * (m[i - 1] - m[i]) - alpha * m[i] ** 2

    return dmdt


def dMdt_not_robust(m, t, N, beta, alpha, dx):
    dmdt = np.zeros(N)

    dmdt[0] = dsdt()

    for i in range(1, N-1):
        dmdt[i] = (beta / dx ** 2) * (m[i+1] - 2*m[i] + m[i-1]) - alpha*m[i]

    dmdt[-1] = (2 * beta / dx ** 2) * (m[i - 1] - m[i]) - alpha * m[i]

    return dmdt


# to solve with own definite integral
def dmdt(ys, t, N, beta, alpha, dx):
    dmdt = np.zeros(N)
    m = ys[t-1,:]

    dmdt[0] = dsdt()

    for i in range(1, N-1):
        dmdt[i] = (beta / dx ** 2) * (m[i+1] - 2*m[i] + m[i-1]) - alpha*m[i] ** 2

    dmdt[-1] = (2 * beta / dx ** 2) * (m[i - 1] - m[i]) - alpha * m[i] ** 2

    return dmdt


def dmdt_not_robust(ys, t, N, beta, alpha, dx):
    dmdt = np.zeros(N)
    m = ys[t - 1, :]


    dmdt[0] = dsdt()

    for i in range(1, N - 1):
        dmdt[i] = (beta / dx ** 2) * (m[i + 1] - 2 * m[i] + m[i - 1]) - alpha * m[i]

    dmdt[-1] = (2 * beta / dx ** 2) * (m[i - 1] - m[i]) - alpha * m[i]

    return dmdt


def defint(y0, dydt, timepoints, dt, args):
    ys = np.zeros((timepoints, len(y0)))
    ys[0] = y0

    for i in range(1, timepoints):
        ys[i] = ys[i - 1] + dydt(ys, i, *args) * dt

    return ys


def get_xi_at_threshold(threshold, ms, time, num_x):
    indices = np.where(ms[time, :] >= threshold)
    # to get the first xi at or above threshold
    xi = indices[0][-1] / float(num_x)

    return xi

# calculate and check steady states
def mst_simple(Mo, x, beta, alpha):
    lam = np.sqrt(beta/alpha)

    return Mo*np.exp(-x/lam)

def mst_power(Mo, x, beta, alpha):
    epsilon = ((alpha * Mo)/ (6*beta))**(-0.5)
    A = 6.0*beta / alpha

    return A*(x+epsilon)**(-2)