import numpy as np

def dsdt():
    return 0




def dIdt(I, N, dI, k, M, alphaI, P, dx):
    dIdt = np.zeros(N)

    # I enters DR at a constant amount
    dIdt[0] = dsdt()
    dIdt[-1] = dsdt()

    for i in range(1, N-1):
        dIdt[i] = (dI/dx**2) * (I[i+1] - 2*I[i] + I[i-1]) - k*I[i] * M[i] - alphaI*P*I[i]
    return dIdt

def dCdt(C, N, dC, k, I, M, alphaC, P, dx):
    dCdt = np.zeros(N)

    ### how to define boundary on both sides properly
    dCdt[0] = (2 * dC / dx ** 2) * (C[1] - C[0]) + k*I[0]*M[0] - alphaC*P[0]*C[0]

    for i in range(1, N - 1):
        dCdt[i] = (dC / dx ** 2) * (C[i + 1] - 2 * C[i] + C[i - 1]) + k*I[i]*M[i] - alphaC*P[i]*C[i]

    dCdt[-1] = (2 * dC / dx ** 2) * (C[-2] - C[-1]) + k*I[-1]*M[-1] - alphaC*P[- 1]*C[- 1]


    return

def dMdt(M, N, dM, k, I, C, alphaC, P, dx):
    dmdt = np.zeros(N)

    dmdt[0] = (2 * dM / dx ** 2) * (M[1] - M[0]) - k*I[0]*M[0] + alphaC*P[0]*C[0]

    for i in range(1, N - 1):
        dmdt[i] = (dM / dx ** 2) * (M[i + 1] - 2 * M[i] + M[i - 1]) - k*I[i]*M[i] + alphaC*P[i]*C[i]

    dmdt[-1] = (2 * dM / dx ** 2) * (M[-2] - M[-1]) - k*I[-1]*M[-1] + alphaC*P[- 1]*C[- 1]

    return

# to solve using odeint


def dYdt(y, t, N, dI, dC, dM, k, alphaI, alphaC, P, dx):

    return [dIdt(y[0], N, dI, k, y[2], alphaI, P, dx),
            dCdt(y[1], N, dC, k, y[0], y[2], alphaC, P, dx),
            dMdt(y[2], N, dM, k, y[0], y[1], alphaC, P, dx)]