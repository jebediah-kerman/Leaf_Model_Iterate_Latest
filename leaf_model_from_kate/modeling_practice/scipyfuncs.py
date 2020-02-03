# compare weighted X1, X2, Y to Ky to determine if threshold is reached
def z_on(X1, X2, Y_o, w1pr, w2pr, w3pr, Kz):
    return (w1pr*X1 + w2pr*X2 + w3pr*Y_o) > Kz

# compare weighted X1, X2 to Ky to determine if threshold is reached
def y_on(X1, X2, w1, w2, Ky):
    return (X1*w1 + X2*w2) > Ky

def dydt_scipy(y, t, X1_pulses, X2_pulses, w1, w2, Ky, alpha, beta):
    X1 = X_pulsing(t, X1_pulses)
    X2 = X_pulsing(t, X2_pulses)
    dydt = beta * y_on(X1, X2, w1, w2, Ky) - alpha*y
    return dydt

def dzdt_scipy(z, t, X1_pulses, X2_pulses, Y, w1pr, w2pr, w3pr, Kz, alphapr, betapr):
    X1 = X_pulsing(t, X1_pulses)
    X2 = X_pulsing(t, X2_pulses)
    dzdt = betapr* z_on(X1, X2, Y, w1pr, w2pr, w3pr, Kz) - alphapr*z
    return dzdt

# functions for x based on timing of pulses
# X_pulses is a list of lists for ranges of pulses
def X_pulsing(t, X_pulses):
    x_t = 0.0

    if any(start <= t <= end for start, end in X_pulses) :
        x_t = 1.0

    return x_t


def dYdt(y, t, X1_pulses, X2_pulses, w1, w2, Ky, alpha, beta, w1pr, w2pr, w3pr, Kz, alphapr, betapr):

    dYdt = [dydt_scipy(y[0], t, X1_pulses, X2_pulses, w1, w2, Ky, alpha, beta), dzdt_scipy(y[1], t, X1_pulses, X2_pulses, y[0], w1pr, w2pr, w3pr, Kz, alphapr, betapr)]


    return dYdt

#X1, X2, w1, w2, Ky, alpha, beta, w1pr, w2pr, w3pr, Kz, alphapr, betapr)

# def f(t, y, c):
#     dydt = [c[0]*np.cos(c[1]*t), c[2]*y[0]+c[3]*t]
#         return dydt