


# C. elegans neuron example from An Intro to Sytems Biology
# rate of Y response from X1 or X2 > Ky
def dydt( ys, t, X1, X2, w1, w2, Ky, alpha, beta):
    Y_t = beta * y_on(X1[t-1], X2[t-1], w1, w2, Ky) - alpha*ys[t-1]
    return Y_t

# compare weighted X1, X2 to Ky to determine if threshold is reached
def y_on(X1, X2, w1, w2, Ky):
    return (X1*w1 + X2*w2) > Ky

# rate of Y response from X1 or X2 and Y > Kz
def dzdt(zs, t, X1, X2, Y, w1pr, w2pr, w3pr, Kz, alphapr, betapr):
    Z_t = betapr* z_on(X1[t-1], X2[t-1], Y[t-1], w1pr, w2pr, w3pr, Kz) - alphapr*zs[t-1]
    return Z_t

# compare weighted X1, X2, Y to Ky to determine if threshold is reached
def z_on(X1, X2, Y_o, w1pr, w2pr, w3pr, Kz):
    return (w1pr*X1 + w2pr*X2 + w3pr*Y_o) > Kz

