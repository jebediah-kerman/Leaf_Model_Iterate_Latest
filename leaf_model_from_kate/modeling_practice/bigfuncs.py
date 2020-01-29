import numpy as np

def coherentFFL(Y, t, b, X1, X2, w1, w2, K, a):
    return b*((w1*X1 + w2*X2) > K) - a*Y

def incoherentFFL(Y, t, b1, b2, X1, X2, w1, w2, K1, K2, a):
    return b1*((w1*X1) > K1) - b2*((w2*X2) > K2) - a*Y

def simpleActivation(Y, t, b, w, X, K, a):
    return b*((w*X) > K) - a*Y

def X1(Xst, S_x_start_time, alphx, t):

    X1 = 0
    # assume saturating expression after turn on of S_X
    if t > S_x_start_time:
        X1 = Xst * (1 - np.exp(-(t - S_x_start_time) * alphx))

    return X1


def dYdt(y, t, betas, ws, Ks, alphas, Xst, S_x_start_time, alphx):

            # dY1 .... y[0]
    return [simpleActivation(y[0], t, betas[0], ws[0], X1(Xst, S_x_start_time, alphx, t), Ks[0], alphas[0]),
            # dX2 .... y[1]
            coherentFFL(y[1], t, betas[1], X1(Xst, S_x_start_time, alphx, t), y[0], ws[1], ws[2], Ks[1], alphas[1]),
            # dY2 .... y[2]
            simpleActivation(y[2], t, betas[2], ws[3], y[1], Ks[2], alphas[2]),
            # dZ1 ..... y[3[
            incoherentFFL(y[3], t, betas[3], betas[4], X1(Xst, S_x_start_time, alphx, t), y[0], ws[4], ws[5], Ks[3], Ks[4], alphas[3]),
            #dZ2 ..... y[4]
            incoherentFFL(y[4], t, betas[5], betas[6], y[1], y[2], ws[6], ws[7], Ks[5], Ks[6], alphas[4]),
            # dZ3 .... y[5]
            coherentFFL(y[5], t, betas[7], y[1], y[2], ws[8], ws[9], Ks[7], alphas[5])
            ]