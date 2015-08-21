import sys
import math
import numpy as np
import cvxopt as cvx
from iopro import genfromtxt

p = genfromtxt('cached_kde/cached_kde.fine6.csv')
N, D = p.shape

I = np.random.choice(N, N)
p = p[I]

s = 1e-8

# Truth proportions:
# 0.004597168626601867
# 0.009011011183443553
# 0.41817119150694954
# 0.07374540140202024
# 0.49265798277959727
def F(x=None, z=None):

    if x is None:
        return 0, cvx.matrix(0.2, (D, 1))

    x_arr = np.array(x.trans()).reshape(-1)
    arg = np.dot(p, x_arr)

    if np.min(arg) <= 0.0:
        return None

    f = -np.sum(np.log(arg))
    f = cvx.matrix(f) * s

    arg_m1 = np.power(arg, -1.0)
    Df = [ -np.dot(arg_m1, p[:,c]) for c in range(D) ]
    Df = cvx.matrix(Df, (1, D)) * s

    if z is None:
        return f, Df

    # avoiding numerical difficulties...
    arg_m2 = -2 * np.log(arg)
    H = np.zeros((D, D))
    for b in range(D):
        for c in range(D):
            H[b, c] = np.sum(np.exp(arg_m2 + np.log(p[:,b]) + np.log(p[:,c]))) * z[0]
    H = cvx.matrix(H) * s

    #arg_m2 = np.power(arg, -2.0)
    #H = np.zeros((D, D))
    #for b in range(D):
    #    for c in range(D):
    #        H[b, c] = np.sum(arg_m2 * p[:,b] * p[:,c]) * z[0]
    #H = cvx.matrix(H) * s

    return f, Df, H


if __name__ == '__main__':

    sys.stdout.write('Setting up the problem...\n\n')
    sys.stdout.flush()

    G = cvx.matrix(np.diag(-np.ones(D)))
    h = cvx.matrix(np.zeros(D))
    A = cvx.matrix(np.ones(D), (1, D))
    b = cvx.matrix(1.0)
    #A = cvx.matrix(np.array([ [1.0,1.0,1.0,1.0,1.0],
    #                          [0.0,0.0,1.0,0.0,0.0],
    #                          [0.0,0.0,0.0,1.0,0.0],
    #                          [0.0,0.0,0.0,0.0,1.0] ]))
    #b = cvx.matrix([1.0, 0.41817119150694954, 0.07374540140202024, 0.49265798277959727 ])
    dims = { 'l': D, 'q': [], 's': [] }

    sys.stdout.write('Solving the problem...\n\n')
    sys.stdout.flush()

    sol = cvx.solvers.cp(F, G, h, dims, A, b)

    print
    print 'Solver status: {0}'.format(sol['status'])
    print 'Minimizers:'
    print np.array(sol['x']).reshape(-1)
