import sys
import math
import numpy as np
import cvxopt as cvx
from iopro import genfromtxt

#p = genfromtxt('cached_kde.csv')
p = genfromtxt('cached_kde.course.csv')
N, D = p.shape

s = 1e-8

def F(x=None, z=None):

    if x is None:
        # return 0, cvx.matrix([0.00458365, 0.00888265, 0.41721128, 0.07378726, 0.49553516], (D,1))
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

    arg_m2 = np.power(arg, -2.0)
    H = np.zeros((D, D))
    for b in range(D):
        for c in range(D):
            H[b, c] = np.sum(arg_m2 * p[:,b] * p[:,c]) * z[0]
    H = cvx.matrix(H) * s

    return f, Df, H


if __name__ == '__main__':

    sys.stdout.write('Setting up the problem...\n\n')
    sys.stdout.flush()

    G = cvx.matrix(np.diag(-np.ones(D)))
    h = cvx.matrix(np.zeros(D))
    A = cvx.matrix(np.array([ [1.0,1.0,1.0,1.0,1.0],
                             [0.0,0.0,0.0,0.0,1.0] ]))
    b = cvx.matrix([1.0, 0.49553516])
    #A = cvx.matrix(np.array([ [1.0,1.0,1.0,1.0,1.0],
    #                          [0.0,0.0,0.0,1.0,0.0],
    #                          [0.0,0.0,0.0,0.0,1.0] ]))
    #b = cvx.matrix([1.0, 0.07378726, 0.49553516])
    #A = cvx.matrix(np.array([ [1.0,1.0,1.0,1.0,1.0],
    #                          [0.0,0.0,1.0,0.0,0.0],
    #                          [0.0,0.0,0.0,1.0,0.0],
    #                          [0.0,0.0,0.0,0.0,1.0] ]))
    #b = cvx.matrix([1.0, 0.41721128, 0.07378726, 0.49553516])
    #A = cvx.matrix(np.array([ [1.0,1.0,1.0,1.0,1.0],
    #                          [0.0,0.0,1.0,0.0,0.0],
    #                          [0.0,0.0,0.0,1.0,0.0] ]))
    #b = cvx.matrix([1.0, 0.41721128, 0.07378726])
    #A = cvx.matrix(np.ones(D), (1, D))
    #b = cvx.matrix(1.0)
    dims = { 'l': D, 'q': [], 's': [] }

    sys.stdout.write('Solving the problem...\n\n')
    sys.stdout.flush()

    sol = cvx.solvers.cp(F, G, h, dims, A, b)

    print
    print 'Solver status: {0}'.format(sol['status'])
    print 'Minimizers:'
    print sol['x']
