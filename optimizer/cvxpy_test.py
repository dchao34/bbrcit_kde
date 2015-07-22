import sys
import math
import numpy as np
import cvxpy as cvx
from iopro import genfromtxt

p = genfromtxt('cached_kde.csv')
N, D = p.shape

#s = 1e-8

if __name__ == '__main__':

    sys.stdout.write('Setting up the problem...\n\n')
    sys.stdout.flush()

    x = cvx.Variable(D)
    s = cvx.Parameter(sign='negative')
    s.value = -1e-8

    obj = cvx.Minimize(s * cvx.sum_entries(cvx.log(p*x)))
    constraints = [ x >= 0, cvx.sum_entries(x) == 1 ]
    prob = cvx.Problem(obj, constraints)

    sys.stdout.write('Solving the problem...\n\n')
    sys.stdout.flush()
    prob.solve(verbose=True)

    print
    print 'Solver status: {0}'.format(prob.status)
    print 'Minimum: {0}'.format(prob.value)
    print 'Minimizers: ', x.value
