import sys
import math
import numpy as np
import scipy.optimize as optimize
from iopro import genfromtxt

p = genfromtxt('cached_kde.csv')
N, D = p.shape

I = np.random.choice(N, N)
p = p[I]

s = 1e-7

def f(x):
    return -np.sum(np.log(np.dot(p,x))) * s

def Df(x):
    arg = np.power(np.dot(p,x), -1.0)
    return np.array([ -np.dot(arg, p[:,c]) for c in range(D) ]) * s

def Hf(x):
    arg = np.power(np.dot(p,x), -2.0)
    H = np.zeros((D, D))
    for b in range(D):
        for c in range(D):
            H[b, c] = np.sum(arg * p[:,b] * p[:,c])

    return H * s


if __name__ == '__main__':

    sys.stdout.write('Setting up the problem...\n\n')
    sys.stdout.flush()

    #x0 = np.random.rand(D)
    #x0 = x0 / np.sum(x0)
    #print x0
    x0 = np.array([ 0.00458365, 0.00888265, 0.41721128, 0.07378726, 0.49553516 ])
    print x0
    cons = ({ 'type': 'eq',
              'fun' : lambda x : np.array(np.sum(x) - 1.0),
              'jac' : lambda x : np.array([1.0, 1.0, 1.0, 1.0, 1.0]) },
            #{ 'type': 'eq',
            #  'fun' : lambda x : np.array(x[4]-0.49553516),
            #  'jac' : lambda x : np.array([0.0, 0.0, 0.0, 0.0, 1.0]) },
            #{ 'type': 'eq',
            #  'fun' : lambda x : np.array(x[3]-0.07378726),
            #  'jac' : lambda x : np.array([0.0, 0.0, 0.0, 1.0, 0.0]) },
            { 'type': 'ineq',
              'fun' : lambda x : np.array(x[0]),
              'jac' : lambda x : np.array([1.0, 0.0, 0.0, 0.0, 0.0]) },
            { 'type': 'ineq',
              'fun' : lambda x : np.array(x[1]),
              'jac' : lambda x : np.array([0.0, 1.0, 0.0, 0.0, 0.0]) },
            { 'type': 'ineq',
              'fun' : lambda x : np.array(x[2]),
              'jac' : lambda x : np.array([0.0, 0.0, 1.0, 0.0, 0.0]) },
            { 'type': 'ineq',
              'fun' : lambda x : np.array(x[3]),
              'jac' : lambda x : np.array([0.0, 0.0, 0.0, 1.0, 0.0]) },
            { 'type': 'ineq',
              'fun' : lambda x : np.array(x[4]),
              'jac' : lambda x : np.array([0.0, 0.0, 0.0, 0.0, 1.0]) },
           )
    bounds = [ (0.0, 1.0) ] * 5
    opts = {'disp': True }

    #x0 = np.array([0.2]*5)
    #print f(x0)
    #print Df(x0)
    #print Hf(x0)


    sys.stdout.write('Solving the problem...\n\n')
    sys.stdout.flush()
    #res = optimize.minimize(f, x0, jac=Df, bounds=bounds, constraints=cons, options=opts)
    res = optimize.minimize(f, x0, jac=Df, constraints=cons, options=opts)


    print
    print 'Solver status: {0}'.format(res.success)
    print 'Minimizers:'
    print res.x
    print res.message
