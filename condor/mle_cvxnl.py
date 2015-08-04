#! /nfs/raid13/babar/software/anaconda/bin/python

import sys
import math
import numpy as np
import cvxopt as cvx
from iopro import genfromtxt

# Truth proportions:
# 0.004597168626601867
# 0.009011011183443553
# 0.41817119150694954
# 0.07374540140202024
# 0.49265798277959727
if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('nbag', type=int,
                        help='Number of trials to bag/fit.')
    parser.add_argument('sample_fname', type=str,
                        help='Path to the data sample cached kde score.')
    parser.add_argument('output_fname', type=str,
                        help='Path to store output.')
    parser.add_argument('--obj_scale', type=float, default=1e-8,
                        help='Scale factor to apply to objective function.')
    args = parser.parse_args()

    # Read in cached KDE evalutions of the data sample to fit
    p_raw = genfromtxt(args.sample_fname)
    N, D = p_raw.shape

    # Open the file to write results
    fout = open(args.output_fname, 'w')

    for i in range(args.nbag):

        sys.stdout.write('Bag iteration {0}:\n\n'.format(i+1))
        sys.stdout.flush()

        # Create a bagged sample.
        p = p_raw[np.random.choice(N, N)]

        # Specify the objective function and its derivatives
        def F(x=None, z=None):

            # Optionally scale the objective to avoid numerical difficulties
            s = args.obj_scale

            # Handle call signature F()
            if x is None:
                return 0, cvx.matrix(0.2, (D, 1))

            # Handle call signature F(x)
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

            # Handle call signature F(x, z)
            x_arr = np.array(x.trans()).reshape(-1)
            arg_m2 = np.power(arg, -2.0)
            H = np.zeros((D, D))
            for b in range(D):
                for c in range(D):
                    H[b, c] = np.sum(arg_m2 * p[:,b] * p[:,c]) * z[0]
            H = cvx.matrix(H) * s

            return f, Df, H

        # Specify the inequality constraints
        G = cvx.matrix(np.diag(-np.ones(D)))
        h = cvx.matrix(np.zeros(D))

        # Specify the equality constraints
        # a. Float all
        #A = cvx.matrix(np.ones(D), (1, D))
        #b = cvx.matrix(1.0)

        # b. Float all but continuum
        A = cvx.matrix(np.array([ [1.0,1.0,1.0,1.0,1.0],
                                [0.0,0.0,0.0,0.0,1.0] ]))
        b = cvx.matrix([1.0, 0.49265798277959727])

        # c. Float all but continuum, and also fix ratio between SL and Had.
        #A = cvx.matrix(np.array([ [1.0,1.0,1.0,1.0,1.0],
        #                          [0.0,0.0,1.0,-5.6705,0.0],
        #                          [0.0,0.0,0.0,0.0,1.0] ]))
        #b = cvx.matrix([1.0, 0.0, 0.49553516])

        # Specify required metadata
        dims = { 'l': D, 'q': [], 's': [] }

        # Solve
        sol = cvx.solvers.cp(F, G, h, dims, A, b)

        # Report results
        print
        print 'Solver status: {0}'.format(sol['status'])
        print 'Minimum: {0}'.format(F(sol['x'])[0][0])
        print 'Minimizers:', np.array(sol['x']).reshape(-1)
        print
        print
        sys.stdout.flush()

        fout.write(' '.join(map(str, np.array(sol['x']).reshape(-1).tolist())) + '\n')


    # Cleanup
    fout.close()
