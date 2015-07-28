import sys
import os
import datetime
import tempfile
import subprocess

b1, e1, s1 = 0.0, 0.1, 0.001
b2, e2, s2 = 0.0, 1.0, 0.01

outdir = 'visual/contour/test'
datadir = 'data'

# apparently over-smoothed?
#opt_h1 = [ 0.00844, 0.00928, 0.00565, 0.00579, 0.000435 ]
#opt_h2 = [ 0.0671, 0.0669, 0.0707, 0.0534, 0.0439 ]

# corrected by guessing...
opt_h1 = [ 0.00844, 0.00928, 0.0012, 0.00579, 0.000435 ]
opt_h2 = [ 0.0671, 0.0669, 0.018, 0.01, 0.01 ]

for i, (h1, h2) in enumerate(zip(opt_h1, opt_h2)):
    sys.stdout.write('Scanning evttype{0}...\n'.format(i+1))
    sys.stdout.flush()
    out2d = outdir + '/evttype{0}.2d.opt.txt'.format(i+1)
    out1d = outdir + '/evttype{0}.1d.opt.txt'.format(i+1)
    in_data = datadir + '/evttype{0}.cv.dat'.format(i+1)

    subprocess.check_call(['./kde_scan', out2d, out1d, in_data,
                           str(h1), str(h2),
                           str(b1), str(e1), str(s1),
                           str(b2), str(e2), str(s2)])
